"""Background jobs with progress, deduplicated by key.

Long operations (cohort aggregates, population matrices, building annotations) run in a small
worker pool. Requests poll the same URL: while a job runs they get ``{"status": "pending",
"progress": ...}``; once done the result is served from the services' caches.
"""
from __future__ import annotations

import threading
import time
import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable, Dict, Optional

JOB_TTL_SECONDS = 600


class JobCancelled(Exception):
    """Raised inside a job when the user (or shutdown) cancelled it."""


class Progress:
    """Progress callback handed to job functions; also exposes cooperative cancellation.

    Calling it after cancellation raises :class:`JobCancelled`; worker loops can check
    ``progress.cancelled`` to skip remaining work quickly.
    """

    def __init__(self, job: "Job", base: float = 0.0, share: float = 1.0, prefix: str = ""):
        self.job = job
        self.base = base
        self.share = share
        self.prefix = prefix

    @property
    def cancelled(self) -> bool:
        return self.job.cancelled

    def __call__(self, progress: float, message: str) -> None:
        if self.job.cancelled:
            raise JobCancelled()
        self.job.update(self.base + self.share * progress, f"{self.prefix}{message}")

    def sub(self, base: float, share: float, prefix: str = "") -> "Progress":
        return Progress(self.job, self.base + self.share * base, self.share * share, self.prefix + prefix)


def is_cancelled(progress: Any) -> bool:
    return bool(getattr(progress, "cancelled", False))


class Job:
    def __init__(self, key: str, title: str):
        self.id = uuid.uuid4().hex[:12]
        self.key = key
        self.title = title
        self.status = "pending"
        self.progress = 0.0
        self.message = "Queued"
        self.result: Any = None
        self.error: Optional[str] = None
        self.created = time.time()
        self.finished: Optional[float] = None
        self.cancelled = False

    def update(self, progress: float, message: str) -> None:
        self.progress = max(0.0, min(1.0, float(progress)))
        self.message = message

    def as_dict(self) -> Dict[str, Any]:
        return {
            "id": self.id,
            "title": self.title,
            "status": self.status,
            "progress": round(self.progress, 4),
            "message": self.message,
            "error": self.error,
            "elapsed": round((self.finished or time.time()) - self.created, 2),
        }


class JobManager:
    def __init__(self, workers: int = 2):
        self._pool = ThreadPoolExecutor(max_workers=workers, thread_name_prefix="visualizer-job")
        self._jobs: Dict[str, Job] = {}
        self._by_id: Dict[str, Job] = {}
        self._lock = threading.Lock()

    def run(self, key: str, title: str, fn: Callable[[Callable[[float, str], None]], Any]) -> Job:
        """Return the job for ``key``, starting it if needed. Finished jobs are consumed once."""
        with self._lock:
            self._expire()
            job = self._jobs.get(key)
            if job is not None and job.status in ("pending", "running"):
                return job
            if job is not None and job.status == "done":
                self._jobs.pop(key, None)
                return job
            if job is not None and job.status == "error" and time.time() - (job.finished or 0) < 5:
                return job
            job = Job(key, title)
            self._jobs[key] = job
            self._by_id[job.id] = job

        def target() -> None:
            job.status = "running"
            job.message = "Starting"
            try:
                job.result = fn(Progress(job))
                if job.cancelled:
                    raise JobCancelled()
                job.status = "done"
                job.progress = 1.0
                job.message = "Done"
            except JobCancelled:
                job.status = "cancelled"
                job.message = "Cancelled"
            except Exception as exc:  # surfaced to the client
                job.status = "error"
                job.error = f"{type(exc).__name__}: {exc}"
                job.message = job.error
                traceback.print_exc()
            finally:
                job.finished = time.time()

        self._pool.submit(target)
        return job

    def get(self, job_id: str) -> Optional[Job]:
        return self._by_id.get(job_id)

    def cancel(self, job_id: str) -> Optional[Job]:
        job = self._by_id.get(job_id)
        if job is not None and job.status in ("pending", "running"):
            job.cancelled = True
            job.message = "Cancelling…"
            with self._lock:
                if self._jobs.get(job.key) is job:
                    self._jobs.pop(job.key, None)
        return job

    def cancel_all(self) -> None:
        for job_id in list(self._by_id):
            self.cancel(job_id)

    def active(self) -> list:
        with self._lock:
            return [job.as_dict() for job in self._by_id.values() if job.status in ("pending", "running")]

    def _expire(self) -> None:
        now = time.time()
        for job_id, job in list(self._by_id.items()):
            if job.finished and now - job.finished > JOB_TTL_SECONDS:
                self._by_id.pop(job_id, None)
                if self._jobs.get(job.key) is job:
                    self._jobs.pop(job.key, None)

    def shutdown(self) -> None:
        self.cancel_all()
        self._pool.shutdown(wait=False)
