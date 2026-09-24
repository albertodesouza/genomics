"""Tests the per-call deadline and retry wrapper around AlphaGenome predictions.

Regression cover for the failure that stalled the specificity-control window
build: `dna_client` issues its streaming RPC with no deadline and
`dna_client.create(timeout=...)` bounds only channel setup, so a half-open
connection never raises and the library's own `@retry_rpc` never fires. Every
worker slept on an open socket for three days at 45/1072 samples.

No `alphagenome` package, API key or network access is needed -- the client is a
stub, and only our wrapper is under test.
"""
import time

import pytest

from genomics.workflows.dataset_builders.non_longevous.build_window_and_predict import (
    ApiCallTimeout,
    _call_with_deadline,
    predict_sequence_resilient,
)


class _StubClient:
    """Records calls and replays a scripted sequence of outcomes.

    Each entry of ``outcomes`` is either an exception instance to raise, the
    string ``"hang"`` to sleep past any plausible deadline, or a value to
    return.
    """

    def __init__(self, outcomes):
        self.outcomes = list(outcomes)
        self.calls = 0

    def predict_sequence(self, sequence, *, requested_outputs, ontology_terms):
        self.calls += 1
        outcome = self.outcomes.pop(0)
        if outcome == "hang":
            time.sleep(30)
            return "unreachable"
        if isinstance(outcome, Exception):
            raise outcome
        return outcome


def _call(client_box, make_client, **kwargs):
    params = dict(
        seq="ACGT",
        requested_outputs=["RNA_SEQ"],
        ontology_terms=["CL:1000458"],
        timeout_s=0.25,
        max_attempts=3,
    )
    params.update(kwargs)
    return predict_sequence_resilient(client_box, make_client, **params)


def test_call_with_deadline_returns_value_when_fast():
    assert _call_with_deadline(lambda: "ok", 5.0) == "ok"


def test_call_with_deadline_raises_on_overrun():
    with pytest.raises(ApiCallTimeout):
        _call_with_deadline(lambda: time.sleep(5), 0.25)


def test_call_with_deadline_disabled_by_zero_timeout():
    # 0 means "no deadline": the call must still run and return normally.
    assert _call_with_deadline(lambda: "ok", 0) == "ok"


def test_call_with_deadline_restores_previous_handler():
    import signal

    before = signal.getsignal(signal.SIGALRM)
    _call_with_deadline(lambda: "ok", 5.0)
    assert signal.getsignal(signal.SIGALRM) is before
    # And the interval timer must be cleared, or a later unrelated call inherits it.
    assert signal.getitimer(signal.ITIMER_REAL)[0] == 0.0


def test_a_hang_is_bounded_and_retried_then_succeeds():
    """The exact production failure: first attempt hangs, retry succeeds."""
    client = _StubClient(["hang", "prediction"])
    rebuilt = []

    def make_client():
        rebuilt.append(True)
        return client

    started = time.time()
    result = _call([client], make_client)
    elapsed = time.time() - started

    assert result == "prediction"
    assert client.calls == 2
    # The deadline fired rather than the 30s sleep completing.
    assert elapsed < 10
    assert len(rebuilt) == 1, "a retry must rebuild the channel, not reuse a broken one"


def test_retry_swaps_in_the_rebuilt_client():
    """A half-open channel stays broken, so the retry must use the new client."""
    broken = _StubClient(["hang"])
    healthy = _StubClient(["prediction"])
    box = [broken]

    result = _call(box, lambda: healthy)

    assert result == "prediction"
    assert broken.calls == 1
    assert healthy.calls == 1
    assert box[0] is healthy


def test_gives_up_and_reraises_after_max_attempts():
    client = _StubClient([RuntimeError("boom")] * 3)

    with pytest.raises(RuntimeError, match="boom"):
        _call([client], lambda: client, max_attempts=3)

    assert client.calls == 3


def test_succeeds_without_retry_when_the_first_call_works():
    client = _StubClient(["prediction"])
    rebuilt = []

    result = _call([client], lambda: rebuilt.append(True) or client)

    assert result == "prediction"
    assert client.calls == 1
    assert rebuilt == [], "no rebuild should happen on the happy path"


def test_rebuild_failure_does_not_mask_the_retry():
    """If rebuilding the channel fails, the previous client is reused, not crashed on."""
    client = _StubClient([RuntimeError("transient"), "prediction"])

    def make_client():
        raise OSError("cannot reach the endpoint")

    assert _call([client], make_client) == "prediction"
    assert client.calls == 2


class _FakeStatus:
    def __init__(self, name):
        self.name = name


class _FakeRpcError(Exception):
    """Mimics grpc.RpcError's `.code()` accessor without importing grpc."""

    def __init__(self, status_name):
        super().__init__(status_name)
        self._status = _FakeStatus(status_name)

    def code(self):
        return self._status


def test_non_retryable_status_fails_on_the_first_attempt():
    """A bad API key must not cost max_attempts x backoff on every sample."""
    client = _StubClient([_FakeRpcError("INVALID_ARGUMENT")])

    with pytest.raises(_FakeRpcError):
        _call([client], lambda: client, max_attempts=5)

    assert client.calls == 1


@pytest.mark.parametrize("status", ["UNAUTHENTICATED", "PERMISSION_DENIED", "NOT_FOUND"])
def test_other_fatal_statuses_also_stop_immediately(status):
    client = _StubClient([_FakeRpcError(status)])

    with pytest.raises(_FakeRpcError):
        _call([client], lambda: client, max_attempts=5)

    assert client.calls == 1


def test_retryable_status_is_still_retried():
    client = _StubClient([_FakeRpcError("UNAVAILABLE"), "prediction"])

    assert _call([client], lambda: client) == "prediction"
    assert client.calls == 2


def test_plain_exceptions_without_a_code_are_retried():
    """`.code` is absent on non-gRPC errors; those must not be mistaken for fatal."""
    client = _StubClient([ValueError("transient"), "prediction"])

    assert _call([client], lambda: client) == "prediction"
    assert client.calls == 2
