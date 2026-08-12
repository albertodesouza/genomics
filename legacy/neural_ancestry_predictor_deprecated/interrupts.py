from __future__ import annotations


class InterruptState:
    """Estado compartilhado para interrupção de treinamento."""

    interrupted = False


interrupt_state = InterruptState()
