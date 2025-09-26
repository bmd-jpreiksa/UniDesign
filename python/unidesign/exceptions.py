"""Custom exception types for the UniDesign Python helpers."""

from __future__ import annotations


class UniDesignError(RuntimeError):
    """Base class for exceptions raised by the UniDesign helpers."""


class BinaryDiscoveryError(UniDesignError):
    """Raised when the compiled UniDesign binary cannot be located or used."""

    def __init__(self, message: str, *, attempted_paths: tuple[str, ...] | None = None) -> None:
        super().__init__(message)
        self.attempted_paths = attempted_paths or ()

    def __str__(self) -> str:  # pragma: no cover - trivial wrapper for repr
        base = super().__str__()
        if self.attempted_paths:
            attempts = ", ".join(str(path) for path in self.attempted_paths)
            return f"{base} (checked: {attempts})"
        return base
