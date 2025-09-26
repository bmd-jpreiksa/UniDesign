"""Utilities for discovering UniDesign assets relative to this package."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable, Tuple

from .exceptions import BinaryDiscoveryError


_PACKAGE_ROOT = Path(__file__).resolve().parent
_PROJECT_ROOT = _PACKAGE_ROOT.parent.parent


def project_root() -> Path:
    """Return the repository root that mirrors ``PROGRAM_PATH`` in ``Main.cpp``."""

    return _PROJECT_ROOT


def _resource_dir(name: str) -> Path:
    candidate = _PROJECT_ROOT / name
    if not candidate.is_dir():
        raise FileNotFoundError(
            f"Expected UniDesign resource directory '{name}' at {candidate!s}."
        )
    return candidate


def library_dir() -> Path:
    """Return the path to the UniDesign ``library`` directory."""

    return _resource_dir("library")


def wread_dir() -> Path:
    """Return the path to the UniDesign ``wread`` directory."""

    return _resource_dir("wread")


def extbin_dir() -> Path:
    """Return the path to the UniDesign ``extbin`` directory."""

    return _resource_dir("extbin")


def _candidate_binary_paths() -> Tuple[Path, ...]:
    search_roots: Iterable[Path] = (
        _PROJECT_ROOT,
        _PROJECT_ROOT / "build",
        _PROJECT_ROOT / "bin",
    )
    binary_names = ("UniDesign", "UniDesign.exe")
    candidates: list[Path] = []
    for root in search_roots:
        for name in binary_names:
            candidate = root / name
            if candidate not in candidates:
                candidates.append(candidate)
    return tuple(candidates)


def discover_binary() -> Path:
    """Locate the compiled UniDesign executable.

    The search mirrors the C++ entry point by resolving paths relative to this
    module. The helper validates executability and raises
    :class:`BinaryDiscoveryError` if discovery fails.
    """

    attempted: list[str] = []
    not_executable: list[str] = []
    for candidate in _candidate_binary_paths():
        attempted.append(str(candidate))
        if candidate.is_file():
            if os.access(candidate, os.X_OK):
                return candidate
            not_executable.append(str(candidate))
    if not_executable:
        raise BinaryDiscoveryError(
            "UniDesign binary exists but is not executable.",
            attempted_paths=tuple(not_executable),
        )
    raise BinaryDiscoveryError(
        "UniDesign binary could not be located. Build the project before proceeding.",
        attempted_paths=tuple(attempted),
    )


__all__ = [
    "project_root",
    "library_dir",
    "wread_dir",
    "extbin_dir",
    "discover_binary",
]
