"""Shared helpers for job wrappers."""

from __future__ import annotations

import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping


@dataclass(slots=True)
class JobArtifact:
    """Representation of a generated UniDesign file."""

    path: Path

    def read_text(self, encoding: str = "utf-8") -> str:
        """Read the artifact as text using the provided encoding."""

        return self.path.read_text(encoding=encoding)

    def read_bytes(self) -> bytes:
        """Read the artifact in binary mode."""

        return self.path.read_bytes()

    def open(self, mode: str = "r", *args, **kwargs):
        """Open the underlying file handle."""

        return self.path.open(mode, *args, **kwargs)


def _existing_artifacts(workdir: Path, candidates: Mapping[str, Path]) -> dict[str, Path]:
    existing: dict[str, Path] = {}
    for name, relative_path in candidates.items():
        candidate = workdir / relative_path
        if candidate.exists():
            existing[name] = candidate
    return existing


def relocate_artifacts(
    workdir: Path,
    candidates: Mapping[str, Path],
    *,
    keep_workspace: bool,
) -> tuple[Path, dict[str, JobArtifact], Callable[[], None]]:
    """Relocate generated files based on caller preferences.

    Parameters
    ----------
    workdir:
        Temporary directory that holds the UniDesign execution artefacts.
    candidates:
        Mapping of artifact names to paths relative to ``workdir``.
    keep_workspace:
        Whether the caller wants to retain the original ``workdir``.

    Returns
    -------
    workspace:
        Directory retained for downstream consumers. When ``keep_workspace`` is
        ``False`` a new directory is created and populated with only the
        generated artefacts.
    artifacts:
        Mapping of artifact names to :class:`JobArtifact` instances.
    cleanup:
        Callable that removes ``workspace`` when invoked.
    """

    existing = _existing_artifacts(workdir, candidates)

    if keep_workspace:
        artifacts = {name: JobArtifact(path=path) for name, path in existing.items()}
        cleanup = lambda: shutil.rmtree(workdir, ignore_errors=True)
        return workdir, artifacts, cleanup

    destination = Path(tempfile.mkdtemp(prefix="unidesign_artifacts_"))
    relocated: dict[str, JobArtifact] = {}

    for name, source in existing.items():
        try:
            relative = source.relative_to(workdir)
        except ValueError:
            relative = Path(source.name)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        relocated[name] = JobArtifact(path=target)

    shutil.rmtree(workdir, ignore_errors=True)
    cleanup = lambda: shutil.rmtree(destination, ignore_errors=True)
    return destination, relocated, cleanup


__all__ = ["JobArtifact", "relocate_artifacts"]
