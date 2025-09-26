"""Shared helpers for job wrappers."""

from __future__ import annotations

import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping

from ..artifacts import UniDesignArtifact


ArtifactFactory = Callable[[Path, str, str], UniDesignArtifact]


@dataclass(slots=True)
class ArtifactSpec:
    """Description of a potential UniDesign artifact."""

    relative_path: Path
    factory: ArtifactFactory

    @classmethod
    def from_type(
        cls, relative_path: Path, artifact_type: type[UniDesignArtifact]
    ) -> ArtifactSpec:
        """Convenience constructor for simple artifact factories."""

        def _factory(path: Path, prefix: str, logical_name: str) -> UniDesignArtifact:
            return artifact_type(path=path, prefix=prefix, logical_name=logical_name)

        return cls(relative_path=relative_path, factory=_factory)


def _existing_artifacts(
    workdir: Path, candidates: Mapping[str, ArtifactSpec]
) -> dict[str, tuple[Path, ArtifactSpec]]:
    existing: dict[str, tuple[Path, ArtifactSpec]] = {}
    for name, spec in candidates.items():
        candidate = workdir / spec.relative_path
        if candidate.exists():
            existing[name] = (candidate, spec)
    return existing


def relocate_artifacts(
    workdir: Path,
    candidates: Mapping[str, ArtifactSpec],
    *,
    keep_workspace: bool,
    prefix: str,
) -> tuple[Path, dict[str, UniDesignArtifact], Callable[[], None]]:
    """Relocate generated files based on caller preferences.

    Parameters
    ----------
    workdir:
        Temporary directory that holds the UniDesign execution artefacts.
    candidates:
        Mapping of artifact names to :class:`ArtifactSpec` entries describing
        potential files.
    keep_workspace:
        Whether the caller wants to retain the original ``workdir``.

    Returns
    -------
    workspace:
        Directory retained for downstream consumers. When ``keep_workspace`` is
        ``False`` a new directory is created and populated with only the
        generated artefacts.
    artifacts:
        Mapping of artifact names to concrete :class:`~unidesign.artifacts.UniDesignArtifact`
        instances.
    cleanup:
        Callable that removes ``workspace`` when invoked.
    """

    existing = _existing_artifacts(workdir, candidates)

    if keep_workspace:
        artifacts = {
            name: spec.factory(path=path, prefix=prefix, logical_name=name)
            for name, (path, spec) in existing.items()
        }
        cleanup = lambda: shutil.rmtree(workdir, ignore_errors=True)
        return workdir, artifacts, cleanup

    destination = Path(tempfile.mkdtemp(prefix="unidesign_artifacts_"))
    relocated: dict[str, UniDesignArtifact] = {}

    for name, (source, spec) in existing.items():
        try:
            relative = source.relative_to(workdir)
        except ValueError:
            relative = Path(source.name)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        relocated[name] = spec.factory(path=target, prefix=prefix, logical_name=name)

    shutil.rmtree(workdir, ignore_errors=True)
    cleanup = lambda: shutil.rmtree(destination, ignore_errors=True)
    return destination, relocated, cleanup


__all__ = ["ArtifactSpec", "relocate_artifacts"]
