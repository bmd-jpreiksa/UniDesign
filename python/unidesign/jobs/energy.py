"""Energy computation job wrappers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping

from ..config import ComputeBindingConfig, ComputeStabilityConfig
from ..runner import UniDesignRunner, UniDesignRunResult
from ._shared import JobArtifact, relocate_artifacts


@dataclass(slots=True)
class StabilityComputationResult:
    """Encapsulates outputs for ``ComputeStability`` runs."""

    run: UniDesignRunResult
    workspace: Path
    rotamer_list: JobArtifact | None
    cleanup: Callable[[], None] | None

    def close(self) -> None:
        if self.cleanup is not None:
            self.cleanup()
            self.cleanup = None


class StabilityComputationJob:
    """Execute the ``ComputeStability`` command."""

    def __init__(self, runner: UniDesignRunner, config: ComputeStabilityConfig) -> None:
        self._runner = runner
        self._config = config

    def run(
        self,
        *,
        keep_workspace: bool = False,
        env: Mapping[str, str] | None = None,
    ) -> StabilityComputationResult:
        run_result = self._runner.run(
            self._config.to_cli_args(), env=env, persist_workdir=True
        )
        candidates = {
            "rotamer_list": Path(f"{run_result.prefix}_rotlist.txt"),
        }
        workspace, artifacts, cleanup = relocate_artifacts(
            run_result.workdir, candidates, keep_workspace=keep_workspace
        )
        run_result.workdir = workspace

        return StabilityComputationResult(
            run=run_result,
            workspace=workspace,
            rotamer_list=artifacts.get("rotamer_list"),
            cleanup=cleanup,
        )


@dataclass(slots=True)
class BindingComputationResult:
    """Encapsulates outputs for ``ComputeBinding`` runs."""

    run: UniDesignRunResult
    workspace: Path
    cleanup: Callable[[], None] | None

    def close(self) -> None:
        if self.cleanup is not None:
            self.cleanup()
            self.cleanup = None


class BindingComputationJob:
    """Execute the ``ComputeBinding`` command."""

    def __init__(self, runner: UniDesignRunner, config: ComputeBindingConfig) -> None:
        self._runner = runner
        self._config = config

    def run(
        self,
        *,
        keep_workspace: bool = False,
        env: Mapping[str, str] | None = None,
    ) -> BindingComputationResult:
        run_result = self._runner.run(
            self._config.to_cli_args(), env=env, persist_workdir=True
        )
        workspace, _, cleanup = relocate_artifacts(
            run_result.workdir, {}, keep_workspace=keep_workspace
        )
        run_result.workdir = workspace
        return BindingComputationResult(
            run=run_result,
            workspace=workspace,
            cleanup=cleanup,
        )


__all__ = [
    "StabilityComputationJob",
    "StabilityComputationResult",
    "BindingComputationJob",
    "BindingComputationResult",
]
