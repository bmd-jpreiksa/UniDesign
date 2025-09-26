"""Ligand preparation job wrappers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping

from ..config import MakeLigParamConfig
from ..runner import UniDesignRunner, UniDesignRunResult
from ._shared import JobArtifact, relocate_artifacts


@dataclass(slots=True)
class LigandParameterizationResult:
    """Outputs for ``MakeLigParamAndTopo`` executions."""

    run: UniDesignRunResult
    workspace: Path
    parameter_file: JobArtifact | None
    topology_file: JobArtifact | None
    cleanup: Callable[[], None] | None

    def close(self) -> None:
        if self.cleanup is not None:
            self.cleanup()
            self.cleanup = None


class LigandParameterizationJob:
    """Execute the ``MakeLigParamAndTopo`` command."""

    def __init__(self, runner: UniDesignRunner, config: MakeLigParamConfig) -> None:
        self._runner = runner
        self._config = config

    def run(
        self,
        *,
        keep_workspace: bool = False,
        env: Mapping[str, str] | None = None,
    ) -> LigandParameterizationResult:
        run_result = self._runner.run(
            self._config.to_cli_args(), env=env, persist_workdir=True
        )
        candidates = {
            "parameter_file": Path(self._config.ligand_parameter_path),
            "topology_file": Path(self._config.ligand_topology_path),
        }
        workspace, artifacts, cleanup = relocate_artifacts(
            run_result.workdir, candidates, keep_workspace=keep_workspace
        )
        run_result.workdir = workspace
        return LigandParameterizationResult(
            run=run_result,
            workspace=workspace,
            parameter_file=artifacts.get("parameter_file"),
            topology_file=artifacts.get("topology_file"),
            cleanup=cleanup,
        )


__all__ = ["LigandParameterizationJob", "LigandParameterizationResult"]
