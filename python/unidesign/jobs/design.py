"""Design-focused job wrappers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping

from ..config import ProteinDesignConfig
from ..runner import UniDesignRunner, UniDesignRunResult
from ._shared import JobArtifact, relocate_artifacts


@dataclass(slots=True)
class ProteinDesignResult:
    """Result bundle returned by :class:`ProteinDesignJob`."""

    run: UniDesignRunResult
    workspace: Path
    self_energy: JobArtifact | None
    rotamer_list: JobArtifact | None
    rotamer_list_secondary: JobArtifact | None
    design_rotamer_indices: JobArtifact | None
    design_sequences: JobArtifact | None
    best_sequences: JobArtifact | None
    best_structure: JobArtifact | None
    best_sites: JobArtifact | None
    best_mutation_sites: JobArtifact | None
    best_ligand_pose: JobArtifact | None
    cleanup: Callable[[], None] | None

    def close(self) -> None:
        """Remove the workspace retained for this result."""

        if self.cleanup is not None:
            self.cleanup()
            self.cleanup = None


class ProteinDesignJob:
    """Execute ``ProteinDesign`` with structured inputs and outputs."""

    def __init__(self, runner: UniDesignRunner, config: ProteinDesignConfig) -> None:
        self._runner = runner
        self._config = config

    def _candidate_files(self, prefix: str) -> Mapping[str, Path]:
        return {
            "self_energy": Path(f"{prefix}_selfenergy.txt"),
            "rotamer_list": Path(f"{prefix}_rotlist.txt"),
            "rotamer_list_secondary": Path(f"{prefix}_rotlistSEC.txt"),
            "design_rotamer_indices": Path(f"{prefix}_desrots"),
            "design_sequences": Path(f"{prefix}_desseqs"),
            "best_sequences": Path(f"{prefix}_bestseqs"),
            "best_structure": Path(f"{prefix}_beststruct"),
            "best_sites": Path(f"{prefix}_bestsites"),
            "best_mutation_sites": Path(f"{prefix}_bestmutsites"),
            "best_ligand_pose": Path(f"{prefix}_bestlig"),
        }

    def run(
        self,
        *,
        keep_workspace: bool = False,
        env: Mapping[str, str] | None = None,
    ) -> ProteinDesignResult:
        """Execute the UniDesign ``ProteinDesign`` command."""

        run_result = self._runner.run(
            self._config.to_cli_args(), env=env, persist_workdir=True
        )
        workspace, artifacts, cleanup = relocate_artifacts(
            run_result.workdir, self._candidate_files(run_result.prefix), keep_workspace=keep_workspace
        )
        run_result.workdir = workspace

        return ProteinDesignResult(
            run=run_result,
            workspace=workspace,
            self_energy=artifacts.get("self_energy"),
            rotamer_list=artifacts.get("rotamer_list"),
            rotamer_list_secondary=artifacts.get("rotamer_list_secondary"),
            design_rotamer_indices=artifacts.get("design_rotamer_indices"),
            design_sequences=artifacts.get("design_sequences"),
            best_sequences=artifacts.get("best_sequences"),
            best_structure=artifacts.get("best_structure"),
            best_sites=artifacts.get("best_sites"),
            best_mutation_sites=artifacts.get("best_mutation_sites"),
            best_ligand_pose=artifacts.get("best_ligand_pose"),
            cleanup=cleanup,
        )


__all__ = ["ProteinDesignJob", "ProteinDesignResult"]
