"""Design-focused job wrappers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping

from ..artifacts import (
    DesignRotamerIndices,
    DesignSequenceSet,
    LigandPoseEnsemble,
    RotamerList,
    SelfEnergyReport,
    SiteSummary,
    StructureModel,
)
from ..config import ProteinDesignConfig
from ..runner import UniDesignRunner, UniDesignRunResult
from ._shared import ArtifactSpec, relocate_artifacts


@dataclass(slots=True)
class ProteinDesignResult:
    """Result bundle returned by :class:`ProteinDesignJob`."""

    run: UniDesignRunResult
    workspace: Path
    self_energy: SelfEnergyReport | None
    rotamer_list: RotamerList | None
    rotamer_list_secondary: RotamerList | None
    design_rotamer_indices: DesignRotamerIndices | None
    design_sequences: DesignSequenceSet | None
    best_sequences: DesignSequenceSet | None
    best_structure: StructureModel | None
    best_sites: SiteSummary | None
    best_mutation_sites: SiteSummary | None
    best_ligand_pose: LigandPoseEnsemble | None
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

    def _candidate_files(self, prefix: str) -> Mapping[str, ArtifactSpec]:
        return {
            "self_energy": ArtifactSpec.from_type(
                Path(f"{prefix}_selfenergy.txt"), SelfEnergyReport
            ),
            "rotamer_list": ArtifactSpec.from_type(
                Path(f"{prefix}_rotlist.txt"), RotamerList
            ),
            "rotamer_list_secondary": ArtifactSpec.from_type(
                Path(f"{prefix}_rotlistSEC.txt"), RotamerList
            ),
            "design_rotamer_indices": ArtifactSpec.from_type(
                Path(f"{prefix}_desrots"), DesignRotamerIndices
            ),
            "design_sequences": ArtifactSpec.from_type(
                Path(f"{prefix}_desseqs"), DesignSequenceSet
            ),
            "best_sequences": ArtifactSpec.from_type(
                Path(f"{prefix}_bestseqs"), DesignSequenceSet
            ),
            "best_structure": ArtifactSpec.from_type(
                Path(f"{prefix}_beststruct"), StructureModel
            ),
            "best_sites": ArtifactSpec.from_type(
                Path(f"{prefix}_bestsites"), SiteSummary
            ),
            "best_mutation_sites": ArtifactSpec.from_type(
                Path(f"{prefix}_bestmutsites"), SiteSummary
            ),
            "best_ligand_pose": ArtifactSpec.from_type(
                Path(f"{prefix}_bestlig"), LigandPoseEnsemble
            ),
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
            run_result.workdir,
            self._candidate_files(run_result.prefix),
            keep_workspace=keep_workspace,
            prefix=run_result.prefix,
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
