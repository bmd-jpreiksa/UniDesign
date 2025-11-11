from __future__ import annotations

import tempfile
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .entities import (
    Structure,
    _AAPP_RELATIVE,
    _ATOM_PARAM_RELATIVE,
    _RAMA_RELATIVE,
    _TOPOLOGY_RELATIVE,
    _WEIGHT_RELATIVE,
    _resolve_data_file,
    _resolve_rotlib_file,
)
from .ligand import Ligand

_REPO_ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class SiteSpec:
    chain: str
    position: int
    allowed: Optional[str] = None


class DesignDomain:
    """Captures design/repack instructions usually stored in a RESFILE."""

    def __init__(
        self,
        design_sites: Optional[Iterable[SiteSpec]] = None,
        repack_sites: Optional[Iterable[SiteSpec]] = None,
        fixed_sites: Optional[Iterable[SiteSpec]] = None,
    ) -> None:
        self.design_sites: List[SiteSpec] = list(design_sites or [])
        self.repack_sites: List[SiteSpec] = list(repack_sites or [])
        self.fixed_sites: List[SiteSpec] = list(fixed_sites or [])

    @staticmethod
    def _format_site(site: SiteSpec) -> str:
        if site.allowed:
            return f"{site.chain:>1s} {site.position:4d}   {site.allowed}"
        return f"{site.chain:>1s} {site.position:4d}"

    def to_resfile(self) -> str:
        lines: List[str] = []
        if self.design_sites:
            lines.append("SITES_DESIGN_START")
            for site in self.design_sites:
                lines.append(self._format_site(site))
            lines.append("SITES_DESIGN_END")
        if self.repack_sites:
            lines.append("SITES_REPACK_START")
            for site in self.repack_sites:
                lines.append(self._format_site(site))
            lines.append("SITES_REPACK_END")
        if self.fixed_sites:
            lines.append("SITES_FIX_START")
            for site in self.fixed_sites:
                lines.append(self._format_site(site))
            lines.append("SITES_FIX_END")
        if not lines:
            return ""
        return "\n".join(lines) + "\n"

    @classmethod
    def from_resfile(cls, path: str | Path) -> "DesignDomain":
        path = Path(path)
        design: List[SiteSpec] = []
        repack: List[SiteSpec] = []
        fixed: List[SiteSpec] = []
        section: Optional[str] = None
        with path.open("r", encoding="utf-8") as handle:
            for raw in handle:
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                if line == "SITES_DESIGN_START":
                    section = "design"
                    continue
                if line == "SITES_DESIGN_END":
                    section = None
                    continue
                if line == "SITES_REPACK_START":
                    section = "repack"
                    continue
                if line == "SITES_REPACK_END":
                    section = None
                    continue
                if line == "SITES_FIX_START":
                    section = "fix"
                    continue
                if line == "SITES_FIX_END":
                    section = None
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                chain = parts[0]
                position = int(parts[1])
                allowed = parts[2] if len(parts) > 2 else None
                spec = SiteSpec(chain=chain, position=position, allowed=allowed)
                if section == "design":
                    design.append(spec)
                elif section == "repack":
                    repack.append(spec)
                elif section == "fix":
                    fixed.append(spec)
        return cls(design_sites=design, repack_sites=repack, fixed_sites=fixed)


class DesignProtein:
    """High-level helper around UniDesign monomer workflow."""

    def __init__(
        self,
        structure: Structure,
        domain: Optional[DesignDomain] = None,
        *,
        interface_only: bool = False,
        ligand: Optional[Ligand] = None,
        rotlib_file: Optional[str | Path] = None,
        quiet: bool = True,
    ) -> None:
        self._structure = structure
        self._domain = domain
        self._interface_only = interface_only
        self._ligand = ligand
        self._rotlib_file = rotlib_file
        self._quiet = quiet

        self.best_sequence: Dict[Tuple[str, int], str] | None = None
        self.best_sequence_string: str | None = None
        self.best_sequence_energy: Dict[str, float] | None = None
        self.best_structure: Optional[Structure] = None
        self.best_sites_structure: Optional[Structure] = None
        self.best_mutable_sites_structure: Optional[Structure] = None
        self.best_residue_energies: Dict[Tuple[str, int], Dict[str, float]] | None = None

    def run(self, *, trajectories: int = 1) -> None:
        resfile_text = ""
        if self._domain is not None:
            resfile_text = self._domain.to_resfile()

        with tempfile.TemporaryDirectory() as tempdir:
            workdir = Path(tempdir)

            with ExitStack() as stack:
                atom_params = _resolve_data_file(_ATOM_PARAM_RELATIVE, stack)
                topology = _resolve_data_file(_TOPOLOGY_RELATIVE, stack)
                weight_file = _resolve_data_file(_WEIGHT_RELATIVE, stack)
                aapp_file = _resolve_data_file(_AAPP_RELATIVE, stack)
                rama_file = _resolve_data_file(_RAMA_RELATIVE, stack)
                rotlib_file = _resolve_rotlib_file(stack, explicit=self._rotlib_file)

            if self._ligand:
                self._structure.handle.read_mol2(
                    str(self._ligand.mol2_path),
                    str(self._ligand.param_path),
                    str(self._ligand.topo_path),
                )
                if self._ligand.conformer_path:
                    self._structure.handle.read_small_mol_rotamers(str(self._ligand.conformer_path))

            options = {
                "program_path": str(_REPO_ROOT),
                "working_directory": tempdir,
                "atom_params": str(atom_params),
                "topology": str(topology),
                "weight_file": str(weight_file),
                "aapp_file": str(aapp_file),
                "rama_file": str(rama_file),
                "rotlib_bin": str(rotlib_file),
                "resfile_text": resfile_text,
                "trajectories": trajectories,
                "interface_only": self._interface_only,
                "use_input_sc": True,
                "rotate_hydroxyl": True,
                "exclude_cys_rotamers": False,
                "wildtype_only": False,
                "quiet": self._quiet,
            }
            if self._ligand:
                options["has_ligand"] = True
                options["ligand_mol2"] = str(self._ligand.mol2_path)
                options["ligand_params"] = str(self._ligand.param_path)
                options["ligand_topology"] = str(self._ligand.topo_path)
                if self._ligand.conformer_path:
                    options["ligand_conformers"] = str(self._ligand.conformer_path)

            payload = self._structure.handle.run_monomer_design(options)
            self._ingest_design_payload(payload)

    def _ingest_design_payload(self, payload: Dict[str, object]) -> None:
        sequence_string = str(payload["sequence_string"])
        trajectory_index = int(payload["trajectory_index"])
        seq_identity = float(payload["sequence_identity"])
        energy_total = float(payload["energy_total"])
        energy_evo = float(payload["energy_evolution"])
        energy_phy = float(payload["energy_physical"])
        energy_bind = float(payload["energy_binding"])
        unsatisfied = float(payload["unsatisfied_constraints"])

        self.best_sequence_string = sequence_string
        self.best_sequence_energy = {
            "trajectory": trajectory_index,
            "sequence_identity": seq_identity,
            "total": energy_total,
            "evolution": energy_evo,
            "physical": energy_phy,
            "binding": energy_bind,
            "unsatisfied_constraints": unsatisfied,
        }
        self.best_sequence = self._map_sequence(sequence_string)
        residue_energy_payload = payload.get("residue_self_energies")
        if residue_energy_payload is None:
            self.best_residue_energies = None
        else:
            energies: Dict[Tuple[str, int], Dict[str, float]] = {}
            for entry in residue_energy_payload:
                chain = str(entry["chain"]).strip()
                position = int(entry["position"])
                energies[(chain, position)] = {
                    "self_energy": float(entry["self_energy"]),
                    "binding_energy": float(entry["binding_energy"]),
                }
            self.best_residue_energies = energies

        self.best_structure = self._wrap_structure(payload.get("best_structure"))
        self.best_sites_structure = self._wrap_structure(payload.get("best_sites_structure"))
        self.best_mutable_sites_structure = self._wrap_structure(payload.get("best_mutable_sites_structure"))

    def _map_sequence(self, sequence_string: str) -> Dict[Tuple[str, int], str]:
        per_chain = sequence_string.split(";")
        sequence_map: Dict[Tuple[str, int], str] = {}
        for chain, seq in zip(self._structure.chains(), per_chain):
            residues = chain.residues()
            if len(seq) != len(residues):
                raise RuntimeError(
                    f"Sequence length {len(seq)} does not match residue count {len(residues)} on chain {chain.name}"
                )
            for residue, aa in zip(residues, seq):
                sequence_map[(chain.name, residue.position)] = aa
        return sequence_map

    @staticmethod
    def _wrap_structure(handle: object) -> Optional[Structure]:
        if handle is None:
            return None
        return Structure(handle=handle)
