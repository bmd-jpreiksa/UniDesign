from __future__ import annotations

import tempfile
from contextlib import ExitStack
from pathlib import Path
from typing import Optional, TYPE_CHECKING

from .entities import (
    Structure,
    _ATOM_PARAM_RELATIVE,
    _TOPOLOGY_RELATIVE,
    _WEIGHT_RELATIVE,
    _resolve_data_file,
    _resolve_rotlib_file,
)
from .ligand import Ligand

if TYPE_CHECKING:
    from .design import DesignDomain

_REPO_ROOT = Path(__file__).resolve().parents[2]


class Minimizer:
    """Convenience wrapper around UniDesign's side-chain minimization workflow."""

    def __init__(
        self,
        structure: Structure,
        *,
        domain: Optional["DesignDomain"] = None,
        ligand: Optional[Ligand] = None,
        rotlib_file: Optional[str | Path] = None,
        atom_params: Optional[str | Path] = None,
        topology: Optional[str | Path] = None,
        weight_file: Optional[str | Path] = None,
        use_input_sc: bool = True,
        rotate_hydroxyl: bool = True,
        quiet: bool = True,
        respect_design_types: bool = False,
    ) -> None:
        self._structure = structure
        self._domain = domain
        self._ligand = ligand
        self._rotlib_file = rotlib_file
        self._atom_params = atom_params
        self._topology = topology
        self._weight_file = weight_file
        self._use_input_sc = use_input_sc
        self._rotate_hydroxyl = rotate_hydroxyl
        self._quiet = quiet
        self._respect_design_types = respect_design_types

        self.minimized_structure: Optional[Structure] = None

    def run(self) -> None:
        with tempfile.TemporaryDirectory() as tempdir:
            workdir = Path(tempdir)

            with ExitStack() as stack:
                atom_param_path = _resolve_data_file(_ATOM_PARAM_RELATIVE, stack, self._atom_params)
                topology_path = _resolve_data_file(_TOPOLOGY_RELATIVE, stack, self._topology)
                rotlib_path = _resolve_rotlib_file(stack, explicit=self._rotlib_file)
                weight_path = _resolve_data_file(_WEIGHT_RELATIVE, stack, self._weight_file)

            if self._ligand:
                self._structure.handle.read_mol2(
                    str(self._ligand.mol2_path),
                    str(self._ligand.param_path),
                    str(self._ligand.topo_path),
                )

            options = {
                "program_path": str(_REPO_ROOT),
                "working_directory": str(workdir),
                "atom_params": str(atom_param_path),
                "topology": str(topology_path),
                "rotlib_bin": str(rotlib_path),
                "weight_file": str(weight_path),
                "use_input_sc": self._use_input_sc,
                "rotate_hydroxyl": self._rotate_hydroxyl,
                "quiet": self._quiet,
                "respect_design_types": self._respect_design_types,
            }
            if self._domain is not None:
                options["design_sites"] = [
                    {"chain": site.chain, "position": site.position, "allowed": site.allowed}
                    for site in self._domain.design_sites
                ]
                options["repack_sites"] = [
                    {"chain": site.chain, "position": site.position, "allowed": site.allowed}
                    for site in self._domain.repack_sites
                ]

            if self._ligand:
                options["has_ligand"] = True
                options["ligand_mol2"] = str(self._ligand.mol2_path)
                options["ligand_params"] = str(self._ligand.param_path)
                options["ligand_topology"] = str(self._ligand.topo_path)

            payload = self._structure.handle.run_minimization(options)

            handle = payload.get("structure")
            if handle is None:
                self.minimized_structure = None
            else:
                self.minimized_structure = Structure(handle=handle)
