from __future__ import annotations

from pathlib import Path


class Ligand:
    def __init__(self, mol2_path: str | Path, param_path: str | Path, topo_path: str | Path, conformer_path: str | Path | None = None):
        self.mol2_path = Path(mol2_path)
        self.param_path = Path(param_path)
        self.topo_path = Path(topo_path)
        self.conformer_path = Path(conformer_path) if conformer_path else None

    @classmethod
    def from_mol2(cls, mol2_path: str | Path, param_path: str | Path, topo_path: str | Path, conformer_path: str | Path | None = None) -> "Ligand":
        return cls(mol2_path, param_path, topo_path, conformer_path)