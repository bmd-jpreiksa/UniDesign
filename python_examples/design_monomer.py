from __future__ import annotations

import os
from pathlib import Path

from unidesign import DesignDomain, DesignProtein, Structure


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    os.environ.setdefault("UNIDESIGN_LIBRARY_PATH", str(repo_root / "library"))
    pdb_path = repo_root / "example" / "MonomerDesign" / "1agy" / "1agy.pdb"
    resfile_path = repo_root / "example" / "MonomerDesign" / "1agy" / "RESFILE.txt"

    structure = Structure.from_pdb(pdb_path)
    domain = DesignDomain.from_resfile(resfile_path)

    designer = DesignProtein(structure, domain)
    designer.run()

    print("Total stability energy:", designer.best_sequence_energy["total"])
    print("Sequence identity:", designer.best_sequence_energy["sequence_identity"])
    print("First 30 residues of designed sequence:")
    seq_items = sorted(designer.best_sequence.items())[:30]
    for (chain_id, position), aa in seq_items:
        print(f"  {chain_id}:{position:4d} -> {aa}")


if __name__ == "__main__":
    main()
