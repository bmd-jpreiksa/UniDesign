#!/usr/bin/env python
"""Minimal example that parses a bundled PDB and prints its chain count."""

from __future__ import annotations

from pathlib import Path

from unidesign import Structure


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    pdb_path = repo_root / "example" / "MonomerDesign" / "1agy" / "1agy.pdb"

    # Data files ship with the package; set UNIDESIGN_LIBRARY_PATH only if you want to override them.
    structure = Structure.from_pdb(pdb_path)

    print(f"Parsed {pdb_path.name} from {pdb_path.parent}")
    print(f"Structure name: {structure.name}")
    print(f"Chain count: {structure.chain_count()}")
    print(f"Chains: {', '.join(structure.chain_ids())}")

    try:
        stability = structure.compute_stability()
    except RuntimeError as error:
        print(f"Stability computation unavailable: {error}")
    else:
        print(f"Total stability energy: {stability.total:.3f}")
        print("Key terms:")
        for term in ("reference_ALA", "intraR_vdwatt", "interS_vdwatt", "aapropensity", "ramachandran"):
            value = stability.terms.get(term)
            if value is not None:
                print(f"  {term:<20} {value:12.3f}")


if __name__ == "__main__":
    main()
