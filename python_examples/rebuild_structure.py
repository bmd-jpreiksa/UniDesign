#!/usr/bin/env python
"""Rebuild a UniDesign structure using Python objects only."""

from __future__ import annotations

from math import isclose
from pathlib import Path

from unidesign import Chain, Residue, Structure

DATA_PDB = Path("example/MonomerDesign/1agy/1agy.pdb")


def rebuild_incremental(source: Structure) -> Structure:
    """
    Rebuild a structure by iterating through chains/residues and copying their native state.

    All heavy lifting (bonds, rotamer metadata, etc.) is handled by the native copy helpers.
    """
    rebuilt = Structure()
    try:
        rebuilt.name = source.name
    except RuntimeError:
        pass
    for chain in source.chains():
        new_chain = Chain()
        new_chain.copy_from(chain)
        rebuilt.add_chain(new_chain)
    return rebuilt


def main() -> None:
    if not hasattr(Residue().handle, "copy_from") or not hasattr(Chain().handle, "copy_from"):
        raise RuntimeError("This example requires the rebuilt UniDesign extension with copy helpers.")

    source = Structure.from_pdb(DATA_PDB)
    source_energy = source.compute_stability()

    # Use fresh copies of the parsed structure so we don't reuse residues whose metadata
    # (e.g. cached propensity scores) was mutated by previous energy evaluations.
    clone_source = Structure.from_pdb(DATA_PDB)
    clone_structure = clone_source.clone()
    clone_energy = clone_structure.compute_stability()

    incremental_source = Structure.from_pdb(DATA_PDB)
    incremental_structure = rebuild_incremental(incremental_source)
    incremental_energy = incremental_structure.compute_stability()

    print("Original total:", source_energy.total)
    print("Clone total   :", clone_energy.total)
    print("Rebuilt total :", incremental_energy.total)

    assert isclose(source_energy.total, clone_energy.total, rel_tol=1e-9, abs_tol=1e-6)
    assert isclose(source_energy.total, incremental_energy.total, rel_tol=1e-9, abs_tol=1e-6)

    print("Structures rebuilt successfully with matching energies.")


if __name__ == "__main__":
    main()
