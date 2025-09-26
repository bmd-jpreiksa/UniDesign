from importlib import resources
import math

import pytest

from unidesign import (
    StructureHandle,
    compute_structure_stability,
    load_energy_weights,
    load_structure_from_pdb,
    set_flag,
    set_path,
)


@pytest.fixture(scope="module")
def fixture_paths():
    base = resources.files("unidesign").joinpath("data/fixtures")
    return {
        "atom_params": base / "param_charmm19_lk.prm",
        "topology": base / "top_polh19.inp",
        "aapropensity": base / "aapropensity.nrg",
        "ramachandran": base / "ramachandran.nrg",
        "weights": base / "weight_all1.wgt",
        "pdb": base / "1igd.pdb",
    }


def test_compute_structure_stability_workflow(fixture_paths):
    set_path("FILE_ATOMPARAM", str(fixture_paths["atom_params"]))
    set_path("FILE_TOPO", str(fixture_paths["topology"]))
    set_path("FILE_AAPROPENSITY", str(fixture_paths["aapropensity"]))
    set_path("FILE_RAMACHANDRAN", str(fixture_paths["ramachandran"]))
    set_path("FILE_WEIGHT_READ", str(fixture_paths["weights"]))
    set_path("PDB", str(fixture_paths["pdb"]))
    set_flag("FLAG_PDB", True)

    structure = StructureHandle()
    load_structure_from_pdb(structure, str(fixture_paths["pdb"]))
    load_energy_weights()

    energies = compute_structure_stability(structure)

    assert len(energies) == 100
    assert all(math.isfinite(value) for value in energies)
    assert any(abs(value) > 1e-6 for value in energies[1:])
