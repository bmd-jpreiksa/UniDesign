from __future__ import annotations

from pathlib import Path

from unidesign.api import Structure

# Example structures to compute stability for.
EXAMPLES = {
    "1igd": (
        Path(__file__).resolve().parents[1]
        / "example"
        / "MonomerDesign"
        / "1igd"
        / "1igd.pdb"
    ),
}


def main() -> None:
    for pdb_id, path in EXAMPLES.items():
        print(f"### Structure {pdb_id} ###")
        structure = Structure.from_pdb(path)
        result = structure.compute_stability()

        print(f"Total stability: {result.total:.6f}")
        print("Selected energy terms:")
        for name in [
            "interD_vdwatt",
            "interD_deslvP",
            "interD_deslvH",
            "interD_hbscbb_dis",
            "interD_hbscbb_the",
            "interD_hbscbb_phi",
            "aapropensity",
            "ramachandran",
            "dunbrack",
        ]:
            value = result.terms.get(name)
            if value is not None:
                print(f"  {name:20s} = {value: .6f}")
        print()


if __name__ == "__main__":
    main()
