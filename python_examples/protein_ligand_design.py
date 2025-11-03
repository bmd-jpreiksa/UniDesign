from __future__ import annotations

from pathlib import Path

from unidesign.api import DesignProtein, Ligand, Structure


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    example_dir = repo_root / "example" / "ProteinLigandInteraction" / "1r091_BS01_JEN"

    receptor = Structure.from_pdb(example_dir / "rec_native.pdb")

    ligand = Ligand.from_mol2(
        example_dir / "lig_charge_fixed.mol2",
        example_dir / "1r091_lig_param.prm",
        example_dir / "1r091_lig_topo.inp",
    )

    designer = DesignProtein(
        receptor,
        interface_only=True,
        ligand=ligand,
        rotlib_file=repo_root / "library" / "rotlib" / "ALLbbdep.bin",
    )
    designer.run(trajectories=1)
    print("Best sequence string:", designer.best_sequence_string)
    if designer.best_sequence_energy:
        for key, value in designer.best_sequence_energy.items():
            print(f"{key:>24s}: {value}")


if __name__ == "__main__":
    main()