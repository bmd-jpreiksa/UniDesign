
from pathlib import Path

from unidesign.api import DesignProtein, Ligand, Structure

BASE_DIR = Path(__file__).resolve().parent
PARAM_DIR = BASE_DIR / "ligand_params"

def main() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    receptor = Structure.from_pdb(BASE_DIR / "protein.pdb")

    ligand = Ligand.from_mol2(
        PARAM_DIR / "G1X_chr.mol2",
        PARAM_DIR / "g1x_param.prm",
        PARAM_DIR / "g1x_topo.inp",
    )

    designer = DesignProtein(
        receptor,
        interface_only=True,
        ligand=ligand,
        rotlib_file=repo_root / "library" / "rotlib" / "ALLbbdep.bin",
    )
    designer.run(trajectories=1)

    if designer.best_sequence_string and designer.best_sequence_energy:
        summary_path = BASE_DIR / "design_summary.txt"
        with summary_path.open("w", encoding="utf-8") as handle:
            handle.write("sequence\ttrajectory\tidentity\ttotal\tphysical\tbinding\n")
            handle.write(
                f"{designer.best_sequence_string}\t"
                f"{designer.best_sequence_energy['trajectory']}\t"
                f"{designer.best_sequence_energy['sequence_identity']:.6f}\t"
                f"{designer.best_sequence_energy['total']:.6f}\t"
                f"{designer.best_sequence_energy['physical']:.6f}\t"
                f"{designer.best_sequence_energy['binding']:.6f}\n"
            )

    if designer.best_structure is not None:
        best_struct_path = BASE_DIR / "design_best_structure.pdb"
        designer.best_structure.handle.write_pdb(str(best_struct_path))

    if designer.best_sites_structure is not None:
        best_sites_path = BASE_DIR / "design_best_sites.pdb"
        designer.best_sites_structure.handle.write_pdb(str(best_sites_path))

    if designer.best_mutable_sites_structure is not None:
        best_mut_path = BASE_DIR / "design_best_mutable_sites.pdb"
        designer.best_mutable_sites_structure.handle.write_pdb(str(best_mut_path))


if __name__ == "__main__":
    main()
