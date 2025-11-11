from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, Set, Tuple

from unidesign.api import Ligand, Minimizer, Structure
from unidesign.api.design import DesignDomain, SiteSpec

BASE_DIR = Path(__file__).resolve().parent
INPUT_DIR = BASE_DIR / "mutation_scan_structures_277A"
OUTPUT_DIR = BASE_DIR / "mutation_scan_structures_277A_minimized"
PARAM_DIR = BASE_DIR / "ligand_params"

# Same domain as the scan: only residue 277A is mutable, rest repack/fixed.
DDOMAIN_WT: Dict[Tuple[int, str], str] = {
    (285, "A"): "R",
    (283, "A"): "E",
    (277, "A"): "K",
    (294, "C"): "R",
    (253, "C"): "E",
}


def _domain_for_mutation(structure: Structure, mut_aa: str, radius: float = 5.0) -> DesignDomain:
    design = [SiteSpec(chain="A", position=277, allowed=mut_aa)]
    repack = {
        (chain, position)
        for (position, chain), aa in DDOMAIN_WT.items()
        if (position, chain) != (277, "A")
    }

    targets = []
    for atom in structure.chains()[0].residues()[276].atoms():
        x, y, z = atom.handle.coords
        targets.append((x, y, z))

    for chain in structure.chains():
        chain_id = chain.name.strip()
        if chain_id != "A":
            continue
        for residue in chain.residues():
            if residue.position == 277:
                continue
            for atom in residue.atoms():
                ax, ay, az = atom.handle.coords
                if any(math.dist((ax, ay, az), ref) <= radius for ref in targets):
                    repack.add((chain_id, residue.position))
                    break

    repack_sites = [SiteSpec(chain=chain, position=pos) for chain, pos in sorted(repack)]
    return DesignDomain(design_sites=design, repack_sites=repack_sites)


def main() -> None:
    if not INPUT_DIR.exists():
        raise RuntimeError(f"Input directory {INPUT_DIR} does not exist. Run run_mutation_scan.py first.")

    OUTPUT_DIR.mkdir(exist_ok=True)

    ligand = Ligand.from_mol2(
        PARAM_DIR / "G1X_chr.mol2",
        PARAM_DIR / "g1x_param.prm",
        PARAM_DIR / "g1x_topo.inp",
    )

    for pdb_file in sorted(INPUT_DIR.glob("*.pdb")):
        label = pdb_file.stem
        out_path = OUTPUT_DIR / f"{label}.pdb"
        if out_path.exists():
            continue

        structure = Structure.from_pdb(pdb_file)
        mut_aa = label.split("_")[-1] if "_" in label else DDOMAIN_WT[(277, "A")]
        domain = _domain_for_mutation(structure, mut_aa)

        minimizer = Minimizer(
            structure,
            domain=domain,
            ligand=ligand,
            quiet=True,
            respect_design_types=False,
        )
        minimizer.run()

        if minimizer.minimized_structure is not None:
            minimizer.minimized_structure.handle.write_pdb(str(out_path))
        else:
            print(f"Failed to minimize {pdb_file.name}")


if __name__ == "__main__":
    main()
