from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Tuple

from unidesign import ResidueDesignType
from unidesign.api import DesignProtein, Ligand, Structure
from unidesign.api.design import DesignDomain, SiteSpec

BASE_DIR = Path(__file__).resolve().parent
PARAM_DIR = BASE_DIR / "ligand_params"

# Canonical interface residues we keep fixed or repacked.
DDOMAIN_WT: Dict[Tuple[int, str], str] = {
    (285, "A"): "R",
    (283, "A"): "E",
    (277, "A"): "K",
    (294, "C"): "R",
    (253, "C"): "E",
}

# All single-letter amino acids we want to scan.
SCAN_AMINO_ACIDS: List[str] = [
    "A",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "K",
    "L",
    "M",
    "N",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "V",
    "W",
    "Y",
]

MUTATION_SITE = (277, "A")

ENERGY_KEYS = ("total", "physical", "binding")


def _find_residue(structure: Structure, chain_id: str, position: int):
    for chain in structure.chains():
        if chain.name.strip() == chain_id:
            for residue in chain.residues():
                if residue.position == position:
                    return residue
    raise ValueError(f"Residue {chain_id}:{position} not found in structure")


def _apply_domain_to_structure(structure: Structure, domain: DesignDomain) -> None:
    for chain in structure.chains():
        for residue in chain.residues():
            residue.design_type = ResidueDesignType.FIXED

    for site in domain.repack_sites:
        residue = _find_residue(structure, site.chain, site.position)
        residue.design_type = ResidueDesignType.REPACKABLE

    for site in domain.design_sites:
        residue = _find_residue(structure, site.chain, site.position)
        residue.design_type = ResidueDesignType.MUTABLE


def _build_design_domain(wild_type: Dict[Tuple[int, str], str],
                         mutant: Dict[Tuple[int, str], str]) -> DesignDomain:
    design_list: List[SiteSpec] = []
    repack_list: List[SiteSpec] = []

    for (position, chain), wt_aa in wild_type.items():
        mt_aa = mutant.get((position, chain), wt_aa)
        if mt_aa != wt_aa:
            design_list.append(SiteSpec(chain=chain, position=position, allowed=mt_aa))
        else:
            repack_list.append(SiteSpec(chain=chain, position=position))

    return DesignDomain(design_sites=design_list, repack_sites=repack_list)


def _restrict_domain_to_structure(structure: Structure, domain: DesignDomain) -> DesignDomain:
    allowed = {
        (site.chain, site.position)
        for site in list(domain.design_sites) + list(domain.repack_sites)
    }
    fixed_sites = []
    for chain in structure.chains():
        chain_id = chain.name.strip()
        for residue in chain.residues():
            key = (chain_id, residue.position)
            if key not in allowed:
                fixed_sites.append(SiteSpec(chain=chain_id, position=residue.position))
    return DesignDomain(
        design_sites=list(domain.design_sites),
        repack_sites=list(domain.repack_sites),
        fixed_sites=fixed_sites,
    )


def _run_design(domain: DesignDomain,
                receptor: Structure,
                ligand: Ligand,
                rotlib: Path,
                *,
                trajectories: int = 1,
                output_dir: Path | None = None,
                label: str | None = None) -> DesignProtein:
    working_structure = receptor.clone()
    _apply_domain_to_structure(working_structure, domain)

    # This is how we can load multiple ligands. However just one of them can be designable:
    working_structure.handle.read_mol2(
        str(ligand.mol2_path),
        str(ligand.param_path),
        str(ligand.topo_path),
    )

    designer = DesignProtein(
        working_structure,
        domain=domain,
        interface_only=False,
        ligand=None,
        rotlib_file=rotlib,
        quiet=True,
    )
    designer.run(trajectories=trajectories)

    if output_dir and label and designer.best_structure is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        designer.best_structure.handle.write_pdb(str(output_dir / f"{label}.pdb"))

    return designer


def _format_energy_line(label: str, energy: Dict[str, float]) -> str:
    return (
        f"{label:>4s}  "
        f"{energy['total']:12.4f}  "
        f"{energy['physical']:12.4f}  "
        f"{energy['binding']:12.4f}"
    )


def run_scan() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    receptor = Structure.from_pdb(BASE_DIR / "protein.pdb")

    ligand = Ligand.from_mol2(
        PARAM_DIR / "G1X_chr.mol2",
        PARAM_DIR / "g1x_param.prm",
        PARAM_DIR / "g1x_topo.inp",
    )
    rotlib = repo_root / "library" / "rotlib" / "ALLbbdep.bin"

    wt_domain = _restrict_domain_to_structure(
        receptor,
        _build_design_domain(DDOMAIN_WT, DDOMAIN_WT),
    )
    output_dir = BASE_DIR / f"mutation_scan_structures_{MUTATION_SITE[0]}{MUTATION_SITE[1]}"
    if (output_dir / "WT.pdb").exists():
        print(f"Structure directory {output_dir} already exists; skipping cleanup to preserve prior scans.")
    else:
        output_dir.mkdir(exist_ok=True)
    wt_designer = _run_design(wt_domain, receptor, ligand, rotlib, output_dir=output_dir, label="WT")
    wt_energy = wt_designer.best_sequence_energy or {}
    if not wt_energy:
        raise RuntimeError("Failed to obtain wild-type energy from DesignProtein")

    print("Wild-type energy terms:")
    print("label        total        physical        binding")
    print(_format_energy_line("WT", wt_energy))

    scan_results = {}
    raw_energies: List[Tuple[str, float]] = []
    for aa in SCAN_AMINO_ACIDS:
        if aa == DDOMAIN_WT[MUTATION_SITE]:
            # Skip redundant WT entry.
            continue
        print(f"\nScanning mutation {MUTATION_SITE} -> {aa}")
        mutant_map = dict(DDOMAIN_WT)
        mutant_map[MUTATION_SITE] = aa
        domain = _restrict_domain_to_structure(
            receptor,
            _build_design_domain(DDOMAIN_WT, mutant_map),
        )
        designer = _run_design(
            domain,
            receptor,
            ligand,
            rotlib,
            output_dir=output_dir,
            label=f"{MUTATION_SITE[0]}{MUTATION_SITE[1]}_{aa}",
        )
        energy = designer.best_sequence_energy or {}
        if not energy:
            print("  (run failed to produce energy; skipping)")
            continue
        delta = {key: energy[key] - wt_energy[key] for key in ENERGY_KEYS}
        scan_results[aa] = {
            "energy": {key: energy[key] for key in ENERGY_KEYS},
            "delta": delta,
        }
        raw_energies.append((aa, delta["total"]))
        print("  energies :", _format_energy_line(aa, energy))
        print(
            "  Δ energies:",
            f"{delta['total']:12.4f}  {delta['physical']:12.4f}  {delta['binding']:12.4f}",
        )

    # Persist results for downstream plotting/analysis.
    output_json = BASE_DIR / f"mutation_scan_{MUTATION_SITE[0]}{MUTATION_SITE[1]}.json"
    with output_json.open("w", encoding="utf-8") as handle:
        json.dump(scan_results, handle, indent=2)

    output_tsv = BASE_DIR / f"mutation_scan_{MUTATION_SITE[0]}{MUTATION_SITE[1]}.tsv"
    with output_tsv.open("w", encoding="utf-8") as handle:
        handle.write("mutation\ttotal\tdelta_total\tphysical\tdelta_physical\tbinding\tdelta_binding\n")
        for aa in SCAN_AMINO_ACIDS:
            if aa not in scan_results:
                continue
            energy = scan_results[aa]["energy"]
            delta = scan_results[aa]["delta"]
            handle.write(
                f"{aa}\t"
                f"{energy['total']:.6f}\t{delta['total']:.6f}\t"
                f"{energy['physical']:.6f}\t{delta['physical']:.6f}\t"
                f"{energy['binding']:.6f}\t{delta['binding']:.6f}\n"
            )

    # Visualize ΔE similar to the Rosetta example.
    if raw_energies:
        try:
            import matplotlib.pyplot as plt  # type: ignore
        except ModuleNotFoundError:
            print("matplotlib not installed; skipping bar-plot generation")
        else:
            labels, deltas = zip(*raw_energies)
            plt.style.use("seaborn-v0_8-darkgrid")
            plt.figure(figsize=(10, 4))
            plt.bar(labels, deltas, color="#2ca25f", edgecolor="black")
            plt.ylabel("ΔE (UniDesign Units)")
            plt.xlabel(f"Mutation at {MUTATION_SITE[0]}{MUTATION_SITE[1]}")
            plt.title("UniDesign Mutation Scan")
            plt.axhline(0.0, color="black", linewidth=0.8, linestyle="--")
            plt.xticks(rotation=45)
            plt.tight_layout()
            plot_path = BASE_DIR / f"mutation_scan_{MUTATION_SITE[0]}{MUTATION_SITE[1]}.png"
            plt.savefig(plot_path, dpi=200)
            print(f"\nPlot saved to {plot_path.name}")

    print(f"Scan complete. Results written to {output_json.name} and {output_tsv.name}")


if __name__ == "__main__":
    run_scan()
