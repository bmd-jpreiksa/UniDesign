from __future__ import annotations

from pathlib import Path
from typing import Dict, Tuple

from unidesign import ResidueDesignType
from unidesign.api import DesignProtein, Ligand, Structure
from unidesign.api.design import DesignDomain, SiteSpec

BASE_DIR = Path(__file__).resolve().parent
PARAM_DIR = BASE_DIR / "ligand_params"

# Wild-type rotamer anchors for the active-site residues we want to keep rigid.
DDOMAIN_WT: Dict[Tuple[int, str], str] = {
    (161, "A"): "W",
    (314, "A"): "D",
    (203, "A"): "E",
    (205, "A"): "K",
    (80, "A"): "H",
    (282, "A"): "D",
}

# Mutant specification switches 161A from W -> S while keeping other sites fixed.
DDOMAIN_MT: Dict[Tuple[int, str], str] = {
    (161, "A"): "W",
    (314, "A"): "D",
    (203, "A"): "E",
    (205, "A"): "K",
    (80, "A"): "R",
    (282, "A"): "D",

}

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
    design_list: list[SiteSpec] = []
    repack_list: list[SiteSpec] = []

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


def _write_structures(label: str, designer: DesignProtein) -> None:
    if designer.best_structure is not None:
        best_struct_path = BASE_DIR / f"design_best_structure_{label}.pdb"
        designer.best_structure.handle.write_pdb(str(best_struct_path))

    if designer.best_sites_structure is not None:
        best_sites_path = BASE_DIR / f"design_best_sites_{label}.pdb"
        designer.best_sites_structure.handle.write_pdb(str(best_sites_path))

    if designer.best_mutable_sites_structure is not None:
        best_mut_path = BASE_DIR / f"design_best_mutable_sites_{label}.pdb"
        designer.best_mutable_sites_structure.handle.write_pdb(str(best_mut_path))


def _summarize_runs(wt: DesignProtein, mt: DesignProtein) -> None:
    if not (wt.best_sequence_energy and mt.best_sequence_energy):
        return

    summary_path = BASE_DIR / "design_summary.txt"
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("label\tsequence\ttrajectory\tidentity\t" +
                     "\t".join(ENERGY_KEYS) + "\n")
        for label, designer in (("WT", wt), ("MT", mt)):
            energy = designer.best_sequence_energy
            handle.write(
                f"{label}\t{designer.best_sequence_string}\t"
                f"{energy['trajectory']}\t"
                f"{energy['sequence_identity']:.6f}\t" +
                "\t".join(f"{energy[key]:.6f}" for key in ENERGY_KEYS) +
                "\n"
            )

        delta = {key: mt.best_sequence_energy[key] - wt.best_sequence_energy[key]
                 for key in ENERGY_KEYS}
        handle.write(
            "DELTA(MT-WT)\t-\t-\t-\t" +
            "\t".join(f"{delta[key]:.6f}" for key in ENERGY_KEYS) +
            "\n"
        )

    print("Delta energies (MT - WT):")
    for key in ENERGY_KEYS:
        print(f"  {key:>8s}: {delta[key]: .6f}")


def _write_residue_energy_summary(wt: DesignProtein, mt: DesignProtein) -> None:
    wt_map = wt.best_residue_energies or {}
    mt_map = mt.best_residue_energies or {}
    if not wt_map and not mt_map:
        print("Per-residue energies unavailable; ensure you're running the rebuilt extension.")
        return

    keys = sorted({*wt_map.keys(), *mt_map.keys()})
    output_path = BASE_DIR / "design_residue_energies.txt"
    def fmt(value: float | None) -> str:
        return f"{value:.6f}" if value is not None else "NA"
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("chain\tposition\twt_self\tmt_self\tdelta_self\twt_binding\tmt_binding\tdelta_binding\n")
        for key in keys:
            wt_entry = wt_map.get(key)
            mt_entry = mt_map.get(key)
            wt_self = wt_entry["self_energy"] if wt_entry else None
            mt_self = mt_entry["self_energy"] if mt_entry else None
            wt_bind = wt_entry["binding_energy"] if wt_entry else None
            mt_bind = mt_entry["binding_energy"] if mt_entry else None
            delta_self = (mt_self - wt_self) if wt_self is not None and mt_self is not None else None
            delta_bind = (mt_bind - wt_bind) if wt_bind is not None and mt_bind is not None else None
            handle.write(
                f"{key[0]}\t{key[1]}\t"
                f"{fmt(wt_self)}\t"
                f"{fmt(mt_self)}\t"
                f"{fmt(delta_self)}\t"
                f"{fmt(wt_bind)}\t"
                f"{fmt(mt_bind)}\t"
                f"{fmt(delta_bind)}\n"
            )

    print("Per-residue self-energy deltas (MT - WT):")
    for key in keys:
        wt_entry = wt_map.get(key)
        mt_entry = mt_map.get(key)
        if not wt_entry or not mt_entry:
            continue
        delta_self = mt_entry["self_energy"] - wt_entry["self_energy"]
        print(f"  {key[0]}:{key[1]:4d}  Δself = {delta_self: .4f}")


def _run_design(label: str, domain: DesignDomain, receptor: Structure, ligand: Ligand, rotlib: Path) -> DesignProtein:
    working_structure = receptor.clone()
    _apply_domain_to_structure(working_structure, domain)
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
    )
    designer.run(trajectories=1)
    _write_structures(label, designer)
    return designer


def run_example() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    receptor = Structure.from_pdb(BASE_DIR / "protein.pdb")

    ligand = Ligand.from_mol2(
        PARAM_DIR / "G1X_chr.mol2",
        PARAM_DIR / "g1x_param.prm",
        PARAM_DIR / "g1x_topo.inp",
    )
    rotlib = repo_root / "library" / "rotlib" / "ALLbbdep.bin"

    wt_domain = _build_design_domain(DDOMAIN_WT, DDOMAIN_WT)
    wt_domain = _restrict_domain_to_structure(receptor, wt_domain)
    mt_domain = _build_design_domain(DDOMAIN_WT, DDOMAIN_MT)
    mt_domain = _restrict_domain_to_structure(receptor, mt_domain)

    wt_designer = _run_design("wt", wt_domain, receptor, ligand, rotlib)
    mt_designer = _run_design("mt", mt_domain, receptor, ligand, rotlib)

    _summarize_runs(wt_designer, mt_designer)
    _write_residue_energy_summary(wt_designer, mt_designer)


if __name__ == "__main__":
    run_example()
