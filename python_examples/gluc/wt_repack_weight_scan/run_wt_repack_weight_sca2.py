from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

from unidesign import ResidueDesignType
from unidesign.api import DesignProtein, Ligand, Structure
from unidesign.api.design import DesignDomain, SiteSpec

BASE_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[3]
TEST_DIR = BASE_DIR / "wt_repack_weight_scan"
WEIGHT_DIR = TEST_DIR / "weights"
STRUCT_DIR = TEST_DIR / "structures"
SUMMARY_JSON = TEST_DIR / "wt_repack_weight_summary.json"
SUMMARY_TSV = TEST_DIR / "wt_repack_weight_summary.tsv"
PARAM_DEBUG_DIR = TEST_DIR / "parameter_debug"

PARAM_DIR = BASE_DIR / "ligand_params"
ROT_LIB = REPO_ROOT / "library" / "rotlib" / "ALLbbdep.bin"
BASE_WEIGHT_FILE = REPO_ROOT / "wread" / "weight_all1.wgt"
ALT_WEIGHT_FILES = {
    "weight_all2": REPO_ROOT / "wread" / "weight_all2.wgt",
    "weight_all3": REPO_ROOT / "wread" / "weight_all3.wgt",
}
ROTAMER_PRUNE_CUTOFF = 0.03  # default UniDesign rotamer probability threshold

# Same canonical domain as run_mutation_scan_close_pocket.py, but strictly repack WT.
DDOMAIN_WT: Dict[Tuple[int, str], str] = {
    (45, "A"): "F",
    (45, "B"): "F",
    (47, "A"): "D",
    (47, "B"): "D",
    (48, "A"): "R",
    (48, "B"): "R",
    (54, "A"): "Y",
    (54, "B"): "Y",
    (316, "A"): "F",
    (316, "B"): "F",
    (318, "A"): "Y",
    (318, "B"): "Y",
    (319, "A"): "R",
    (319, "B"): "R",
}

# Variant label -> configuration for producing weight files.
WEIGHT_VARIANTS: Dict[str, Dict[str, object]] = {
    "baseline": {"source": "weight_all1", "overrides": {}},
    "repel_d": {
        "source": "weight_all1",
        "overrides": {
            "interS_vdwrep": 1.500,
            "interD_vdwrep": 1.500,
            "intraR_vdwrep": 0.060,
        },
    },
}


def _ensure_dirs() -> None:
    WEIGHT_DIR.mkdir(parents=True, exist_ok=True)
    STRUCT_DIR.mkdir(parents=True, exist_ok=True)


def _read_weight_lines(source: str) -> List[str]:
    if source == "weight_all1":
        path = BASE_WEIGHT_FILE
    else:
        path = ALT_WEIGHT_FILES.get(source)
        if path is None:
            raise ValueError(f"Unknown weight source '{source}'")
    return path.read_text().splitlines()


def _write_weight_file(label: str, overrides: Dict[str, float], base_lines: Iterable[str]) -> Path:
    output_lines: List[str] = []
    override_keys = {key.lower(): value for key, value in overrides.items()}
    for line in base_lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            output_lines.append(line)
            continue
        term = stripped.split()[0]
        key = term.lower()
        if key in override_keys:
            value = override_keys[key]
            output_lines.append(f"{term:<25s} {value:>10.3f}")
        else:
            output_lines.append(line)
    out_path = WEIGHT_DIR / f"weight_all1_{label}.wgt"
    out_path.write_text("\n".join(output_lines) + "\n")
    return out_path


def _find_residue(structure: Structure, chain_id: str, position: int):
    for chain in structure.chains():
        if chain.name.strip() == chain_id:
            for residue in chain.residues():
                if residue.position == position:
                    return residue
    raise ValueError(f"Residue {chain_id}:{position} not found in structure")


def _apply_design_domain(structure: Structure, domain: DesignDomain) -> None:
    for chain in structure.chains():
        for residue in chain.residues():
            residue.design_type = ResidueDesignType.FIXED

    for site in domain.repack_sites:
        residue = _find_residue(structure, site.chain, site.position)
        residue.design_type = ResidueDesignType.REPACKABLE

    for site in domain.design_sites:
        residue = _find_residue(structure, site.chain, site.position)
        residue.design_type = ResidueDesignType.MUTABLE


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


def _build_wt_fixed_domain(structure: Structure) -> DesignDomain:
    design_sites = [
        SiteSpec(chain=chain, position=position, allowed=amino)
        for (position, chain), amino in DDOMAIN_WT.items()
    ]
    return _restrict_domain_to_structure(structure, DesignDomain(design_sites=design_sites))


def run_weight_scan() -> None:
    _ensure_dirs()
    structure = Structure.from_pdb(BASE_DIR / "protein.pdb")
    domain = _build_wt_fixed_domain(structure)

    ligand = Ligand.from_mol2(
        PARAM_DIR / "G1X_chr.mol2",
        PARAM_DIR / "g1x_param.prm",
        PARAM_DIR / "g1x_topo.inp",
    )

    cached_bases: Dict[str, List[str]] = {}
    summary: Dict[str, Dict[str, float]] = {}
    parameter_debug: Dict[str, Dict[str, object]] = {}

    for label, config in WEIGHT_VARIANTS.items():
        print(f"Running WT repack with weight variant '{label}'")
        source = str(config.get("source", "weight_all1"))
        overrides = dict(config.get("overrides", {}))
        if source not in cached_bases:
            cached_bases[source] = _read_weight_lines(source)
        base_lines = cached_bases[source]
        weight_path = _write_weight_file(label, overrides, base_lines)

        working_structure = structure.clone()
        _apply_design_domain(working_structure, domain)
        working_structure.handle.read_mol2(
            str(ligand.mol2_path),
            str(ligand.param_path),
            str(ligand.topo_path),
        )

        designer = DesignProtein(
            working_structure,
            domain=domain,
            interface_only=False,
            ligand=ligand,
            rotlib_file=ROT_LIB,
            weight_file=weight_path,
            rotamer_probability_cutoff=ROTAMER_PRUNE_CUTOFF,
            quiet=True,
        )
        designer.run(trajectories=1)

        if designer.best_structure is not None:
            out_path = STRUCT_DIR / f"best_structure_{label}.pdb"
            designer.best_structure.handle.write_pdb(str(out_path))

        if designer.best_sequence_energy:
            summary[label] = {
                key: designer.best_sequence_energy[key]
                for key in ("total", "physical", "binding", "trajectory")
            }
            parameter_debug[label] = {
                "before": designer.best_sequence_energy.get("energy_terms_initial") or {},
                "after": designer.best_sequence_energy.get("energy_terms_final") or {},
                "residue_energies": designer.best_residue_energies or {},
            }

    if summary:
        SUMMARY_JSON.write_text(json.dumps(summary, indent=2))
        with SUMMARY_TSV.open("w", encoding="utf-8") as handle:
            handle.write("label\ttotal\tphysical\tbinding\ttrajectory\n")
            for label in WEIGHT_VARIANTS:
                data = summary.get(label)
                if not data:
                    continue
                handle.write(
                    f"{label}\t{data['total']:.6f}\t{data['physical']:.6f}"
                    f"\t{data['binding']:.6f}\t{int(data['trajectory'])}\n"
                )

    if parameter_debug:
        PARAM_DEBUG_DIR.mkdir(exist_ok=True)
        (PARAM_DEBUG_DIR / "parameter_debug.json").write_text(json.dumps(parameter_debug, indent=2))
        for label in WEIGHT_VARIANTS:
            entries = parameter_debug.get(label)
            if not entries:
                continue
            before = dict(entries.get("before") or {})
            after = dict(entries.get("after") or {})
            terms = sorted(set(before) | set(after))
            output = PARAM_DEBUG_DIR / f"parameter_terms_{label}.tsv"
            with output.open("w", encoding="utf-8") as handle:
                handle.write("term\tbefore\tafter\tdelta\n")
                for term in terms:
                    b = before.get(term)
                    a = after.get(term)
                    delta = (a - b) if (a is not None and b is not None) else None
                    handle.write(
                        f"{term}\t"
                        f"{'' if b is None else f'{b:.6f}'}\t"
                        f"{'' if a is None else f'{a:.6f}'}\t"
                        f"{'' if delta is None else f'{delta:.6f}'}\n"
                    )

            residues = entries.get("residue_energies") or {}
            if residues:
                res_out = PARAM_DEBUG_DIR / f"residue_energies_{label}.tsv"
                with res_out.open("w", encoding="utf-8") as handle:
                    handle.write("chain\tposition\tself_energy\tbinding_energy\n")
                    for (chain, position), values in sorted(residues.items()):
                        handle.write(
                            f"{chain}\t{position}\t{values['self_energy']:.6f}\t{values['binding_energy']:.6f}\n"
                        )


if __name__ == "__main__":
    run_weight_scan()
