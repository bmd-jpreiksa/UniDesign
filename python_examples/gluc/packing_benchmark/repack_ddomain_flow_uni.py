#!/usr/bin/env python3
"""Repack Rosetta design domains using UniDesign but seed from FlowPacker poses."""
from __future__ import annotations

import argparse
import csv
import pickle
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import numpy as np

from unidesign import ResidueDesignType
from unidesign.api import DesignProtein, Structure
from unidesign.api.design import DesignDomain, SiteSpec


REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_ROTLIB = REPO_ROOT / "library" / "rotlib" / "ALLbbdep.bin"

DDomainEntry = Dict[Tuple[int, str], str]
DDomain = Dict[Tuple[int, str], DDomainEntry]


@dataclass(frozen=True)
class HeavyAtomBundle:
    positions: np.ndarray
    count: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Repack design domains via UniDesign while using FlowPacker structures "
            "as the starting pose for each domain."
        ),
    )
    parser.add_argument(
        "--flow-dir",
        default="flowpacker_init/wt_domain_repacking_flowpacker",
        help="Directory containing FlowPacker outputs named flow_repacked_<site>.pdb.",
    )
    parser.add_argument(
        "--ddomain",
        default="ddomain.pickle",
        help="Pickled Rosetta design-domain dictionary (default: %(default)s).",
    )
    parser.add_argument(
        "--structure",
        default="3m0m_minimized_Hs_only.pdb",
        help="Reference WT structure used for RMSD calculations (default: %(default)s).",
    )
    parser.add_argument(
        "--rotlib",
        default=str(DEFAULT_ROTLIB),
        help="Rotamer library passed to UniDesign (default: %(default)s).",
    )
    parser.add_argument(
        "--trajectories",
        type=int,
        default=1,
        help="Number of UniDesign trajectories per domain (default: %(default)s).",
    )
    parser.add_argument(
        "--sites",
        nargs="*",
        help="Optional subset of design domains to process (e.g. 344A 390A).",
    )
    parser.add_argument(
        "--output-dir",
        default="wt_domain_repacking_flow_uni",
        help="Directory storing repacked PDBs and RMSD CSV (default: %(default)s).",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Silence UniDesign logging (stderr still shows per-domain summary).",
    )
    return parser.parse_args()


def load_ddomain(path: Path) -> DDomain:
    with path.open("rb") as handle:
        ddomain = pickle.load(handle)
    if not isinstance(ddomain, dict):
        raise ValueError(f"Expected ddomain dict, got {type(ddomain)}")
    return ddomain  # type: ignore[return-value]


def _resolve_with_fallback(path_arg: str) -> Path:
    path = Path(path_arg).expanduser()
    if path.is_absolute() and path.exists():
        return path
    if path.exists():
        return path.resolve()
    candidate = SCRIPT_DIR / path
    if candidate.exists():
        return candidate.resolve()
    raise FileNotFoundError(path if path.is_absolute() else candidate)


def _to_site_specs(entry: Mapping[Tuple[int, str], str]) -> List[SiteSpec]:
    return [SiteSpec(chain=chain, position=position) for (position, chain) in entry.keys()]


def _find_residue(structure: Structure, chain_id: str, position: int):
    for chain in structure.chains():
        name = chain.name.strip() or "_"
        if name == chain_id:
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


def _restrict_domain_to_structure(structure: Structure, domain: DesignDomain) -> DesignDomain:
    allowed = {
        (site.chain, site.position)
        for site in list(domain.design_sites) + list(domain.repack_sites) + list(domain.fixed_sites)
    }
    fixed_sites = []
    for chain in structure.chains():
        chain_id = chain.name.strip() or "_"
        for residue in chain.residues():
            key = (chain_id, residue.position)
            if key not in allowed:
                fixed_sites.append(SiteSpec(chain=chain_id, position=residue.position))
    return DesignDomain(
        design_sites=list(domain.design_sites),
        repack_sites=list(domain.repack_sites),
        fixed_sites=fixed_sites,
    )


def _run_repack(
    base_structure: Structure,
    domain: DesignDomain,
    rotlib: str | Path,
    trajectories: int,
    quiet: bool,
) -> DesignProtein:
    working_structure = base_structure.clone()
    _apply_domain_to_structure(working_structure, domain)
    designer = DesignProtein(
        working_structure,
        domain=domain,
        interface_only=False,
        rotlib_file=rotlib,
        quiet=quiet,
    )
    designer.run(trajectories=trajectories)
    if designer.best_structure is None:
        raise RuntimeError("UniDesign run did not produce a best_structure payload")
    return designer


def _is_hydrogen(atom_name: str, element: str) -> bool:
    if element:
        return element.upper() == "H"
    return atom_name.strip().upper().startswith("H")


def _heavy_atom_lookup(pdb_path: Path) -> Dict[Tuple[str, int], Dict[str, np.ndarray]]:
    residues: Dict[Tuple[str, int], Dict[str, np.ndarray]] = {}
    with pdb_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_name = line[12:16].strip()
            element = line[76:78].strip()
            if _is_hydrogen(atom_name, element):
                continue
            chain = line[21].strip() or "_"
            try:
                position = int(line[22:26])
            except ValueError:
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            residues.setdefault((chain, position), {})[atom_name] = np.array([x, y, z], dtype=float)
    return residues


def _collect_heavy_atoms(
    atoms: Mapping[Tuple[str, int], Mapping[str, np.ndarray]],
    residues: Sequence[Tuple[int, str]],
) -> HeavyAtomBundle:
    coords: List[np.ndarray] = []
    for position, chain in residues:
        chain_id = chain.strip() or "_"
        residue_atoms = atoms.get((chain_id, position))
        if residue_atoms is None:
            raise KeyError(f"Residue {chain_id}:{position} missing from PDB data")
        for atom_name in sorted(residue_atoms):
            coords.append(residue_atoms[atom_name])
    if not coords:
        return HeavyAtomBundle(positions=np.empty((0, 3)), count=0)
    stacked = np.vstack(coords)
    return HeavyAtomBundle(positions=stacked, count=stacked.shape[0])


def _kabsch_rmsd(reference: np.ndarray, mobile: np.ndarray) -> float:
    if reference.shape != mobile.shape or reference.size == 0:
        return 0.0
    ref_center = reference.mean(axis=0)
    mob_center = mobile.mean(axis=0)
    ref_shifted = reference - ref_center
    mob_shifted = mobile - mob_center
    covariance = mob_shifted.T @ ref_shifted
    v, _, w = np.linalg.svd(covariance)
    d = np.sign(np.linalg.det(v @ w))
    rotation = v @ np.diag([1.0, 1.0, d]) @ w
    aligned = mob_shifted @ rotation
    diff = ref_shifted - aligned
    return float(np.sqrt((diff * diff).sum() / reference.shape[0]))


def write_results_csv(results: List[Dict[str, object]], path: Path) -> None:
    if not results:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "site",
                "residue_count",
                "heavy_atom_count",
                "heavy_atom_rmsd",
                "structure_path",
                "flowpacker_input",
                "total_energy",
            ],
        )
        writer.writeheader()
        writer.writerows(results)


def main() -> None:
    args = parse_args()
    start = time.perf_counter()
    flow_dir = _resolve_with_fallback(args.flow_dir)
    if not flow_dir.is_dir():
        raise NotADirectoryError(flow_dir)
    ddomain = load_ddomain(_resolve_with_fallback(args.ddomain))
    reference_path = _resolve_with_fallback(args.structure)
    reference_atoms = _heavy_atom_lookup(reference_path)
    all_sites = sorted(ddomain.keys())
    if args.sites:
        requested = {site.upper() for site in args.sites}
        all_sites = [site for site in all_sites if f"{site[0]}{site[1].upper()}" in requested]
        if not all_sites:
            raise ValueError("No requested sites matched the ddomain keys.")

    rotlib = Path(args.rotlib).expanduser().resolve() if args.rotlib else None
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results: List[Dict[str, object]] = []
    for site in all_sites:
        entry = ddomain[site]
        residues = sorted(entry.keys())
        site_label = f"{site[0]}{site[1]}"
        flow_path = flow_dir / f"flow_repacked_{site_label}.pdb"
        if not flow_path.is_file():
            raise FileNotFoundError(f"Missing FlowPacker seed for site {site_label}: {flow_path}")

        template_structure = Structure.from_pdb(flow_path)
        base_domain = DesignDomain(repack_sites=_to_site_specs(entry))
        domain = _restrict_domain_to_structure(template_structure, base_domain)

        designer = _run_repack(
            template_structure,
            domain,
            rotlib if rotlib is not None else DEFAULT_ROTLIB,
            args.trajectories,
            args.quiet,
        )
        best_structure = designer.best_structure
        if best_structure is None:
            raise RuntimeError("Best structure missing despite successful run.")
        best_energy = designer.best_sequence_energy

        pdb_out = output_dir / f"wt_repacked_{site_label}.pdb"
        best_structure.handle.write_pdb(str(pdb_out))

        repacked_atoms = _heavy_atom_lookup(pdb_out)
        ref = _collect_heavy_atoms(reference_atoms, residues)
        repacked = _collect_heavy_atoms(repacked_atoms, residues)
        if ref.count != repacked.count:
            raise RuntimeError(
                f"Heavy-atom mismatch for site {site_label}: {ref.count} vs {repacked.count}"
            )
        rmsd = _kabsch_rmsd(ref.positions, repacked.positions)

        print(
            f"Site {site_label}: residues={len(residues):2d}, "
            f"heavy-atom RMSD={rmsd:.4f} Å, input={flow_path}, output={pdb_out}"
        )
        results.append(
            {
                "site": site_label,
                "residue_count": len(residues),
                "heavy_atom_count": repacked.count,
                "heavy_atom_rmsd": rmsd,
                "structure_path": str(pdb_out),
                "flowpacker_input": str(flow_path),
                "total_energy": None if best_energy is None else best_energy.get("total"),
            }
        )

    csv_path = output_dir / "wt_domain_repacking_flow_uni_rmsd.csv"
    write_results_csv(results, csv_path)
    elapsed = time.perf_counter() - start
    print(f"\nSaved {len(results)} structures to {output_dir.resolve()}")
    print(f"RMSD summary: {csv_path}")
    print(f"Total runtime: {elapsed:.2f} seconds")


if __name__ == "__main__":
    main()
