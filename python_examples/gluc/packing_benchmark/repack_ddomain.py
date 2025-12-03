#!/usr/bin/env python3
"""
Repack each design domain defined in ddomain.pickle using the same Rosetta helpers
as the site_48 notebook and report heavy-atom RMSDs (domain-only) versus the input pose.
"""
from __future__ import annotations

import argparse
import csv
import pickle
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import pyrosetta as pr
from pyrosetta import Pose, rosetta

from bmdairs.utilities.rosetta import movers as rst_mv
from bmdairs.utilities.rosetta import scores as rst_sc
from bmdairs.utilities.rosetta.init import rosetta_custom_init


DDomainEntry = Dict[Tuple[int, str], str]
DDomain = Dict[Tuple[int, str], DDomainEntry]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="WT repacking and RMSD calculation for Rosetta design domains."
    )
    parser.add_argument(
        "--structure",
        default="3m0m_prepared_openmm.pdb",
        help="Starting WT structure (default: %(default)s).",
    )
    parser.add_argument(
        "--ddomain",
        default="ddomain.pickle",
        help="Design-domain pickle generated previously (default: %(default)s).",
    )
    parser.add_argument(
        "--params-dir",
        default="params/rosetta",
        help="Directory containing Rosetta params to mirror the notebook setup.",
    )
    parser.add_argument(
        "--extra-options",
        default="-in:auto_setup_metals",
        help="Extra command-line options passed to Rosetta (default: %(default)s).",
    )
    parser.add_argument(
        "--sites",
        nargs="*",
        help="Optional subset of sites (e.g. 344A 390A). Process all domains if omitted.",
    )
    parser.add_argument(
        "--output-dir",
        default="wt_domain_repacking",
        help="Directory to store repacked structures and RMSD CSV (default: %(default)s).",
    )
    return parser.parse_args()


def init_rosetta(params_dir: Path, extra_options: str) -> None:
    param_files = sorted(params_dir.glob("*.params"))
    rosetta_custom_init(
        residue_parameters=[str(p.resolve()) for p in param_files],
        extra_options=extra_options,
    )


def load_ddomain(path: Path) -> DDomain:
    with path.open("rb") as handle:
        ddomain = pickle.load(handle)
    if not isinstance(ddomain, dict):
        raise ValueError(f"Expected ddomain dict, got {type(ddomain)}")
    return ddomain


def to_one_letter(name: str) -> str:
    aa_enum = rosetta.core.chemical.aa_from_one_or_three(name.upper())
    if aa_enum == rosetta.core.chemical.aa_unk:
        raise ValueError(f"Unrecognized residue name '{name}' for PIKAA conversion.")
    return rosetta.core.chemical.oneletter_code_from_aa(aa_enum)


def domain_list(entry: DDomainEntry) -> List[Tuple[int, str, str]]:
    return [(resnum, chain, to_one_letter(aa)) for (resnum, chain), aa in entry.items()]


def residue_indices(pose: Pose, residue_ids: Iterable[Tuple[int, str]]) -> List[int]:
    pdb_info = pose.pdb_info()
    indices: List[int] = []
    for resnum, chain in residue_ids:
        idx = pdb_info.pdb2pose(chain, resnum)
        if idx < 1:
            raise ValueError(f"Residue {resnum}{chain} not present in pose numbering.")
        indices.append(idx)
    return indices


def residue_subset(indices: Sequence[int]) -> rosetta.std.list_unsigned_long_t:
    subset = rosetta.std.list_unsigned_long_t()
    for idx in indices:
        subset.append(idx)
    return subset


def count_heavy_atoms(pose: Pose, indices: Sequence[int]) -> int:
    count = 0
    for idx in indices:
        residue = pose.residue(idx)
        for atom_idx in range(1, residue.natoms() + 1):
            if residue.atom_type(atom_idx).is_hydrogen():
                continue
            count += 1
    return count


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
            ],
        )
        writer.writeheader()
        writer.writerows(results)


def main() -> None:
    args = parse_args()
    init_rosetta(Path(args.params_dir), args.extra_options)

    ddomain = load_ddomain(Path(args.ddomain))
    all_sites = sorted(ddomain.keys())
    if args.sites:
        requested = {site.upper() for site in args.sites}
        all_sites = [
            site for site in all_sites if f"{site[0]}{site[1].upper()}" in requested
        ]
        if not all_sites:
            raise ValueError("No requested sites matched the ddomain keys.")

    pose = pr.pose_from_file(args.structure)
    scorefxn = rst_sc.get_rosetta_full_atom_score_function()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for site in all_sites:
        entry = ddomain[site]
        ddomain_list = domain_list(entry)
        residues = residue_indices(pose, entry.keys())
        subset = residue_subset(residues)

        wt_pose = Pose()
        wt_pose.assign(pose)
        repacked = Pose()
        repacked.assign(pose)

        tf = rst_mv.get_task_factory(ddomain=ddomain_list)
        pack = rst_mv.PackRotamersMover(scorefxn)
        pack.task_factory(tf)
        pack.apply(repacked)

        rmsd = rosetta.core.scoring.all_atom_rmsd(wt_pose, repacked, subset)
        heavy_atoms = count_heavy_atoms(repacked, residues)

        site_label = f"{site[0]}{site[1]}"
        pdb_out = output_dir / f"wt_repacked_{site_label}.pdb"
        repacked.dump_pdb(str(pdb_out))

        print(
            f"Site {site_label}: residues={len(residues):2d}, heavy-atom RMSD={rmsd:.4f} Å, output={pdb_out}"
        )
        results.append(
            {
                "site": site_label,
                "residue_count": len(residues),
                "heavy_atom_count": heavy_atoms,
                "heavy_atom_rmsd": rmsd,
                "structure_path": str(pdb_out),
            }
        )

    csv_path = output_dir / "wt_domain_repacking_rmsd.csv"
    write_results_csv(results, csv_path)
    print(f"\nSaved {len(results)} structures to {output_dir.resolve()}")
    print(f"RMSD summary: {csv_path}")


if __name__ == "__main__":
    main()
