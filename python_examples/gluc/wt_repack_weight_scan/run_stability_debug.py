from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List

from unidesign.api import Structure

BASE_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[3]
WEIGHT_DIR = BASE_DIR / "wt_repack_weight_scan" / "weights"
STRUCT_DIR = BASE_DIR / "wt_repack_weight_scan" / "structures"
OUTPUT_DIR = BASE_DIR / "wt_repack_weight_scan" / "stability_debug"

ROT_LIB = REPO_ROOT / "library" / "rotlib" / "ALLbbdep.bin"

VARIANTS = [
    "baseline",
    "repel_d",
    "repel_d_plus",
    "repel_lower_S",
    "repel_h_S",
    "repel_ultra",
]


def _run_compute_stability(pdb: Path, weight: Path) -> Dict[str, float]:
    structure = Structure.from_pdb(pdb)
    result = structure.compute_stability(weight_file=weight, rotlib_file=ROT_LIB)
    return result.terms


def run_debug() -> None:
    OUTPUT_DIR.mkdir(exist_ok=True)
    summary_lines: List[str] = ["label\tstructure\ttotal"]
    for label in VARIANTS:
        print(f"Computing stability diagnostics for variant '{label}'")
        weight_path = WEIGHT_DIR / f"weight_all1_{label}.wgt"
        if not weight_path.exists():
            print(f"  Weight file {weight_path} missing; skipping")
            continue

        initial_terms = _run_compute_stability(BASE_DIR / "protein.pdb", weight_path)
        best_path = STRUCT_DIR / f"best_structure_{label}.pdb"
        best_terms = _run_compute_stability(best_path, weight_path) if best_path.exists() else {}

        out_json = OUTPUT_DIR / f"stability_{label}.json"
        out_json.write_text(json.dumps({"initial": initial_terms, "best": best_terms}, indent=2))

        out_tsv = OUTPUT_DIR / f"stability_terms_{label}.tsv"
        terms = sorted(set(initial_terms) | set(best_terms))
        with out_tsv.open("w", encoding="utf-8") as handle:
            handle.write("term\tinitial\tbest\tdelta\n")
            for term in terms:
                before = initial_terms.get(term)
                after = best_terms.get(term)
                delta = (after - before) if (after is not None and before is not None) else None
                handle.write(
                    f"{term}\t"
                    f"{'' if before is None else f'{before:.6f}'}\t"
                    f"{'' if after is None else f'{after:.6f}'}\t"
                    f"{'' if delta is None else f'{delta:.6f}'}\n"
                )

        total_line = lambda terms, tag: f"{label}\t{tag}\t{terms.get('Total', float('nan')):.6f}"
        summary_lines.append(total_line(initial_terms, "initial"))
        if best_terms:
            summary_lines.append(total_line(best_terms, "best"))

    (OUTPUT_DIR / "stability_summary.tsv").write_text("\n".join(summary_lines) + "\n")


if __name__ == "__main__":
    run_debug()
