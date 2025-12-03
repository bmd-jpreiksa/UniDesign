from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description="Plot per-term stability deltas")
    parser.add_argument('--variant', default='baseline', help='weight variant label (default: baseline)')
    args = parser.parse_args()

    terms_path = Path(f"python_examples/gluc/wt_repack_weight_scan/stability_debug/stability_terms_{args.variant}.tsv")
    output_path = Path(f"python_examples/gluc/wt_repack_weight_scan/stability_debug/{args.variant}_delta.png")

    terms = []
    with terms_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                delta = float(row["delta"])
            except (KeyError, ValueError):
                continue
            if abs(delta) >= 1.0:
                terms.append((row["term"], delta))

    if not terms:
        print("No terms exceed |Δ| >= 1.0; nothing to plot")
        return

    terms.sort(key=lambda item: item[1])
    labels, deltas = zip(*terms)
    colors = ["#2c7bb6" if value < 0 else "#d7191c" for value in deltas]

    plt.figure(figsize=(8, 5))
    plt.barh(labels, deltas, color=colors)
    plt.axvline(0.0, color="black", linewidth=0.8)
    plt.xlabel("Δ energy (after - before)")
    plt.title(f"{args.variant} Stability Δ (|Δ| ≥ 1)")
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=200)
    print(f"Plot written to {output_path}")


if __name__ == "__main__":
    main()
