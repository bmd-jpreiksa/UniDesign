from __future__ import annotations

from pathlib import Path

from unidesign.api import Structure

# Example complexes for binding-energy evaluation.
EXAMPLES = {
    "1e44": (
        Path(__file__).resolve().parents[1]
        / "example"
        / "ProteinProteinInteractionDesign"
        / "1e44"
        / "1e44.pdb",
        [("AB", "C")],
    ),
    "1ay7": (
        Path(__file__).resolve().parents[1]
        / "example"
        / "ProteinProteinInteractionDesign"
        / "1ay7"
        / "1ay7.pdb",
        [],
    ),
}


def main() -> None:
    for pdb_id, (path, splits) in EXAMPLES.items():
        print(f"### Complex {pdb_id} ###")
        structure = Structure.from_pdb(path)

        print("=== Binding energy between chains A and B ===", flush=True)
        structure.compute_binding()

        for split1, split2 in splits:
            print()
            print(
                f"=== Binding energy after chain split ({split1} vs {split2}) ===",
                flush=True,
            )
            try:
                structure.compute_binding(split1=split1, split2=split2)
            except ValueError as exc:
                print(f"Skipping split: {exc}")

        print()


if __name__ == "__main__":
    main()
