"""Run a ProteinDesign monomer optimization using the Python wrapper."""

from __future__ import annotations

from pathlib import Path

from unidesign import (
    ProteinDesignConfig,
    ProteinDesignJob,
    UniDesignRunner,
    discover_binary,
)


def main() -> None:
    binary = discover_binary()
    runner = UniDesignRunner(binary)
    config = ProteinDesignConfig(
        pdb_path=Path("example/MonomerDesign/1agy/1agy.pdb"),
        design_chains="A",
        resfile_path=Path("example/MonomerDesign/1agy/RESFILE.txt"),
    )

    job = ProteinDesignJob(runner, config)
    result = job.run()
    try:
        print("stdout:\n", result.run.stdout)
        if result.best_sequences:
            print("Best sequences:\n", result.best_sequences.read_text())
        if result.best_structure:
            output = result.best_structure.persist(
                "outputs/monomer_design", filename="best_structure.pdb"
            )
            print("Persisted best structure to", output)
    finally:
        result.close()


if __name__ == "__main__":
    main()
