"""Example script running a monomer-design trajectory via the Python wrapper."""

from __future__ import annotations

from unidesign import (
    BinaryDiscoveryError,
    ProteinDesignConfig,
    UniDesignRunner,
    discover_binary,
)
from unidesign.jobs import ProteinDesignJob
from unidesign.paths import project_root


def main() -> None:
    try:
        binary = discover_binary()
    except BinaryDiscoveryError as exc:  # pragma: no cover - runtime safeguard
        raise SystemExit(f"UniDesign binary could not be located: {exc}") from exc

    runner = UniDesignRunner(binary)

    pdb_path = project_root() / "example/MonomerDesign/1igd/1igd.pdb"
    if not pdb_path.is_file():
        raise SystemExit(f"Example PDB file not found: {pdb_path}")

    config = ProteinDesignConfig(
        pdb_path=pdb_path,
        design_chains="A",
        n_trajectories=1,
        rotate_hydroxyl=True,
    )

    job = ProteinDesignJob(runner, config)
    result = job.run(keep_workspace=False)

    try:
        print(f"Return code: {result.run.returncode}")
        print(f"Prefix: {result.prefix}")
        if result.best_structure:
            print("Best structure:", result.best_structure.path)
        if result.best_sequences:
            print("Best sequences:", result.best_sequences.path)
    finally:
        result.close()


if __name__ == "__main__":  # pragma: no cover - manual execution hook
    main()
