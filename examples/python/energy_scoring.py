"""Example script demonstrating energy-scoring helpers."""

from __future__ import annotations

from unidesign import (
    BinaryDiscoveryError,
    ComputeBindingConfig,
    ComputeStabilityConfig,
    UniDesignRunner,
    discover_binary,
)
from unidesign.jobs import BindingComputationJob, StabilityComputationJob
from unidesign.paths import project_root


def main() -> None:
    try:
        binary = discover_binary()
    except BinaryDiscoveryError as exc:  # pragma: no cover - runtime safeguard
        raise SystemExit(f"UniDesign binary could not be located: {exc}") from exc

    runner = UniDesignRunner(binary)
    pdb_path = project_root() / "example/ProteinProteinInteractionDesign/1ay7/1ay7.pdb"
    if not pdb_path.is_file():
        raise SystemExit(f"Example PDB file not found: {pdb_path}")

    stability_job = StabilityComputationJob(
        runner,
        ComputeStabilityConfig(pdb_path=pdb_path),
    )
    stability_result = stability_job.run()
    try:
        print("Stability return code:", stability_result.run.returncode)
        if stability_result.rotamer_list:
            print("Rotamer list:", stability_result.rotamer_list.path)
    finally:
        stability_result.close()

    binding_job = BindingComputationJob(
        runner,
        ComputeBindingConfig(pdb_path=pdb_path),
    )
    binding_result = binding_job.run()
    try:
        print("Binding return code:", binding_result.run.returncode)
    finally:
        binding_result.close()


if __name__ == "__main__":  # pragma: no cover - manual execution hook
    main()
