"""Score a structure with the ComputeStability command."""

from __future__ import annotations

from pathlib import Path

from unidesign import (
    ComputeStabilityConfig,
    StabilityComputationJob,
    UniDesignRunner,
    discover_binary,
)


def main() -> None:
    binary = discover_binary()
    runner = UniDesignRunner(binary)
    config = ComputeStabilityConfig(
        pdb_path=Path("example/MonomerDesign/1igd/1igd.pdb"),
    )

    job = StabilityComputationJob(runner, config)
    result = job.run()
    try:
        print("stdout:\n", result.run.stdout)
        print("stderr:\n", result.run.stderr)
        if result.rotamer_list:
            report = result.rotamer_list.persist("outputs/stability")
            print("Stored rotamer list at", report)
    finally:
        result.close()


if __name__ == "__main__":
    main()
