"""Generate ligand parameter and topology files via UniDesign."""

from __future__ import annotations

from pathlib import Path

from unidesign import (
    LigandParameterizationJob,
    MakeLigParamConfig,
    UniDesignRunner,
    discover_binary,
)


def main() -> None:
    binary = discover_binary()
    runner = UniDesignRunner(binary)
    config = MakeLigParamConfig(
        mol2_path=Path(
            "example/ProteinLigandInteraction/1r091_BS01_JEN/lig_charge.mol2"
        ),
        ligand_parameter_path="lig_charge.prm",
        ligand_topology_path="lig_charge.top",
    )

    job = LigandParameterizationJob(runner, config)
    result = job.run()
    try:
        print("stdout:\n", result.run.stdout)
        if result.parameter_file:
            prm_target = result.parameter_file.persist("outputs/ligand")
            print("Parameter file copied to", prm_target)
        if result.topology_file:
            top_target = result.topology_file.persist("outputs/ligand")
            print("Topology file copied to", top_target)
    finally:
        result.close()


if __name__ == "__main__":
    main()
