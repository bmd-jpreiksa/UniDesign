"""Generate ligand parameters using the UniDesign Python wrapper."""

from __future__ import annotations

from unidesign import (
    BinaryDiscoveryError,
    MakeLigParamConfig,
    UniDesignRunner,
    discover_binary,
)
from unidesign.jobs import LigandParameterizationJob
from unidesign.paths import project_root


def main() -> None:
    try:
        binary = discover_binary()
    except BinaryDiscoveryError as exc:  # pragma: no cover - runtime safeguard
        raise SystemExit(f"UniDesign binary could not be located: {exc}") from exc

    runner = UniDesignRunner(binary)
    mol2_path = (
        project_root()
        / "example/ProteinLigandInteraction/1r091_BS01_JEN/lig_charge.mol2"
    )
    if not mol2_path.is_file():
        raise SystemExit(f"Example MOL2 file not found: {mol2_path}")

    config = MakeLigParamConfig(
        mol2_path=mol2_path,
        ligand_parameter_path="example_LIG_PARAM.prm",
        ligand_topology_path="example_LIG_TOPO.inp",
    )

    job = LigandParameterizationJob(runner, config)
    result = job.run(keep_workspace=True)
    try:
        print("Parameter file:", result.parameter_file.path if result.parameter_file else "<missing>")
        print("Topology file:", result.topology_file.path if result.topology_file else "<missing>")
    finally:
        result.close()


if __name__ == "__main__":  # pragma: no cover - manual execution hook
    main()
