"""Example script running BuildMutant on an existing PDB file."""

from __future__ import annotations

from unidesign import (
    BinaryDiscoveryError,
    BuildMutantConfig,
    UniDesignRunner,
    discover_binary,
)
from unidesign.jobs import MutantStructureBuildJob
from unidesign.paths import project_root


def main() -> None:
    try:
        binary = discover_binary()
    except BinaryDiscoveryError as exc:
        raise SystemExit(f"UniDesign binary could not be located: {exc}") from exc

    runner = UniDesignRunner(binary)

    # Example PDB bundled with UniDesign
    pdb_path = project_root() / "example/ProteinProteinInteractionDesign/1ay7/1ay7.pdb"
    if not pdb_path.is_file():
        raise SystemExit(f"Example PDB file not found: {pdb_path}")

    # Configuration for BuildMutant
    config = BuildMutantConfig(pdb_path=pdb_path)

    # Define mutations: dict[chain_id -> mutation(s)]
    # Example: mutate residue 50 in chain A from Valine to Alanine (V50A)
    mutants = [
        {"A": "I70V"},
        {"A": ["I70F", "T72R"]},  # multiple substitutions in chain A
    ]

    job = MutantStructureBuildJob(runner, config)
    result = job.run(mutants, keep_workspace=True)

    try:
        print("Return code:", result.run.returncode)
        for label, model in result.mutant_models.items():
            print(f"Mutant {label}: {model.path}")
    finally:
        result.close()


if __name__ == "__main__":
    main()
