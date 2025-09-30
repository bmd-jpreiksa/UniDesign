"""Integration-style smoke tests for job helpers."""

from __future__ import annotations

from pathlib import Path
import subprocess

import pytest

from unidesign import (
    BuildMutantConfig,
    ComputeBindingConfig,
    ComputeStabilityConfig,
    MakeLigParamConfig,
    ProteinDesignConfig,
)
from unidesign.jobs import (
    BindingComputationJob,
    LigandParameterizationJob,
    MutantStructureBuildJob,
    ProteinDesignJob,
    StabilityComputationJob,
)
from unidesign.runner import UniDesignRunner


@pytest.fixture
def runner_with_fake_binary(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Create a runner wired to a fake subprocess implementation."""

    binary = tmp_path / "UniDesign"
    binary.write_text("#!/bin/sh\n")
    binary.chmod(0o755)

    resources: list[tuple[str, Path]] = []
    for name in ("library", "wread", "extbin"):
        resource = tmp_path / f"resource_{name}"
        resource.mkdir()
        (resource / "placeholder.txt").write_text(name)
        resources.append((name, resource))

    monkeypatch.setattr(UniDesignRunner, "_STATIC_RESOURCES", tuple(resources))
    runner = UniDesignRunner(binary)
    calls: list[tuple[str | None, Path, dict[str, object]]] = []

    def fake_run(argv, cwd, env, check, capture_output, text):
        workdir = Path(cwd)
        prefix = argv[argv.index("--prefix") + 1]
        command = None
        if "--command" in argv:
            command = argv[argv.index("--command") + 1]

        metadata: dict[str, object] = {}

        if command == "ProteinDesign":
            (workdir / f"{prefix}_selfenergy.txt").write_text("energy")
            (workdir / f"{prefix}_bestseqs").write_text("SEQ\n")
            (workdir / f"{prefix}_beststruct").write_text("MODEL\n")
        elif command == "ComputeStability":
            (workdir / f"{prefix}_rotlist.txt").write_text("rotamer")
        elif command == "MakeLigParamAndTopo":
            param_path = Path(argv[argv.index("--lig_param") + 1])
            topo_path = Path(argv[argv.index("--lig_topo") + 1])
            (workdir / param_path).write_text("PARAMS")
            (workdir / topo_path).write_text("TOPO")
        elif command == "BuildMutant":
            pdb_path = Path(argv[argv.index("--pdb") + 1])
            mutant_file = Path(argv[argv.index("--mutant_file") + 1])
            lines = mutant_file.read_text(encoding="utf-8").splitlines()
            metadata["mutant_lines"] = lines
            generated: list[str] = []
            for index, _ in enumerate(lines, start=1):
                #filename = f"{pdb_path.stem}_Model_{index:04d}.pdb"
                filename = f"pdb_Model_{index:04d}.pdb"
                (workdir / filename).write_text(f"MODEL {index}\n")
                generated.append(filename)
            metadata["generated_files"] = generated
        # ComputeBinding does not write additional artefacts in this smoke test.

        calls.append((command, workdir, metadata))
        return subprocess.CompletedProcess(argv, 0, stdout=f"{command} stdout", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    return runner, calls


def _create_minimal_pdb(path: Path) -> None:
    path.write_text(
        "HEADER    TEST PDB\n"
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "TER\n"
        "END\n"
    )


def _create_minimal_mol2(path: Path) -> None:
    path.write_text(
        "@<TRIPOS>MOLECULE\n"
        "LIG\n"
        " 1 0 0 0 0\n"
        "SMALL\n"
        "NO_CHARGES\n"
        "@<TRIPOS>ATOM\n"
        "      1 C1         0.0    0.0    0.0 C.3       1 <0>        0.0000\n"
        "@<TRIPOS>BOND\n"
    )


def test_protein_design_job_end_to_end(runner_with_fake_binary, tmp_path: Path):
    runner, calls = runner_with_fake_binary

    pdb_path = tmp_path / "input.pdb"
    _create_minimal_pdb(pdb_path)
    config = ProteinDesignConfig(pdb_path=pdb_path, design_chains="A", n_trajectories=1)

    job = ProteinDesignJob(runner, config)
    result = job.run()

    try:
        assert result.run.returncode == 0
        assert result.self_energy is not None
        assert result.self_energy.read_text() == "energy"
        assert result.best_sequences is not None
        assert "SEQ" in result.best_sequences.read_text()
        assert result.best_structure is not None
        assert result.best_structure.read_text().startswith("MODEL")
        assert calls and calls[0][0] == "ProteinDesign"
        original_workdir = calls[0][1]
        assert not original_workdir.exists(), "Original workspace should be relocated"
        assert result.workspace.exists()
    finally:
        result.close()
        assert not result.workspace.exists()


def test_energy_jobs_cover_rotamer_and_binding(runner_with_fake_binary, tmp_path: Path):
    runner, calls = runner_with_fake_binary

    pdb_path = tmp_path / "energy_input.pdb"
    _create_minimal_pdb(pdb_path)

    stability_job = StabilityComputationJob(runner, ComputeStabilityConfig(pdb_path=pdb_path))
    stability_result = stability_job.run(keep_workspace=True)
    try:
        assert stability_result.rotamer_list is not None
        assert stability_result.rotamer_list.read_text() == "rotamer"
        assert calls and any(entry[0] == "ComputeStability" for entry in calls)
    finally:
        stability_result.close()

    binding_job = BindingComputationJob(runner, ComputeBindingConfig(pdb_path=pdb_path))
    binding_result = binding_job.run()
    try:
        assert binding_result.run.returncode == 0
        assert any(entry[0] == "ComputeBinding" for entry in calls)
    finally:
        binding_result.close()


def test_mutant_structure_job_generates_named_models(
    runner_with_fake_binary, tmp_path: Path
):
    runner, calls = runner_with_fake_binary

    pdb_path = tmp_path / "mutant_input.pdb"
    _create_minimal_pdb(pdb_path)

    config = BuildMutantConfig(pdb_path=pdb_path)
    job = MutantStructureBuildJob(runner, config)
    mutants = [
        {"A": "Q22D"},
        {"A": ["H18F", "Q22D"]},
        {"A": ["H18F"], "B": ["M20A"]},
    ]

    result = job.run(mutants)
    try:
        expected_labels = ["QA22D", "HA18F,QA22D", "HA18F,MB20A"]
        assert list(result.mutant_models.keys()) == expected_labels
        expected_files = {
            "QA22D": "QA22D.pdb",
            "HA18F,QA22D": "HA18F_QA22D.pdb",
            "HA18F,MB20A": "HA18F_MB20A.pdb",
        }
        assert {
            label: artifact.path.name for label, artifact in result.mutant_models.items()
        } == expected_files

        command, _, metadata = next(entry for entry in calls if entry[0] == "BuildMutant")
        assert command == "BuildMutant"
        assert metadata.get("mutant_lines") == [f"{label};" for label in expected_labels]
        assert metadata.get("generated_files") == [
            f"pdb_Model_{index:04d}.pdb" for index in range(1, 4)
        ]

        for artifact in result.mutant_models.values():
            assert artifact.path.exists()
            assert artifact.path.parent == result.workspace
    finally:
        result.close()
        assert not result.workspace.exists()


def test_ligand_parameter_job_handles_outputs(runner_with_fake_binary, tmp_path: Path):
    runner, calls = runner_with_fake_binary

    mol2_path = tmp_path / "ligand.mol2"
    _create_minimal_mol2(mol2_path)
    config = MakeLigParamConfig(mol2_path=mol2_path)

    job = LigandParameterizationJob(runner, config)
    result = job.run()
    try:
        assert result.parameter_file is not None
        assert result.parameter_file.read_text() == "PARAMS"
        assert result.topology_file is not None
        assert result.topology_file.read_text() == "TOPO"
        assert any(entry[0] == "MakeLigParamAndTopo" for entry in calls)
    finally:
        result.close()
