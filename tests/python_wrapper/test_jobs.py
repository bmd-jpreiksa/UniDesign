from __future__ import annotations

from pathlib import Path
from unittest.mock import Mock

from unidesign import (
    ComputeStabilityConfig,
    LigandParameterizationJob,
    MakeLigParamConfig,
    ProteinDesignConfig,
    ProteinDesignJob,
    StabilityComputationJob,
)
from unidesign.jobs import JobArtifact
from unidesign.runner import UniDesignRunResult, UniDesignRunner


def _make_run_result(workdir: Path, prefix: str = "testprefix") -> UniDesignRunResult:
    return UniDesignRunResult(
        args=("--command", "ProteinDesign"),
        returncode=0,
        stdout="ok",
        stderr="",
        workdir=workdir,
        prefix=prefix,
    )


def test_protein_design_job_collects_artifacts(tmp_path: Path) -> None:
    pdb_file = tmp_path / "input.pdb"
    pdb_file.write_text("PDB")

    runner = Mock(spec=UniDesignRunner)
    workdir = tmp_path / "workdir"
    workdir.mkdir()

    prefix = "designcase"
    (workdir / f"{prefix}_bestseqs").write_text("SEQ")
    (workdir / f"{prefix}_beststruct").write_text("ATOM")

    run_result = _make_run_result(workdir=workdir, prefix=prefix)
    runner.run.return_value = run_result

    job = ProteinDesignJob(runner, ProteinDesignConfig(pdb_path=pdb_file))
    result = job.run()

    runner.run.assert_called_once()
    call_args, call_kwargs = runner.run.call_args
    assert call_args == (job._config.to_cli_args(),)
    assert call_kwargs["env"] is None
    assert call_kwargs["persist_workdir"] is True

    assert result.best_sequences is not None
    assert result.best_sequences.read_text() == "SEQ"
    assert result.best_structure is not None
    assert result.best_structure.read_text() == "ATOM"
    assert result.workspace != workdir
    assert not workdir.exists()
    assert result.run.workdir == result.workspace

    export = result.best_sequences.persist(tmp_path / "exports")
    assert export.read_text() == "SEQ"

    result.close()
    assert not result.workspace.exists()


def test_ligand_job_can_keep_workspace(tmp_path: Path) -> None:
    mol2_file = tmp_path / "ligand.mol2"
    mol2_file.write_text("mol2")

    runner = Mock(spec=UniDesignRunner)
    workdir = tmp_path / "ligand_work"
    workdir.mkdir()
    (workdir / "lig_charge.prm").write_text("param")
    (workdir / "lig_charge.top").write_text("topo")

    run_result = _make_run_result(workdir=workdir)
    runner.run.return_value = run_result

    config = MakeLigParamConfig(
        mol2_path=mol2_file,
        ligand_parameter_path="lig_charge.prm",
        ligand_topology_path="lig_charge.top",
    )
    job = LigandParameterizationJob(runner, config)
    result = job.run(keep_workspace=True)

    assert result.workspace == workdir
    assert result.parameter_file is not None
    assert result.parameter_file.read_text() == "param"
    assert result.topology_file is not None
    assert result.topology_file.read_text() == "topo"

    result.close()
    assert not workdir.exists()


def test_stability_job_relocates_rotamer_list(tmp_path: Path) -> None:
    pdb_file = tmp_path / "input.pdb"
    pdb_file.write_text("PDB")

    runner = Mock(spec=UniDesignRunner)
    workdir = tmp_path / "stability_work"
    workdir.mkdir()
    (workdir / "testprefix_rotlist.txt").write_text("rotamers")

    run_result = _make_run_result(workdir=workdir)
    runner.run.return_value = run_result

    job = StabilityComputationJob(
        runner,
        ComputeStabilityConfig(pdb_path=pdb_file),
    )
    result = job.run()

    assert result.rotamer_list is not None
    assert isinstance(result.rotamer_list, JobArtifact)
    assert result.rotamer_list.read_text() == "rotamers"
    assert result.workspace != workdir
    assert not workdir.exists()

    result.close()
