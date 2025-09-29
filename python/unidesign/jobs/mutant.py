"""Mutant modelling job wrappers."""

from __future__ import annotations

import shutil
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Callable, Mapping

from ..artifacts import MutantStructureModel
from ..config import BuildMutantConfig
from ..runner import UniDesignRunner, UniDesignRunResult


def _discover_mutant_models(workdir: Path) -> list[Path]:
    """Return potential mutant PDB files emitted by UniDesign."""

    candidates = [
        path
        for path in workdir.rglob("*.pdb")
        if "Model_" in path.name and not path.name.endswith("_WT.pdb")
    ]
    candidates.sort()
    return candidates


def _relocate_mutant_models(
    workdir: Path,
    *,
    prefix: str,
    keep_workspace: bool,
) -> tuple[Path, list[MutantStructureModel], Callable[[], None]]:
    """Relocate mutant models respecting caller workspace preferences."""

    models = _discover_mutant_models(workdir)

    if keep_workspace:
        artifacts = [
            MutantStructureModel(
                path=path,
                prefix=prefix,
                logical_name=f"mutant_model_{index}",
                mutant_index=index,
            )
            for index, path in enumerate(models, start=1)
        ]
        cleanup = lambda: shutil.rmtree(workdir, ignore_errors=True)
        return workdir, artifacts, cleanup

    destination = Path(tempfile.mkdtemp(prefix="unidesign_artifacts_"))

    relocated: list[MutantStructureModel] = []
    for index, source in enumerate(models, start=1):
        try:
            relative = source.relative_to(workdir)
        except ValueError:
            relative = Path(source.name)
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        relocated.append(
            MutantStructureModel(
                path=target,
                prefix=prefix,
                logical_name=f"mutant_model_{index}",
                mutant_index=index,
            )
        )

    shutil.rmtree(workdir, ignore_errors=True)

    def _cleanup() -> None:
        shutil.rmtree(destination, ignore_errors=True)

    return destination, relocated, _cleanup


@dataclass(slots=True)
class MutantModelingResult:
    """Encapsulates outputs for ``BuildMutant`` runs."""

    run: UniDesignRunResult
    workspace: Path
    mutant_models: tuple[MutantStructureModel, ...]
    cleanup: Callable[[], None] | None

    def close(self) -> None:
        if self.cleanup is not None:
            self.cleanup()
            self.cleanup = None


class MutantModelingJob:
    """Execute the ``BuildMutant`` command with structured inputs."""

    def __init__(self, runner: UniDesignRunner, config: BuildMutantConfig) -> None:
        self._runner = runner
        self._config = config

    def run(
        self,
        *,
        keep_workspace: bool = False,
        env: Mapping[str, str] | None = None,
    ) -> MutantModelingResult:
        with TemporaryDirectory(prefix="unidesign_mutants_") as tmp_dir:
            mutant_file = Path(tmp_dir) / "mutants.txt"
            mutant_file.write_text(self._config.mutant_file_contents())
            config = replace(self._config, mutant_file_path=mutant_file)
            run_result = self._runner.run(
                config.to_cli_args(), env=env, persist_workdir=True
            )

        workspace, models, cleanup = _relocate_mutant_models(
            run_result.workdir,
            prefix=run_result.prefix,
            keep_workspace=keep_workspace,
        )
        run_result.workdir = workspace

        return MutantModelingResult(
            run=run_result,
            workspace=workspace,
            mutant_models=tuple(models),
            cleanup=cleanup,
        )


__all__ = ["MutantModelingJob", "MutantModelingResult"]
