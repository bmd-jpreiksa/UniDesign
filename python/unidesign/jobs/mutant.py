"""Mutant structure modelling job wrappers."""

from __future__ import annotations

import re
import tempfile
from collections.abc import Mapping as MappingCollection
from collections.abc import Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Callable, Mapping

from ..artifacts import StructureModel
from ..config import BuildMutantConfig
from ..runner import UniDesignRunResult, UniDesignRunner
from ._shared import ArtifactSpec, relocate_artifacts

MutationValue = str | Sequence[str]
"""Accepted value types describing a mutation for a specific chain."""

MutantSpecification = MappingCollection[str, MutationValue]
"""Mapping between chain identifiers and mutation descriptors."""

_WITH_CHAIN_PATTERN = re.compile(
    r"^(?P<native>[A-Za-z])(?P<chain>[A-Za-z])(?P<position>\d+)(?P<mutant>[A-Za-z])$"
)
_WITHOUT_CHAIN_PATTERN = re.compile(
    r"^(?P<native>[A-Za-z])(?P<position>\d+)(?P<mutant>[A-Za-z])$"
)
_SANITISE_PATTERN = re.compile(r"[^0-9A-Za-z]+")


def _normalise_chain_id(chain_id: str) -> str:
    chain = chain_id.strip()
    if len(chain) != 1 or not chain.isalnum():
        raise ValueError(f"Invalid chain identifier: {chain_id!r}")
    return chain.upper()


def _coerce_mutation_values(value: MutationValue) -> list[str]:
    if isinstance(value, str):
        return [value]
    if isinstance(value, Sequence):
        return [str(item) for item in value]
    raise TypeError(f"Unsupported mutation descriptor type: {type(value)!r}")


def _parse_mutation_descriptor(chain_id: str, descriptor: str) -> str:
    cleaned = descriptor.strip()
    if not cleaned:
        raise ValueError("Mutation descriptor cannot be empty")

    match = _WITH_CHAIN_PATTERN.match(cleaned)
    if match:
        native = match.group("native").upper()
        chain = match.group("chain").upper()
        position = match.group("position")
        mutant = match.group("mutant").upper()
        if chain != chain_id:
            raise ValueError(
                "Chain identifier mismatch between dictionary key and mutation descriptor"
            )
        return f"{native}{chain}{position}{mutant}"

    match = _WITHOUT_CHAIN_PATTERN.match(cleaned)
    if match:
        native = match.group("native").upper()
        position = match.group("position")
        mutant = match.group("mutant").upper()
        return f"{native}{chain_id}{position}{mutant}"

    raise ValueError(f"Invalid mutation descriptor: {descriptor!r}")


def _serialise_mutants(
    mutants: Sequence[MutantSpecification],
) -> tuple[list[str], list[str]]:
    if not mutants:
        raise ValueError("At least one mutant must be specified")

    manifest_lines: list[str] = []
    labels: list[str] = []

    for specification in mutants:
        if not specification:
            raise ValueError("Mutation specification entries cannot be empty")

        components: list[str] = []
        for chain_id, descriptor in specification.items():
            normalised_chain = _normalise_chain_id(str(chain_id))
            for value in _coerce_mutation_values(descriptor):
                components.append(_parse_mutation_descriptor(normalised_chain, value))

        if not components:
            raise ValueError("Mutation specification must include at least one substitution")

        label = ",".join(components)
        manifest_lines.append(f"{label};")
        labels.append(label)

    return manifest_lines, labels


def _write_mutant_manifest(mutants: Sequence[MutantSpecification]) -> tuple[Path, list[str]]:
    manifest_lines, labels = _serialise_mutants(mutants)
    with tempfile.NamedTemporaryFile(
        "w", prefix="unidesign_mutants_", suffix=".txt", delete=False, encoding="utf-8"
    ) as handle:
        handle.write("\n".join(manifest_lines))
        manifest_path = Path(handle.name)
    return manifest_path, labels


def _sanitise_label(label: str) -> str:
    sanitised = _SANITISE_PATTERN.sub("_", label).strip("_")
    return sanitised or "mutant"


def _unique_label(base: str, used: set[str], index: int) -> str:
    candidate = base or f"mutant_{index}"
    if candidate not in used:
        return candidate
    suffix = 2
    while f"{candidate}_{suffix}" in used:
        suffix += 1
    return f"{candidate}_{suffix}"


@dataclass(slots=True)
class MutantStructureBuildResult:
    """Bundle representing outputs from :class:`MutantStructureBuildJob`."""

    run: UniDesignRunResult
    workspace: Path
    mutant_models: dict[str, StructureModel]
    cleanup: Callable[[], None] | None

    def close(self) -> None:
        if self.cleanup is not None:
            self.cleanup()
            self.cleanup = None


class MutantStructureBuildJob:
    """Execute the ``BuildMutant`` command with structured mutation inputs."""

    def __init__(self, runner: UniDesignRunner, config: BuildMutantConfig) -> None:
        self._runner = runner
        self._config = config

    def run(
        self,
        mutants: Sequence[MutantSpecification],
        *,
        keep_workspace: bool = False,
        env: Mapping[str, str] | None = None,
    ) -> MutantStructureBuildResult:
        manifest_path, labels = _write_mutant_manifest(mutants)
        config = replace(self._config, mutant_file=manifest_path)

        try:
            run_result = self._runner.run(
                config.to_cli_args(), env=env, persist_workdir=True
            )
        finally:
            manifest_path.unlink(missing_ok=True)

        pdb_stem = Path(config.pdb_path).stem
        candidates: dict[str, ArtifactSpec] = {}
        for index, label in enumerate(labels, start=1):
            expected_name = f"pdb_Model_{index:04d}.pdb"
            sanitised = _sanitise_label(label)

            def factory(path: Path, prefix: str, logical_name: str, *, _sanitised=sanitised) -> StructureModel:
                return StructureModel(path=path, prefix=prefix, logical_name=_sanitised)

            candidates[label] = ArtifactSpec(relative_path=Path(expected_name), factory=factory)

        workspace, artifacts, cleanup = relocate_artifacts(
            run_result.workdir,
            candidates,
            keep_workspace=keep_workspace,
            prefix=run_result.prefix,
        )
        run_result.workdir = workspace

        missing = [label for label in labels if label not in artifacts]
        if missing:
            if not keep_workspace:
                cleanup()
            raise FileNotFoundError(
                "Expected mutant model(s) not produced: " + ", ".join(missing)
            )

        used_names: set[str] = set()
        mutant_models: dict[str, StructureModel] = {}

        for index, label in enumerate(labels, start=1):
            artifact = artifacts[label]
            base_name = _unique_label(_sanitise_label(label), used_names, index)
            used_names.add(base_name)
            target_name = f"{base_name}.pdb"
            target_path = artifact.path.with_name(target_name)
            if artifact.path.name != target_name:
                artifact.path.rename(target_path)
                artifact.path = target_path
            artifact.logical_name = base_name
            mutant_models[label] = artifact

        return MutantStructureBuildResult(
            run=run_result,
            workspace=workspace,
            mutant_models=mutant_models,
            cleanup=cleanup,
        )


__all__ = ["MutantStructureBuildJob", "MutantStructureBuildResult"]

