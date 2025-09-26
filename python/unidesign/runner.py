"""Execution helpers for invoking the UniDesign binary."""

from __future__ import annotations

import os
import shutil
import subprocess
import uuid
from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Mapping, MutableMapping, Sequence

from . import paths


@dataclass(slots=True)
class UniDesignRunResult:
    """Structured response describing a UniDesign execution."""

    args: tuple[str, ...]
    """The argument vector used to invoke the binary (excluding the executable)."""

    returncode: int
    """Process return code reported by :func:`subprocess.run`."""

    stdout: str
    """Captured standard output from the UniDesign binary."""

    stderr: str
    """Captured standard error from the UniDesign binary."""

    workdir: Path
    """Working directory that held input/output files for the execution."""

    prefix: str
    """Prefix automatically injected into CLI arguments for unique output names."""


class UniDesignRunner:
    """Convenience wrapper around the UniDesign command line binary."""

    _STATIC_RESOURCES: tuple[tuple[str, Path], ...] = (
        ("library", paths.library_dir()),
        ("wread", paths.wread_dir()),
        ("extbin", paths.extbin_dir()),
    )

    def __init__(
        self,
        binary_path: os.PathLike[str] | str,
        *,
        default_env: Mapping[str, str] | None = None,
        base_working_dir: os.PathLike[str] | str | None = None,
    ) -> None:
        self._binary_path = Path(binary_path)
        self._default_env = dict(default_env or {})
        self._base_working_dir = Path(base_working_dir) if base_working_dir else None

    @property
    def binary_path(self) -> Path:
        return self._binary_path

    def _prepare_environment(
        self, overrides: Mapping[str, str] | None = None
    ) -> MutableMapping[str, str]:
        env: MutableMapping[str, str] = os.environ.copy()
        env.update(self._default_env)
        if overrides:
            env.update(overrides)
        return env

    def _ensure_resource(self, workdir: Path, name: str, source: Path) -> None:
        target = workdir / name
        if target.exists():
            return
        try:
            os.symlink(source, target, target_is_directory=source.is_dir())
        except (OSError, NotImplementedError):
            if source.is_dir():
                shutil.copytree(source, target, dirs_exist_ok=True)
            else:
                shutil.copy2(source, target)

    def _prepare_workdir(self, persist: bool) -> tuple[TemporaryDirectory, Path]:
        tmp_dir = TemporaryDirectory(
            prefix="unidesign_",
            dir=str(self._base_working_dir) if self._base_working_dir else None,
            delete=not persist,
        )
        workdir = Path(tmp_dir.name)
        for name, source in self._STATIC_RESOURCES:
            self._ensure_resource(workdir, name, source)
        return tmp_dir, workdir

    def run(
        self,
        args: Sequence[str] | None = None,
        *,
        env: Mapping[str, str] | None = None,
        persist_workdir: bool = False,
    ) -> UniDesignRunResult:
        """Invoke the UniDesign binary and capture its output."""

        extra_args = tuple(args or ())
        if any(arg.startswith("--prefix") for arg in extra_args):
            raise ValueError("UniDesignRunner manages the --prefix argument automatically.")

        prefix = f"unidesign_{uuid.uuid4().hex}"
        tmp_mgr, workdir = self._prepare_workdir(persist_workdir)
        prepared_env = self._prepare_environment(env)

        try:
            argv = (str(self._binary_path), "--prefix", prefix, *extra_args)
            completed = subprocess.run(
                argv,
                cwd=str(workdir),
                env=prepared_env,
                check=False,
                capture_output=True,
                text=True,
            )
            return UniDesignRunResult(
                args=argv[1:],
                returncode=completed.returncode,
                stdout=completed.stdout,
                stderr=completed.stderr,
                workdir=workdir,
                prefix=prefix,
            )
        finally:
            if not persist_workdir:
                tmp_mgr.cleanup()


__all__ = ["UniDesignRunner", "UniDesignRunResult"]
