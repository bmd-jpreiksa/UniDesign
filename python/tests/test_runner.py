"""Tests covering :mod:`unidesign.runner`."""

from __future__ import annotations

from pathlib import Path
import subprocess

import pytest

from unidesign.runner import UniDesignRunner


@pytest.fixture
def dummy_binary(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """Return a dummy binary path and patch runner resources with lightweight data."""

    binary = tmp_path / "UniDesign"
    binary.write_text("#!/bin/sh\n")
    binary.chmod(0o755)

    resources: list[tuple[str, Path]] = []
    for name in ("library", "wread", "extbin"):
        resource_path = tmp_path / f"resource_{name}"
        resource_path.mkdir()
        (resource_path / "placeholder.txt").write_text(name)
        resources.append((name, resource_path))

    monkeypatch.setattr(UniDesignRunner, "_STATIC_RESOURCES", tuple(resources))
    return binary, resources


def test_runner_prepares_workspace(dummy_binary, monkeypatch: pytest.MonkeyPatch):
    binary, resources = dummy_binary
    captured: dict[str, object] = {}

    def fake_run(argv, cwd, env, check, capture_output, text):
        workdir = Path(cwd)
        prefix = argv[argv.index("--prefix") + 1]
        (workdir / f"{prefix}_marker.txt").write_text("marker")
        captured.update({
            "argv": tuple(argv),
            "cwd": workdir,
            "env": dict(env),
            "prefix": prefix,
        })
        return subprocess.CompletedProcess(argv, 0, stdout="ok", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    runner = UniDesignRunner(binary, default_env={"DEFAULT": "1"})
    result = runner.run(["--command", "Smoke"], env={"EXTRA": "2"}, persist_workdir=True)

    assert captured["argv"][:3] == (str(binary), "--prefix", result.prefix)
    assert result.args[0] == "--prefix"
    assert result.args[1] == captured["prefix"]
    assert result.stdout == "ok"
    assert result.stderr == ""
    assert result.workdir == captured["cwd"]
    assert (result.workdir / f"{result.prefix}_marker.txt").read_text() == "marker"

    for name, _ in resources:
        assert (result.workdir / name).exists()

    assert captured["env"]["DEFAULT"] == "1"
    assert captured["env"]["EXTRA"] == "2"


def test_runner_cleans_workspace_by_default(dummy_binary, monkeypatch: pytest.MonkeyPatch):
    binary, _ = dummy_binary

    def fake_run(argv, cwd, env, check, capture_output, text):
        return subprocess.CompletedProcess(argv, 0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    runner = UniDesignRunner(binary)
    result = runner.run(["--command", "Smoke"])

    assert result.returncode == 0
    assert not result.workdir.exists(), "Temporary workspace should be removed when not persisted"
