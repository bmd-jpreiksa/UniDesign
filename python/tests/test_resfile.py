"""Unit tests for the RESFILE helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

from unidesign import Resfile, ResfileEntry


def test_resfile_entry_validation() -> None:
    with pytest.raises(ValueError):
        ResfileEntry("", 10)
    with pytest.raises(ValueError):
        ResfileEntry("A", 0)


def test_resfile_to_text_round_trip(tmp_path: Path) -> None:
    resfile = Resfile(
        design=(
            ResfileEntry("A", 782, "ACDEFGHIKLMNPQRSTVWY"),
            ResfileEntry("A", 785, ("A", "C", "F")),
        ),
        repack=(ResfileEntry("B", 12),),
        catalytic=(ResfileEntry("C", 5, ["S", "T"]),),
    )

    text = resfile.to_text()
    assert text.startswith("SITES_DESIGN_START\n")
    assert "A  782" in text
    assert "A  785   ACF" in text
    assert "SITES_REPACK_START" in text
    assert "B   12" in text
    assert "SITES_CATALYTIC_END" in text

    output_path = tmp_path / "custom.resfile"
    written_path = resfile.write(output_path)
    assert written_path == output_path
    assert output_path.read_text(encoding="utf-8") == text


def test_resfile_empty_renders_blank() -> None:
    assert Resfile().to_text() == ""
