"""Utilities for describing and writing UniDesign RESFILE inputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


@dataclass(slots=True)
class ResfileEntry:
    """A single residue entry within a RESFILE section."""

    chain: str
    residue: int
    allowed: str | tuple[str, ...] | None = None

    def __post_init__(self) -> None:
        if not self.chain:
            raise ValueError("chain must contain at least one character")
        if self.residue <= 0:
            raise ValueError("residue indices must be positive")
        allowed = self.allowed
        if allowed is not None and not isinstance(allowed, str):
            self.allowed = tuple(allowed)

    def format_line(self) -> str:
        """Render the entry in UniDesign RESFILE format."""

        line = f"{self.chain} {self.residue:4d}"
        allowed = self.allowed
        if allowed:
            allowed_str = "".join(allowed) if isinstance(allowed, tuple) else allowed
            line = f"{line}   {allowed_str}"
        return line


@dataclass(slots=True)
class Resfile:
    """Structured representation of a UniDesign RESFILE."""

    design: Sequence[ResfileEntry] | None = None
    repack: Sequence[ResfileEntry] | None = None
    catalytic: Sequence[ResfileEntry] | None = None

    def _section_lines(
        self, header: str, footer: str, entries: Iterable[ResfileEntry] | None
    ) -> list[str]:
        if not entries:
            return []
        lines = [header]
        lines.extend(entry.format_line() for entry in entries)
        lines.append(footer)
        return lines

    def to_text(self) -> str:
        """Render the RESFILE as a string."""

        lines: list[str] = []
        lines.extend(
            self._section_lines("SITES_DESIGN_START", "SITES_DESIGN_END", self.design)
        )
        lines.extend(
            self._section_lines("SITES_REPACK_START", "SITES_REPACK_END", self.repack)
        )
        lines.extend(
            self._section_lines(
                "SITES_CATALYTIC_START", "SITES_CATALYTIC_END", self.catalytic
            )
        )
        if not lines:
            return ""
        return "\n".join(lines) + "\n"

    def write(self, path: str | Path) -> Path:
        """Write the RESFILE to ``path`` and return the resulting :class:`Path`."""

        target = Path(path)
        if target.parent and not target.parent.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(self.to_text(), encoding="utf-8")
        return target


__all__ = ["Resfile", "ResfileEntry"]

