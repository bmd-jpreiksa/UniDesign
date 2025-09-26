"""Sphinx configuration for the UniDesign Python wrapper documentation."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
import sys

# Ensure the bundled Python package can be discovered when using autodoc
PROJECT_ROOT = Path(__file__).resolve().parents[2]
PYTHON_DIR = PROJECT_ROOT / "python"
if PYTHON_DIR.exists():
    sys.path.insert(0, str(PYTHON_DIR))

project = "UniDesign Python Wrapper"
copyright = f"{datetime.utcnow():%Y}, UniDesign developers"
author = "UniDesign developers"

extensions: list[str] = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

templates_path: list[str] = []
exclude_patterns: list[str] = []

html_theme = "alabaster"
html_static_path: list[str] = []
