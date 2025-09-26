"""Pytest fixtures for the UniDesign Python wrapper."""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure the ``unidesign`` package is importable when tests run from the
# repository root without an installed wheel.
PACKAGE_ROOT = Path(__file__).resolve().parents[1]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))
