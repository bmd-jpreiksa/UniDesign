"""Python bindings and high-level API for the UniDesign engine."""

from . import _core
from .api import Atom, Chain, Residue, StabilityResult, Structure
from .api.design import DesignDomain, DesignProtein

ChainType = _core.ChainType
try:
    ResidueDesignType = _core.ResidueDesignType
except AttributeError:  # older extension build; provide placeholder
    ResidueDesignType = None
print_version = _core.print_version

__all__ = [
    "Atom",
    "Residue",
    "Chain",
    "Structure",
    "DesignDomain",
    "DesignProtein",
    "StabilityResult",
    "ChainType",
    "ResidueDesignType",
    "print_version",
    "_core",
]
__version__ = "0.1.0"
