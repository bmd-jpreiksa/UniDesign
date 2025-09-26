"""Convenience exports for the UniDesign Python bindings."""

from importlib import metadata as _metadata

try:  # pragma: no cover - fallback for editable installs
    __version__ = _metadata.version(__name__)
except _metadata.PackageNotFoundError:  # pragma: no cover
    __version__ = "0.0.0"

from ._core import (
    StructureHandle,
    SequenceHandle,
    compute_binding,
    compute_structure_stability,
    get_cutoff,
    get_flag,
    get_integer,
    get_path,
    list_cutoffs,
    list_flags,
    list_integers,
    list_paths,
    load_energy_weights,
    load_structure_from_pdb,
    protein_design,
    set_cutoff,
    set_flag,
    set_integer,
    set_path,
)

__all__ = [
    "StructureHandle",
    "SequenceHandle",
    "compute_binding",
    "compute_structure_stability",
    "get_cutoff",
    "get_flag",
    "get_integer",
    "get_path",
    "list_cutoffs",
    "list_flags",
    "list_integers",
    "list_paths",
    "load_energy_weights",
    "load_structure_from_pdb",
    "protein_design",
    "set_cutoff",
    "set_flag",
    "set_integer",
    "set_path",
    "__version__",
]
