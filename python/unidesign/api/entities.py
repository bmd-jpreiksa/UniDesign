from __future__ import annotations

import os
from contextlib import ExitStack
from dataclasses import dataclass
import importlib.resources as pkg_resources
from pathlib import Path
from typing import Dict, Iterable, List, Sequence

from .. import _core

_ATOM_PARAM_RELATIVE = Path("toppar") / "param_charmm19_lk.prm"
_TOPOLOGY_RELATIVE = Path("toppar") / "top_polh19.inp"
_WEIGHT_RELATIVE = Path("wread") / "weight_all1.wgt"
_AAPP_RELATIVE = Path("eterms") / "aapropensity.nrg"
_RAMA_RELATIVE = Path("eterms") / "ramachandran.nrg"
_ROTLIB_RELATIVE = Path("rotlib") / "ALLbbdep.bin"

ENERGY_TERM_NAMES: Dict[int, str] = {
    0: "total",
    1: "reference_ALA",
    2: "reference_CYS",
    3: "reference_ASP",
    4: "reference_GLU",
    5: "reference_PHE",
    6: "reference_GLY",
    7: "reference_HIS",
    8: "reference_ILE",
    9: "reference_LYS",
    10: "reference_LEU",
    11: "reference_MET",
    12: "reference_ASN",
    13: "reference_PRO",
    14: "reference_GLN",
    15: "reference_ARG",
    16: "reference_SER",
    17: "reference_THR",
    18: "reference_VAL",
    19: "reference_TRP",
    20: "reference_TYR",
    21: "intraR_vdwatt",
    22: "intraR_vdwrep",
    23: "intraR_electr",
    24: "intraR_deslvP",
    25: "intraR_deslvH",
    26: "intraR_hbscbb_dis",
    27: "intraR_hbscbb_the",
    28: "intraR_hbscbb_phi",
    31: "interS_vdwatt",
    32: "interS_vdwrep",
    33: "interS_electr",
    34: "interS_deslvP",
    35: "interS_deslvH",
    36: "interS_ssbond",
    41: "interS_hbbbbb_dis",
    42: "interS_hbbbbb_the",
    43: "interS_hbbbbb_phi",
    44: "interS_hbscbb_dis",
    45: "interS_hbscbb_the",
    46: "interS_hbscbb_phi",
    47: "interS_hbscsc_dis",
    48: "interS_hbscsc_the",
    49: "interS_hbscsc_phi",
    51: "interD_vdwatt",
    52: "interD_vdwrep",
    53: "interD_electr",
    54: "interD_deslvP",
    55: "interD_deslvH",
    56: "interD_ssbond",
    61: "interD_hbbbbb_dis",
    62: "interD_hbbbbb_the",
    63: "interD_hbbbbb_phi",
    64: "interD_hbscbb_dis",
    65: "interD_hbscbb_the",
    66: "interD_hbscbb_phi",
    67: "interD_hbscsc_dis",
    68: "interD_hbscsc_the",
    69: "interD_hbscsc_phi",
    71: "prolig_vdwatt",
    72: "prolig_vdwrep",
    73: "prolig_electr",
    74: "prolig_deslvP",
    75: "prolig_deslvH",
    81: "prolig_hbscbb_dis",
    82: "prolig_hbscbb_the",
    83: "prolig_hbscbb_phi",
    84: "prolig_hbscsc_dis",
    85: "prolig_hbscsc_the",
    86: "prolig_hbscsc_phi",
    91: "aapropensity",
    92: "ramachandran",
    93: "dunbrack",
}


def map_energy_terms(values: Sequence[float]) -> Dict[str, float]:
    """Map indexed energy term list into a named dictionary."""
    terms: Dict[str, float] = {}
    for idx, name in ENERGY_TERM_NAMES.items():
        try:
            value = values[idx]  # type: ignore[index]
        except (IndexError, TypeError):
            continue
        terms[name] = float(value)
    return terms
def _env_override(relative: Path) -> Path | None:
    env = os.environ.get("UNIDESIGN_LIBRARY_PATH")
    if not env:
        return None
    candidate = Path(env).expanduser().resolve() / relative
    if candidate.exists():
        return candidate
    return None


def _resource_path(relative: Path, stack: ExitStack) -> Path:
    traversable = pkg_resources.files("unidesign.data")
    for part in relative.parts:
        traversable = traversable.joinpath(part)
    return Path(stack.enter_context(pkg_resources.as_file(traversable)))


def _resolve_data_file(
    relative: Path,
    stack: ExitStack,
    explicit: str | Path | None = None,
) -> Path:
    if explicit is not None:
        path = Path(explicit).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        return path

    env_path = _env_override(relative)
    if env_path:
        return env_path

    try:
        return _resource_path(relative, stack)
    except FileNotFoundError:
        repo_root = Path(__file__).resolve().parents[3]
        candidate = (repo_root / relative).resolve()
        if candidate.exists():
            return candidate
        raise


def _resolve_rotlib_file(stack: ExitStack, explicit: str | Path | None = None) -> Path:
    if explicit is not None:
        path = Path(explicit).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(path)
        return path

    env = os.environ.get("UNIDESIGN_LIBRARY_PATH")
    if env:
        base = Path(env).expanduser().resolve()
        for rel in (_ROTLIB_RELATIVE, Path("library") / _ROTLIB_RELATIVE):
            candidate = (base / rel).resolve()
            if candidate.exists():
                return candidate

    resource_candidates = (
        Path("rotlib") / _ROTLIB_RELATIVE.name,
        Path("library") / _ROTLIB_RELATIVE,
    )
    for rel in resource_candidates:
        try:
            candidate = _resource_path(rel, stack)
        except FileNotFoundError:
            candidate = None
        if candidate is not None and candidate.exists():
            return candidate

    repo_root = Path(__file__).resolve().parents[3]
    for rel in (
        _ROTLIB_RELATIVE,
        Path("library") / _ROTLIB_RELATIVE,
    ):
        candidate = (repo_root / rel).resolve()
        if candidate.exists():
            return candidate
    fallback = (repo_root / "library" / "rotlib" / _ROTLIB_RELATIVE.name).resolve()
    if fallback.exists():
        return fallback

    raise FileNotFoundError("Unable to locate ALLbbdep.bin rotamer library")


class NativeBacked:
    """Small utility base to expose the underlying pybind11 handle."""

    __slots__ = ("_handle",)

    def __init__(self, handle) -> None:
        self._handle = handle

    @property
    def handle(self):
        return self._handle


class Atom(NativeBacked):
    def __init__(
        self,
        name: str | None = None,
        chain: str | None = None,
        position: int | None = None,
        handle: _core.Atom | None = None,
    ) -> None:
        super().__init__(handle if handle is not None else _core.Atom())
        if name is not None:
            self.name = name
        if chain is not None:
            self.chain = chain
        if position is not None:
            self.position = position

    @property
    def name(self) -> str:
        return self._handle.name

    @name.setter
    def name(self, value: str) -> None:
        if len(value) > 5:
            value = value[:5]
        self._handle.name = value

    @property
    def chain(self) -> str:
        return self._handle.chain

    @chain.setter
    def chain(self, value: str) -> None:
        self._handle.chain = value

    @property
    def position(self) -> int:
        return self._handle.position

    @position.setter
    def position(self, value: int) -> None:
        self._handle.position = value


class Residue(NativeBacked):
    def __init__(
        self,
        name: str | None = None,
        chain: str | None = None,
        position: int | None = None,
        handle: _core.Residue | None = None,
    ) -> None:
        super().__init__(handle if handle is not None else _core.Residue())
        if name is not None:
            self.name = name
        if chain is not None:
            self.chain = chain
        if position is not None:
            self.position = position

    @property
    def name(self) -> str:
        return self._handle.name

    @name.setter
    def name(self, value: str) -> None:
        self._handle.name = value

    @property
    def chain(self) -> str:
        return self._handle.chain

    @chain.setter
    def chain(self, value: str) -> None:
        self._handle.chain = value

    @property
    def position(self) -> int:
        return self._handle.position

    @position.setter
    def position(self, value: int) -> None:
        self._handle.position = value

    def atom_count(self) -> int:
        return self._handle.atom_count()

    def atom_names(self) -> List[str]:
        return list(self._handle.atom_names())

    def _assert_support(self, attr: str) -> None:
        if not hasattr(self._handle, attr):
            raise RuntimeError(
                f"Residue handle missing '{attr}'. Rebuild the UniDesign extension to enable programmatic builders."
            )

    @property
    def design_type(self) -> _core.ResidueDesignType:
        self._assert_support("design_type")
        return self._handle.design_type

    @design_type.setter
    def design_type(self, value: _core.ResidueDesignType) -> None:
        self._assert_support("design_type")
        self._handle.design_type = value

    def copy_from(self, other: "Residue") -> None:
        self._assert_support("copy_from")
        self._handle.copy_from(other.handle)

    def add_atom(self, atom: "Atom") -> None:
        self._assert_support("add_atom")
        self._handle.add_atom(atom.handle)

    def atoms(self) -> List["Atom"]:
        self._assert_support("atoms")
        return [Atom(handle=atom) for atom in self._handle.atoms()]


class Chain(NativeBacked):
    def __init__(
        self,
        name: str | None = None,
        chain_type: _core.ChainType | None = None,
        handle: _core.Chain | None = None,
    ) -> None:
        super().__init__(handle if handle is not None else _core.Chain())
        if name is not None:
            self.name = name
        if chain_type is not None:
            self.chain_type = chain_type

    @property
    def name(self) -> str:
        return self._handle.name

    @name.setter
    def name(self, value: str) -> None:
        self._handle.name = value

    @property
    def chain_type(self) -> _core.ChainType:
        return self._handle.type

    @chain_type.setter
    def chain_type(self, value: _core.ChainType) -> None:
        self._handle.type = value

    def append_residue(self, residue: Residue) -> None:
        self._handle.append_residue(residue.handle)

    def residue_count(self) -> int:
        return self._handle.residue_count()

    def residues(self) -> List[Residue]:
        return [Residue(handle=self._handle.residue(i)) for i in range(self.residue_count())]

    def copy_from(self, other: "Chain") -> None:
        if not hasattr(self._handle, "copy_from"):
            raise RuntimeError(
                "Chain handle missing 'copy_from'. Rebuild the UniDesign extension to enable programmatic builders."
            )
        self._handle.copy_from(other.handle)


class Structure(NativeBacked):
    def __init__(self, name: str | None = None, handle: _core.Structure | None = None) -> None:
        super().__init__(handle if handle is not None else _core.Structure())
        if name is not None:
            try:
                self.name = name
            except RuntimeError:
                # Native setter enforces 5-character cap; ignore if exceeded.
                pass

    @property
    def name(self) -> str:
        return self._handle.name

    @name.setter
    def name(self, value: str) -> None:
        self._handle.name = value

    def chain_count(self) -> int:
        return self._handle.chain_count()

    def chains(self) -> List[Chain]:
        return [Chain(handle=self._handle.chain(i)) for i in range(self.chain_count())]

    def add_chain(self, chain: Chain) -> None:
        self._handle.add_chain(chain.handle)

    def recalc_phi_psi(self) -> None:
        if not hasattr(self._handle, "calc_phi_psi"):
            raise RuntimeError("calc_phi_psi requires the compiled UniDesign extension")
        self._handle.calc_phi_psi()

    def reset_energy_terms(self) -> None:
        if not hasattr(self._handle, "reset_energy_terms"):
            raise RuntimeError("reset_energy_terms requires the compiled UniDesign extension")
        self._handle.reset_energy_terms()

    def prepare_propensity(
        self,
        *,
        aapp_file: str | Path | None = None,
        rama_file: str | Path | None = None,
        reset: bool = False,
    ) -> None:
        if not hasattr(self._handle, "calc_propensity"):
            raise RuntimeError("calc_propensity requires the compiled UniDesign extension")
        if reset:
            self.reset_energy_terms()
        with ExitStack() as stack:
            aapp_path = _resolve_data_file(_AAPP_RELATIVE, stack, aapp_file)
            rama_path = _resolve_data_file(_RAMA_RELATIVE, stack, rama_file)
            self._handle.calc_propensity(str(aapp_path), str(rama_path))

    def prepare_dunbrack(
        self,
        *,
        rotlib_file: str | Path | None = None,
        reset: bool = False,
    ) -> None:
        if not hasattr(self._handle, "calc_dunbrack"):
            raise RuntimeError("calc_dunbrack requires the compiled UniDesign extension")
        if reset:
            self.reset_energy_terms()
        with ExitStack() as stack:
            rotlib_path = _resolve_rotlib_file(stack, rotlib_file)
            self._handle.calc_dunbrack(str(rotlib_path))

    def clone(self) -> "Structure":
        clone = Structure()
        if hasattr(clone._handle, "copy_from"):
            clone._handle.copy_from(self._handle)
        else:
            try:
                clone.name = self.name
            except RuntimeError:
                pass
            for chain in self.chains():
                new_chain = Chain()
                new_chain.copy_from(chain)
                clone.add_chain(new_chain)
        return clone

    @classmethod
    def from_pdb(
        cls,
        pdb_path: str | Path,
        *,
        atom_params: str | Path | None = None,
        topology: str | Path | None = None,
    ) -> Structure:
        """
        Load a Structure from a PDB using UniDesign's native parser.

        When atom/topology files are not supplied, the loader searches for the packaged
        UniDesign library data (or the directory pointed at by UNIDESIGN_LIBRARY_PATH).
        """

        pdb_file = Path(pdb_path).expanduser().resolve()
        if not pdb_file.exists():
            raise FileNotFoundError(pdb_file)

        structure = cls()
        if not hasattr(structure.handle, "read_pdb"):
            return cls._fallback_from_pdb(pdb_file)

        with ExitStack() as stack:
            atom_param_file = _resolve_data_file(_ATOM_PARAM_RELATIVE, stack, atom_params)
            topology_file = _resolve_data_file(_TOPOLOGY_RELATIVE, stack, topology)
            structure.handle.read_pdb(str(pdb_file), str(atom_param_file), str(topology_file))

        stem = pdb_file.stem[:5]
        try:
            structure.name = stem
        except RuntimeError:
            pass

        return structure

    def chain_ids(self) -> Iterable[str]:
        for chain in self.chains():
            yield chain.name

    @classmethod
    def _fallback_from_pdb(cls, pdb_file: Path) -> Structure:
        structure = cls(name=pdb_file.stem)
        chains: dict[str, Chain] = {}

        with pdb_file.open("r", encoding="utf-8") as handle:
            for raw_line in handle:
                if not raw_line.startswith(("ATOM", "HETATM")):
                    continue
                if len(raw_line) < 22:
                    continue
                chain_id = raw_line[21].strip() or "_"
                if chain_id not in chains:
                    chain = Chain(name=chain_id, chain_type=_core.ChainType.UNKNOWN)
                    chains[chain_id] = chain

        for chain in chains.values():
            structure.add_chain(chain)
        return structure

    def compute_stability(
        self,
        *,
        weight_file: str | Path | None = None,
        aapp_file: str | Path | None = None,
        rama_file: str | Path | None = None,
        rotlib_file: str | Path | None = None,
    ) -> "StabilityResult":
        if not hasattr(self._handle, "compute_stability"):
            raise RuntimeError("compute_stability requires the compiled UniDesign extension")

        with ExitStack() as stack:
            weight_path = _resolve_data_file(_WEIGHT_RELATIVE, stack, weight_file)
            aapp_path = _resolve_data_file(_AAPP_RELATIVE, stack, aapp_file)
            rama_path = _resolve_data_file(_RAMA_RELATIVE, stack, rama_file)
            rotlib_path = _resolve_rotlib_file(stack, rotlib_file)

            raw_terms = list(
                self._handle.compute_stability(
                    str(weight_path), str(aapp_path), str(rama_path), str(rotlib_path)
                )
            )

        per_term = {name: raw_terms[index] for index, name in ENERGY_TERM_NAMES.items()}
        return StabilityResult(total=raw_terms[0], terms=per_term, raw=raw_terms)

    def compute_binding(
        self,
        *,
        weight_file: str | Path | None = None,
        split1: str | None = None,
        split2: str | None = None,
    ) -> None:
        if not hasattr(self._handle, "compute_binding"):
            raise RuntimeError("compute_binding requires the compiled UniDesign extension")
        if (split1 is None) ^ (split2 is None):
            raise ValueError("split1 and split2 must both be provided or both omitted")
        with ExitStack() as stack:
            weight_path = _resolve_data_file(_WEIGHT_RELATIVE, stack, weight_file)
            try:
                self._handle.compute_binding(str(weight_path), split1, split2)
            except ValueError as exc:
                # Recast C++ invalid_argument into a Python ValueError with context.
                raise ValueError(str(exc)) from None


@dataclass(frozen=True)
class StabilityResult:
    total: float
    terms: Dict[str, float]
    raw: List[float]
