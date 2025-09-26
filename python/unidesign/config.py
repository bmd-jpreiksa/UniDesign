"""Structured configuration objects for UniDesign CLI commands."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Literal, Protocol, Sequence


def _as_path(value: str | Path) -> str:
    """Return a filesystem argument as a string path."""

    return str(Path(value))


def _format_bool(value: bool, *, true: str = "yes", false: str = "no") -> str:
    """Translate booleans into the UniDesign CLI "yes"/"no" vocabulary."""

    return true if value else false


class CommandConfig(Protocol):
    """Protocol implemented by configuration objects that render CLI arguments."""

    def to_cli_args(self) -> list[str]:  # pragma: no cover - structural protocol
        """Return the argument vector fragment for this command (excluding the binary)."""


@dataclass(slots=True)
class ProteinDesignConfig:
    """Configuration for the ``ProteinDesign`` command.

    The defaults mirror the native initialization performed in ``src/Main.cpp`` where the
    binary defaults to monomer design with backbone-dependent rotamers, the first chain
    selected for design, and energy weights sourced from ``wread/weight_all1.wgt``.
    """

    pdb_path: str | Path
    """Path to the input structure passed via ``--pdb``."""

    design_chains: str | None = None
    """Chain identifiers forwarded to ``--design_chains`` (defaults to ``"A"`` when omitted)."""

    mode: Literal["monomer", "ppi", "protlig", "enzyme"] = "monomer"
    """Design mode toggled through ``--monomer``, ``--ppint``, ``--protlig``, or ``--enzyme``."""

    use_bbdep_rotlib: bool | None = None
    """Override for ``--bbdep`` (defaults to ``yes``)."""

    use_input_sidechains: bool | None = None
    """Controls ``--use_input_sc``; defaults to ``yes`` to reuse input rotamers."""

    rotate_hydroxyl: bool | None = None
    """Toggles ``--rotate_hydroxyl``; defaults to ``yes`` for serine/threonine/tyrosine."""

    enable_evolution: bool = False
    """Emit ``--evolution`` to incorporate evolutionary restraints (default disabled)."""

    sequence_profile: str | Path | None = None
    """Optional profile supplied through ``--seq`` for sequence biasing."""

    profile_weight: float | None = None
    """When provided, forwarded via ``--wprof`` (native default ``1.0``)."""

    weight_file: str | Path | None = None
    """Energy weight file path for ``--wread`` (defaults to ``wread/weight_all1.wgt``)."""

    rotamer_library: str | None = None
    """Named Dunbrack/Honig library passed to ``--rotlib`` when supplied."""

    n_trajectories: int | None = None
    """Number of Monte Carlo trajectories (``--ntraj``); defaults to ``1``."""

    n_trajectory_start_index: int | None = None
    """Starting index for trajectory enumeration via ``--ntraj_start_ndx`` (default ``1``)."""

    exclude_low_prob_rotamers_cutoff: float | None = None
    """Cut-off for ``--excl_low_prob``; the binary defaults to ``0.03``."""

    ppi_shell1: float | None = None
    """Distance for ``--ppi_shell1`` (first interface shell, default ``5.0`` Å)."""

    ppi_shell2: float | None = None
    """Distance for ``--ppi_shell2`` (second interface shell, default ``8.0`` Å)."""

    pli_shell1: float | None = None
    """Distance for ``--pli_shell1`` (ligand shell 1, default ``5.0`` Å)."""

    pli_shell2: float | None = None
    """Distance for ``--pli_shell2`` (ligand shell 2, default ``8.0`` Å)."""

    clash_ratio: float | None = None
    """Steric clash ratio threshold supplied through ``--clash_ratio`` (default ``0.6``)."""

    design_type: Literal["natro", "nataa", "allaa", "allaaxc"] | None = None
    """Default rotamer type selection for ``--init_rotype`` (native default ``natro``)."""

    resfile_path: str | Path | None = None
    """Path to a resfile consumed via ``--resfile``."""

    wildtype_only: bool = False
    """Emit ``--wildtype_only`` to restrict positions to native amino acids."""

    interface_only: bool = False
    """Emit ``--interface_only`` to limit design to interface residues."""

    seed_from_native_sequence: bool = False
    """Emit ``--seed_from_nat_seq`` to initialize trajectories from the input sequence."""

    exclude_cysteine_rotamers: bool = False
    """Emit ``--excl_cys_rots`` to drop cysteine rotamers."""

    write_hydrogen: bool | None = None
    """Explicitly control hydrogen output via ``--show_hydrogen`` (default ``yes``)."""

    rotate_reference_residues: str | None = None
    """Residues fed to ``--within_residues`` for distance-based selection."""

    reference_distance: float | None = None
    """Maximum Å radius passed through ``--within_range`` (default ``9.0`` Å)."""

    binding_weight: float | None = None
    """Scalar for binding term weighting using ``--wbind`` (default ``1.0``)."""

    ligand_parameter_path: str | Path | None = None
    """Ligand parameter override forwarded via ``--lig_param`` (default ``LIG_PARAM.prm``)."""

    ligand_topology_path: str | Path | None = None
    """Ligand topology override supplied through ``--lig_topo`` (default ``LIG_TOPO.inp``)."""

    ligand_constraint_path: str | Path | None = None
    """Constraint definition path for ``--lig_catacons`` (default ``LIG_CATACONS.txt``)."""

    ligand_placement_path: str | Path | None = None
    """Placement configuration for ``--lig_placing`` (default ``LIG_PLACING.txt``)."""

    ligand_pose_input: str | Path | None = None
    """Ligand pose input file consumed via ``--read_lig_poses``."""

    ligand_pose_output: str | Path | None = None
    """Ligand pose output file produced with ``--write_lig_poses``."""

    ligand_orientation_screen: str | Path | None = None
    """Orientation screen manifest forwarded to ``--scrn_by_orien``."""

    ligand_vdw_percentile: float | None = None
    """Percentile cutoff for ``--scrn_by_vdw_pctl`` (default ``0.5`` corresponding to 50%)."""

    ligand_rmsd_cutoff: float | None = None
    """Ångström cutoff for ``--scrn_by_rmsd`` (default ``0.5`` Å)."""

    prefix: str | None = None
    """Output prefix forwarded through ``--prefix`` (defaults to ``UniDesign``)."""

    def _mode_args(self) -> Iterable[str]:
        if self.mode == "monomer":
            return ()
        if self.mode == "ppi":
            return ("--ppint",)
        if self.mode == "protlig":
            return ("--protlig",)
        return ("--enzyme",)

    def _append_yes_no(self, args: list[str], flag: str, value: bool | None) -> None:
        if value is None:
            return
        args.extend((flag, _format_bool(value)))

    def to_cli_args(self) -> list[str]:
        args: list[str] = ["--command", "ProteinDesign", "--pdb", _as_path(self.pdb_path)]
        args.extend(self._mode_args())

        if self.design_chains:
            args.extend(("--design_chains", self.design_chains))

        self._append_yes_no(args, "--bbdep", self.use_bbdep_rotlib)
        self._append_yes_no(args, "--use_input_sc", self.use_input_sidechains)
        self._append_yes_no(args, "--rotate_hydroxyl", self.rotate_hydroxyl)

        if self.enable_evolution:
            args.append("--evolution")
        if self.sequence_profile is not None:
            args.extend(("--seq", _as_path(self.sequence_profile)))
        if self.profile_weight is not None:
            args.extend(("--wprof", str(self.profile_weight)))
        if self.weight_file is not None:
            args.extend(("--wread", _as_path(self.weight_file)))
        if self.rotamer_library is not None:
            args.extend(("--rotlib", self.rotamer_library))
        if self.n_trajectories is not None:
            if self.n_trajectories <= 0:
                raise ValueError("n_trajectories must be positive")
            args.extend(("--ntraj", str(self.n_trajectories)))
        if self.n_trajectory_start_index is not None:
            if self.n_trajectory_start_index <= 0:
                raise ValueError("n_trajectory_start_index must be positive")
            args.extend(("--ntraj_start_ndx", str(self.n_trajectory_start_index)))
        if self.exclude_low_prob_rotamers_cutoff is not None:
            args.extend(("--excl_low_prob", str(self.exclude_low_prob_rotamers_cutoff)))
        if self.ppi_shell1 is not None:
            args.extend(("--ppi_shell1", str(self.ppi_shell1)))
        if self.ppi_shell2 is not None:
            args.extend(("--ppi_shell2", str(self.ppi_shell2)))
        if self.pli_shell1 is not None:
            args.extend(("--pli_shell1", str(self.pli_shell1)))
        if self.pli_shell2 is not None:
            args.extend(("--pli_shell2", str(self.pli_shell2)))
        if self.clash_ratio is not None:
            args.extend(("--clash_ratio", str(self.clash_ratio)))
        if self.design_type is not None:
            args.extend(("--init_rotype", self.design_type))
        if self.resfile_path is not None:
            args.extend(("--resfile", _as_path(self.resfile_path)))
        if self.wildtype_only:
            args.append("--wildtype_only")
        if self.interface_only:
            args.append("--interface_only")
        if self.seed_from_native_sequence:
            args.append("--seed_from_nat_seq")
        if self.exclude_cysteine_rotamers:
            args.append("--excl_cys_rots")
        self._append_yes_no(args, "--show_hydrogen", self.write_hydrogen)
        if self.rotate_reference_residues is not None:
            args.extend(("--within_residues", self.rotate_reference_residues))
        if self.reference_distance is not None:
            args.extend(("--within_range", str(self.reference_distance)))
        if self.binding_weight is not None:
            args.extend(("--wbind", str(self.binding_weight)))
        if self.ligand_parameter_path is not None:
            args.extend(("--lig_param", _as_path(self.ligand_parameter_path)))
        if self.ligand_topology_path is not None:
            args.extend(("--lig_topo", _as_path(self.ligand_topology_path)))
        if self.ligand_constraint_path is not None:
            args.extend(("--lig_catacons", _as_path(self.ligand_constraint_path)))
        if self.ligand_placement_path is not None:
            args.extend(("--lig_placing", _as_path(self.ligand_placement_path)))
        if self.ligand_pose_input is not None:
            args.extend(("--read_lig_poses", _as_path(self.ligand_pose_input)))
        if self.ligand_pose_output is not None:
            args.extend(("--write_lig_poses", _as_path(self.ligand_pose_output)))
        if self.ligand_orientation_screen is not None:
            args.extend(("--scrn_by_orien", _as_path(self.ligand_orientation_screen)))
        if self.ligand_vdw_percentile is not None:
            args.extend(("--scrn_by_vdw_pctl", str(self.ligand_vdw_percentile)))
        if self.ligand_rmsd_cutoff is not None:
            args.extend(("--scrn_by_rmsd", str(self.ligand_rmsd_cutoff)))
        if self.prefix is not None:
            args.extend(("--prefix", self.prefix))

        return args


@dataclass(slots=True)
class ComputeStabilityConfig:
    """Configuration for the ``ComputeStability`` command.

    Native defaults enable the backbone-dependent Dunbrack library (``--bbdep=yes``) and
    read weights from ``wread/weight_all1.wgt`` while expecting ``model_1.pdb`` when a PDB
    path is not supplied.
    """

    pdb_path: str | Path
    """Structure to analyse with ``--pdb``."""

    use_bbdep_rotlib: bool | None = None
    """Optional override for ``--bbdep`` (defaults to ``yes``)."""

    rotamer_library: str | None = None
    """Named text rotamer library passed through ``--rotlib`` if provided."""

    weight_file: str | Path | None = None
    """Alternate weights via ``--wread`` (default ``wread/weight_all1.wgt``)."""

    def to_cli_args(self) -> list[str]:
        args: list[str] = ["--command", "ComputeStability", "--pdb", _as_path(self.pdb_path)]
        if self.use_bbdep_rotlib is not None:
            args.extend(("--bbdep", _format_bool(self.use_bbdep_rotlib)))
        if self.rotamer_library is not None:
            args.extend(("--rotlib", self.rotamer_library))
        if self.weight_file is not None:
            args.extend(("--wread", _as_path(self.weight_file)))
        return args


def _validate_split_parts(part1: str, part2: str) -> None:
    overlaps = set(part1) & set(part2)
    if overlaps:
        overlap = ", ".join(sorted(overlaps))
        raise ValueError(
            f"split_chains groups must be disjoint; overlapping chain(s): {overlap}"
        )


@dataclass(slots=True)
class ComputeBindingConfig:
    """Configuration for the ``ComputeBinding`` command.

    The binary defaults to splitting the complex into ``"AB"`` and ``"C"`` when chain
    splitting is engaged. The validation below mirrors the runtime guard that rejects
    overlapping chain identifiers in ``src/Main.cpp``.
    """

    pdb_path: str | Path
    """Structure analysed via ``--pdb``."""

    split_part1: str = "AB"
    """First chain group for ``--split_chains`` (default ``"AB"``)."""

    split_part2: str = "C"
    """Second chain group for ``--split_chains`` (default ``"C"``)."""

    def __post_init__(self) -> None:
        _validate_split_parts(self.split_part1, self.split_part2)

    def to_cli_args(self) -> list[str]:
        args = ["--command", "ComputeBinding", "--pdb", _as_path(self.pdb_path)]
        split_value = f"{self.split_part1},{self.split_part2}"
        args.extend(("--split_chains", split_value))
        return args


def _normalise_atom_triplet(atoms: Sequence[str]) -> tuple[str, str, str]:
    if len(atoms) != 3:
        raise ValueError("Exactly three atom names are required for init_3atoms")
    return tuple(atoms)  # type: ignore[return-value]


@dataclass(slots=True)
class MakeLigParamConfig:
    """Configuration for the ``MakeLigParamAndTopo`` command.

    By default the binary expects a ligand ``mol2.mol2`` file and seeds topology generation
    with the ``C1``, ``C2``, ``C3`` atom triplet while writing ``LIG_PARAM.prm`` and
    ``LIG_TOPO.inp`` outputs.
    """

    mol2_path: str | Path
    """Ligand MOL2 file supplied via ``--mol2``."""

    ligand_parameter_path: str | Path = "LIG_PARAM.prm"
    """Output parameter file for ``--lig_param`` (native default ``LIG_PARAM.prm``)."""

    ligand_topology_path: str | Path = "LIG_TOPO.inp"
    """Output topology file for ``--lig_topo`` (native default ``LIG_TOPO.inp``)."""

    initial_atoms: Sequence[str] = ("C1", "C2", "C3")
    """Ordered atom triplet forwarded via ``--init_3atoms`` (default ``C1,C2,C3``)."""

    def __post_init__(self) -> None:
        _normalise_atom_triplet(self.initial_atoms)

    def to_cli_args(self) -> list[str]:
        atom1, atom2, atom3 = _normalise_atom_triplet(self.initial_atoms)
        atom_arg = ",".join((atom1, atom2, atom3))
        return [
            "--command",
            "MakeLigParamAndTopo",
            "--mol2",
            _as_path(self.mol2_path),
            "--lig_param",
            _as_path(self.ligand_parameter_path),
            "--lig_topo",
            _as_path(self.ligand_topology_path),
            "--init_3atoms",
            atom_arg,
        ]


__all__ = [
    "CommandConfig",
    "ProteinDesignConfig",
    "ComputeStabilityConfig",
    "ComputeBindingConfig",
    "MakeLigParamConfig",
]

