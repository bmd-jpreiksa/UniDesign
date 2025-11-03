# Agent Notes – PyBindings

## Current pybind11 surface
- `Atom`, `Residue`, `Chain`, and `Structure` handles are live and mirror the native structs; enums (`ChainType`, `ResidueDesignType`) and the `run_monomer_design` orchestration are exposed.
- High-level façade lives in `python/unidesign/api/`, with `DesignProtein` wiring the monomer workflow and optional ligand payloads; packaged data (`toppar/`, `eterms/`, `wread/`, `rotlib/`) back the defaults.
- Ligand helpers exist as a thin container (`Ligand` class) and the binding currently forwards `Structure.read_mol2`, `Structure.read_small_mol_rotamers`, and the ligand options map.

## Protein–ligand parity gap
- Running the CLI example (`UniDesign --command=ProteinDesign --protlig …`) produces `UniDesign_bestseqs.txt` with total energy ≈ -648, physical ≈ -595, binding ≈ -53 and, crucially, the `StructureShowDesignSites` log lists the ligand as `site 44 … smallmol`.
- The Python example (`python_examples/protein_ligand_design.py`) yields total energy ≈ -669, physical ≈ -646, binding ≈ -23 and its design-site dump never mentions a `smallmol` entry. This indicates the ligand never becomes an active design site in the binding-driven workflow, so protein residues are repacked around a static ligand pose, leading to a different energy decomposition.
- Root cause: `RunMonomerDesignWorkflow` does not mirror the CLI fallback that writes `LIG_POSES2.pdb` and re-reads it through `StructureReadSmallMolRotamers`. We only pre-load the ligand via `Structure.read_mol2`, and unless the caller supplies an explicit conformer file, the small-molecule site is never promoted to `Type_DesType_SmallMol`. As a result, ligand rotamers are absent, binding interactions are under-sampled, and the reported binding term diverges.
- Evidence trail:
  - Python run log (`/tmp/python_design.log`) lacks the `site … smallmol` line that appears in `/tmp/cli_interface.log`.
  - Energy totals reported by `DesignProtein.best_sequence_energy` versus `example/ProteinLigandInteraction/1r091_BS01_JEN/UniDesign_bestseqs.txt`.
- Status: Mirroring the CLI pose fallback and flipping `FLAG_MONOMER` to `FALSE` for ligand jobs restores parity (`python_examples/protein_ligand_design.py` now reports total ≈ -647.15, physical ≈ -595.27, binding ≈ -51.88 to match the CLI run).

## Follow-up hints
- Inside `RunMonomerDesignWorkflow` replicate the CLI steps around ligand preparation: set `FILE_LIG_POSES_IN/OUT` within the scratch directory, force `FLAG_LIG_POSES = TRUE`, write the default single-pose file when no conformers are supplied, and invoke `StructureReadSmallMolRotamers` so the ligand residue is promoted to `Type_DesType_SmallMol`.
- After bringing the ligand rotamers online, rerun both the CLI example and `python_examples/protein_ligand_design.py` with the shortened MOL2 to confirm the sequence/energy columns match.
- Longer term, add coverage around ligand attach/detach in the Python test suite to catch regressions in the binding layer.
