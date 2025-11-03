# UniDesign Repository Overview

## Project scope
UniDesign is a protein design engine written in C++ (`src/`). The downstream roadmap is to surface the complete functionality as a Python package backed by `pybind11` so that every computational building block is discoverable from a well-structured Python API.

## Native code layout
The C++ source files are organised by domain concepts:

| Module | Responsibility highlights | Key public types/functions |
| --- | --- | --- |
| `Atom.*`, `AtomParamsSet.*` | Atom topology, parameter loading, array utilities | `Atom`, `AtomArray`, `AtomParamsSetLoadFromFile` |
| `Residue.*`, `ResidueTopology.*` | Residue representation and topology templates | `Residue`, `ResiTopo`, `ResidueCalcXYZ` |
| `Chain.*`, `Structure.*` | Protein/nucleic acid chains and whole assemblies | `Chain`, `Structure`, `StructureReadPDB`, `StructureShowInPDBFormat` |
| `Rotamer*.{cpp,h}` | Rotamer libraries, optimisation, builder utilities | `Rotamer`, `BBdepRotamerLib`, `RotamerOptimizerRun` |
| `Energy*.{cpp,h}` | Energy evaluation core and supporting tables | `EnergyFunction`, `EnergyMatrix`, `ComputeStructureStability` |
| `ProgramFunction.*`, `ProgramPreprocess.*` | High-level orchestration used by CLI binaries | Workflow functions such as `BuildMutant`, `ProteinSiteAddDesignSite` |
| `Evolution*`, `Evo*` | Evolutionary sequence analysis helpers | `Evolution`, `EvoUtility`, `EvoSeqAlign`, etc. |
| `Utility.*`, `ErrorTracker.*`, `Getopt.*` | Shared utilities, logging, CLI parsing | `StringArray`, `ErrorTracker`, command-line helpers |
| `SmallMol*.{cpp,h}` | Small molecule support and parameter loading | `SmallMol`, `SmallMolEEF1`, `SmallMolParAndTopoLoad` |

Every header follows a procedural style (structs + free functions), so the binding layer will promote them into Python-friendly classes without altering native code.

## Binding strategy
The `pybind11` layer will live under `python/unidesign/_bindings/`. Each native module receives a dedicated binding translation unit to keep compile times manageable. Example layout:

```
python/unidesign/_bindings/
  atom.cpp        -> wraps Atom/AtomArray APIs
  residue.cpp     -> wraps residue/topology helpers
  chain.cpp       -> wraps Chain and related functions
  structure.cpp   -> exposes Structure, Protein site management, IO helpers
  energy.cpp      -> exposes scoring & optimisation routines
  evolution.cpp   -> exposes Evo* analytics
  smallmol.cpp    -> exposes ligand utilities
  module.cpp      -> aggregates the per-domain submodules
```

Within each translation unit:

- `py::class_` is used for struct wrappers so Python owns the C++ memory and calls the existing `*Create`/`*Destroy` helpers in constructors/destructors.
- Enumerations (`Type_AtomPolarity`, `Type_Chain`, etc.) are mapped with `py::enum_`.
- Functions that act on pointers (e.g. `StructureReadPDB`) become instance methods on the owning Python class to keep the API object-oriented.
- `Structure.read_pdb(...)` underpins the high-level `Structure.from_pdb` helper; the binding now recomputes backbone torsions (`StructureCalcPhiPsi`) so energy tables match CLI behaviour.
- `Structure.compute_stability(...)` loads packaged reference tables and the Dunbrack BB-dependent library (when available) before delegating to the native scorer; results now match `UniDesign --command=ComputeStability` bit-for-bit.
- `Structure.run_monomer_design(...)` executes the monomer workflow completely in-memory (no CLI shell-out), wiring the native rotamer builders, simulated annealing loop, and best-decoy extraction into a single call that returns the final sequence/energy summary plus optional structure snapshots.
- Chains and residues expose `copy_from`, `design_type`, and `atoms()/add_atom()` utilities so new structures can be assembled programmatically; `Structure.clone()` builds deep copies from Python.
- Heavy I/O helpers returning status codes surface as methods returning `(ok, message)` tuples to keep error handling pythonic while preserving original semantics.
- `ChainType` and `ResidueDesignType` enums are exported alongside the structural handles for convenient use in Python.
- Additional native helpers (e.g., `StructureCalcAminoAcidDunbrackEnergy`) are wired internally to keep the Python façade aligned with CLI workflows.
- Python examples:
  - `python_examples/stability.py` runs `Structure.compute_stability()` on 1igd using the packaged data files and reports the same totals as `UniDesign --command=ComputeStability`.
  - `python_examples/binding_energy.py` demonstrates calling `Structure.compute_binding` on the 1e44 and 1ay7 complexes (with optional chain splitting). Invalid chain identifiers now raise a `ValueError`, which the example reports before continuing.
  - `python_examples/design_monomer.py` exercises the binding-driven `DesignProtein` façade (which delegates to `Structure.run_monomer_design`) to reproduce the `MonomerDesign/1agy` reference trajectory without spawning the CLI binary.
  - `python_examples/design_monomer_domain.py` shows how to build a `DesignDomain` in Python, mixing fully flexible sites with amino-acid-restricted positions (e.g., limiting chain A position 35 to `AVIL`) before running the in-memory workflow.
- Ligand bindings are **not** exposed yet. The current Python façade cannot create or parameterise small-molecule residues, and `DesignProtein` rejects ligand payloads. The native workflow still expects pre-generated parameter/topology snippets on disk.

### Ligand workflow status

- ✅ Native C++ now deep-copies design sites correctly (`StructureCopy` rebinds design-site residue pointers), eliminating crashes when the CLI handles ligand jobs.
- ⚠️ Missing pybind surface: there is no `Ligand` handle/class in `_core`, so Python users cannot supply MOL2/parameter/topology data programmatically.
- ⚠️ `Structure.run_monomer_design` ignores ligand-related options because `PyMonomerDesignOptions` does not carry ligand payloads into `RunMonomerDesignWorkflow`.
- ⚠️ Example gap: `python_examples/protein_ligand_design.py` cannot be executed; it imports a non-existent `Ligand` helper and therefore fails before reaching the workflow.
- 📋 To finish ligand exposure we need:
  1. A `LigandHandle` binding (MOL2 ingestion, parameter/topology setters, conformer loading).
  2. Extended `PyMonomerDesignOptions`/`RunMonomerDesignWorkflow` to persist ligand data into the temporary working directory and enable the `FLAG_PROT_LIG` code path.
  3. A high-level `Ligand` façade in `python/unidesign/api/` plus documentation/tests, after which `python_examples/protein_ligand_design.py` can be re-enabled.

All bindings importable through `unidesign._core`, while a thin, user-friendly Python façade will live inside `python/unidesign/api/` for higher-level workflows (`ProteinDesigner`, `EnergyScorer`, etc.).

### Current pybind11 coverage

| Component | Python exposure | Status |
| --- | --- | --- |
| `AtomHandle` | Name/chain/position setters, cartesian coordinates, B-factor accessors | ✅ implemented |
| `ResidueHandle` | Name/chain/position setters, design type enum, `copy_from`, `add_atom`, `atoms()`, atom listing | ✅ implemented |
| `ChainHandle` | Name/type setters, `append_residue`, residue lookup, `copy_from` | ✅ implemented |
| `StructureHandle` | Name setter, chain management, `read_pdb` (phi/psi recomputation), `compute_stability` (weights, tables, rotlib), `compute_binding` (optional weight/splitting), `run_monomer_design` (full native workflow with optional resfile/domains), `calc_phi_psi`, `calc_propensity`, `calc_dunbrack`, `reset_energy_terms`, `copy_from` | ✅ implemented |
| Enums | `ChainType`, `ResidueDesignType` | ✅ implemented |
| Functions | `print_version` passthrough | ✅ implemented |
| Energy helpers | `ComputeStructureStabilitySilent`, `ComputeStructureStabilityByBBdepRotLib2` (invoked internally) | ✅ leveraged |
| Remaining domains | Rotamers, energy matrix construction, CLI workflow orchestration, evolution utilities, small-molecule support, broader IO | ⏳ pending |

### Ligand workflow notes

- Ligands live inside the main `Structure` as residues on chains marked `Type_Chain_SmallMol`. The design site for the ligand uses `Type_DesType_SmallMol`, so rotamer packing and energy evaluation reuse the same machinery as protein residues.
- CHARMM-style atom parameters (`param` file) and topology (`top` file) must be supplied for every ligand. `GenerateSmallMolParameterAndTopologyFromMol2` builds these from a MOL2 by assigning EEF1 atom types, charges, Lennard–Jones radii/epsilons, solvation parameters (ΔG_free, volume, λ), hydrogen-bond roles, and generating IC records.
- Protein residues within configurable shells (5 Å mutable, 8 Å repackable by default) are automatically flagged and receive rotamer sets so side-chain design happens around the ligand. Protein–ligand scores use `EnergyResidueAndLigandResidue`, which relies on the ligand atom parameters to populate vdW, electrostatic, desolvation, and H-bond terms (energy indices 71–86).
- Ligand conformers enter the design site either by reading pose files (`StructureReadSmallMolRotamers`) or by in-silico placement. `StructureGenerateSmallMolRotamers` executes user-defined placing rules and catalytic constraints, sampling ligand torsions/translations against a truncated backbone and enforcing distance/angle/torsion checks between ligand atoms and protein residues.
- Catalytic constraint files support arbitrary atom names (including side chains and pseudo atoms). During placement UniDesign inspects the relevant rotamers to locate the requested atoms—constraints are not limited to backbone Cα atoms.
- Desolvation contributions stem from per-atom EEF1 parameters (volume, ΔG_free, λ) attached to the ligand atoms; these are exported in the ligand `.prm` and consumed by `LKDesolvationEnergyAtomAndAtom` during scoring.

## Python package layout

```
python/
  unidesign/
    __init__.py        # re-export high level API
    _core.pyi          # optional stub file for IDEs (generated post-build)
    api/
      __init__.py
      structure.py     # object-oriented façade around Structure bindings
      energy.py        # energy evaluation helpers
      evolution.py
    data/              # packaged reference tables (toppar, eterms, wread, ...)
```

The package directory doubles as the source tree for `setuptools`. `pyproject.toml` declares the build backend and compiles a single extension `unidesign._core` from every file in `src/`. The build uses `pybind11.setup_helpers.Pybind11Extension` with sensible optimisation flags (`-O3`, `-ffast-math`) mirroring `build.sh`.

## Build & installation pipeline
1. **Conda environment** (`environment.yml`) supplies compilers and Python dependencies (`python>=3.9`, `pybind11`, `numpy`, `cmake`, `ninja`).
2. `pip install .` triggers the `pybind11` build, compiling all C++ sources into the `unidesign._core` extension.
3. Console entry points defined in `pyproject.toml` expose existing CLI flows (`unidesign-protein-design`, etc.) while Pythonists interact with the object-oriented API.

## Testing hooks
The `python/tests/` package will target both the raw bindings and the higher-level façade. Smoke tests will validate:

- Packaging metadata (import succeeds, version available).
- Structure parsing via bindings (`Structure.read_pdb("example/...")`).
- Energy evaluation on sample data.

Native regression tests (if any) can be wired in later via `CTest` or a custom harness.

---

This plan allows incremental binding delivery: start with the structural core (`Atom`, `Residue`, `Chain`, `Structure`), expand to energy and evolutionary helpers, and keep documentation in this file in sync as modules land.
