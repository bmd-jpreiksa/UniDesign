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
- `Structure.read_pdb(...)` underpins the high-level `Structure.from_pdb` helper now used in Python examples.
- `Structure.compute_stability(...)` uses the packaged reference tables (`unidesign.data`) by default, evaluates energies with the silent native routine, and returns a structured summary (env overrides still honoured).
- Heavy I/O helpers returning status codes surface as methods returning `(ok, message)` tuples to keep error handling pythonic while preserving original semantics.

All bindings importable through `unidesign._core`, while a thin, user-friendly Python façade will live inside `python/unidesign/api/` for higher-level workflows (`ProteinDesigner`, `EnergyScorer`, etc.).

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
