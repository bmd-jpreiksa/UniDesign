UniDesign Python Wrapper Guide
==============================

.. contents:: Table of contents
   :local:
   :depth: 2

Overview
--------

The :mod:`unidesign` Python package wraps the UniDesign command line binary with
high-level helpers that simplify the preparation of command arguments, the
management of temporary working directories, and the collection of generated
artefacts. This guide explains how to set up the wrapper, documents the
available commands, clarifies how temporary files are handled, and highlights
options for persisting results.

Setup
-----

Binary discovery
~~~~~~~~~~~~~~~~

#. Build the UniDesign executable using the provided :mod:`build.sh` helper or
   by invoking your preferred C++ toolchain manually. The wrapper will look for
   ``UniDesign`` or ``UniDesign.exe`` in the repository root, ``build/`` and
   ``bin/`` directories.
#. When the executable is not present or lacks the executable bit, the helper
   :func:`unidesign.discover_binary` raises
   :class:`unidesign.BinaryDiscoveryError` describing the attempted search
   locations. The error can be surfaced to users of CLI tools to instruct them
   to build the project first.

Python environment
~~~~~~~~~~~~~~~~~~

The wrapper ships as a plain Python package under ``python/unidesign`` with no
external dependencies. Add the ``python`` directory to :envvar:`PYTHONPATH` or
install the package into a virtual environment using ``pip install -e
python``. Tests and examples in this repository run with Python 3.9+, although
newer interpreters are supported.

Initialising a runner
~~~~~~~~~~~~~~~~~~~~~

The :class:`unidesign.UniDesignRunner` class is responsible for invoking the
binary while ensuring required resource directories (``library/``, ``wread/``
and ``extbin/``) are available in the working directory. Typical initialisation
looks like:

.. code-block:: python

   from unidesign import UniDesignRunner, discover_binary

   binary_path = discover_binary()
   runner = UniDesignRunner(
       binary_path,
       default_env={"OMP_NUM_THREADS": "8"},
   )

The optional ``default_env`` mapping is merged with the current process
environment for every invocation. Per-call overrides can be supplied via the
``env=`` keyword argument of :meth:`UniDesignRunner.run`.

Command usage
-------------

Monomer and interface design
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use :class:`unidesign.ProteinDesignConfig` to describe protein-design settings
and :class:`unidesign.jobs.ProteinDesignJob` to execute the command. For
example, the script at :file:`examples/python/monomer_design.py` runs a
single-trajectory monomer design against the bundled
:file:`example/MonomerDesign/1igd/1igd.pdb` structure:

.. code-block:: python

   from pathlib import Path

   from unidesign import ProteinDesignConfig, UniDesignRunner, discover_binary
   from unidesign.jobs import ProteinDesignJob
   from unidesign.paths import project_root

   pdb_path = project_root() / "example/MonomerDesign/1igd/1igd.pdb"
   config = ProteinDesignConfig(pdb_path=pdb_path, design_chains="A", n_trajectories=1)

   runner = UniDesignRunner(discover_binary())
   job = ProteinDesignJob(runner, config)
   result = job.run()
   try:
       if result.best_structure:
           print("Best structure saved to", result.best_structure.path)
   finally:
       result.close()

The job wrapper returns a :class:`unidesign.jobs.ProteinDesignResult` bundling
the captured stdout/stderr, relocated artefacts (self-energy report, sequence
sets, structures, ligand poses), and a ``close`` method that removes the
workspace on demand.

Energy scoring
~~~~~~~~~~~~~~

The wrapper offers dedicated helpers for the ``ComputeStability`` and
``ComputeBinding`` commands. The script at
:file:`examples/python/energy_scoring.py` demonstrates both jobs:

.. code-block:: python

   from unidesign import ComputeBindingConfig, ComputeStabilityConfig, UniDesignRunner, discover_binary
   from unidesign.jobs import BindingComputationJob, StabilityComputationJob
   from unidesign.paths import project_root

   runner = UniDesignRunner(discover_binary())
   pdb_path = project_root() / "example/ProteinProteinInteractionDesign/1ay7/1ay7.pdb"

   stability = StabilityComputationJob(runner, ComputeStabilityConfig(pdb_path=pdb_path))
   stability_result = stability.run()
   try:
       if stability_result.rotamer_list:
           print("Rotamer list saved to", stability_result.rotamer_list.path)
   finally:
       stability_result.close()

   binding = BindingComputationJob(runner, ComputeBindingConfig(pdb_path=pdb_path))
   binding_result = binding.run()
   binding_result.close()

Ligand parameter generation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ligand parameterisation is handled by
:class:`unidesign.MakeLigParamConfig` and
:class:`unidesign.jobs.LigandParameterizationJob`. The script at
:file:`examples/python/ligand_parameters.py` turns the bundled
:file:`example/ProteinLigandInteraction/1r091_BS01_JEN/lig_charge.mol2` file
into parameter and topology artefacts, persisting them to the chosen output
location.

Temporary workspace behaviour
------------------------------

Each :meth:`UniDesignRunner.run` invocation creates a dedicated temporary
workspace. By default these directories are deleted as soon as the subprocess
completes, mirroring the behaviour of :class:`tempfile.TemporaryDirectory`.
Callers can set ``persist_workdir=True`` to retain the original directory for
inspection or additional post-processing. The generated result object includes
these fields:

``args``
    Argument vector used to run the binary (excluding the executable path).
``stdout`` / ``stderr``
    Captured process output streams.
``workdir``
    Path to the directory containing inputs and outputs for the execution.
``prefix``
    Unique prefix injected automatically into CLI arguments to namespace
    generated files.

Persistence options
-------------------

Job helpers always invoke the runner with ``persist_workdir=True`` and then
relocate artefacts into a caller-visible directory:

* When ``keep_workspace=False`` (default) a new directory under
  ``unidesign_artifacts_<uuid>`` is created, populated only with discovered
  artefacts, and returned via ``result.workspace``.
* Setting ``keep_workspace=True`` preserves the original workspace produced by
  the runner. This is useful when debugging or when additional intermediate
  files are required.

All result containers expose a ``close()`` method that removes the workspace and
invalidates further file access. Individual artefacts also provide a
:meth:`unidesign.artifacts.UniDesignArtifact.persist` helper to copy the file
into a user-controlled location for long-term storage.

Testing and CI
--------------

Pytest-based smoke tests under :file:`python/tests/` simulate each command and
verify that artefacts are captured correctly. The GitHub Actions workflow
:file:`.github/workflows/python-wrapper.yml` runs these tests and builds this
Sphinx documentation to guard against regressions in the Python wrapper.

Additional resources
--------------------

* :file:`examples/python/` contains runnable scripts demonstrating the wrapper
  APIs.
* :mod:`unidesign.config` documents the full set of command-line options
  exposed by the wrapper. Autodoc stubs can be added to this documentation set
  to include API reference material as the wrapper evolves.
