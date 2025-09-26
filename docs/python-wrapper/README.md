# UniDesign Python Wrapper Guide

The Python helper package under `python/unidesign` wraps the UniDesign CLI so
that pipelines can assemble strongly-typed configurations, execute the binary in
isolated workspaces, and capture the generated artifacts in a structured way.
This guide walks through environment setup, core concepts, and persistence
options for working with the wrapper.

## Environment Setup

1. **Build or obtain the UniDesign binary** – compile the project with
   `./build.sh` or download the prebuilt executable. Place the binary in the
   repository root so that :func:`unidesign.paths.discover_binary` can find it.
2. **Install Python dependencies** – the wrapper only requires the standard
   library at runtime, but development and testing use `pytest`.
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   pip install -r requirements-test.txt  # optional convenience file
   ```
3. **Add the repository to `PYTHONPATH`** – run scripts from the repository
   root or install the package in editable mode: `pip install -e python`.

## Core Abstractions

The wrapper is organized into a few layers:

- **Configuration dataclasses** (:mod:`unidesign.config`) map Python-friendly
  fields to CLI arguments. For example,
  :class:`unidesign.config.ProteinDesignConfig` renders all flags required to
  run the `ProteinDesign` command.
- **Runner** (:class:`unidesign.runner.UniDesignRunner`) locates the binary,
  prepares a private workspace populated with `library/`, `wread/`, and
  `extbin/`, injects a unique `--prefix`, and captures stdout/stderr.
- **Jobs** (modules under :mod:`unidesign.jobs`) tie the configuration and the
  runner together, returning typed result objects that expose generated files as
  :class:`unidesign.jobs.JobArtifact` instances.

Using these building blocks keeps orchestration logic declarative while
preserving full access to the underlying CLI.

## Executing Commands

The typical flow for running a command is:

```python
from pathlib import Path

from unidesign import (
    ProteinDesignConfig,
    ProteinDesignJob,
    UniDesignRunner,
    discover_binary,
)

binary = discover_binary()
runner = UniDesignRunner(binary)
config = ProteinDesignConfig(
    pdb_path=Path("example/MonomerDesign/1agy/1agy.pdb"),
    design_chains="A",
)

job = ProteinDesignJob(runner, config)
result = job.run()
print(result.run.stdout)
print(result.best_sequences.read_text())
```

Every job accepts an optional ``keep_workspace`` flag to retain the original
scratch directory, which is useful during debugging.

## Temporary Workspace Handling

Each call to :meth:`UniDesignRunner.run <unidesign.runner.UniDesignRunner.run>`
uses ``tempfile.TemporaryDirectory`` to provide an isolated workspace. Static
resources are symlinked or copied into the directory so that concurrent runs do
not interfere with one another. Jobs call
:func:`unidesign.jobs.relocate_artifacts` to either keep the workspace intact or
copy only the generated outputs into a second managed directory. The
:class:`JobArtifact` objects returned by result bundles always reference the
retained directory, making it safe to read even after the original workspace was
cleaned up.

## Persisting Outputs

Artifacts expose helpers for common persistence workflows:

```python
design_result = job.run()
# Copy the best-scoring structure into a reports directory with a custom name.
if design_result.best_structure:
    target = design_result.best_structure.persist(
        "reports/design-run-001", filename="best_structure.pdb"
    )
    print("Wrote", target)
```

Passing ``keep_workspace=True`` to ``job.run`` is another option when the entire
scratch tree—including intermediate files—is required downstream. Remember to
call ``result.close()`` when you are done so the temporary directories can be
removed.

## Parallel Execution

Because every job invocation receives its own temporary directory and prefix,
you can safely run multiple tasks in parallel using ``concurrent.futures``.
For CPU-bound workloads, prefer :class:`concurrent.futures.ProcessPoolExecutor`;
for I/O-bound orchestration that mostly shells out to UniDesign, threads are
usually sufficient. Each job's result object includes ``stdout``, ``stderr``,
and the return code, enabling callers to detect and handle failures gracefully.

## Next Steps

- Review the example scripts in `examples/python/` for end-to-end workflows.
- Inspect the automated tests under `tests/python_wrapper/` to see how the
  runner and job classes can be mocked for integration testing.
- Extend the package with additional job classes as new UniDesign CLI commands
  are needed.
