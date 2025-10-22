from __future__ import annotations

from pathlib import Path
import sys

from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext


def build_extensions() -> list[Pybind11Extension]:
    src_dir = Path("src")
    binding_dir = Path("python") / "unidesign" / "_bindings"
    sources = [str(path) for path in src_dir.glob("*.cpp")]
    sources += [str(path) for path in binding_dir.glob("*.cpp")]
    include_dirs = [str(src_dir)]

    extra_compile_args = []
    if sys.platform != "win32":
        extra_compile_args.extend(["-O3", "-ffast-math"])

    return [
        Pybind11Extension(
            "unidesign._core",
            sources=sources,
            include_dirs=include_dirs,
            cxx_std=17,
            extra_compile_args=extra_compile_args,
        )
    ]


setup(
    ext_modules=build_extensions(),
    cmdclass={"build_ext": build_ext},
    package_data={"unidesign": ["data/**/*", "py.typed"]},
)
