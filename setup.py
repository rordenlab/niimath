"""Compile source code and setup Python 3 package"""
import re
from pathlib import Path

from setuptools_scm import get_version
from skbuild import setup

# Get version from SCM (Git tags)
__version__ = get_version(root=".", relative_to=__file__)
build_ver = ".".join(__version__.split(".")[:3]).split(".dev")[0]

# Clean problematic pip env entries in CMakeCache
for i in (Path(__file__).resolve().parent / "_skbuild").rglob("CMakeCache.txt"):
    i.write_text(re.sub("^//.*$\n^[^#].*pip-build-env.*$", "", i.read_text(), flags=re.M))

# Build the package with scikit-build and CMake
setup(
    name="niimath",
    use_scm_version=True,
    packages=["niimath"],
    cmake_languages=("C",),
    cmake_minimum_required_version="3.5",
    entry_points={
        "console_scripts": [
            "niimath=niimath:main",
        ],
    }
)
