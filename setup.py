"""Compile source code and setup Python 3 package"""
import re
from pathlib import Path

from setuptools_scm import get_version
from skbuild import setup

# Get version from SCM (Git tags)
__version__ = get_version(root=".", relative_to=__file__)
build_ver = ".".join(__version__.split(".")[:3]).split(".dev")[0]

# Clean problematic pip env entries in CMakeCache
cache_dir = Path(__file__).resolve().parent / "_skbuild"
for cmake_file in cache_dir.rglob("CMakeCache.txt"):
    cmake_file.write_text(re.sub(
        r"^//.*$\n^[^#].*pip-build-env.*$", "",
        cmake_file.read_text(),
        flags=re.M
    ))

# Build the package with scikit-build and CMake
setup(
    name="niimath",
    use_scm_version=True,
    packages=["niimath"],
    cmake_languages=("C",),  # "C" since niimath is a C program
    cmake_minimum_required_version="3.5",
    # Optional: explicitly state where the source is
    # package_dir={"": "src"},
)
