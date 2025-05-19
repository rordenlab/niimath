"""Compile source code and setup Python 3 package"""
import re
from pathlib import Path

from setuptools_scm import get_version
from skbuild import setup

# Get version directly from most recent git tag
import subprocess
try:
    # Get the latest tag without commits since then
    tag_version = subprocess.check_output(['git', 'describe', '--tags', '--abbrev=0'], 
                                         universal_newlines=True).strip()
    # Remove leading 'v' if present
    if tag_version.startswith('v'):
        tag_version = tag_version[1:]
except (subprocess.SubprocessError, OSError):
    # Fallback version if git command fails
    tag_version = "1.0.20250518"

# Clean problematic pip env entries in CMakeCache
for i in (Path(__file__).resolve().parent / "_skbuild").rglob("CMakeCache.txt"):
    i.write_text(re.sub("^//.*$\n^[^#].*pip-build-env.*$", "", i.read_text(), flags=re.M))

# Build the package with scikit-build and CMake
setup(
    name="niimath",
    version=tag_version,
    packages=["niimath"],
    cmake_languages=("C",),
    cmake_minimum_required_version="3.5",
    entry_points={
        "console_scripts": [
            "niimath=niimath:main",
        ],
    }
)
