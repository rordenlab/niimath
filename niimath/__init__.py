"""Thin wrapper around niimath binary"""
__author__ = "Casper da Costa-Luis <https://github.com/casperdcl>"
__date__ = "2025"
# version detector. Precedence: installed dist, git, 'UNKNOWN'
try:
    from ._dist_ver import __version__
except ImportError: # pragma: nocover
    try:
        from setuptools_scm import get_version

        __version__ = get_version(root="../..", relative_to=__file__)
    except (ImportError, LookupError):
        __version__ = "UNKNOWN"
__all__ = ['bin', 'bin_path', 'main']

import os
import platform
from pathlib import Path

# Handle platform-specific binary naming
system = platform.system()
if system == "Windows":
    bin_name = "niimath.exe"
else:
    bin_name = "niimath"

# Check for binary in bin directory (standard installation)
bin_path = Path(__file__).resolve().parent / "bin" / bin_name

# Fallback to direct directory for backward compatibility
if not bin_path.exists():
    bin_path = Path(__file__).resolve().parent / bin_name

bin = str(bin_path)

# Make the binary executable on Unix systems
if system != "Windows" and bin_path.exists():
    os.chmod(bin, os.stat(bin).st_mode | 0o111)


def main(args=None, **run_kwargs):
    """
    Arguments:
      args: defaults to `sys.argv[1:]`.
      **run_kwargs: passed to `subprocess.run`.
    Returns: `int` exit code from niimath execution.
    """
    if args is None:
        import sys
        args = sys.argv[1:]
    from subprocess import run
    return run([bin] + args, **run_kwargs).returncode
