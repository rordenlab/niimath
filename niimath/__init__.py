"""Thin wrapper around niimath binary"""
__author__ = "Casper da Costa-Luis <https://github.com/casperdcl>"
__date__ = "2025"

try:
    from ._dist_ver import version as __version__
except ImportError:  # pragma: nocover
    __version__ = "UNKNOWN"

__all__ = ['bin', 'bin_path', 'main']

from pathlib import Path

bin_path = Path(__file__).resolve().parent / "niimath"
bin = str(bin_path)


def main(args=None, **run_kwargs):
    """
    Arguments:
      args: defaults to `sys.argv[1:]`.
      **run_kwargs: passed to `subprocess.run`.
    Returns: `int` exit code from dcm2niix execution.
    """
    if args is None:
        import sys
        args = sys.argv[1:]
    from subprocess import run
    return run([bin] + args, **run_kwargs).returncode

