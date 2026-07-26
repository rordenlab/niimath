#!/usr/bin/env python3
"""Turn a BIDS ``SliceTiming`` sidecar into slice times for niimath's ``-stc``.

    python3 stc_slicetiming.py bold.nii.gz              # comma list for niimath
    python3 stc_slicetiming.py bold.nii.gz -f afni      # whitespace list for 3dTshift @file
    python3 stc_slicetiming.py bold.nii.gz -o times.1D

    niimath bold.nii.gz -stc --slicetiming "$(python3 stc_slicetiming.py bold.nii.gz)" out.nii.gz

Standard library only: it reads the NIfTI-1 header (352 bytes) without touching image data, so it
runs anywhere Python does and adds no dependency to niimath.

What it checks, and why each check exists:

* ``SliceTiming`` must have exactly one value per slice along the NIfTI storage axis ``k``. A
  count mismatch almost always means the sidecar belongs to a different acquisition.
* ``SliceEncodingDirection`` is honoured: absent or ``k`` passes the array through unchanged,
  ``k-`` reverses it into storage slice-index order. ``i``, ``i-``, ``j`` and ``j-`` are rejected,
  because niimath's ``-stc`` v1 corrects along ``k`` only and silently emitting a ``k``-ordered
  list for ``j`` data would corrupt the result without any error.
* ``RepetitionTime`` (seconds, per BIDS) is compared against the NIfTI ``pixdim[4]`` normalised
  through ``xyzt_units``. A mismatch is an error rather than an implicit override: one of the two
  files is wrong, and guessing which produces plausible, silently misaligned output.

These checks deliberately DUPLICATE limits that ``src/stc.c`` also enforces -- the minimum of 5
volumes (``STC_MIN_NT``), the ``[0, TR]`` slice-time range, and the seconds/milliseconds/
microseconds ``xyzt_units`` table (``stc_time_scale``). The duplication is the point: the helper
exists to fail before you build a command line, with a message naming the sidecar. If you change
one of those rules in ``stc.c``, change it here too.

NIfTI-2 and paired ``.hdr``/``.img`` inputs are out of scope and are reported as such.
"""

import argparse
import gzip
import json
import os
import struct
import sys

HDR_BYTES = 348
DIM_OFF = 40          # short dim[8]
PIXDIM_OFF = 76       # float pixdim[8]
XYZT_UNITS_OFF = 123  # unsigned char
MAGIC_OFF = 344

UNITS_SEC = 8
UNITS_MSEC = 16
UNITS_USEC = 24
TIME_SCALE = {UNITS_SEC: 1.0, UNITS_MSEC: 1.0e-3, UNITS_USEC: 1.0e-6}


class Bad(Exception):
    pass


def read_header(path):
    """Return (nz, nt, tr_seconds) from a NIfTI-1 header, without reading image data."""
    opener = gzip.open if path.endswith(".gz") else open
    try:
        with opener(path, "rb") as f:
            raw = f.read(HDR_BYTES)
    except OSError as e:
        raise Bad("cannot read '%s': %s" % (path, e))
    if len(raw) < HDR_BYTES:
        raise Bad("'%s' is shorter than a NIfTI-1 header" % path)
    (native,) = struct.unpack("<i", raw[:4])
    (swapped,) = struct.unpack(">i", raw[:4])
    if native == 348:
        end = "<"
    elif swapped == 348:
        end = ">"
    elif 540 in (native, swapped):
        raise Bad("'%s' is NIfTI-2; this helper reads NIfTI-1 only" % path)
    else:
        raise Bad("'%s' is not a NIfTI-1 file (sizeof_hdr = %d)" % (path, native))
    magic = raw[MAGIC_OFF:MAGIC_OFF + 4]
    if magic not in (b"n+1\x00", b"ni1\x00"):
        raise Bad("'%s' has an unrecognised NIfTI magic %r" % (path, magic))
    if magic == b"ni1\x00":
        raise Bad("'%s' is a paired .hdr/.img dataset; this helper reads single-file NIfTI" % path)
    dim = struct.unpack(end + "8h", raw[DIM_OFF:DIM_OFF + 16])
    pixdim = struct.unpack(end + "8f", raw[PIXDIM_OFF:PIXDIM_OFF + 32])
    units = struct.unpack(end + "B", raw[XYZT_UNITS_OFF:XYZT_UNITS_OFF + 1])[0]
    if dim[0] < 3:
        raise Bad("'%s' has dim[0] = %d; a 4D series is required" % (path, dim[0]))
    nz = dim[3]
    nt = dim[4] if dim[0] >= 4 else 1
    tunit = units & 0x38
    scale = TIME_SCALE.get(tunit)
    if scale is None:
        raise Bad("'%s' has no usable temporal unit (xyzt_units time field = %d)" % (path, tunit))
    tr = pixdim[4] * scale
    if not (tr > 0.0) or tr != tr or tr in (float("inf"), float("-inf")):
        raise Bad("'%s' has a non-positive or non-finite TR (pixdim[4] = %r)" % (path, pixdim[4]))
    return nz, nt, tr


def sidecar_for(nifti_path):
    base = nifti_path
    for ext in (".nii.gz", ".nii"):
        if base.endswith(ext):
            return base[: -len(ext)] + ".json"
    return base + ".json"


def main(argv=None):
    ap = argparse.ArgumentParser(description="BIDS SliceTiming -> niimath -stc slice times")
    ap.add_argument("nifti", help="4D NIfTI-1 file (.nii or .nii.gz)")
    ap.add_argument("--json", help="BIDS sidecar (default: the .json beside the NIfTI)")
    ap.add_argument("-f", "--format", choices=("niimath", "afni"), default="niimath",
                    help="niimath: one comma-separated line; afni: whitespace, for 3dTshift @file")
    ap.add_argument("-o", "--out", help="write to this file instead of stdout")
    ap.add_argument("--tr-tol", type=float, default=1e-4,
                    help="allowed |JSON RepetitionTime - header TR| in seconds (default 1e-4)")
    args = ap.parse_args(argv)
    # A NaN or negative tolerance would make every comparison below vacuously pass.
    if not (args.tr_tol >= 0.0) or args.tr_tol == float("inf"):
        sys.stderr.write("stc_slicetiming: --tr-tol must be a finite, non-negative number\n")
        return 1

    try:
        nz, nt, tr = read_header(args.nifti)
        jpath = args.json or sidecar_for(args.nifti)
        try:
            with open(jpath, "r") as f:
                meta = json.load(f)
        except OSError as e:
            raise Bad("cannot read sidecar '%s': %s" % (jpath, e))
        except ValueError as e:
            raise Bad("'%s' is not valid JSON: %s" % (jpath, e))
        if not isinstance(meta, dict):
            raise Bad("'%s' does not hold a JSON object" % jpath)

        times = meta.get("SliceTiming")
        if times is None:
            raise Bad("'%s' has no SliceTiming" % jpath)
        if not isinstance(times, list) or not times:
            raise Bad("'%s' SliceTiming is not a non-empty array" % jpath)
        vals = []
        for i, t in enumerate(times):
            if isinstance(t, bool) or not isinstance(t, (int, float)):
                raise Bad("SliceTiming[%d] is not a number (%r)" % (i, t))
            t = float(t)
            if t != t or t in (float("inf"), float("-inf")):
                raise Bad("SliceTiming[%d] is not finite" % i)
            vals.append(t)
        if len(vals) != nz:
            raise Bad("SliceTiming has %d values but '%s' has %d slices"
                      % (len(vals), args.nifti, nz))

        direction = meta.get("SliceEncodingDirection")
        if direction is None or direction == "k":
            pass
        elif direction == "k-":
            vals.reverse()
        elif direction in ("i", "i-", "j", "j-"):
            raise Bad("SliceEncodingDirection '%s' is not supported: niimath -stc corrects "
                      "along the storage k axis only" % direction)
        else:
            raise Bad("unrecognised SliceEncodingDirection %r" % (direction,))

        rt = meta.get("RepetitionTime")
        if rt is not None:
            if isinstance(rt, bool) or not isinstance(rt, (int, float)):
                raise Bad("RepetitionTime is not a number (%r)" % (rt,))
            rt = float(rt)
            # Python's json accepts the non-standard literals NaN/Infinity, and every comparison
            # against NaN is False -- so an unchecked NaN would sail through the tolerance test
            # below and silently disable the very check it is meant to fail.
            if rt != rt or rt in (float("inf"), float("-inf")) or rt <= 0.0:
                raise Bad("RepetitionTime is not a finite positive number (%r)" % (rt,))
            if abs(rt - tr) > args.tr_tol:
                raise Bad("RepetitionTime %g s in '%s' disagrees with the header TR %g s in '%s'"
                          % (rt, jpath, tr, args.nifti))

        lo, hi = min(vals), max(vals)
        if lo < 0.0 or hi > tr:
            raise Bad("SliceTiming spans [%g, %g] s, outside [0, TR] with TR = %g s" % (lo, hi, tr))
        if nt < 5:
            raise Bad("'%s' has nt = %d; -stc needs at least 5 volumes" % (args.nifti, nt))
    except Bad as e:
        sys.stderr.write("stc_slicetiming: %s\n" % e)
        return 1

    sep = "," if args.format == "niimath" else " "
    text = sep.join(repr(v) for v in vals) + "\n"
    if args.out:
        try:
            with open(args.out, "w") as f:
                f.write(text)
        except OSError as e:
            sys.stderr.write("stc_slicetiming: cannot write '%s': %s\n" % (args.out, e))
            return 1
    else:
        sys.stdout.write(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
