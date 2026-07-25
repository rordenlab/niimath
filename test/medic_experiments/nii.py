# Tiny NIfTI-1 read/write helpers for the MEDIC black-box experiments.
# Analysis-only: NOT shipped, NOT run in CI. Uses numpy (warpkit venv python).
# ponytail: no nibabel dependency -- one-file NIfTI-1 is 348 bytes of struct.
import gzip
import struct

import numpy as np

DT = {2: np.uint8, 4: np.int16, 8: np.int32, 16: np.float32, 64: np.float64, 512: np.uint16, 256: np.int8}
DT_INV = {v: k for k, v in DT.items()}


def _open(path, mode):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def read(path, scaled=True):
    """Return (data[x,y,z,t] float64 if scaled else raw, header dict)."""
    with _open(path, "rb") as f:
        raw = f.read()
    little = struct.unpack("<i", raw[:4])[0] == 348
    e = "<" if little else ">"
    h = {}
    h["dim"] = struct.unpack(e + "8h", raw[40:56])
    h["intent_code"] = struct.unpack(e + "h", raw[68:70])[0]
    h["datatype"] = struct.unpack(e + "h", raw[70:72])[0]
    h["bitpix"] = struct.unpack(e + "h", raw[72:74])[0]
    h["pixdim"] = struct.unpack(e + "8f", raw[76:108])
    h["vox_offset"] = struct.unpack(e + "f", raw[108:112])[0]
    h["scl_slope"], h["scl_inter"] = struct.unpack(e + "2f", raw[112:120])
    h["cal_max"], h["cal_min"] = struct.unpack(e + "2f", raw[124:132])
    h["xyzt_units"] = raw[123]
    h["descrip"] = raw[148:228].split(b"\0")[0].decode("ascii", "replace")
    h["qform_code"], h["sform_code"] = struct.unpack(e + "2h", raw[252:256])
    h["quatern"] = struct.unpack(e + "6f", raw[256:280])
    h["srow"] = np.array(struct.unpack(e + "12f", raw[280:328])).reshape(3, 4)
    n = h["dim"][0]
    shape = tuple(h["dim"][1 : n + 1])
    off = int(h["vox_offset"])
    a = np.frombuffer(raw, dtype=np.dtype(DT[h["datatype"]]).newbyteorder(e), count=int(np.prod(shape)), offset=off)
    a = a.reshape(shape, order="F")
    if scaled and h["scl_slope"] not in (0.0, 1.0) or (scaled and h["scl_inter"] != 0.0):
        a = a.astype(np.float64) * (h["scl_slope"] or 1.0) + h["scl_inter"]
    return np.asarray(a), h


def write(path, data, ref=None, dtype=np.float32, scl=(1.0, 0.0), descrip=b"medic-experiment"):
    """Write a one-file NIfTI-1. `ref` supplies pixdim/srow/qform when given."""
    data = np.asarray(data, dtype=dtype, order="F")
    shape = data.shape
    dim = [len(shape)] + list(shape) + [1] * (7 - len(shape))
    hdr = bytearray(352)
    e = "<"
    struct.pack_into(e + "i", hdr, 0, 348)
    struct.pack_into(e + "8h", hdr, 40, *dim)
    struct.pack_into(e + "h", hdr, 70, DT_INV[np.dtype(dtype).type])
    struct.pack_into(e + "h", hdr, 72, np.dtype(dtype).itemsize * 8)
    pix = list(ref["pixdim"]) if ref else [1.0] * 8
    struct.pack_into(e + "8f", hdr, 76, *pix)
    struct.pack_into(e + "f", hdr, 108, 352.0)
    struct.pack_into(e + "2f", hdr, 112, *scl)
    hdr[123] = ref["xyzt_units"] if ref else 10
    hdr[148 : 148 + len(descrip)] = descrip
    hdr[344:348] = b"n+1\0"
    if ref is not None:
        struct.pack_into(e + "2h", hdr, 252, ref["qform_code"], ref["sform_code"])
        struct.pack_into(e + "6f", hdr, 256, *ref["quatern"])
        struct.pack_into(e + "12f", hdr, 280, *ref["srow"].ravel())
    else:
        struct.pack_into(e + "2h", hdr, 252, 0, 1)
        srow = np.array([[pix[1], 0, 0, 0], [0, pix[2], 0, 0], [0, 0, pix[3], 0]], dtype=np.float32)
        struct.pack_into(e + "12f", hdr, 280, *srow.ravel())
    with _open(path, "wb") as f:
        f.write(bytes(hdr))
        f.write(data.tobytes(order="F"))
    return path


def stats(a, name=""):
    a = np.asarray(a, dtype=np.float64)
    fin = np.isfinite(a)
    return (
        f"{name} shape={a.shape} finite={fin.sum()}/{a.size} "
        f"min={np.nanmin(a[fin]) if fin.any() else float('nan'):.6g} "
        f"max={np.nanmax(a[fin]) if fin.any() else float('nan'):.6g} "
        f"mean={a[fin].mean() if fin.any() else float('nan'):.6g}"
    )
