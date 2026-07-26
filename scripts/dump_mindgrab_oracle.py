#!/usr/bin/env python3
"""Dump per-layer MindGrab intermediates from the tinygrad reference, for test/mindgrab_parity.c.

    pip install brainchop            # brings tinygrad; needs `niimath` on PATH for -conform
    MG_BACKEND=METAL python3 scripts/dump_mindgrab_oracle.py t1.nii.gz oracle_metal
    MG_BACKEND=CPU   MG_FULL=-1 python3 scripts/dump_mindgrab_oracle.py t1.nii.gz oracle_cpu

then build and run the comparison. Its compile line lives in test/mindgrab_parity.c's header --
deliberately in ONE place: a second copy here drifted from it (wrong -I, missing -lomp) and was
removed in audit.

MG_BACKEND selects the tinygrad device (default CPU). MG_FULL lists which layers get a full
1 GB dump (default 0,1,5,24; pass -1 for none) -- per-layer statistics are always written.
The full set is about 4 GB and mindgrab_parity.c EXITS NONZERO without it.
Comparing a METAL run against a CPU run is worth doing first: it measures how much numerical
slack the argmax actually has (measured: logits differ by up to 3.1e-5, argmax by 0 voxels).

Writes raw little-endian float32 (or uint8) blobs in the RAW NIfTI axis order
(z slowest, x fastest) so the C port can memcmp/compare without transposing:

  conformed.u8      256^3 uint8   (what the model actually sees, pre-normalise)
  normalized.f32    256^3 float32 (qnormalize output)
  layerNN.f32       256^3*15 float32, channel-LAST (voxel-major), post-GELU
  logits.f32        256^3*2  float32, channel-LAST
  argmax.u8         256^3 uint8
"""
import os, sys, json, subprocess, struct
import numpy as np

# tinygrad picks its device from an env var named after the device, e.g. METAL=1 / CPU=1.
backend = os.environ.pop("MG_BACKEND", "CPU")
os.environ.setdefault(backend, "1")

from tinygrad.tensor import Tensor
from tinygrad import nn
from brainchop.tiny_meshnet import load_meshnet, qnormalize

MODEL_DIR = os.path.expanduser("~/.cache/brainchop/models/mindgrab")


def conform(path):
    cmd = ["niimath", path, "-conform", "-gz", "0", "-", "-odt", "char"]
    out = subprocess.run(cmd, capture_output=True, check=True).stdout
    vox_offset = int(struct.unpack("<f", out[108:112])[0])
    hdr, data = out[:vox_offset], out[vox_offset:]
    return np.frombuffer(data, dtype=np.uint8).reshape((256, 256, 256)), hdr


def main():
    t1, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    volume, hdr = conform(t1)
    volume.tofile(os.path.join(outdir, "conformed.u8"))
    # the 352-byte NIfTI header of the conformed volume; mindgrab_parity.c does not read it,
    # but it is what lets you wrap any of these raw blobs back into a viewable .nii by hand
    open(os.path.join(outdir, "conformed.hdr"), "wb").write(hdr)

    model = load_meshnet(os.path.join(MODEL_DIR, "model.json"),
                         os.path.join(MODEL_DIR, "model.pth"))

    # exactly what main.py does
    image = Tensor(volume.transpose((2, 1, 0)).astype(np.float32)).rearrange("... -> 1 1 ...")
    image = model.normalize(image)
    # normalized, back to raw (z,y,x) order
    np.ascontiguousarray(image.numpy()[0, 0].transpose((2, 1, 0))).tofile(
        os.path.join(outdir, "normalized.f32"))

    # full volumes are 1 GB each; keep only the ones the plan pins, stats for all
    FULL = {int(v) for v in os.environ.get("MG_FULL", "0,1,5,24").split(",")}
    stats = {}
    x = image
    layer_idx = 0
    for i, layer in enumerate(model.model):
        x = layer(x)
        # the 3-op groups are (Conv, GroupNorm, gelu); dump after each gelu
        is_gelu = not isinstance(layer, (nn.Conv2d, nn.GroupNorm))
        if is_gelu:
            arr = x.numpy()[0]  # (C, X, Y, Z)
            # -> raw (z,y,x) with channel last
            arr = np.ascontiguousarray(arr.transpose((3, 2, 1, 0)))
            if layer_idx in FULL:
                arr.tofile(os.path.join(outdir, "layer%02d.f32" % layer_idx))
            f = arr.reshape(-1, arr.shape[-1]).astype(np.float64)
            stats["layer%02d" % layer_idx] = {
                "min": f.min(0).tolist(), "max": f.max(0).tolist(),
                "mean": f.mean(0).tolist(), "sumsq": (f * f).sum(0).tolist()}
            del f
            layer_idx += 1
    logits = x.numpy()[0]  # (2, X, Y, Z)
    logits = np.ascontiguousarray(logits.transpose((3, 2, 1, 0)))
    logits.tofile(os.path.join(outdir, "logits.f32"))
    am = np.argmax(logits, axis=-1).astype(np.uint8)
    am.tofile(os.path.join(outdir, "argmax.u8"))
    stats["ties"] = int((logits[..., 0] == logits[..., 1]).sum())
    with open(os.path.join(outdir, "stats.json"), "w") as fh:
        json.dump(stats, fh)
    print("layers dumped:", layer_idx, "argmax fg voxels:", int(am.sum()),
          "exact logit ties:", stats["ties"])


if __name__ == "__main__":
    main()
