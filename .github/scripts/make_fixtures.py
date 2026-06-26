#!/usr/bin/env python3
"""Generate tiny synthetic NIfTI fixtures for the GPL coregistration CI smoke test.

Writes /tmp/{src,ref,mask}.nii.gz: small 3D blobs with a usable sform/qform so
spm_coreg can run (it only needs to execute, not converge well). `ref` is a
slightly shifted copy of `src`; `mask` is a binary sphere on the same grid.
"""
import numpy as np
import nibabel as nib


def blob(shape=(24, 24, 24), center=(12, 12, 12), r=6.0):
    z, y, x = np.ogrid[: shape[0], : shape[1], : shape[2]]
    d = np.sqrt((x - center[2]) ** 2 + (y - center[1]) ** 2 + (z - center[0]) ** 2)
    return np.clip(1.0 - d / r, 0.0, 1.0).astype(np.float32)


aff = np.diag([2.0, 2.0, 2.0, 1.0])
aff[:3, 3] = [-24.0, -24.0, -24.0]

vols = {
    "/tmp/src.nii.gz": blob(center=(12, 12, 12)),
    "/tmp/ref.nii.gz": blob(center=(13, 12, 11)),                  # slightly offset
    "/tmp/mask.nii.gz": (blob(center=(12, 12, 12), r=8.0) > 0.1).astype(np.float32),
}

for path, data in vols.items():
    img = nib.Nifti1Image(data, aff)
    img.header.set_sform(aff, code=1)   # spm_vox_mat needs a usable sform/qform
    img.header.set_qform(aff, code=1)
    nib.save(img, path)
    print("wrote", path, data.shape)

# 4D fixture for the dimension-rejection regression: registration/defacing are
# 3D-only, and a 4D input must fail closed (not silently use volume 0 or write a
# corrupt header). See the "4D inputs are rejected" CI step.
src3d = vols["/tmp/src.nii.gz"]
src4d = np.stack([src3d, src3d], axis=-1)            # (24, 24, 24, 2)
img4d = nib.Nifti1Image(src4d, aff)
img4d.header.set_sform(aff, code=1)
img4d.header.set_qform(aff, code=1)
nib.save(img4d, "/tmp/src4d.nii.gz")
print("wrote /tmp/src4d.nii.gz", src4d.shape)

# Larger fixture pair for the OpenMP thread-parity step. The parallel histogram
# branch in coreg_hist2_cached only triggers when base samples exceed 100000;
# 64^3 = 262144 voxels clears that at the 2mm sampling step (the 24^3 blobs do
# not, so they would only ever exercise the serial path). The phantom is
# ANISOTROPIC + multi-lobe (not a symmetric sphere): a symmetric blob hides
# rotational/trajectory differences, so it would not expose a thread-dependent
# divergence in the optimizer. bigref is a translated copy (a real registration
# target); the asymmetry makes the cost landscape rotation-sensitive.
def phantom(shape=(64, 64, 64), shift=(0, 0, 0)):
    z, y, x = np.ogrid[: shape[0], : shape[1], : shape[2]]
    sz, sy, sx = shift
    # anisotropic main ellipsoid (distinct per-axis extents -> rotation-sensitive)
    e = ((x - (30 + sx)) / 22.0) ** 2 + ((y - (32 + sy)) / 16.0) ** 2 + ((z - (34 + sz)) / 12.0) ** 2
    v = np.clip(1.0 - e, 0.0, 1.0)
    # two off-center lobes break the remaining symmetry
    for cx, cy, cz, r in [(44, 40, 28, 7.0), (20, 26, 42, 5.0)]:
        d = np.sqrt((x - (cx + sx)) ** 2 + (y - (cy + sy)) ** 2 + (z - (cz + sz)) ** 2)
        v = np.maximum(v, np.clip(1.0 - d / r, 0.0, 1.0) * 0.8)
    return v.astype(np.float32)

affbig = np.diag([2.0, 2.0, 2.0, 1.0])
affbig[:3, 3] = [-64.0, -64.0, -64.0]
big = {
    "/tmp/big.nii.gz": phantom(shift=(0, 0, 0)),
    "/tmp/bigref.nii.gz": phantom(shift=(1, 0, -1)),   # translated target
}
for path, data in big.items():
    img = nib.Nifti1Image(data, affbig)
    img.header.set_sform(affbig, code=1)
    img.header.set_qform(affbig, code=1)
    nib.save(img, path)
    print("wrote", path, data.shape)
