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
