#!/usr/bin/env python3
import numpy as np
import nibabel as nib

# Define volume dimensions
dim = 64
mx = 1.0
fnm = "borg.nii.gz"

# Create grid coordinates scaled to [-π, π]
lin = np.linspace(-np.pi, np.pi, dim)
X, Y, Z = np.meshgrid(lin, lin, lin, indexing='ij')

# Evaluate the Borg function
val = np.sin(X * Y) + np.sin(Y * Z) + np.sin(Z * X)

# Clamp negative values to 0
val[val < 0] = 0

# Normalize to range 0..255
val -= val.min()
val /= val.max()
val *= mx

# Convert to float32 for saving
volume = val.astype(np.float32)

# Identity affine (1mm voxel spacing)
affine = np.eye(4)

# Save as NIfTI
nii = nib.Nifti1Image(volume, affine)
nib.save(nii, fnm)

print(f"Saved NIfTI image as {fnm}")
