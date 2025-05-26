#!/usr/bin/env python3
import numpy as np
import nibabel as nib
import subprocess

def make_binary_sphere(dim=63, radius=20):
    """Create a binary sphere in a NIfTI image."""
    X, Y, Z = np.meshgrid(
        np.arange(dim), np.arange(dim), np.arange(dim), indexing='ij'
    )
    cx = cy = cz = dim // 2
    dist = np.sqrt((X - cx) ** 2 + (Y - cy) ** 2 + (Z - cz) ** 2)
    data = (dist <= radius).astype(np.float32)
    affine = np.eye(4)
    return nib.Nifti1Image(data, affine)

def run_niimath_stream(nifti_img, niimath_args):
    """Send NIfTI image to niimath via stdin and return output as nibabel image."""
    nii_bytes = nifti_img.to_bytes()
    cmd = ['niimath', '-', *niimath_args, '-']

    try:
        result = subprocess.run(
            cmd,
            input=nii_bytes,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True
        )
    except subprocess.CalledProcessError as e:
        print("niimath failed:")
        print(e.stderr.decode())
        raise

    return nib.Nifti1Image.from_bytes(result.stdout)

if __name__ == "__main__":
    sphere_img = make_binary_sphere(dim=63, radius=20)
    blurred_img = run_niimath_stream(sphere_img, ['-s', '3'])

    nib.save(blurred_img, 'blurred.nii.gz')
    print("Saved blurred image to blurred.nii.gz")
