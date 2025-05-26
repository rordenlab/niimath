#!/usr/bin/env python3
import subprocess
import numpy as np
import struct
import io

input_file = "borg.nii.gz"

# Step 1: Run niimath and get output to stdout
try:
    result = subprocess.run(
        ['niimath', input_file, '-add', '1', '-.nii'],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True
    )
except subprocess.CalledProcessError as e:
    print("niimath failed:")
    print(e.stderr.decode())
    exit(1)

# Step 2: Load binary output and parse NIfTI header
nii_bytes = result.stdout
buffer = io.BytesIO(nii_bytes)

# Validate NIfTI header size
sizeof_hdr = struct.unpack('<I', buffer.read(4))[0]
if sizeof_hdr != 348:
    raise ValueError(f"Invalid NIfTI header size: {sizeof_hdr} (expected 348)")

# Seek to dim[1..3] at offset 42
buffer.seek(40)  # dim array starts at byte 40
dim = struct.unpack('<8H', buffer.read(16))  # dim[0..7]
nx, ny, nz = dim[1], dim[2], dim[3]
print(f"Image dimensions: {nx} x {ny} x {nz}")

# Read datatype (UINT16) at offset 70
buffer.seek(70)
datatype = struct.unpack('<H', buffer.read(2))[0]
if datatype != 16:
    raise ValueError(f"Unsupported datatype {datatype} (expected 16 for float32)")

# Read vox_offset (FLOAT32) at offset 108
buffer.seek(108)
vox_offset = struct.unpack('<f', buffer.read(4))[0]
vox_offset = int(vox_offset)
print(f"vox_offset: {vox_offset}")

# Step 3: Read and reshape image data
buffer.seek(vox_offset)
nvox = nx * ny * nz
data = np.frombuffer(buffer.read(nvox * 4), dtype=np.float32).reshape((nx, ny, nz))

# Step 4: Report max value
print(f"Max voxel value: {np.max(data):.3f}")
