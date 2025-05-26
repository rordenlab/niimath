#!/usr/bin/env python3
import subprocess
import gzip
import struct
import io
import numpy as np

input_file = "borg.nii.gz"

def is_gzipped(path):
    with open(path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

# Step 1: Read and decompress image if gzipped
if is_gzipped(input_file):
    with gzip.open(input_file, 'rb') as f:
        image_data = f.read()
else:
    with open(input_file, 'rb') as f:
        image_data = f.read()

# Step 2: Run niimath with piped input and receive piped output
try:
    result = subprocess.run(
        ['niimath', '-', '-add', '1', '-.nii'],
        input=image_data,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True
    )
except subprocess.CalledProcessError as e:
    print("niimath failed:")
    print(e.stderr.decode())
    exit(1)

# Step 3: Parse the NIfTI output stream
nii_bytes = result.stdout
buffer = io.BytesIO(nii_bytes)

# Read and validate NIfTI header
sizeof_hdr = struct.unpack('<I', buffer.read(4))[0]
if sizeof_hdr != 348:
    raise ValueError(f"Invalid NIfTI header size: {sizeof_hdr} (expected 348)")

# Read dimensions from offset 40 (dim[1..3])
buffer.seek(40)
dim = struct.unpack('<8H', buffer.read(16))  # dim[0..7]
nx, ny, nz = dim[1], dim[2], dim[3]
print(f"Image dimensions: {nx} x {ny} x {nz}")

# Read datatype from offset 70
buffer.seek(70)
datatype = struct.unpack('<H', buffer.read(2))[0]
if datatype != 16:  # float32
    raise ValueError(f"Unsupported datatype {datatype} (expected 16 for float32)")

# Read vox_offset (FLOAT32) from offset 108
buffer.seek(108)
vox_offset = struct.unpack('<f', buffer.read(4))[0]
vox_offset = int(vox_offset)
print(f"vox_offset: {vox_offset}")

# Step 4: Load and reshape image data
buffer.seek(vox_offset)
nvox = nx * ny * nz
data = np.frombuffer(buffer.read(nvox * 4), dtype=np.float32).reshape((nx, ny, nz))

# Step 5: Display max value in the output volume
print(f"Max voxel value: {np.max(data):.3f}")
