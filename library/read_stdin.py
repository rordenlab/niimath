#!/usr/bin/env python3
import subprocess
import gzip

input_file = "borg.nii.gz"

def is_gzipped(path):
    with open(path, 'rb') as f:
        magic = f.read(2)
    return magic == b'\x1f\x8b'

# Load and decompress if needed
if is_gzipped(input_file):
    with gzip.open(input_file, 'rb') as f:
        image_data = f.read()
else:
    with open(input_file, 'rb') as f:
        image_data = f.read()

# Run niimath with piped input
try:
    result = subprocess.run(
        ['niimath', '-', '-add', '1', 'borgPlus1'],
        input=image_data,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True
    )
    print("niimath completed successfully.")
except subprocess.CalledProcessError as e:
    print("niimath failed:")
    print(e.stderr.decode())
