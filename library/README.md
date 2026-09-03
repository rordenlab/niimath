# Streaming NIfTI images through niimath

## About

`niimath` is a high-performance command-line tool for NIfTI images. Because `niimath` is so fast, disk I/O is the main bottleneck in many workflows. To avoid it, niimath can read from standard input and write to standard output. Use the special filename `-` for either stream. Both the input and the output must be single-file NIfTI-1 images (for example `img.nii`). Decompress compressed files (for example `img.nii.gz`) before you pipe them in.

Streaming lets you chain tools without writing intermediate results to disk:

```bash
niimath - -add 1 -     # read from stdin, add 1, write to stdout
```

Only one image can be piped at a time. A second image is read from disk:

```bash
niimath - -add img.nii -   # the input image is piped; the second image is read from disk
```

## Demo scripts

This folder contains minimal Python scripts that show these features. They assume that `niimath` is on your system path. Some scripts use nibabel. Others show direct interaction, for when you do not want nibabel as a dependency.

```bash
# Generate a synthetic NIfTI volume with a 3D pattern
python generate_borg.py

# Process a NIfTI file and return the result via stdout (no disk writes)
python write_stdout.py

# Pipe a NIfTI image to niimath via stdin, and save the output to disk
python read_stdin.py

# Pipe a NIfTI image to niimath and capture the output via stdout, fully in memory
python read_write_stream.py

# Use nibabel to_bytes() and from_bytes() to call niimath through pipes
python nibabel_niimath.py
```
