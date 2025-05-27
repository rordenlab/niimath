## About

`niimath` is a high-performance command-line tool for manipulating NIfTI images. Because `niimath` is so fast, the main bottleneck in many workflows is disk I/O. To enable efficient workflows, niimath supports reading from and writing to standard input/output streams using the special filename `-`. It is important to note that both the input and output must be single-file NIfTI-1 images (e.g. `img.nii`). Compressed files (e.g. `img.nii.gz`) must be decompressed prior to using the memory stream input.

This streaming capability is especially useful for chaining tools together without writing intermediate results to disk. For example:

```bash
niimath - add 1 -     # read from stdin, add 1, write to stdout
```

Note that only one image can be piped at a time. For example:

```bash
niimath - add img.nii -   # input image is piped; second image is read from disk
```

## Demo Scripts

This folder contains minimal Python scripts demonstrating how to use these features. These assume niimath is available in your system path. Note that some scripts use nibabel while others demonstrate direct interaction if you do not wish to include nibabel as a dependency.

```bash
# Generate a synthetic NIfTI volume with a 3D pattern
python generate_borg.py

# Use niimath to process a NIfTI file and return the result via stdout (no disk writes)
python write_stdout.py

# Pipe a NIfTI image to niimath via stdin, output saved to disk
python read_stdin.py

# Pipe a NIfTI image to niimath and capture the output via stdout, fully in memory
python read_write_stream.py

# Use nibabel to_bytes() and from_bytes() to call niimath using pipes
python nibabel_niimath.py
```
