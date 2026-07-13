#!/usr/bin/env python3
"""Huge-image (> INT_MAX voxel) smoke test for niimath (issue #67).

Opt-in / manual — NOT part of `make test` or CI (like the WASI suite): a real run needs
~9 GB RAM to hold one INT_MAX-voxel float32 image. Input files are created SPARSE (a valid
352-byte NIfTI-1 header + a hole extended to the full logical data size), so they cost a few
KB on disk but read back as an all-zero > INT_MAX buffer.

Cases (see the "Huge image support" section of AGENTS.md):
  A. huge 4D `-p 1 -add 1 -ing 100 -Tmean -gz 0` -> every output voxel 100.0 (exercises
     harmless execution/output modifiers, the > INT_MAX rescale and global-normalization
     counts, and temporal-collapse source addressing at inPos > INT_MAX).
  B. huge 4D `-Tmean` on the zero image -> every output voxel 0.0.
  C. huge 4D `-mul 3 -sub 3 -Tmean` -> -3.0 (index-wrap sanity across chained > INT_MAX ops).
  D. out-of-scope op (`-fmean`) on a huge image -> rejected, exit 2 (fail-closed whitelist).
  E. double-dash out-of-scope op (`--compare`) -> rejected too (exact matching, not a
     leading-dash heuristic, classifies long options).
  F. a skewed huge image (one 1, >INT_MAX zeros) through `-thrp 50 -Tsum` -> one 1 and
     otherwise zero (exercises a robust-range histogram bin count above INT_MAX).
  G. a disguised token (`x-otsu`) -> rejected (exact admission and dispatch grammar).
  H. `-Tstd` on a huge image -> rejected (not promoted without a dedicated huge regression).
  I. `-restart` to a huge image, then `-otsu` -> rejected by replacement-header preflight
     before the huge payload is loaded.
  J. `small3d -add huge4d`, then `-otsu` -> rejected (re-admission after a binary op adopts a
     huge operand). D/E/G/H/I/J together cover both bypass classes and the parser grammar.

Usage: python3 test/huge_smoke.py /path/to/niimath  [--keep]
"""
import os, sys, struct, subprocess, tempfile, shutil

FLOAT32 = 16

def write_sparse_nifti(path, nx, ny, nz, nt, first_value=0.0):
    """Write a valid NIfTI-1 (.nii) header, then extend the file to the full logical size."""
    ndim = 4 if nt > 1 else 3
    hdr = bytearray(352)
    struct.pack_into('<i', hdr, 0, 348)                 # sizeof_hdr
    struct.pack_into('<h', hdr, 40, ndim)               # dim[0]
    struct.pack_into('<h', hdr, 42, nx)                 # dim[1]
    struct.pack_into('<h', hdr, 44, ny)                 # dim[2]
    struct.pack_into('<h', hdr, 46, nz)                 # dim[3]
    struct.pack_into('<h', hdr, 48, nt)                 # dim[4]
    struct.pack_into('<h', hdr, 50, 1)                  # dim[5]
    struct.pack_into('<h', hdr, 52, 1)                  # dim[6]
    struct.pack_into('<h', hdr, 54, 1)                  # dim[7]
    struct.pack_into('<h', hdr, 70, FLOAT32)            # datatype
    struct.pack_into('<h', hdr, 72, 32)                 # bitpix
    # pixdim[0..4] = 1.0
    for i, off in enumerate((76, 80, 84, 88, 92)):
        struct.pack_into('<f', hdr, off, 1.0)
    struct.pack_into('<f', hdr, 108, 352.0)             # vox_offset
    struct.pack_into('<f', hdr, 112, 1.0)               # scl_slope
    hdr[344:348] = b'n+1\x00'                           # magic
    nvox = nx * ny * nz * nt
    total = 352 + nvox * 4
    with open(path, 'wb') as f:
        f.write(hdr)
        f.truncate(total)                               # sparse hole -> reads as zeros
        if first_value != 0.0:
            f.seek(352)
            f.write(struct.pack('<f', first_value))
    return nvox

def read_nifti_data(path):
    with open(path, 'rb') as f:
        h = f.read(352)
        dim = struct.unpack_from('<8h', h, 40)
        dt = struct.unpack_from('<h', h, 70)[0]
        vox_offset = int(struct.unpack_from('<f', h, 108)[0])
        nvox = 1
        for d in range(1, dim[0] + 1):
            nvox *= dim[d]
        f.seek(vox_offset)
        raw = f.read(nvox * 4)
    assert dt == FLOAT32, f"unexpected output datatype {dt}"
    return struct.unpack('<%df' % nvox, raw), dim

def run(niimath, args, workdir):
    env = dict(os.environ, FSLOUTPUTTYPE='NIFTI')        # force plain .nii (no .gz)
    p = subprocess.run([niimath] + args, cwd=workdir, capture_output=True, text=True, env=env)
    return p.returncode, (p.stdout + p.stderr)

def main():
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(2)
    niimath = os.path.abspath(sys.argv[1])
    keep = '--keep' in sys.argv[2:]
    d = tempfile.mkdtemp(prefix='niimath_huge_')
    fails = []
    try:
        # nifti-1 dim[] is int16 (<= 32767). nvox3D = 256*256*2 = 131072; nt = 16385 ->
        # total = 2,147,614,720 = INT_MAX + 131073 (just over, and a genuine 4D collapse).
        NX, NY, NZ, NT = 256, 256, 2, 16385
        INT_MAX = 2147483647
        huge = os.path.join(d, 'huge4d.nii')
        nvox = write_sparse_nifti(huge, NX, NY, NZ, NT)
        assert nvox > INT_MAX, f"test image nvox {nvox} must exceed INT_MAX"
        phys = os.stat(huge).st_blocks * 512
        print(f"# huge4d.nii: {nvox:,} voxels ({nvox*4/1e9:.2f} GB logical), "
              f"{phys/1e3:.0f} KB on disk; {NX}x{NY}x{NZ}x{NT}")

        def check_all(outpath, expect, label):
            data, dim = read_nifti_data(outpath)
            n = len(data)
            bad = [(i, v) for i, v in enumerate(data) if abs(v - expect) > 1e-4]
            # spot-check first, middle, last too (already covered by full scan)
            if bad:
                fails.append(f"{label}: {len(bad)}/{n} voxels != {expect} (e.g. idx {bad[0][0]}={bad[0][1]})")
                print(f"  FAIL {label}: {len(bad)}/{n} wrong")
            else:
                print(f"  ok   {label}: all {n} output voxels == {expect}")

        # A: safe modifiers, >INT_MAX positive count, and Tmean -> 100.0 everywhere
        outA = os.path.join(d, 'A.nii')
        rc, log = run(niimath, [huge, '-p', '1', '-add', '1', '-ing', '100',
                               '-Tmean', '-gz', '0', outA], d)
        if rc != 0:
            fails.append(f"A exit {rc}: {log.strip()}"); print(f"  FAIL A exit {rc}: {log.strip()[:200]}")
        else:
            check_all(outA, 100.0, "A (-p 1 -add 1 -ing 100 -Tmean -gz 0)")

        # B: -Tmean of zeros -> 0.0
        outB = os.path.join(d, 'B.nii')
        rc, log = run(niimath, [huge, '-Tmean', outB], d)
        if rc != 0:
            fails.append(f"B exit {rc}: {log.strip()}"); print(f"  FAIL B exit {rc}: {log.strip()[:200]}")
        else:
            check_all(outB, 0.0, "B (-Tmean of zeros)")

        # C: chained > INT_MAX ops then reduce -> 0.0
        outC = os.path.join(d, 'C.nii')
        rc, log = run(niimath, [huge, '-mul', '3', '-sub', '3', '-Tmean', outC], d)
        if rc != 0:
            fails.append(f"C exit {rc}: {log.strip()}"); print(f"  FAIL C exit {rc}: {log.strip()[:200]}")
        else:
            check_all(outC, -3.0, "C (-mul 3 -sub 3 -Tmean)")   # 0*3-3 = -3

        # D: out-of-scope op must be rejected (exit 2), not silently processed
        outD = os.path.join(d, 'D.nii')
        rc, log = run(niimath, [huge, '-fmean', outD], d)
        if rc == 2 and 'huge-image capable' in log:
            print(f"  ok   D (-fmean rejected exit 2)")
        else:
            fails.append(f"D not rejected: exit {rc}: {log.strip()}")
            print(f"  FAIL D: expected exit 2 reject, got {rc}: {log.strip()[:200]}")

        # E: double-dash unsafe aliases must also be recognized as flags and rejected.
        outE = os.path.join(d, 'E.nii')
        rc, log = run(niimath, [huge, '--compare', huge, outE], d)
        if rc == 2 and 'huge-image capable' in log:
            print(f"  ok   E (--compare rejected exit 2)")
        else:
            fails.append(f"E not rejected: exit {rc}: {log.strip()}")
            print(f"  FAIL E: expected exit 2 reject, got {rc}: {log.strip()[:200]}")

        # F: one histogram bin contains >INT_MAX zeros; its counter must remain wide.
        skewed = os.path.join(d, 'skewed4d.nii')
        write_sparse_nifti(skewed, NX, NY, NZ, NT, first_value=1.0)
        outF = os.path.join(d, 'F.nii')
        rc, log = run(niimath, [skewed, '-thrp', '50', '-Tsum', outF], d)
        if rc != 0:
            fails.append(f"F exit {rc}: {log.strip()}")
            print(f"  FAIL F exit {rc}: {log.strip()[:200]}")
        else:
            data, _ = read_nifti_data(outF)
            if abs(data[0] - 1.0) > 1e-4 or any(abs(v) > 1e-4 for v in data[1:]):
                fails.append("F robust histogram produced unexpected output")
                print("  FAIL F: expected one leading 1.0 and otherwise zero")
            else:
                print("  ok   F (-thrp huge-bin count): one leading 1.0, otherwise zero")

        # Fail-closed admission is re-checked per-op against the CURRENT image, with EXACT
        # op matching. These must all exit 2 (rejected), never silently corrupt.
        small = os.path.join(d, 'small3d.nii')
        write_sparse_nifti(small, NX, NY, NZ, 1)          # matches huge's 3D dims (nvox3D)
        outX = os.path.join(d, 'X.nii')

        def check_reject(args, label):
            rc, log = run(niimath, args, d)
            if rc == 2 and 'not huge-image capable' in log:
                print(f"  ok   {label} (rejected exit 2)")
            else:
                fails.append(f"{label}: expected exit 2 reject, got {rc}: {log.strip()[:160]}")
                print(f"  FAIL {label}: expected exit 2, got {rc}")

        # G: disguised operation token is rejected by exact admission/dispatch matching.
        check_reject([huge, 'x-otsu', outX], "G (disguised 'x-otsu' on huge)")
        # H: scratch is hoisted, but this temporal op is not promoted without its own huge oracle.
        check_reject([huge, '-Tstd', outX], "H (-Tstd not huge-safe)")
        # I: restart replacement is header-preflighted against its next operation.
        check_reject([small, '-restart', huge, '-otsu', outX], "I (-restart-to-huge then -otsu)")
        # J: image becomes huge via an image-image binary adopting a huge 4D operand.
        check_reject([small, '-add', huge, '-otsu', outX], "J (small -add huge4d then -otsu)")

    finally:
        if keep:
            print(f"# kept scratch dir {d}")
        else:
            shutil.rmtree(d, ignore_errors=True)

    if fails:
        print("\nHUGE SMOKE FAILED:")
        for m in fails:
            print("  - " + m)
        sys.exit(1)
    print("\nhuge smoke passed")

if __name__ == '__main__':
    main()
