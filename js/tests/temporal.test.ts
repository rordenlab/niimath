// Numerical gate for the two 4D temporal operations in the WASM bundle: `-stc` (slice-time
// correction) and `-moco` (motion correction).
//
// These are NOT smoke tests. Both commands are gated OFF on any platform that has not run a
// numerical check, and this file IS that check for wasm32 — so it asserts values, not exit codes.
// The `-stc` case has a closed form: an INTEGER sample shift makes the whole pipeline analytic
// (least-squares detrend, circular shift of the zero-padded residual, clip to the residual range,
// retrend, clip to the original range), so the expected output can be computed here in
// TypeScript with no FFT and no reference tool. That is the same strategy
// .github/scripts/release_smoke.py uses for the native build, and it is what makes the two
// platforms comparable.
//
// wasm32 differs from the native build in ways that could plausibly break either command:
// -DFORCE_INT32_MAX, a 32-bit size_t, no OpenMP (so `-stc` runs its batches serially and `-moco`
// its volumes serially), and a different code generator under the same -ffast-math contract.
import { describe, test, expect } from 'bun:test';
import { loadModule, run, BSD_MODULE } from './helpers';

// Presence is asserted functionally rather than by grepping the help banner: if either command
// were missing from the bundle, every test below would fail with a nonzero exit and no output.

/** 4D NIfTI-1 with a real TR and a temporal unit — `-stc` rejects a header without one. */
function make4d(
  nx: number, ny: number, nz: number, nt: number,
  fill: (x: number, y: number, z: number, t: number) => number,
  tr = 2.0,
): Uint8Array {
  const nvox = nx * ny * nz * nt;
  const buf = new ArrayBuffer(352 + nvox * 4);
  const dv = new DataView(buf);
  dv.setInt32(0, 348, true);
  const dim = [4, nx, ny, nz, nt, 1, 1, 1];
  for (let i = 0; i < 8; i++) dv.setInt16(40 + i * 2, dim[i], true);
  dv.setInt16(70, 16, true); // DT_FLOAT32
  dv.setInt16(72, 32, true);
  const pixdim = [1, 1, 1, 1, tr, 0, 0, 0];
  for (let i = 0; i < 8; i++) dv.setFloat32(76 + i * 4, pixdim[i], true);
  dv.setFloat32(108, 352, true);
  dv.setUint8(123, 2 | 8); // mm | sec  -- the temporal unit is load-bearing for -stc
  dv.setInt16(254, 1, true); // sform_code
  dv.setFloat32(280, 1, true);
  dv.setFloat32(296 + 4, 1, true);
  dv.setFloat32(312 + 8, 1, true);
  for (let i = 0; i < 4; i++) dv.setUint8(344 + i, 'n+1\0'.charCodeAt(i));
  const f = new Float32Array(buf, 352, nvox);
  const n3 = nx * ny * nz;
  for (let t = 0; t < nt; t++)
    for (let z = 0; z < nz; z++)
      for (let y = 0; y < ny; y++)
        for (let x = 0; x < nx; x++) f[t * n3 + x + y * nx + z * nx * ny] = fill(x, y, z, t);
  return new Uint8Array(buf);
}

function series(nii: Uint8Array, v: number, n3: number, nt: number): number[] {
  const dv = new DataView(nii.buffer, nii.byteOffset, nii.byteLength);
  const out: number[] = [];
  for (let t = 0; t < nt; t++) out.push(dv.getFloat32(352 + (t * n3 + v) * 4, true));
  return out;
}

/** Closed-form `-stc` output for an integer sample shift. Mirrors release_smoke.py. */
function stcReference(x: number[], shift: number): number[] {
  const nt = x.length;
  const half = (nt - 1) / 2;
  const mean = x.reduce((a, b) => a + b, 0) / nt;
  const sdd = (nt * (nt * nt - 1)) / 12;
  let sxy = 0;
  for (let i = 0; i < nt; i++) sxy += (i - half) * x[i];
  const slope = sxy / sdd;
  const trend = x.map((_, i) => mean + slope * (i - half));
  const xd = x.map((v, i) => v - trend[i]);
  const lo = Math.min(...xd), hi = Math.max(...xd);
  const xlo = Math.min(...x), xhi = Math.max(...x);
  return x.map((_, i) => {
    const j = i - shift;
    // outside [0, nt) the padded array is zero
    let y = j >= 0 && j < nt ? xd[j] : 0;
    y = Math.min(Math.max(y, lo), hi);
    return Math.min(Math.max(y + trend[i], xlo), xhi);
  });
}

const NX = 3, NY = 3, NZ = 2, NT = 16, TR = 2.0;
const N3 = NX * NY * NZ;
// Spiky, with a strong trend, so both clips have a chance to engage.
const sig = (v: number) => (t: number) =>
  100 + 6 * t + 30 * Math.sin(2.7 * t + v) + (t === 5 + (v % 3) ? 25 : 0);

describe('-stc in the WASM build', () => {
  test('matches the closed form for an exact +1 sample shift, and copies the zero-shift slice', async () => {
    const { mod, log } = await loadModule(BSD_MODULE);
    const input = make4d(NX, NY, NZ, NT, (x, y, z, t) => sig(x + y * NX + z * NX * NY)(t), TR);
    // times [0, TR] with -tzero 0: slice 0 shifts by 0 (skipped), slice 1 by exactly +1 sample.
    const r = await run(mod, log, [{ name: 'in.nii', data: input }],
      ['in.nii', '-stc', '--slicetiming', `0,${TR}`, '-tzero', '0', '-gz', '0', 'out.nii'], 'out.nii');
    expect(r.exitCode).toBe(0);
    expect(r.output).not.toBeNull();
    const out = r.output!;
    let clipSeen = false;
    for (let z = 0; z < NZ; z++) {
      for (let v0 = 0; v0 < NX * NY; v0++) {
        const v = v0 + z * NX * NY;
        const have = series(input, v, N3, NT);
        const got = series(out, v, N3, NT);
        if (z === 0) {
          expect(got).toEqual(have); // verbatim copy, bit-for-bit
          continue;
        }
        const want = stcReference(have, 1);
        for (let t = 0; t < NT; t++)
          expect(Math.abs(got[t] - want[t])).toBeLessThan(2e-3 * Math.max(1, Math.abs(want[t])));
        if (Math.min(...want) <= Math.min(...have) || Math.max(...want) >= Math.max(...have)) clipSeen = true;
      }
    }
    expect(clipSeen).toBe(true); // the fixture must actually exercise the output clip
  });

  test('reverses sign with -tzero at the far end', async () => {
    const { mod, log } = await loadModule(BSD_MODULE);
    const input = make4d(NX, NY, NZ, NT, (x, y, z, t) => sig(x + y * NX + z * NX * NY)(t), TR);
    const r = await run(mod, log, [{ name: 'in.nii', data: input }],
      ['in.nii', '-stc', '--slicetiming', `0,${TR}`, '-tzero', `${TR}`, '-gz', '0', 'out.nii'], 'out.nii');
    expect(r.exitCode).toBe(0);
    for (let v = 0; v < NX * NY; v++) {   // slice 0 now shifts by exactly -1
      const want = stcReference(series(input, v, N3, NT), -1);
      const got = series(r.output!, v, N3, NT);
      for (let t = 0; t < NT; t++)
        expect(Math.abs(got[t] - want[t])).toBeLessThan(2e-3 * Math.max(1, Math.abs(want[t])));
    }
  });

  test('leaves a pure linear time trend untouched, for any shift', async () => {
    const { mod, log } = await loadModule(BSD_MODULE);
    // the detrend annihilates a line and the retrend restores it, so the output must equal the input
    const input = make4d(NX, NY, NZ, NT, (x, y, z, t) => 5 + 3 * t + 0.25 * (x + y * NX + z * NX * NY), TR);
    const r = await run(mod, log, [{ name: 'in.nii', data: input }],
      ['in.nii', '-stc', '--slicetiming', '0,1.3', '-gz', '0', 'out.nii'], 'out.nii');
    expect(r.exitCode).toBe(0);
    for (let v = 0; v < N3; v++) {
      const have = series(input, v, N3, NT), got = series(r.output!, v, N3, NT);
      for (let t = 0; t < NT; t++) expect(Math.abs(got[t] - have[t])).toBeLessThan(1e-3);
    }
  });

  test('all-equal slice times are a bit-exact no-op', async () => {
    const { mod, log } = await loadModule(BSD_MODULE);
    const input = make4d(NX, NY, NZ, NT, (x, y, z, t) => sig(x + y * NX + z * NX * NY)(t), TR);
    const r = await run(mod, log, [{ name: 'in.nii', data: input }],
      ['in.nii', '-stc', '--slicetiming', '0.8,0.8', '-gz', '0', 'out.nii'], 'out.nii');
    expect(r.exitCode).toBe(0);
    for (let v = 0; v < N3; v++)
      expect(series(r.output!, v, N3, NT)).toEqual(series(input, v, N3, NT));
  });

  test('rejects a bad option group rather than computing something', async () => {
    const { mod, log } = await loadModule(BSD_MODULE);
    const input = make4d(NX, NY, NZ, NT, () => 1, TR);
    for (const args of [
      ['in.nii', '-stc', 'out.nii'],                                  // no --slicetiming
      ['in.nii', '-stc', '--slicetiming', '0', 'out.nii'],            // wrong count
      ['in.nii', '-stc', '--slicetiming', `0,${TR}`, '-tzero', '9', 'out.nii'], // tzero out of range
      ['in.nii', '-stc', '--slicetiming', '0,9', 'out.nii'],          // slice time above TR
    ]) {
      const r = await run(mod, log, [{ name: 'in.nii', data: input }], args, 'out.nii');
      expect(r.exitCode).not.toBe(0);
      expect(r.output).toBeNull();
    }
  });
});

describe('-moco in the WASM build', () => {
  const MX = 24, MY = 24, MZ = 20, MT = 3, SHIFT = 2;
  const cell = (x: number, y: number, z: number) =>
    140 * Math.exp(-(((x - 10) ** 2 + (y - 11) ** 2 + (z - 9) ** 2) / 18)) +
    90 * Math.exp(-(((x - 15) ** 2 + (y - 8) ** 2 + (z - 12) ** 2) / 9.7)) +
    70 * Math.exp(-(((x - 9) ** 2 + (y - 15) ** 2 + (z - 11) ** 2) / 8));

  test('recovers a known whole-voxel shift and passes volume 0 through unchanged', async () => {
    const { mod, log } = await loadModule(BSD_MODULE);
    // volume 1 is volume 0 displaced by an exact whole number of voxels along k
    const input = make4d(MX, MY, MZ, MT, (x, y, z, t) => {
      const src = t === 1 ? z - SHIFT : z;
      return src >= 0 && src < MZ ? cell(x, y, src) : 0;
    }, 2.0);
    const r = await run(mod, log, [{ name: 'in.nii', data: input }],
      ['in.nii', '-moco', '-gz', '0', 'out.nii'], 'out.nii');
    expect(r.exitCode).toBe(0);
    expect(r.output).not.toBeNull();
    const out = r.output!, n3 = MX * MY * MZ;

    // volume 0 is the base: copied through bit-for-bit
    expect(series(out, 0, n3, MT).slice(0, 1)).toEqual(series(input, 0, n3, MT).slice(0, 1));
    for (let v = 0; v < n3; v += 97) {
      const dv = new DataView(out.buffer, out.byteOffset, out.byteLength);
      const iv = new DataView(input.buffer, input.byteOffset, input.byteLength);
      expect(dv.getFloat32(352 + v * 4, true)).toBe(iv.getFloat32(352 + v * 4, true));
    }

    // and volume 1 now lands on volume 0 far better than it did before
    const interiorRms = (buf: Uint8Array, offA: number, offB: number) => {
      const dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
      let tot = 0, cnt = 0;
      for (let z = 4; z < MZ - 4; z++)
        for (let y = 4; y < MY - 4; y++)
          for (let x = 4; x < MX - 4; x++) {
            const i = x + y * MX + z * MX * MY;
            const d = dv.getFloat32(352 + (offA + i) * 4, true) - dv.getFloat32(352 + (offB + i) * 4, true);
            tot += d * d; cnt++;
          }
      return Math.sqrt(tot / cnt);
    };
    const before = interiorRms(input, 0, n3);
    const after = interiorRms(out, 0, n3);
    expect(after).toBeLessThan(before * 0.1);
  });

  test('rejects a 3D image', async () => {
    const { mod, log } = await loadModule(BSD_MODULE);
    const vol3d = make4d(8, 8, 8, 1, () => 1, 2.0);
    const r = await run(mod, log, [{ name: 'in.nii', data: vol3d }],
      ['in.nii', '-moco', '-gz', '0', 'out.nii'], 'out.nii');
    expect(r.exitCode).not.toBe(0);
  });

  // -ref reads a SECOND file, which on wasm32 means an emscripten virtual-FS read rather than a
  // real one, and it must be told apart from a volume number without touching that FS at all.
  test('-ref accepts a volume number and an external reference image', async () => {
    const { mod, log } = await loadModule(BSD_MODULE);
    // volume 0 is displaced by +SHIFT along k; volumes 1 and 2 are the undisplaced object, so
    // registering onto EITHER of them has the same known answer for volume 0.
    const input = make4d(MX, MY, MZ, MT, (x, y, z, t) => {
      const src = t === 0 ? z - SHIFT : z;
      return src >= 0 && src < MZ ? cell(x, y, src) : 0;
    }, 2.0);
    const ref = make4d(MX, MY, MZ, 1, (x, y, z) => cell(x, y, z), 2.0);
    const n3 = MX * MY * MZ;
    const at = (buf: Uint8Array, i: number) =>
      new DataView(buf.buffer, buf.byteOffset, buf.byteLength).getFloat32(352 + i * 4, true);
    const interiorRms = (a: Uint8Array, offA: number, b: Uint8Array, offB: number) => {
      let tot = 0, cnt = 0;
      for (let z = 4; z < MZ - 4; z++)
        for (let y = 4; y < MY - 4; y++)
          for (let x = 4; x < MX - 4; x++) {
            const i = x + y * MX + z * MX * MY;
            const d = at(a, offA + i) - at(b, offB + i);
            tot += d * d; cnt++;
          }
      return Math.sqrt(tot / cnt);
    };
    const before = interiorRms(input, 0, input, n3);

    // "-ref 1": volume 1 is the base, so IT is the sub-brick copied through -- not volume 0
    const byIndex = await run(mod, log, [{ name: 'in.nii', data: input }],
      ['in.nii', '-moco', '-ref', '1', '-gz', '0', 'out.nii'], 'out.nii');
    expect(byIndex.exitCode).toBe(0);
    expect(byIndex.output).not.toBeNull();
    for (let v = 0; v < n3; v += 97)
      expect(at(byIndex.output!, n3 + v)).toBe(at(input, n3 + v));
    expect(interiorRms(byIndex.output!, 0, byIndex.output!, n3)).toBeLessThan(before * 0.1);

    // "-ref ref.nii": the same base supplied as a SECOND file, read from the emscripten virtual
    // filesystem. Same registration, so the same corrected volume 0.
    const byFile = await run(mod, log,
      [{ name: 'in.nii', data: input }, { name: 'ref.nii', data: ref }],
      ['in.nii', '-moco', '-ref', 'ref.nii', '-gz', '0', 'out.nii'], 'out.nii');
    expect(byFile.exitCode).toBe(0);
    expect(byFile.output).not.toBeNull();
    expect(interiorRms(byFile.output!, 0, ref, 0)).toBeLessThan(before * 0.1);
    for (let v = 0; v < n3; v += 97)
      expect(at(byFile.output!, v)).toBeCloseTo(at(byIndex.output!, v), 2);

    // a volume number is decided by its digits, never by a filesystem probe, so an index past the
    // end is refused even though a file called "9" does not exist either
    const bad = await run(mod, log, [{ name: 'in.nii', data: input }],
      ['in.nii', '-moco', '-ref', '9', '-gz', '0', 'out.nii'], 'out.nii');
    expect(bad.exitCode).not.toBe(0);

    // a reference on a different grid is rejected, not silently resliced
    const small = make4d(MX / 2, MY, MZ, 1, () => 1, 2.0);
    const mismatch = await run(mod, log,
      [{ name: 'in.nii', data: input }, { name: 'small.nii', data: small }],
      ['in.nii', '-moco', '-ref', 'small.nii', '-gz', '0', 'out.nii'], 'out.nii');
    expect(mismatch.exitCode).not.toBe(0);
  });
});
