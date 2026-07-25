/* featureParity.test.ts — feature-completeness gate for the WASI-C backend.
 *
 * Exercises each supported op against the REAL reactor (js/src/niimath-wasi.wasm) and, where a
 * numeric oracle applies, compares to the native binary on the same input. Floating paths use
 * explicit per-op tolerances (recorded here); integer/simple paths are exact.
 *
 * `bun run test:wasi` verifies/regenerates every fixture and requires both the native
 * `src/niimath` reference binary and `make -C src wasm-wasi` reactor.
 */
import { test, expect, beforeAll, describe } from "bun:test";
import { readFileSync, writeFileSync, existsSync, mkdirSync, rmSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { fileURLToPath } from "node:url";
import { WasiRunner } from "./WasiRunner";

const ROOT = fileURLToPath(new URL("../../../", import.meta.url)); // repo root
const FX = `${ROOT}js/src/wasi/fixtures`;
const NII = `${ROOT}src/niimath${process.platform === "win32" ? ".exe" : ""}`;
const WASM = `${ROOT}js/src/niimath-wasi.wasm`;
const TMP = join(tmpdir(), "niimath-wasi-tests");
mkdirSync(TMP, { recursive: true });
const tmp = (name: string) => join(TMP, name);

const rd = (p: string) => new Uint8Array(readFileSync(p));
function nativeRaw(args: string[]): void {
  execFileSync(NII, args, { env: { ...process.env, FSLOUTPUTTYPE: "NIFTI" } });
}
function payloadFloat(buf: Uint8Array): Float32Array {
  const dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
  const off = dv.getFloat32(108, true);
  return new Float32Array(buf.buffer, buf.byteOffset + off, (buf.byteLength - off) / 4);
}
function maxAbsDiff(a: Float32Array, b: Float32Array): number {
  expect(a.length).toBe(b.length);
  let m = 0;
  for (let i = 0; i < a.length; i++) {
    const d = Math.abs(a[i] - b[i]);
    if (!Number.isFinite(d)) return Infinity;
    if (d > m) m = d;
  }
  return m;
}

function withCenterValue(raw: Uint8Array, value: number): Uint8Array {
  const out = raw.slice();
  const dv = new DataView(out.buffer, out.byteOffset, out.byteLength);
  const nx = dv.getInt16(42, true);
  const ny = dv.getInt16(44, true);
  const nz = dv.getInt16(46, true);
  const off = dv.getFloat32(108, true);
  const i = Math.floor(nx / 2) + Math.floor(ny / 2) * nx + Math.floor(nz / 2) * nx * ny;
  dv.setFloat32(off + i * 4, value, true);
  return out;
}

function nonFiniteCounts(buf: Uint8Array) {
  let nan = 0, pos = 0, neg = 0;
  for (const v of payloadFloat(buf)) {
    if (Number.isNaN(v)) nan++;
    else if (v === Infinity) pos++;
    else if (v === -Infinity) neg++;
  }
  return { nan, pos, neg };
}
function corr(a: Float32Array, b: Float32Array): number {
  let sa = 0, sb = 0, n = a.length;
  for (let i = 0; i < n; i++) { sa += a[i]; sb += b[i]; }
  const ma = sa / n, mb = sb / n;
  let num = 0, da = 0, db = 0;
  for (let i = 0; i < n; i++) { const x = a[i] - ma, y = b[i] - mb; num += x * y; da += x * x; db += y * y; }
  return num / Math.sqrt(da * db);
}

let runner: WasiRunner;
const fixture = () => rd(`${FX}/t2w_256x256x192_int16.nii`);

beforeAll(async () => {
  if (!existsSync(WASM)) throw new Error(`missing ${WASM} — run: make -C src wasm-wasi`);
  if (!existsSync(`${FX}/t2w_256x256x192_int16.nii`)) throw new Error("run gen_fixture.py + gen_test_fixtures.py");
  runner = await WasiRunner.create(rd(WASM));
});

describe("reactor", () => {
  test("two sequential runs on one instance; reset removes all files", async () => {
    let r = await runner.runFiles({ argv: ["in.nii", "-add", "1", "a.nii"], inputs: { "in.nii": fixture() }, outputs: ["a.nii"] });
    expect(r.exitCode).toBe(0);
    expect(r.files["a.nii"]).toBeTruthy();
    r = await runner.runFiles({ argv: ["in.nii", "-mul", "2", "b.nii"], inputs: { "in.nii": fixture() }, outputs: ["b.nii"] });
    expect(r.exitCode).toBe(0);
    expect(r.files["b.nii"]).toBeTruthy();
    await runner.reset();
    expect(runner.listFiles()).toEqual([]);
  });
});

describe("errors", () => {
  test("malformed input returns nonzero; next run succeeds on a fresh instance", async () => {
    let r = await runner.runFiles({ argv: ["nope.nii", "-add", "1", "x.nii"], inputs: {}, outputs: ["x.nii"] });
    expect(r.exitCode).not.toBe(0);
    expect(r.files["x.nii"]).toBeUndefined();
    r = await runner.runFiles({ argv: ["in.nii", "-add", "1", "ok.nii"], inputs: { "in.nii": fixture() }, outputs: ["ok.nii"] });
    expect(r.exitCode).toBe(0);
    expect(r.files["ok.nii"]).toBeTruthy();
  });
  test("unsupported option returns a clear nonzero error, not a silent no-op", async () => {
    const r = await runner.runFiles({ argv: ["in.nii", "-nosuchop", "y.nii"], inputs: { "in.nii": fixture() }, outputs: ["y.nii"] });
    expect(r.exitCode).not.toBe(0);
  });
});

describe("core ops vs native (exact for integer/simple; tolerance for floating)", () => {
  // op label -> [argv, tolerance]. tol 0 = byte-exact payload.
  const cases: [string, string[], number][] = [
    ["scalar add+mul+sqrt", ["-add", "10", "-mul", "2", "-sqrt"], 0],
    ["binary op (-thr)", ["-thr", "500"], 0],
    ["datatype convert (-bin)", ["-bin"], 0],
    ["temporal no-op (-Tmean 3D)", ["-Tmean"], 0],
    ["filter gaussian (-s 3)", ["-s", "3"], 2e-3],
    ["filter DoG (-dog 2 3.2)", ["-dog", "2", "3.2"], 1.5],
    // float box-gather sums in a different order on native (SSE) vs wasm (scalar/simd128), so
    // -fmean diverges ~1.7e-3 abs (~1e-5 relative) across builds. The exact-match double-intermediate
    // separable path that once justified 1e-3 was reverted (see noncompliant.md), so this tolerance
    // reflects the restored gather's real cross-build rounding margin.
    ["filter box mean (-fmean)", ["-kernel", "boxv", "5", "-fmean"], 5e-3],
    ["edge (-edge)", ["-edge"], 0],
  ];
  for (const [label, ops, tol] of cases) {
    test(label, async () => {
      const ref = tmp(`parity_${label.replace(/\W+/g, "_")}.nii`);
      nativeRaw([`${FX}/t2w_256x256x192_int16.nii`, ...ops, ref]);
      const r = await runner.runFiles({ argv: ["in.nii", ...ops, "out.nii"], inputs: { "in.nii": fixture() }, outputs: ["out.nii"] });
      expect(r.exitCode).toBe(0);
      const d = maxAbsDiff(payloadFloat(r.files["out.nii"]), payloadFloat(rd(ref)));
      expect(d).toBeLessThanOrEqual(tol);
    }, 30000); // heavy: 25 MB volume + native reference spawn
  }
});

describe("local non-finite filter semantics", () => {
  for (const [label, value, expected] of [
    ["NaN", NaN, { nan: 27, pos: 0, neg: 0 }],
    ["+Inf", Infinity, { nan: 0, pos: 27, neg: 0 }],
    ["-Inf", -Infinity, { nan: 0, pos: 0, neg: 27 }],
  ] as const) {
    for (const op of ["-fmean", "-fmeanu"]) {
      test(`${op} boxv3 keeps one ${label} within its 3x3x3 neighborhood`, async () => {
        const input = withCenterValue(rd(`${FX}/t1_small.nii`), value);
        const r = await runner.runFiles({
          argv: ["in.nii", "-kernel", "boxv", "3", op, "out.nii"],
          inputs: { "in.nii": input },
          outputs: ["out.nii"],
        });
        expect(r.exitCode).toBe(0);
        expect(nonFiniteCounts(r.files["out.nii"])).toEqual(expected);
      });
    }
  }

  for (const [op, label, value, expected] of [
    ["-dilF", "NaN", NaN, { nan: 1, pos: 0, neg: 0 }],
    ["-dilF", "+Inf", Infinity, { nan: 0, pos: 27, neg: 0 }],
    ["-dilF", "-Inf", -Infinity, { nan: 0, pos: 0, neg: 0 }],
    ["-eroF", "NaN", NaN, { nan: 1, pos: 0, neg: 0 }],
    ["-eroF", "+Inf", Infinity, { nan: 0, pos: 0, neg: 0 }],
    ["-eroF", "-Inf", -Infinity, { nan: 0, pos: 0, neg: 27 }],
    ["-ero", "NaN", NaN, { nan: 1, pos: 0, neg: 0 }],
    ["-ero", "+Inf", Infinity, { nan: 0, pos: 1, neg: 0 }],
    ["-ero", "-Inf", -Infinity, { nan: 0, pos: 0, neg: 1 }],
  ] as const) {
    test(`${op} boxv3 preserves the historical ${label} footprint`, async () => {
      const input = withCenterValue(rd(`${FX}/t1_small.nii`), value);
      const r = await runner.runFiles({
        argv: ["in.nii", "-kernel", "boxv", "3", op, "out.nii"],
        inputs: { "in.nii": input },
        outputs: ["out.nii"],
      });
      expect(r.exitCode).toBe(0);
      expect(nonFiniteCounts(r.files["out.nii"])).toEqual(expected);
    });
  }
});

describe("temporal filter (bandpass)", () => {
  test("-bandpass runs on a 4D input and matches native within tolerance", async () => {
    const bold = rd(`${FX}/bold_4x4x4x20.nii`); // 20 timepoints (>= 12 the filter needs)
    const ref = tmp("parity_bandpass.nii");
    nativeRaw([`${FX}/bold_4x4x4x20.nii`, "-bandpass", "0.01", "0.1", "2", ref]);
    const r = await runner.runFiles({ argv: ["in.nii", "-bandpass", "0.01", "0.1", "2", "out.nii"], inputs: { "in.nii": bold }, outputs: ["out.nii"] });
    expect(r.exitCode).toBe(0);
    const d = maxAbsDiff(payloadFloat(r.files["out.nii"]), payloadFloat(rd(ref)));
    expect(d).toBeLessThanOrEqual(1e-2);
  });
});

describe("NIfTI formats", () => {
  test("NIfTI-1 single-file (covered above)", () => { expect(true).toBe(true); });
  test("NIfTI-2 single-file reads and processes", async () => {
    const ref = tmp("parity_nifti2.nii");
    nativeRaw([`${FX}/t2w_small_nifti2.nii`, "-add", "5", ref]);
    const r = await runner.runFiles({ argv: ["in.nii", "-add", "5", "out.nii"], inputs: { "in.nii": rd(`${FX}/t2w_small_nifti2.nii`) }, outputs: ["out.nii"] });
    expect(r.exitCode).toBe(0);
    const d = maxAbsDiff(payloadFloat(r.files["out.nii"]), payloadFloat(rd(ref)));
    expect(d).toBe(0);
  });
});

describe("multi-input: allineate/deface", () => {
  test("allineate runs the multi-input reactor path and reslices onto the base grid", async () => {
    // Functional only: the default cost is now the fast SPM/FLIRT-inspired engine, which is tuned
    // for real brain resolution/contrast and is NOT meaningful on this 24^3 synthetic fixture (it
    // can converge worse than the AFNI engine here). Registration QUALITY is validated on real data
    // (the gpl-build CI + benchmark/register pairs), not on a synthetic phantom. This test only
    // proves the reactor loads two inputs, runs -allineate, and returns a valid resliced volume.
    const r = await runner.runFiles({
      argv: ["mov.nii", "-allineate", "tmpl.nii", "out.nii"],
      inputs: { "mov.nii": rd(`${FX}/mov_small.nii`), "tmpl.nii": rd(`${FX}/tmpl_small.nii`) },
      outputs: ["out.nii"],
    });
    expect(r.exitCode).toBe(0);
    const out = payloadFloat(r.files["out.nii"]);
    const tmpl = payloadFloat(rd(`${FX}/tmpl_small.nii`));
    expect(out.length).toBe(tmpl.length);                    // resliced onto the base (template) grid
    expect(out.every((v) => Number.isFinite(v))).toBe(true); // valid output, no NaN/garbage
  });
  // Deface native/WASI parity, WITHIN A TOLERANCE. Registration (both the default fast engine
  // and -cost hel) is not byte-reproducible across builds — NEWUOA/pyramid + -ffast-math land on
  // a neighboring, equally-valid optimum under different codegen — so the resliced mask flips a
  // thin shell of BOUNDARY voxels between native and WASI. We assert the two agree everywhere
  // except that small shell, plus that the mask actually removed the identifiable NON-BRAIN head
  // tissue (face/scalp/skull), rather than byte identity. The realistic 24-object simbrain head
  // (js/src/wasi/simbrain.py) is what lets fast converge to essentially the same defacing as
  // native (fast-vs-hel KEPT Dice ~0.999); a lone Gaussian blob did not. NOTE on the removal
  // metric: the brain mask keeps only ~15% of the volume and zeros the rest, but ~67% of the
  // volume is air (already 0), so voxels VISIBLY changed by defacing (~19%) are exactly the
  // removed non-air head tissue — the privacy-critical part. We threshold on that, not on a
  // ">50% of voxels" count that would merely be re-zeroing background.
  for (const [label, extra] of [["fast (default)", []], ["-cost hel", ["-cost", "hel"]]] as const) {
    test(`deface ${label} native/WASI agree within tolerance (privacy parity)`, async () => {
      const ref = tmp(`parity_deface_${label.includes("hel") ? "hel" : "fast"}.nii`);
      nativeRaw([`${FX}/mov_small.nii`, "-deface", `${FX}/tmpl_small.nii`, `${FX}/mask_small.nii`, ...extra, ref]);
      const r = await runner.runFiles({
        argv: ["mov.nii", "-deface", "tmpl.nii", "mask.nii", ...extra, "out.nii"],
        inputs: {
          "mov.nii": rd(`${FX}/mov_small.nii`),
          "tmpl.nii": rd(`${FX}/tmpl_small.nii`),
          "mask.nii": rd(`${FX}/mask_small.nii`),
        },
        outputs: ["out.nii"],
      });
      expect(r.exitCode).toBe(0);
      const out = payloadFloat(r.files["out.nii"]);
      const nat = payloadFloat(rd(ref));
      const mov = payloadFloat(rd(`${FX}/mov_small.nii`));
      expect(out.every((v) => Number.isFinite(v))).toBe(true);
      let removed = 0, disagree = 0;
      for (let i = 0; i < out.length; i++) {
        if (out[i] !== mov[i]) removed++;                 // non-air head tissue zeroed by defacing
        if (Math.abs(out[i] - nat[i]) > 1e-3) disagree++; // native vs WASI mismatch
      }
      console.log(`[deface ${label}] removed(non-air)=${(100 * removed / out.length).toFixed(1)}% native/WASI-disagree=${(100 * disagree / out.length).toFixed(2)}%`);
      // Identifiable non-brain head tissue (face/scalp/skull) was actually removed (~19% of the
      // volume on this phantom; a no-op would be ~0), not just background re-zeroed.
      expect(removed / out.length).toBeGreaterThan(0.10);
      // Native/WASI agree except a thin registration-sensitive boundary shell. Measured on this
      // corrected phantom: 0.07% (fast) / 0.09% (hel) disagreement; the 1% bound keeps ~10x margin
      // for platform codegen variation while still catching a gross reslice/transform regression.
      expect(disagree / out.length).toBeLessThan(0.01);
    }, 30000); // -cost hel runs a full NEWUOA fit single-threaded in WASI; allow headroom
  }
});

describe("dtifit (text bvec/bval inputs + all expected outputs)", () => {
  const outs = ["FA", "MD", "L1", "L2", "L3", "V1", "V2", "V3", "S0", "MO", "tensor"].map((s) => `dti_${s}.nii`);
  test("produces every expected output file and matches native FA", async () => {
    for (const o of outs) rmSync(tmp(o), { force: true });
    execFileSync(NII, ["--dtifit", "-k", `${FX}/dwi_6x6x6x7.nii`, "-r", `${FX}/dwi.bvec`, "-b", `${FX}/dwi.bval`, "-m", `${FX}/dwi_mask.nii`, "-o", tmp("dti")], { env: { ...process.env, FSLOUTPUTTYPE: "NIFTI" } });
    const r = await runner.runFiles({
      argv: ["--dtifit", "-k", "dwi.nii", "-r", "dwi.bvec", "-b", "dwi.bval", "-m", "mask.nii", "-o", "dti"],
      inputs: {
        "dwi.nii": rd(`${FX}/dwi_6x6x6x7.nii`),
        "dwi.bvec": rd(`${FX}/dwi.bvec`),
        "dwi.bval": rd(`${FX}/dwi.bval`),
        "mask.nii": rd(`${FX}/dwi_mask.nii`),
      },
      outputs: outs,
    });
    expect(r.exitCode).toBe(0);
    for (const o of outs) expect(r.files[o], `missing ${o}`).toBeTruthy();
    const d = maxAbsDiff(payloadFloat(r.files["dti_FA.nii"]), payloadFloat(rd(tmp("dti_FA.nii"))));
    expect(d).toBeLessThanOrEqual(1e-3);
  });
});

describe("QC (TSV schema + numeric parity)", () => {
  test("writes a TSV whose schema + values match native", async () => {
    const qcRef = tmp("qc_native.tsv");
    execFileSync(NII, ["--qc", `${FX}/t1_small.nii`, "--seg", `${FX}/seg_small.nii`, "--csf", "1", "--wm", "3", "--out", qcRef]);
    const r = await runner.runFiles({
      argv: ["--qc", "t1.nii", "--seg", "seg.nii", "--csf", "1", "--wm", "3", "--out", "qc.tsv"],
      inputs: { "t1.nii": rd(`${FX}/t1_small.nii`), "seg.nii": rd(`${FX}/seg_small.nii`) },
      outputs: ["qc.tsv"],
    });
    expect(r.exitCode).toBe(0);
    const tsv = new TextDecoder().decode(r.files["qc.tsv"]).trim().split("\n");
    const nat = readFileSync(qcRef, "utf8").trim().split("\n");
    expect(tsv.length).toBe(2); // header + one value row
    expect(tsv[0]).toBe(nat[0]); // identical schema
    // numeric parity per column, small tolerance for floating metrics
    const hv = tsv[0].split("\t"), wv = tsv[1].split("\t"), nv = nat[1].split("\t");
    for (let i = 0; i < hv.length; i++) {
      const a = parseFloat(wv[i]), b = parseFloat(nv[i]);
      if (Number.isFinite(a) && Number.isFinite(b)) {
        const rel = Math.abs(a - b) / (Math.abs(b) + 1e-9);
        expect(rel, `col ${hv[i]}: wasi=${wv[i]} native=${nv[i]}`).toBeLessThan(1e-3);
      }
    }
  });
});

/* Minimal NIfTI-1 (n+1) float32 writer, so the ROMEO fixtures below are synthesized in-test
   rather than committed: -romeo needs a WRAPPED phase, which none of the existing fixtures are. */
function writeF32Nifti(nx: number, ny: number, nz: number, nt: number, data: Float32Array): Uint8Array {
  const buf = new Uint8Array(352 + data.length * 4);
  const dv = new DataView(buf.buffer);
  dv.setInt32(0, 348, true);
  const dims = [nt > 1 ? 4 : 3, nx, ny, nz, nt, 1, 1, 1];
  for (let i = 0; i < 8; i++) dv.setInt16(40 + i * 2, dims[i], true);
  dv.setInt16(70, 16, true);  // DT_FLOAT32
  dv.setInt16(72, 32, true);  // bitpix
  const pix = [1, 1, 1, 1, 0, 0, 0, 0];
  for (let i = 0; i < 8; i++) dv.setFloat32(76 + i * 4, pix[i], true);
  dv.setFloat32(108, 352, true); // vox_offset
  dv.setFloat32(112, 1, true);   // scl_slope
  buf[123] = 2;                  // xyz_units = mm
  dv.setInt16(252, 1, true);     // qform_code
  dv.setInt16(254, 1, true);     // sform_code
  dv.setFloat32(280, 1, true); dv.setFloat32(300, 1, true); dv.setFloat32(320, 1, true);
  buf.set(new Uint8Array([0x6e, 0x2b, 0x31, 0x00]), 344); // "n+1\0"
  new Float32Array(buf.buffer, 352, data.length).set(data);
  return buf;
}

describe("ROMEO phase unwrapping (-romeo)", () => {
  const NX = 12, NY = 10, NZ = 8, N = NX * NY * NZ, TWO_PI = 2 * Math.PI;
  const truth = new Float32Array(N);
  const wrapped = new Float32Array(N);
  const mag = new Float32Array(N);
  for (let k = 0, i = 0; k < NZ; k++) for (let j = 0; j < NY; j++) for (let x = 0; x < NX; x++, i++) {
    const t = 0.9 * x + 0.4 * j + 0.25 * k;
    truth[i] = t;
    wrapped[i] = t - TWO_PI * Math.round(t / TWO_PI);
    const core = x >= 2 && x < NX - 2 && j >= 2 && j < NY - 2 && k >= 1 && k < NZ - 1;
    mag[i] = core ? 900 + ((x + j + k) % 7) : 10;
  }
  const phaseNii = writeF32Nifti(NX, NY, NZ, 1, wrapped);
  const magNii = writeF32Nifti(NX, NY, NZ, 1, mag);

  test("wasm32 reactor unwraps identically to the native binary (and to ground truth)", async () => {
    const pPath = tmp("romeo_phase.nii"), mPath = tmp("romeo_mag.nii"), nPath = tmp("romeo_native.nii");
    writeFileSync(pPath, phaseNii);
    writeFileSync(mPath, magNii);
    nativeRaw([pPath, "-romeo", mPath, "-t", "5.0", "-k", "nomask", "-no-phase-rescale", nPath]);

    const r = await runner.runFiles({
      argv: ["p.nii", "-romeo", "m.nii", "-t", "5.0", "-k", "nomask", "-no-phase-rescale", "out.nii"],
      inputs: { "p.nii": phaseNii, "m.nii": magNii },
      outputs: ["out.nii"],
    });
    expect(r.exitCode).toBe(0);
    const wasi = payloadFloat(r.files["out.nii"]);
    const native = payloadFloat(rd(nPath));
    // wasm32 vs native arm64/x86: same strict-FP source, so require exact agreement
    expect(maxAbsDiff(wasi, native)).toBe(0);
    // and the unwrap itself must be right: ONE constant 2*pi offset from the ground truth
    const offsets = new Set<number>();
    for (let i = 0; i < N; i++) {
      const d = wasi[i] - truth[i];
      const w = Math.round(d / TWO_PI);
      expect(Math.abs(d - w * TWO_PI)).toBeLessThan(1e-3);
      offsets.add(w);
    }
    expect(offsets.size).toBe(1);
  });

  test("robustmask writes a binary <out>_mask side output", async () => {
    const r = await runner.runFiles({
      argv: ["p.nii", "-romeo", "m.nii", "-t", "5.0", "-no-phase-rescale", "out.nii"],
      inputs: { "p.nii": phaseNii, "m.nii": magNii },
      outputs: ["out.nii", "out_mask.nii"],
    });
    expect(r.exitCode).toBe(0);
    const mask = payloadFloat(r.files["out_mask.nii"]);
    let inside = 0;
    for (const v of mask) { expect(v === 0 || v === 1).toBe(true); if (v === 1) inside++; }
    expect(inside).toBeGreaterThan(0);
    expect(inside).toBeLessThan(N);
  });

  test("an unported ROMEO option fails with a specific message, not a silent no-op", async () => {
    const r = await runner.runFiles({
      argv: ["p.nii", "-romeo", "m.nii", "-t", "5.0", "-w", "bestpath", "o.nii"],
      inputs: { "p.nii": phaseNii, "m.nii": magNii },
      outputs: [],
    });
    expect(r.exitCode).not.toBe(0);
    expect(r.stderr).toContain("not implemented");
  });
});

describe("reach goals excluded from initial scope fail clearly", () => {
  test("-bitmap is unsupported (not compiled) and errors", async () => {
    const r = await runner.runFiles({ argv: ["in.nii", "-bitmap", "o.png"], inputs: { "in.nii": fixture() }, outputs: [] });
    expect(r.exitCode).not.toBe(0);
    expect(r.stderr.toLowerCase()).toContain("unsupported");
  });
  test("-mesh is unsupported (not compiled) and errors", async () => {
    const r = await runner.runFiles({ argv: ["in.nii", "-mesh", "o.mz3"], inputs: { "in.nii": fixture() }, outputs: [] });
    expect(r.exitCode).not.toBe(0);
    expect(r.stderr.toLowerCase()).toContain("not compiled");
  });
});
