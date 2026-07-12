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
import { readFileSync, existsSync, mkdirSync, rmSync } from "node:fs";
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
  test("deface -cost hel matches native byte-stable output (privacy parity)", async () => {
    // Byte-exact cross-build parity requires the ORDINARY engine: -deface -cost hel converges
    // identically here and reslices via nii_reslice_affine, which is byte-stable across builds
    // (CLAUDE.md). The DEFAULT fast engine's registration is NOT byte-reproducible across builds
    // (independent SIMD/codegen -> a neighboring optimum), so it is checked for agreement, not
    // identity, in the next test. This must also actually modify the input (mask applied).
    const ref = tmp("parity_deface_hel.nii");
    nativeRaw([`${FX}/mov_small.nii`, "-deface", `${FX}/tmpl_small.nii`, `${FX}/mask_small.nii`, "-cost", "hel", ref]);
    const r = await runner.runFiles({
      argv: ["mov.nii", "-deface", "tmpl.nii", "mask.nii", "-cost", "hel", "out.nii"],
      inputs: {
        "mov.nii": rd(`${FX}/mov_small.nii`),
        "tmpl.nii": rd(`${FX}/tmpl_small.nii`),
        "mask.nii": rd(`${FX}/mask_small.nii`),
      },
      outputs: ["out.nii"],
    });
    expect(r.exitCode).toBe(0);
    const out = payloadFloat(r.files["out.nii"]);
    const mov = payloadFloat(rd(`${FX}/mov_small.nii`));
    let changed = 0;
    for (let i = 0; i < out.length; i++) if (out[i] !== mov[i]) changed++;
    expect(changed).toBeGreaterThan(0); // the mask actually removed/altered content
    expect(maxAbsDiff(out, payloadFloat(rd(ref)))).toBe(0); // byte-stable parity with native (ordinary engine)
  });
  test("default (fast) deface removes the same region as native (privacy agreement)", async () => {
    // Default -deface uses the fast engine (not byte-reproducible across builds). Assert a
    // privacy-MEANINGFUL agreement rather than byte identity: WASI and native remove nearly the
    // same voxels (high overlap of the removed set), not merely "some voxels changed".
    const ref = tmp("parity_deface_fast.nii");
    nativeRaw([`${FX}/mov_small.nii`, "-deface", `${FX}/tmpl_small.nii`, `${FX}/mask_small.nii`, ref]);
    const r = await runner.runFiles({
      argv: ["mov.nii", "-deface", "tmpl.nii", "mask.nii", "out.nii"],
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
    // Removed set = voxels the deface changed from the original. Compare WASI vs native.
    let inter = 0, uni = 0, wasiRemoved = 0;
    for (let i = 0; i < out.length; i++) {
      const a = out[i] !== mov[i];
      const b = nat[i] !== mov[i];
      if (a) wasiRemoved++;
      if (a && b) inter++;
      if (a || b) uni++;
    }
    expect(wasiRemoved).toBeGreaterThan(0);       // the mask was actually applied
    expect(inter / uni).toBeGreaterThan(0.9);     // same region removed as native (Jaccard > 0.9)
  });
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
