/* compression.test.ts — M3 "Compression" feature row + M4 gzip wrapper.
 * gzip input and gzip output through the host-side CompressionStream path (module stays zlib-free).
 * Runs under Bun; the same assertions are re-driven by scripts/wasi-node-smoke.mjs under Node 20+.
 */
import { test, expect, beforeAll } from "bun:test";
import { readFileSync, existsSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { tmpdir } from "node:os";
import { join } from "node:path";
import { fileURLToPath } from "node:url";
import { packArgv, WasiRunner } from "./WasiRunner";
import { isGzip, gzip, gunzip } from "./gzip";
import { ERRNO, WasiVfs } from "./WasiVfs";

const ROOT = fileURLToPath(new URL("../../../", import.meta.url));
const FX = `${ROOT}js/src/wasi/fixtures`;
const NII = `${ROOT}src/niimath${process.platform === "win32" ? ".exe" : ""}`;
const WASM = `${ROOT}js/src/niimath-wasi.wasm`;
const rd = (p: string) => new Uint8Array(readFileSync(p));

function payloadFloat(buf: Uint8Array): Float32Array {
  const dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
  const off = dv.getFloat32(108, true);
  return new Float32Array(buf.buffer, buf.byteOffset + off, (buf.byteLength - off) / 4);
}

let runner: WasiRunner;
beforeAll(async () => {
  if (!existsSync(WASM)) throw new Error("run: make -C src wasm-wasi");
  runner = await WasiRunner.create(rd(WASM));
});

test("gzip round-trips through CompressionStream", async () => {
  const raw = rd(`${FX}/t2w_small_nifti2.nii`);
  const gz = await gzip(raw);
  expect(isGzip(gz)).toBe(true);
  expect(isGzip(raw)).toBe(false);
  const back = await gunzip(gz);
  expect(Buffer.compare(Buffer.from(back), Buffer.from(raw))).toBe(0);
});

test("gzip input is decompressed + argv rewritten; gzip output is recompressed; matches native raw", async () => {
  // native raw reference
  const ref = join(tmpdir(), "niimath-wasi-comp-ref.nii");
  execFileSync(NII, [`${FX}/t2w_small_nifti2.nii`, "-add", "3", ref], { env: { ...process.env, FSLOUTPUTTYPE: "NIFTI" } });

  const gzIn = await gzip(rd(`${FX}/t2w_small_nifti2.nii`));
  const r = await runner.runFilesGz({
    argv: ["in.nii.gz", "-add", "3", "out.nii.gz"],
    inputs: { "in.nii.gz": gzIn },
    outputs: ["out.nii.gz"],
  });
  expect(r.exitCode).toBe(0);
  const outGz = r.files["out.nii.gz"];
  expect(outGz).toBeTruthy();
  expect(isGzip(outGz)).toBe(true); // output is gzip-framed

  const outRaw = await gunzip(outGz);
  const d = maxAbs(payloadFloat(outRaw), payloadFloat(rd(ref)));
  expect(d).toBe(0);
});

test("mixed: gzip input with a raw output", async () => {
  const gzIn = await gzip(rd(`${FX}/t2w_small_nifti2.nii`));
  const r = await runner.runFilesGz({
    argv: ["in.nii.gz", "-mul", "2", "out.nii"],
    inputs: { "in.nii.gz": gzIn },
    outputs: ["out.nii"],
  });
  expect(r.exitCode).toBe(0);
  expect(isGzip(r.files["out.nii"])).toBe(false); // raw output stays raw
});

test("gzip input and output may use the same logical path", async () => {
  const raw = rd(`${FX}/t2w_small_nifti2.nii`);
  const gz = await gzip(raw);
  const r = await runner.runFilesGz({
    argv: ["same.nii.gz", "-add", "1", "same.nii.gz"],
    inputs: { "same.nii.gz": gz },
    outputs: ["same.nii.gz"],
  });
  expect(r.exitCode).toBe(0);
  expect(isGzip(r.files["same.nii.gz"])).toBe(true);
  const back = await gunzip(r.files["same.nii.gz"]);
  const separate = await runner.runFilesGz({
    argv: ["same.nii.gz", "-add", "1", "out.nii.gz"],
    inputs: { "same.nii.gz": gz },
    outputs: ["out.nii.gz"],
  });
  expect(separate.exitCode).toBe(0);
  const expected = await gunzip(separate.files["out.nii.gz"]);
  expect(Buffer.compare(Buffer.from(back), Buffer.from(expected))).toBe(0);
  expect(runner.listFiles()).toEqual([]); // gzip adapter releases raw VFS state before returning
});

test("multiple decompressed inputs share one aggregate quota", async () => {
  const capped = await WasiRunner.create(rd(WASM), { env: { FSLOUTPUTTYPE: "NIFTI" } }, { maxBytes: 12 });
  const gz = await gzip(new Uint8Array(8));
  await expect(capped.runFilesGz({
    argv: ["a.nii.gz", "b.nii.gz"],
    inputs: { "a.nii.gz": gz, "b.nii.gz": gz },
  })).rejects.toThrow("exceeds cap");
});

test("file bytes and captured logs share one VFS quota", () => {
  const vfs = new WasiVfs({ maxBytes: 10, maxFiles: 2 });
  vfs.addFile("in.nii", new Uint8Array(8));
  vfs.captureStd(1, new Uint8Array(5));
  expect(vfs.takeStdout().byteLength).toBe(2);
  const opened = vfs.open("out.nii", 1);
  expect(opened.fd).toBeNumber();
  expect(vfs.writeFile(opened.fd!, [new Uint8Array(1)])).toBe(-ERRNO.NOSPC);
});

test("runner is single-flight: an overlapping run is rejected, not interleaved", async () => {
  const raw = rd(`${FX}/t2w_small_nifti2.nii`);
  const first = runner.runFiles({ argv: ["in.nii", "-add", "1", "a.nii"], inputs: { "in.nii": raw }, outputs: ["a.nii"] });
  await expect(
    runner.runFiles({ argv: ["in.nii", "-add", "2", "b.nii"], inputs: { "in.nii": raw }, outputs: ["b.nii"] }),
  ).rejects.toThrow("single-flight");
  const r1 = await first; // the in-flight run still completes cleanly
  expect(r1.exitCode).toBe(0);
});

test("staging is single-flight: a reset/addFile scheduled mid-run is rejected, run still completes", async () => {
  const raw = rd(`${FX}/t2w_small_nifti2.nii`);
  const run = runner.runFiles({ argv: ["in.nii", "-add", "1", "out.nii"], inputs: { "in.nii": raw }, outputs: ["out.nii"] });
  // inFlight is set synchronously by the lock, so these mid-run mutations (the audit repro: an
  // external reset landing in the run's first await gap, which would clear its staged input) must
  // be rejected here, before we await the run.
  expect(() => runner.addFile("x.nii", raw)).toThrow("single-flight");
  await expect(runner.reset()).rejects.toThrow("single-flight");
  const r = await run; // the in-flight run completes with its input intact
  expect(r.exitCode).toBe(0);
  expect(r.files["out.nii"]).toBeTruthy();
});

test("public readFile returns a detached copy", async () => {
  await runner.reset();
  runner.addFile("owned-by-runner.nii", new Uint8Array([1, 2, 3]));
  const copy = runner.readFile("owned-by-runner.nii")!;
  copy[0] = 99;
  expect(Array.from(runner.readFile("owned-by-runner.nii")!)).toEqual([1, 2, 3]);
});

test("argv packing rejects NUL injection and the reactor argument limit", () => {
  expect(() => packArgv(["niimath", "in.nii\0-add", "out.nii"])).toThrow("NUL");
  expect(() => packArgv(Array.from({ length: 256 }, () => "x"))).toThrow("255 arguments");
});

function maxAbs(a: Float32Array, b: Float32Array): number {
  expect(a.length).toBe(b.length);
  let m = 0;
  for (let i = 0; i < a.length; i++) { const d = Math.abs(a[i] - b[i]); if (d > m) m = d; }
  return m;
}
