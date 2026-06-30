// Shared test helpers. These tests load the *built* Emscripten WASM modules
// (dist/niimath.js and dist/niimath-gpl.js) directly and drive them through the
// in-memory filesystem + callMain — the same entry points the worker uses, but
// without a browser Worker, so they run under `bun test` in Node/Bun.
//
// Run `bun run build` (or at least `bun run makeWasm`) before these tests so the
// dist artifacts exist.
import { gunzipSync } from 'node:zlib';
import { existsSync } from 'node:fs';

export interface EmscriptenModule {
  FS_createDataFile(parent: string, name: string, data: Uint8Array, canRead: boolean, canWrite: boolean): void;
  callMain(args: string[]): number;
  FS_readFile(filename: string): Uint8Array;
  FS_unlink(filename: string): void;
}
type Factory = (overrides?: Record<string, unknown>) => Promise<EmscriptenModule>;

export const BSD_MODULE = new URL('../dist/niimath.js', import.meta.url).pathname;
export const GPL_MODULE = new URL('../dist/niimath-gpl.js', import.meta.url).pathname;

export const gplBuilt = existsSync(GPL_MODULE);

/** Load an Emscripten module factory and capture its stdout/stderr into `log`. */
export async function loadModule(modulePath: string): Promise<{ mod: EmscriptenModule; log: string[] }> {
  const log: string[] = [];
  const factory = (await import(modulePath)).default as Factory;
  const mod = await factory({ print: (s: string) => log.push(s), printErr: (s: string) => log.push(s) });
  return { mod, log };
}

/**
 * Build a minimal, valid little-endian NIfTI-1 (.nii) volume of the given cubic
 * size, filled by `fill(index)` (defaults to the voxel's linear index).
 */
export function makeNifti(side = 8, fill: (i: number) => number = (i) => i): Uint8Array {
  const nvox = side * side * side;
  const buf = new ArrayBuffer(352 + nvox * 4);
  const dv = new DataView(buf);
  dv.setInt32(0, 348, true); // sizeof_hdr
  const dim = [3, side, side, side, 1, 1, 1, 1];
  for (let i = 0; i < 8; i++) dv.setInt16(40 + i * 2, dim[i], true);
  dv.setInt16(70, 16, true); // datatype = DT_FLOAT32
  dv.setInt16(72, 32, true); // bitpix
  for (let i = 0; i < 8; i++) dv.setFloat32(76 + i * 4, 1, true); // pixdim = 1mm iso
  dv.setFloat32(108, 352, true); // vox_offset
  // sform: an affine is required by some ops (e.g. -spm_coreg). Identity rotation,
  // 1mm spacing on the diagonal (srow_x/y/z at offsets 280/296/312).
  dv.setInt16(254, 1, true); // sform_code = NIFTI_XFORM_SCANNER_ANAT
  dv.setFloat32(280, 1, true); // srow_x[0]
  dv.setFloat32(296 + 4, 1, true); // srow_y[1]
  dv.setFloat32(312 + 8, 1, true); // srow_z[2]
  const magic = 'n+1\0';
  for (let i = 0; i < 4; i++) dv.setUint8(344 + i, magic.charCodeAt(i));
  const f = new Float32Array(buf, 352, nvox);
  for (let i = 0; i < nvox; i++) f[i] = fill(i);
  return new Uint8Array(buf);
}

export interface RunResult {
  exitCode: number;
  log: string[];
  /** Output bytes, transparently gunzipped if niimath wrote a .gz (FSL default). */
  output: Uint8Array | null;
}

/**
 * Stage `inputs` into MEMFS, run `args`, and read back `outName` (trying the
 * literal name and the .gz variant niimath writes by default). Always cleans up.
 */
export async function run(
  mod: EmscriptenModule,
  log: string[],
  inputs: { name: string; data: Uint8Array }[],
  args: string[],
  outName: string,
): Promise<RunResult> {
  log.length = 0;
  const staged: string[] = [];
  try {
    for (const f of inputs) {
      mod.FS_createDataFile('.', f.name, f.data, true, true);
      staged.push(f.name);
    }
    const exitCode = mod.callMain(args);
    let output: Uint8Array | null = null;
    for (const candidate of [outName, `${outName}.gz`]) {
      try {
        const raw = mod.FS_readFile(candidate);
        output = candidate.endsWith('.gz') ? new Uint8Array(gunzipSync(raw)) : new Uint8Array(raw);
        staged.push(candidate);
        break;
      } catch {
        /* not this name */
      }
    }
    return { exitCode, log: [...log], output };
  } finally {
    for (const name of staged) {
      try { mod.FS_unlink(name); } catch { /* already gone */ }
    }
  }
}

/** Read float32 voxel `i` from a decoded NIfTI-1 buffer (352-byte header). */
export function voxel(nii: Uint8Array, i: number): number {
  const dv = new DataView(nii.buffer, nii.byteOffset, nii.byteLength);
  return dv.getFloat32(352 + i * 4, true);
}
