# WASI-C niimath backend

A lean, zlib-free, single-threaded WASI reactor of niimath's **computational BSD** feature set
(core math, multi-image ops, allineate/deface, dtifit, QC, conform, butterworth). It is the
experimental WASI-C backend; the shipped default backend remains the Emscripten build
(`../niimath.js`/`.wasm`).

## What it is (and isn't)

- **Is:** the whole scoped compute engine compiled with `zig cc -target wasm32-wasi
  -mexec-model=reactor` (Zig bundles wasi-libc — no wasi-sdk install needed). File I/O flows
  through real wasi-libc stdio, serviced by a small host-side in-memory VFS. NIfTI gzip is done
  host-side (`CompressionStream`/`node:zlib`), so the module only sees raw `.nii` bytes.
- **Isn't:** a drop-in replacement yet. **Bitmap (`-bitmap`) and mesh (`-mesh`) are not compiled**
  — they fail clearly (nonzero exit), never silently no-op. No GPL, zstd, or OpenMP.

## Files

| file | role |
|---|---|
| `../../../src/wasi_shim.c` | reactor entry (`nii_run`) — the only C added; no algorithm changes |
| `WasiVfs.ts` | narrow host-side in-memory VFS (flat relative namespace, quota, stdout/stderr capture) |
| `WasiImports.ts` | the 15 frozen WASI Preview-1 imports over the VFS; `proc_exit` → `WasiExit` |
| `WasiRunner.ts` | reactor lifecycle, argv marshaling, instance recycling, `runFiles`/`runFilesGz` |
| `gzip.ts` | host-side magic-byte gzip (CompressionStream, falls back to `node:zlib` in Bun) |
| `featureParity.test.ts`, `compression.test.ts` | Bun conformance, lifecycle, quota, and gzip tests |

## Build

```bash
make -C ../../../src wasm-wasi        # -> ../niimath-wasi.wasm (+ import-surface check)
make -C ../../../src wasm-emcc-core   # feature-matched Emscripten comparison artifact

cd ../..
bun run test:wasi                     # fixtures + Bun conformance/lifecycle tests
bun run wasiNodeSmoke                 # plain-Node stream/gzip/quota smoke
bun run wasiChromium                  # optional real-browser correctness/agreement check
```

`wasm-wasi` verifies the toolchain (Zig ≥ 0.16), builds the reactor `-O3 -flto -ffast-math
--strip-all`, and fails if the WASI import surface grows beyond
[`import-manifest.json`](./import-manifest.json).

## Use

```ts
import { WasiRunner } from "./WasiRunner";
const runner = await WasiRunner.create(wasmBytes);

// low-level multi-file API
const r = await runner.runFiles({
  argv: ["in.nii", "-add", "1", "out.nii"],
  inputs: { "in.nii": rawNiftiBytes },
  outputs: ["out.nii"],
});
// r.exitCode, r.files["out.nii"], r.stdout, r.stderr

// gzip-aware: gzip inputs decompressed + argv rewritten; .gz outputs recompressed
const g = await runner.runFilesGz({
  argv: ["in.nii.gz", "-s", "3", "out.nii.gz"],
  inputs: { "in.nii.gz": gzippedBytes },
  outputs: ["out.nii.gz"],
});
```

`create()` calls `_initialize` once. A normal run reuses the instance after `reset()`; a
`proc_exit` (error path) recreates it automatically. Multi-output commands (dtifit → 11 files)
list every output name in `outputs`. One runner is single-flight: await a run before calling another
run, `reset()`, or `addFile()`, and use one runner per concurrent job. Returned file buffers are
detached copies and may be mutated by the caller. `runFilesGz()` clears its raw staging state before
compression/return; retain its returned buffers rather than expecting `readFile()` to expose them.

## Lifecycle contract

1. `_initialize` exactly once per instance.
2. Normal `nii_run` return → reuse after `reset()` (clears the VFS).
3. `proc_exit` → `WasiExit` thrown, exit code captured, instance discarded and recreated (C
   cleanup was skipped, so linear memory / libc state is untrusted).

## Status

Functionally complete and verified: feature parity + payload/header conformance vs the native
binary pass (`bun run test:wasi`). It is a correct, portable, glue-free backend kept for
audit/portability. Performance evidence is intentionally deferred until the current source is
committed; the benchmark lives in `niimath_tests/wasm_benchmark` and is not a release-quality claim.
