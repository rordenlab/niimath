/* WasiRunner.ts — drives the zlib-free WASI-C niimath reactor.
 *
 * Lifecycle contract:
 *   - call _initialize exactly once per instance, before any allocator/reactor export;
 *   - a normal nii_run return reuses the instance after vfs.reset();
 *   - a proc_exit (error path) throws WasiExit -> we capture the code, then DISCARD and recreate
 *     the instance, because C cleanup was skipped and linear memory / libc state may be stale.
 *
 * Exposes the shared backend-adapter contract (reset/addFile/run/readFile/listFiles) plus the
 * low-level runFiles({argv, inputs, outputs}) API the public wrapper builds on. The VFS is the
 * host-side source of truth: inputs are staged directly into it and outputs read directly from
 * it (the production zero-copy path — no host->WASM->host result copy).
 */
import { QuotaError, WasiVfs } from "./WasiVfs";
import { makeWasiImports, WasiExit, WasiHostOptions } from "./WasiImports";
import { isGzip, gunzip, gzip } from "./gzip";

export interface RunResult {
  exitCode: number;
  stdout: string;
  stderr: string;
}

export interface RunFilesArgs {
  argv: string[]; // WITHOUT argv[0]; "niimath" is prepended
  inputs?: Record<string, Uint8Array>;
  outputs?: string[]; // names to extract after the run
}

export interface RunFilesResult extends RunResult {
  files: Record<string, Uint8Array>;
}

interface ReactorExports {
  memory: WebAssembly.Memory;
  _initialize: () => void;
  nii_run: (argvPtr: number, argvLen: number) => number;
  malloc: (len: number) => number;
  free: (ptr: number) => void;
}

export class WasiRunner {
  private module: WebAssembly.Module;
  private vfs: WasiVfs;
  private hostOpts: WasiHostOptions;
  private instance: WebAssembly.Instance | null = null;
  private exports: ReactorExports | null = null;
  private dead = false;

  private constructor(module: WebAssembly.Module, hostOpts: WasiHostOptions, vfs: WasiVfs) {
    this.module = module;
    this.hostOpts = hostOpts;
    this.vfs = vfs;
  }

  static async create(
    wasm: BufferSource | WebAssembly.Module,
    hostOpts: WasiHostOptions = { env: { FSLOUTPUTTYPE: "NIFTI" } },
    vfsLimits?: { maxBytes?: number; maxFiles?: number },
  ): Promise<WasiRunner> {
    const module =
      wasm instanceof WebAssembly.Module ? wasm : await WebAssembly.compile(wasm);
    const vfs = new WasiVfs(vfsLimits ?? {});
    const runner = new WasiRunner(module, hostOpts, vfs);
    await runner.instantiate();
    return runner;
  }

  private async instantiate(): Promise<void> {
    const imports = makeWasiImports(this.vfs, () => this.exports!.memory, this.hostOpts);
    const instance = await WebAssembly.instantiate(this.module, {
      wasi_snapshot_preview1: imports,
    });
    this.instance = instance;
    this.exports = instance.exports as unknown as ReactorExports;
    this.exports._initialize(); // exactly once, before any other export
    this.dead = false;
  }

  private async ensureLive(): Promise<void> {
    if (this.dead || !this.instance) await this.instantiate();
  }

  // ---- shared backend-adapter contract ----

  /** Clear the VFS. If the last run poisoned the instance (proc_exit), recreate it. Single-flight:
   *  rejects if a run is in progress — an external reset mid-run would clear its staged inputs, and
   *  two concurrent resets would race instance recreation. */
  reset(): Promise<void> {
    return this.lock("reset", () => this._reset());
  }
  private async _reset(): Promise<void> {
    this.vfs.reset();
    if (this.dead) await this.instantiate();
  }

  /** Stage a file. Throws if a run is in progress — staging mid-run corrupts its VFS. */
  addFile(name: string, data: Uint8Array): void {
    if (this.inFlight)
      throw new Error("WasiRunner is single-flight: cannot addFile() while a run is in progress.");
    this.vfs.addFile(name, data);
  }

  readFile(name: string): Uint8Array | null {
    const data = this.vfs.readFile(name);
    return data ? data.slice() : null;
  }

  listFiles(): string[] {
    return this.vfs.listFiles();
  }

  // ---- single-flight guard ----
  // This runner wraps ONE stateful reactor instance + VFS, so two overlapping runs would interleave
  // (reset() clears the VFS mid-run; linear memory / libc state is shared). The public run methods
  // are therefore single-flight: a call made while another is in progress throws rather than
  // corrupting state — use one runner per concurrent job. Internal calls use the un-guarded
  // _-prefixed variants so the nesting (runFilesGz -> runFiles -> run) cannot self-deadlock.
  private inFlight = false;
  private async lock<T>(label: string, fn: () => Promise<T>): Promise<T> {
    if (this.inFlight)
      throw new Error(
        `WasiRunner is single-flight: a run is already in progress — await it before calling ${label}(). Use one runner per concurrent job.`,
      );
    this.inFlight = true;
    try {
      return await fn();
    } finally {
      this.inFlight = false;
    }
  }

  /** Run niimath with argv (argv[0] "niimath" is prepended). Files must already be staged. */
  run(argv: string[]): Promise<RunResult> {
    return this.lock("run", () => this._run(argv));
  }
  private async _run(argv: string[]): Promise<RunResult> {
    await this.ensureLive();
    const ex = this.exports!;
    const packed = packArgv(["niimath", ...argv]);
    const ptr = ex.malloc(packed.length);
    if (!ptr) throw new Error("malloc failed");
    let exitCode: number;
    try {
      new Uint8Array(ex.memory.buffer).set(packed, ptr);
      exitCode = ex.nii_run(ptr, packed.length);
      ex.free(ptr);
    } catch (e) {
      if (e instanceof WasiExit) {
        // Error path took exit(): instance state is now untrustworthy.
        exitCode = e.code | 0;
        this.dead = true;
      } else {
        this.dead = true;
        throw e;
      }
    }
    const stdout = new TextDecoder().decode(this.vfs.takeStdout());
    const stderr = new TextDecoder().decode(this.vfs.takeStderr());
    return { exitCode, stdout, stderr };
  }

  /** Low-level multi-file API: reset, stage inputs, run, then extract output copies. */
  runFiles(args: RunFilesArgs): Promise<RunFilesResult> {
    return this.lock("runFiles", () => this._runFiles(args));
  }
  private async _runFiles(args: RunFilesArgs, owned?: Set<string>): Promise<RunFilesResult> {
    await this._reset();
    if (args.inputs) {
      for (const [name, data] of Object.entries(args.inputs)) {
        // owned names transfer ownership (no defensive copy) — see _runFilesGz for gunzip buffers.
        this.vfs.addFile(name, data, owned?.has(name)); // throws QuotaError on quota exceed
      }
    }
    const res = await this._run(args.argv);
    const files: Record<string, Uint8Array> = {};
    for (const name of args.outputs ?? []) {
      const data = this.vfs.readFile(name);
      if (data) files[name] = data.slice();
    }
    return { ...res, files };
  }

  /** gzip-aware runFiles: gzip-framed inputs are decompressed and staged under raw names (argv
   * tokens rewritten); outputs whose requested name ends in .gz are gzipped on the way out. The
   * module itself only ever sees raw .nii bytes (it has no zlib). Detection is by magic bytes.
   *
   * Staged names are allocated UNIQUELY across every VFS name that will exist during the run —
   * fixed names (non-gz inputs and non-gz outputs) plus every original argv path — so a second
   * decompressed input or rewritten .gz output can never silently clobber a first, and two argv
   * operands can never resolve to one file. A gz-derived name that would collide is renamed (the
   * matching argv token is rewritten with it); if the namespace is somehow exhausted we throw a
   * clear error rather than overwrite. Decompression is bounded by the VFS byte quota (gzip-bomb
   * safety) — see gzip.ts. */
  runFilesGz(args: RunFilesArgs): Promise<RunFilesResult> {
    return this.lock("runFilesGz", () => this._runFilesGz(args));
  }
  private async _runFilesGz(args: RunFilesArgs): Promise<RunFilesResult> {
    // Release the preceding session before allocating decompression buffers. _runFiles resets
    // again immediately before staging, but doing it here prevents old VFS data and a new gunzip
    // result from occupying the quota simultaneously.
    await this._reset();
    const cap = this.vfs.limits.maxBytes;
    let remaining = cap;
    let argv = [...args.argv];
    const rawInputs: Record<string, Uint8Array> = {};
    const stagedInputs = new Map<string, string>();
    const ownedRaw = new Set<string>(); // gunzip buffers WE created — stage by ownership transfer (no re-copy)
    // Reserve every name that is fixed (cannot be renamed): original argv paths, non-gz input
    // names, and non-gz output names. gz-derived staged names are then allocated around these.
    const reserved = new Set<string>(argv);
    for (const [name, data] of Object.entries(args.inputs ?? {})) {
      if (!isGzip(data)) reserved.add(name);
    }
    for (const name of args.outputs ?? []) {
      if (!name.endsWith(".gz")) reserved.add(name);
    }

    // Non-gz inputs keep their literal name (argv already references it).
    for (const [name, data] of Object.entries(args.inputs ?? {})) {
      if (isGzip(data)) continue;
      if (data.byteLength > remaining)
        throw new QuotaError(`byte quota exceeded (${cap - remaining + data.byteLength} > ${cap})`);
      remaining -= data.byteLength;
      rawInputs[name] = data;
      stagedInputs.set(name, name);
    }
    // gz inputs: decompress (bounded by the byte quota) and stage under a unique raw name,
    // rewriting the matching argv token(s).
    for (const [name, data] of Object.entries(args.inputs ?? {})) {
      if (!isGzip(data)) continue;
      const rawName = uniqueStagedName(stagedRawName(name), reserved);
      reserved.add(rawName);
      const raw = await gunzip(data, remaining);
      remaining -= raw.byteLength;
      rawInputs[rawName] = raw;
      ownedRaw.add(rawName);
      stagedInputs.set(name, rawName);
      argv = argv.map((tok) => (tok === name ? rawName : tok));
    }

    // For .gz outputs, ask niimath for a unique raw name, then compress after the run.
    const seenOutputs = new Set<string>();
    const outSpec = (args.outputs ?? []).map((name) => {
      if (seenOutputs.has(name))
        throw new Error(`WASI staging: duplicate requested output '${name}'`);
      seenOutputs.add(name);
      if (!name.endsWith(".gz")) return { requested: name, raw: name, compress: false };
      // An in-place gzip operation uses the already-staged raw input for its raw output. This
      // preserves the single logical pathname instead of rewriting the input and output tokens to
      // different files (which would either clobber the input or make the requested output vanish).
      const raw = stagedInputs.get(name) ?? uniqueStagedName(name.slice(0, -3), reserved);
      reserved.add(raw);
      return { requested: name, raw, compress: true };
    });
    argv = argv.map((tok) => {
      const hit = outSpec.find((o) => o.compress && o.requested === tok);
      return hit ? hit.raw : tok;
    });
    let res: RunFilesResult;
    try {
      res = await this._runFiles({ argv, inputs: rawInputs, outputs: outSpec.map((o) => o.raw) }, ownedRaw);
    } finally {
      // Drop both references to internally-owned decompression buffers, including on a trap.
      for (const name of Object.keys(rawInputs)) delete rawInputs[name];
      this.vfs.reset();
    }
    const files: Record<string, Uint8Array> = {};
    for (const o of outSpec) {
      const raw = res.files[o.raw];
      if (!raw) continue;
      files[o.requested] = o.compress ? await gzip(raw) : raw;
      delete res.files[o.raw];
    }
    return { exitCode: res.exitCode, stdout: res.stdout, stderr: res.stderr, files };
  }
}

/** Raw name for a staged (decompressed) input: strip a trailing .gz, else add .nii. */
function stagedRawName(name: string): string {
  if (name.endsWith(".gz")) return name.slice(0, -3);
  if (name.endsWith(".nii")) return name;
  return name + ".nii";
}

/** Return `desired` if free, else the first `<stem>-wasiN<ext>` not already reserved (extension
 *  preserved so niimath still infers the format). Throws if the namespace is exhausted. */
function uniqueStagedName(desired: string, reserved: Set<string>): string {
  if (!reserved.has(desired)) return desired;
  const dot = desired.lastIndexOf(".");
  const stem = dot > 0 ? desired.slice(0, dot) : desired;
  const ext = dot > 0 ? desired.slice(dot) : "";
  for (let i = 1; i < 100000; i++) {
    const cand = `${stem}-wasi${i}${ext}`;
    if (!reserved.has(cand)) return cand;
  }
  throw new Error(`WASI staging: could not allocate a unique staged name for '${desired}'`);
}

/** Pack argv into a single NUL-separated, NUL-terminated buffer (matches wasi_shim.c nii_run). */
export function packArgv(argv: string[]): Uint8Array {
  if (argv.length >= 256)
    throw new RangeError("WASI reactor supports at most 255 arguments including argv[0]");
  for (const arg of argv) {
    if (arg.includes("\0")) throw new Error("WASI reactor arguments must not contain NUL bytes");
  }
  const enc = new TextEncoder();
  const parts = argv.map((a) => enc.encode(a));
  let n = 0;
  for (const p of parts) n += p.byteLength + 1;
  const out = new Uint8Array(n);
  let o = 0;
  for (const p of parts) {
    out.set(p, o);
    o += p.byteLength;
    out[o++] = 0;
  }
  return out;
}
