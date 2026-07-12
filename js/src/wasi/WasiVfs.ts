/* WasiVfs.ts — a narrow, host-backed, in-memory VFS for the WASI-C niimath reactor.
 *
 * This is intentionally NOT Emscripten MEMFS. It is a flat relative namespace with exactly the
 * open/create/truncate, read/write/seek/tell, close, and stat semantics that the computational
 * feature suite exercises, plus stdout/stderr capture. One session byte-limit and file-count
 * limit live here (the quota is enforced once, at the host boundary — C never re-checks). The same
 * byte-limit also bounds captured stdout/stderr (excess is dropped, not grown unbounded) and the
 * reachable seek range. It bounds retained VFS memory; caller-owned inputs and returned copies are
 * outside this accounting.
 *
 * The VFS is the source of truth for both wasi-libc stdio (via WasiImports) and reactor staging.
 * A leading "./" is normalized away; absolute paths, "..", NUL, and directory creation are
 * rejected at the boundary.
 */

export const PREOPEN_FD = 3; // the single preopened "." directory
export const FIRST_FILE_FD = 4;

export interface VfsFile {
  name: string;
  data: Uint8Array; // logical bytes [0, size)
  size: number;
}

type FdEntry =
  | { kind: "dir"; name: string }
  | { kind: "stdout" }
  | { kind: "stderr" }
  | { kind: "stdin" }
  | { kind: "file"; name: string; offset: number; append: boolean };

export interface VfsLimits {
  maxBytes: number; // total logical bytes across all files plus captured stdout/stderr
  maxFiles: number;
}

const DEFAULT_LIMITS: VfsLimits = { maxBytes: 512 * 1024 * 1024, maxFiles: 256 };

export class QuotaError extends Error {}

function normalizeName(raw: string): string {
  if (raw.length === 0) throw new QuotaError("empty path");
  if (raw.indexOf("\0") >= 0) throw new QuotaError("NUL in path");
  let p = raw;
  while (p.startsWith("./")) p = p.slice(2);
  if (p === "." || p === "") throw new QuotaError("path resolves to directory");
  if (p.startsWith("/")) throw new QuotaError(`absolute path rejected: ${raw}`);
  if (p.split("/").some((seg) => seg === "..")) throw new QuotaError(`'..' rejected: ${raw}`);
  return p;
}

export class WasiVfs {
  private files = new Map<string, VfsFile>();
  private fds = new Map<number, FdEntry>();
  private nextFd = FIRST_FILE_FD;
  private stdoutChunks: Uint8Array[] = [];
  private stderrChunks: Uint8Array[] = [];
  private stdBytes = 0; // total captured stdout+stderr bytes this session (bounded by maxBytes)
  readonly limits: VfsLimits;

  constructor(limits: Partial<VfsLimits> = {}) {
    const merged = { ...DEFAULT_LIMITS, ...limits };
    if (!Number.isSafeInteger(merged.maxBytes) || merged.maxBytes < 0)
      throw new RangeError("maxBytes must be a non-negative safe integer");
    if (!Number.isSafeInteger(merged.maxFiles) || merged.maxFiles < 0)
      throw new RangeError("maxFiles must be a non-negative safe integer");
    this.limits = Object.freeze(merged);
    this.installStd();
  }

  private installStd() {
    this.fds.set(0, { kind: "stdin" });
    this.fds.set(1, { kind: "stdout" });
    this.fds.set(2, { kind: "stderr" });
    this.fds.set(PREOPEN_FD, { kind: "dir", name: "." });
  }

  /** Drop all files, open fds, and captured output; reinstall std fds + preopen. */
  reset() {
    this.files.clear();
    this.fds.clear();
    this.nextFd = FIRST_FILE_FD;
    this.stdoutChunks = [];
    this.stderrChunks = [];
    this.stdBytes = 0;
    this.installStd();
  }

  // ---- host-side file table API (used by the runner/adapter, not by WASM directly) ----

  private totalBytes(exclude?: string): number {
    let t = 0;
    for (const [n, f] of this.files) if (n !== exclude) t += f.size;
    return t;
  }

  /** Stage a file. By default `data` is defensively copied (the caller may still own/mutate it).
   *  Pass `owned=true` to TRANSFER ownership: the buffer is stored as-is with no copy — used by the
   *  runner for buffers it just created and will not touch again (e.g. gunzip output), which avoids
   *  a redundant full copy of the decompressed input at the host boundary. */
  addFile(name: string, data: Uint8Array, owned = false): void {
    const n = normalizeName(name);
    const existing = this.files.get(n);
    const projected = this.totalBytes(n) + this.stdBytes + data.byteLength;
    if (projected > this.limits.maxBytes)
      throw new QuotaError(`byte quota exceeded (${projected} > ${this.limits.maxBytes})`);
    if (!existing && this.files.size >= this.limits.maxFiles)
      throw new QuotaError(`file-count quota exceeded (${this.limits.maxFiles})`);
    this.files.set(n, { name: n, data: owned ? data : data.slice(), size: data.byteLength });
  }

  readFile(name: string): Uint8Array | null {
    const f = this.files.get(normalizeName(name));
    return f ? f.data.subarray(0, f.size) : null;
  }

  hasFile(name: string): boolean {
    try {
      return this.files.has(normalizeName(name));
    } catch {
      return false;
    }
  }

  listFiles(): string[] {
    return [...this.files.keys()].sort();
  }

  takeStdout(): Uint8Array {
    return concat(this.stdoutChunks);
  }
  takeStderr(): Uint8Array {
    return concat(this.stderrChunks);
  }

  // ---- syscall-level operations invoked by WasiImports (errno-style, 0 = ok) ----

  /** path_open. Returns {fd} or {errno}. */
  open(rawName: string, oflags: number): { fd?: number; errno?: number } {
    let name: string;
    try {
      name = normalizeName(rawName);
    } catch {
      return { errno: ERRNO.INVAL };
    }
    const OFLAG_CREAT = 1,
      OFLAG_EXCL = 4,
      OFLAG_TRUNC = 8;
    let f = this.files.get(name);
    if (!f) {
      if (!(oflags & OFLAG_CREAT)) return { errno: ERRNO.NOENT };
      if (this.files.size >= this.limits.maxFiles) return { errno: ERRNO.NOSPC };
      f = { name, data: new Uint8Array(0), size: 0 };
      this.files.set(name, f);
    } else if (oflags & OFLAG_EXCL) {
      return { errno: ERRNO.EXIST };
    } else if (oflags & OFLAG_TRUNC) {
      f.data = new Uint8Array(0);
      f.size = 0;
    }
    const fd = this.nextFd++;
    this.fds.set(fd, { kind: "file", name, offset: 0, append: false });
    return { fd };
  }

  close(fd: number): number {
    if (fd <= PREOPEN_FD) return ERRNO.SUCCESS; // std/preopen: no-op close
    if (!this.fds.has(fd)) return ERRNO.BADF;
    this.fds.delete(fd);
    return ERRNO.SUCCESS;
  }

  private ensureCapacity(f: VfsFile, needed: number): boolean {
    if (needed <= f.data.byteLength) return true;
    const over = this.totalBytes(f.name) + this.stdBytes + needed;
    if (over > this.limits.maxBytes) return false;
    // Grow exactly to the requested logical extent. Geometric growth can allocate beyond the
    // configured quota even when `needed` itself fits, defeating the VFS memory boundary.
    const grown = new Uint8Array(needed);
    grown.set(f.data.subarray(0, f.size));
    f.data = grown;
    return true;
  }

  /** fd_write for a file fd. Returns bytes written or -errno. stdout/stderr handled separately. */
  writeFile(fd: number, chunks: Uint8Array[]): number {
    const e = this.fds.get(fd);
    if (!e || e.kind !== "file") return -ERRNO.BADF;
    const f = this.files.get(e.name)!;
    let written = 0;
    for (const c of chunks) {
      const end = e.offset + c.byteLength;
      if (!this.ensureCapacity(f, end)) return written > 0 ? written : -ERRNO.NOSPC;
      f.data.set(c, e.offset);
      e.offset = end;
      if (end > f.size) f.size = end;
      written += c.byteLength;
    }
    return written;
  }

  captureStd(fd: number, chunk: Uint8Array): void {
    // Bound total captured output by the session byte-limit; drop the excess (WASM is still told
    // the whole write succeeded, like /dev/null) so a runaway log can't grow host memory unbounded.
    const available = this.limits.maxBytes - this.totalBytes() - this.stdBytes;
    if (available <= 0) return;
    const c = chunk.byteLength > available ? chunk.subarray(0, available) : chunk;
    this.stdBytes += c.byteLength;
    (fd === 2 ? this.stderrChunks : this.stdoutChunks).push(c.slice());
  }

  /** fd_read. Fills iov targets (as arrays) from file at offset; returns bytes read or -errno. */
  readInto(fd: number, targets: Uint8Array[]): number {
    const e = this.fds.get(fd);
    if (!e) return -ERRNO.BADF;
    if (e.kind === "stdin") return 0; // no stdin
    if (e.kind !== "file") return -ERRNO.BADF;
    const f = this.files.get(e.name)!;
    let read = 0;
    for (const t of targets) {
      if (e.offset >= f.size) break;
      const n = Math.min(t.byteLength, f.size - e.offset);
      t.set(f.data.subarray(e.offset, e.offset + n));
      e.offset += n;
      read += n;
      if (n < t.byteLength) break;
    }
    return read;
  }

  /** fd_seek. whence 0=SET,1=CUR,2=END. Returns {offset} or {errno}. Rejects a non-integer delta
   *  and any resulting offset outside [0, maxBytes] with EINVAL, so a crafted i64 offset yields a
   *  clean WASI errno instead of an imprecise offset or a later throw. */
  seek(fd: number, delta: number, whence: number): { offset?: number; errno?: number } {
    const e = this.fds.get(fd);
    if (!e || e.kind !== "file") return { errno: ERRNO.BADF };
    if (!Number.isSafeInteger(delta)) return { errno: ERRNO.INVAL };
    const f = this.files.get(e.name)!;
    let base = 0;
    if (whence === 0) base = 0;
    else if (whence === 1) base = e.offset;
    else if (whence === 2) base = f.size;
    else return { errno: ERRNO.INVAL };
    const off = base + delta;
    if (off < 0 || off > this.limits.maxBytes || !Number.isSafeInteger(off))
      return { errno: ERRNO.INVAL };
    e.offset = off;
    return { offset: off };
  }

  fdIsFile(fd: number): boolean {
    return this.fds.get(fd)?.kind === "file";
  }
  fdIsDir(fd: number): boolean {
    return this.fds.get(fd)?.kind === "dir";
  }

  /** path_filestat_get / stat by name. Returns {size, filetype} or {errno}. */
  statByName(rawName: string): { size?: number; filetype?: number; errno?: number } {
    let name: string;
    try {
      name = normalizeName(rawName);
    } catch {
      return { errno: ERRNO.INVAL };
    }
    const f = this.files.get(name);
    if (!f) return { errno: ERRNO.NOENT };
    return { size: f.size, filetype: FILETYPE.REGULAR_FILE };
  }
}

function concat(chunks: Uint8Array[]): Uint8Array {
  let n = 0;
  for (const c of chunks) n += c.byteLength;
  const out = new Uint8Array(n);
  let o = 0;
  for (const c of chunks) {
    out.set(c, o);
    o += c.byteLength;
  }
  return out;
}

export const ERRNO = {
  SUCCESS: 0,
  BADF: 8,
  EXIST: 20,
  INVAL: 28,
  IO: 29,
  ISDIR: 31,
  NOENT: 44,
  NOSPC: 51,
  NOSYS: 52,
  NOTDIR: 54,
} as const;

export const FILETYPE = {
  UNKNOWN: 0,
  CHARACTER_DEVICE: 2,
  DIRECTORY: 3,
  REGULAR_FILE: 4,
} as const;
