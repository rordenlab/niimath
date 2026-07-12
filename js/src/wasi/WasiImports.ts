/* WasiImports.ts — the 15 WASI Preview 1 imports the niimath reactor requests, serviced from a
 * host-side WasiVfs. This is the frozen import surface (see ./import-manifest.json);
 * CI fails if the module asks for anything not implemented here.
 *
 * Memory growth safety: ALLOW_MEMORY_GROWTH detaches the ArrayBuffer on grow, so every import
 * re-reads `memory.buffer` through fresh DataView/typed-array views — no cached views.
 */
import { WasiVfs, ERRNO, FILETYPE, PREOPEN_FD } from "./WasiVfs";

const MAX_SAFE_BIGINT = BigInt(Number.MAX_SAFE_INTEGER);

/** Thrown by proc_exit; the runner converts it to a normal result and discards the instance. */
export class WasiExit extends Error {
  constructor(public readonly code: number) {
    super(`wasi proc_exit(${code})`);
  }
}

export interface WasiHostOptions {
  env?: Record<string, string>;
  /** deterministic clock start (ns); increments by a fixed step per call for reproducibility */
  clockStartNs?: bigint;
}

export function makeWasiImports(
  vfs: WasiVfs,
  getMemory: () => WebAssembly.Memory,
  opts: WasiHostOptions = {},
) {
  const enc = new TextEncoder();
  const dec = new TextDecoder();
  const env = opts.env ?? {};
  const envPairs = Object.entries(env).map(([k, v]) => `${k}=${v}`);
  let clockNs = opts.clockStartNs ?? 0n;

  const dv = () => new DataView(getMemory().buffer);
  const u8 = () => new Uint8Array(getMemory().buffer);
  const readStr = (ptr: number, len: number) => dec.decode(u8().subarray(ptr, ptr + len));

  // Read an iovec array [(bufPtr u32, bufLen u32) * count] into a list of subarray views.
  function iovs(ptr: number, count: number): Uint8Array[] {
    const view = dv();
    const mem = u8();
    const out: Uint8Array[] = [];
    for (let i = 0; i < count; i++) {
      const p = view.getUint32(ptr + i * 8, true);
      const l = view.getUint32(ptr + i * 8 + 4, true);
      out.push(mem.subarray(p, p + l));
    }
    return out;
  }

  return {
    proc_exit(code: number): never {
      throw new WasiExit(code >>> 0 ? code : code);
    },

    fd_write(fd: number, iovsPtr: number, iovsLen: number, nwrittenPtr: number): number {
      const chunks = iovs(iovsPtr, iovsLen);
      if (fd === 1 || fd === 2) {
        let total = 0;
        for (const c of chunks) {
          vfs.captureStd(fd, c);
          total += c.byteLength;
        }
        dv().setUint32(nwrittenPtr, total, true);
        return ERRNO.SUCCESS;
      }
      // file write: pass the linear-memory subarrays straight to the VFS, which copies them into
      // its file buffer immediately (no memory growth occurs within a single fd_write). Slicing
      // first would be a redundant second copy of the whole payload — a real cost on large writes.
      const n = vfs.writeFile(fd, chunks);
      if (n < 0) return -n;
      dv().setUint32(nwrittenPtr, n, true);
      return ERRNO.SUCCESS;
    },

    fd_read(fd: number, iovsPtr: number, iovsLen: number, nreadPtr: number): number {
      const targets = iovs(iovsPtr, iovsLen);
      const n = vfs.readInto(fd, targets);
      if (n < 0) return -n;
      dv().setUint32(nreadPtr, n, true);
      return ERRNO.SUCCESS;
    },

    fd_seek(fd: number, offset: bigint, whence: number, newOffsetPtr: number): number {
      // Reject an i64 delta that can't round-trip as an exact JS integer before it reaches the VFS,
      // so a crafted offset yields EINVAL rather than a silently-narrowed (imprecise) seek.
      if (offset > MAX_SAFE_BIGINT || offset < -MAX_SAFE_BIGINT) return ERRNO.INVAL;
      const r = vfs.seek(fd, Number(offset), whence);
      if (r.errno !== undefined) return r.errno;
      dv().setBigUint64(newOffsetPtr, BigInt(r.offset!), true);
      return ERRNO.SUCCESS;
    },

    fd_close(fd: number): number {
      return vfs.close(fd);
    },

    fd_fdstat_get(fd: number, retPtr: number): number {
      const view = dv();
      let filetype: number = FILETYPE.UNKNOWN;
      if (fd === 0 || fd === 1 || fd === 2) filetype = FILETYPE.CHARACTER_DEVICE;
      else if (vfs.fdIsDir(fd)) filetype = FILETYPE.DIRECTORY;
      else if (vfs.fdIsFile(fd)) filetype = FILETYPE.REGULAR_FILE;
      else return ERRNO.BADF;
      // fdstat: fs_filetype u8@0, fs_flags u16@2, fs_rights_base u64@8, fs_rights_inheriting u64@16
      view.setUint8(retPtr, filetype);
      view.setUint16(retPtr + 2, 0, true);
      view.setBigUint64(retPtr + 8, 0xffffffffffffffffn, true);
      view.setBigUint64(retPtr + 16, 0xffffffffffffffffn, true);
      return ERRNO.SUCCESS;
    },

    fd_fdstat_set_flags(_fd: number, _flags: number): number {
      return ERRNO.SUCCESS;
    },

    fd_prestat_get(fd: number, retPtr: number): number {
      if (fd !== PREOPEN_FD) return ERRNO.BADF;
      // prestat: tag u8@0 (0=dir), pr_name_len u32@4
      const view = dv();
      view.setUint8(retPtr, 0);
      view.setUint32(retPtr + 4, 1, true); // "." length
      return ERRNO.SUCCESS;
    },

    fd_prestat_dir_name(fd: number, pathPtr: number, pathLen: number): number {
      if (fd !== PREOPEN_FD) return ERRNO.BADF;
      const bytes = enc.encode(".");
      if (pathLen < bytes.length) return ERRNO.INVAL;
      u8().set(bytes, pathPtr);
      return ERRNO.SUCCESS;
    },

    path_open(
      dirFd: number,
      _dirFlags: number,
      pathPtr: number,
      pathLen: number,
      oflags: number,
      _rightsBase: bigint,
      _rightsInheriting: bigint,
      _fdFlags: number,
      openedFdPtr: number,
    ): number {
      if (dirFd !== PREOPEN_FD) return ERRNO.BADF;
      const name = readStr(pathPtr, pathLen);
      const r = vfs.open(name, oflags);
      if (r.errno !== undefined) return r.errno;
      dv().setUint32(openedFdPtr, r.fd!, true);
      return ERRNO.SUCCESS;
    },

    path_filestat_get(
      dirFd: number,
      _flags: number,
      pathPtr: number,
      pathLen: number,
      retPtr: number,
    ): number {
      if (dirFd !== PREOPEN_FD) return ERRNO.BADF;
      const name = readStr(pathPtr, pathLen);
      const r = vfs.statByName(name);
      if (r.errno !== undefined) return r.errno;
      writeFilestat(dv(), retPtr, r.size!, r.filetype!);
      return ERRNO.SUCCESS;
    },

    clock_time_get(_id: number, _precision: bigint, retPtr: number): number {
      clockNs += 1_000_000n; // deterministic monotonic step (1ms); no host entropy/time
      dv().setBigUint64(retPtr, clockNs, true);
      return ERRNO.SUCCESS;
    },

    environ_sizes_get(countPtr: number, bufSizePtr: number): number {
      const view = dv();
      let bufSize = 0;
      for (const p of envPairs) bufSize += enc.encode(p).length + 1;
      view.setUint32(countPtr, envPairs.length, true);
      view.setUint32(bufSizePtr, bufSize, true);
      return ERRNO.SUCCESS;
    },

    environ_get(environPtr: number, bufPtr: number): number {
      const view = dv();
      const mem = u8();
      let ptr = bufPtr;
      let arr = environPtr;
      for (const p of envPairs) {
        view.setUint32(arr, ptr, true);
        arr += 4;
        const bytes = enc.encode(p);
        mem.set(bytes, ptr);
        mem[ptr + bytes.length] = 0;
        ptr += bytes.length + 1;
      }
      return ERRNO.SUCCESS;
    },

    poll_oneoff(_in: number, _out: number, _nsubs: number, neventsPtr: number): number {
      dv().setUint32(neventsPtr, 0, true);
      return ERRNO.SUCCESS;
    },
  };
}

function writeFilestat(view: DataView, ptr: number, size: number, filetype: number) {
  // filestat: dev u64@0, ino u64@8, filetype u8@16, nlink u64@24, size u64@32, atim@40, mtim@48, ctim@56
  view.setBigUint64(ptr, 0n, true);
  view.setBigUint64(ptr + 8, 0n, true);
  view.setUint8(ptr + 16, filetype);
  view.setBigUint64(ptr + 24, 1n, true);
  view.setBigUint64(ptr + 32, BigInt(size), true);
  view.setBigUint64(ptr + 40, 0n, true);
  view.setBigUint64(ptr + 48, 0n, true);
  view.setBigUint64(ptr + 56, 0n, true);
}
