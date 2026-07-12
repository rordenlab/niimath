/* gzip.ts — host-side gzip for the zlib-free WASI-C backend.
 *
 * The module sees only raw .nii bytes; compression happens here. Detection is by magic bytes
 * (0x1f 0x8b), not extension.
 *
 * Runtime portability: browsers and Node expose the WHATWG CompressionStream/DecompressionStream
 * globals; Bun (1.3.x) does not, but implements node:zlib. We prefer the streaming globals when
 * present (the plan's canonical path) and fall back to node:zlib otherwise — so the identical
 * source runs in Bun, Node 20+, and Chromium.
 *
 * Gzip-bomb safety: gunzip() takes an optional `maxBytes` cap (the caller passes the VFS session
 * byte quota). The streaming path fails fast the instant the accumulated output crosses the cap;
 * the node:zlib fallback enforces it via `maxOutputLength`. This bounds the LOGICAL output size.
 *
 * Peak-memory accounting (browser stream path): decompression itself is the residual. To return one
 * contiguous buffer, `viaStream` transiently holds the retained chunks AND the final result during
 * the copy — up to ~2× the logical cap (plus the caller-owned compressed input). The runner keeps it
 * to that: it stages the decompressed buffer into the VFS by OWNERSHIP TRANSFER (no extra copy, see
 * WasiVfs.addFile(owned)), clears prior-session data before decompression, and drops VFS/local input
 * references before compressing detached output copies. A raw output and its gzip result necessarily
 * coexist while that output is compressed. For untrusted large inputs, lower the VFS `maxBytes`
 * accordingly. Removing even the in-`viaStream` doubling requires streaming chunks straight
 * into a chunk-aware VFS file (no contiguous concat) — a productization item, tracked, not yet done.
 */

export function isGzip(bytes: Uint8Array): boolean {
  return bytes.length >= 2 && bytes[0] === 0x1f && bytes[1] === 0x8b;
}

const hasStreams =
  typeof (globalThis as any).CompressionStream !== "undefined" &&
  typeof (globalThis as any).DecompressionStream !== "undefined";

/** Wrap bytes as an ArrayBuffer-backed view for the DOM Blob API (never a SharedArrayBuffer);
 *  zero-copy in the common case, copying only a SharedArrayBuffer-backed view. */
function blobPart(bytes: Uint8Array): Uint8Array<ArrayBuffer> {
  if (bytes.buffer instanceof ArrayBuffer) return bytes as Uint8Array<ArrayBuffer>;
  const copy = new Uint8Array(bytes.byteLength);
  copy.set(bytes);
  return copy;
}

async function viaStream(bytes: Uint8Array, dir: "c" | "d", maxBytes: number): Promise<Uint8Array> {
  const Ctor = dir === "c" ? (globalThis as any).CompressionStream : (globalThis as any).DecompressionStream;
  const stream = new Blob([blobPart(bytes)]).stream().pipeThrough(new Ctor("gzip"));
  const chunks: Uint8Array[] = [];
  let n = 0;
  for await (const chunk of stream as any) {
    const c = chunk as Uint8Array;
    n += c.byteLength;
    if (n > maxBytes)
      throw new Error(`gunzip: decompressed output exceeds cap (${n} > ${maxBytes} bytes)`);
    chunks.push(c);
  }
  const out = new Uint8Array(n);
  let o = 0;
  for (const c of chunks) { out.set(c, o); o += c.byteLength; }
  return out;
}

export async function gzip(bytes: Uint8Array): Promise<Uint8Array> {
  if (hasStreams) return viaStream(bytes, "c", Infinity);
  const { gzipSync } = await import("node:zlib");
  return new Uint8Array(gzipSync(bytes));
}

export async function gunzip(bytes: Uint8Array, maxBytes = Infinity): Promise<Uint8Array> {
  if (hasStreams) return viaStream(bytes, "d", maxBytes);
  const { gunzipSync } = await import("node:zlib");
  const opts = Number.isFinite(maxBytes) ? { maxOutputLength: maxBytes } : undefined;
  try {
    return new Uint8Array(gunzipSync(bytes, opts));
  } catch (error) {
    const code = (error as { code?: string }).code;
    const message = error instanceof Error ? error.message : "";
    if (Number.isFinite(maxBytes) &&
        (code === "ERR_BUFFER_TOO_LARGE" || message.includes("Buffer larger than")))
      throw new Error(`gunzip: decompressed output exceeds cap (${maxBytes} bytes)`);
    throw error;
  }
}

/** Decompress if the bytes are gzip-framed; otherwise return them unchanged. Honors the same cap. */
export async function maybeGunzip(bytes: Uint8Array, maxBytes = Infinity): Promise<Uint8Array> {
  return isGzip(bytes) ? gunzip(bytes, maxBytes) : bytes;
}
