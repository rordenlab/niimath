// Threaded (OpenMP) build, run under Node: Bun's workers cannot host Emscripten pthreads.
// Catches a missing pool cap (in a browser a -p above the pool deadlocks; see niimath.c main) and
// -allineate output depending on the thread count.
import assert from 'node:assert/strict';

const side = 24;
function blob(cx) {
  const buf = new ArrayBuffer(352 + side ** 3 * 4);
  const dv = new DataView(buf);
  dv.setInt32(0, 348, true);
  [3, side, side, side, 1, 1, 1, 1].forEach((d, i) => dv.setInt16(40 + i * 2, d, true));
  dv.setInt16(70, 16, true); // float32
  dv.setInt16(72, 32, true);
  for (let i = 0; i < 8; i++) dv.setFloat32(76 + i * 4, 1, true);
  dv.setFloat32(108, 352, true);
  dv.setUint8(123, 2);
  [...'n+1\0'].forEach((c, i) => dv.setUint8(344 + i, c.charCodeAt(0)));
  const f = new Float32Array(buf, 352);
  for (let i = 0; i < f.length; i++) {
    const x = i % side, y = Math.floor(i / side) % side, z = Math.floor(i / side ** 2);
    f[i] = 100 * Math.exp(-((x - cx) ** 2 + (y - 12) ** 2 + (z - 12) ** 2) / 20);
  }
  return new Uint8Array(buf);
}

const factory = (await import('../dist/niimath-mt.js')).default;
const log = [];
const mod = await factory({ print() {}, printErr: (s) => log.push(s) });
mod.FS_createDataFile('/', 'a.nii', blob(10), true, true);
mod.FS_createDataFile('/', 'b.nii', blob(12), true, true);
const reg = (p) => {
  assert.equal(mod.callMain(['a.nii', '-p', p, '-allineate', 'b.nii', '-gz', '0', 'out.nii', '-odt', 'float']), 0);
  return Buffer.from(mod.FS_readFile('out.nii'));
};
assert.ok(reg('1').equals(reg('64')), '-allineate differs between -p 1 and -p 64');
// Node can spawn workers on demand, so only this proves the browser-deadlock cap is applied.
assert.match(log.join('\n'), /Cannot form a team with 64 threads/);
console.log('mt.node.mjs: ok');
process.exit(0); // pool workers keep Node alive
