// -n4 in the WASM build: recover a known multiplicative bias field (single-threaded there).
import { describe, test, expect } from 'bun:test';
import { BSD_MODULE, loadModule, makeNifti, run } from './helpers';

describe('-n4 in the WASM build', () => {
  test('recovers a known exp(0.3u + 0.2w) field and rejects a bad option', async () => {
    const { mod, log } = await loadModule(BSD_MODULE);
    const side = 24, c = (side - 1) / 2;
    const at = (i: number) => [(i % side - c) / 10, (Math.floor(i / side) % side - c) / 10, (Math.floor(i / (side * side)) - c) / 10];
    const truth = (i: number) => { const [u, , w] = at(i); return 0.3 * u + 0.2 * w; };
    const r2 = (i: number) => at(i).reduce((s, v) => s + v * v, 0);
    const base = (i: number) => (r2(i) < 0.35 ? 100 : r2(i) < 1 ? 70 : 0);
    mod.FS_createDataFile('.', 'in.nii', makeNifti(side, (i) => base(i) * Math.exp(truth(i))), true, true);
    mod.FS_createDataFile('.', 'w.nii', makeNifti(side, (i) => (r2(i) < 1 ? 1 : 0)), true, true);
    expect(mod.callMain(['in.nii', '-gz', '0', '-n4', '-weight', 'w.nii', '-shrink', '1', '-iter', '20x20',
      '-field', 'f.nii', 'out.nii'])).toBe(0);
    const voxels = (f: string) => new Float32Array(mod.FS_readFile(f).slice(352).buffer);
    const field = voxels('f.nii');
    const inside = [...Array(side ** 3).keys()].filter((i) => r2(i) < 1);
    const est = inside.map((i) => Math.log(field[i])), truth_vals = inside.map(truth);
    const mean = (a: number[]) => a.reduce((s, v) => s + v, 0) / a.length;
    const me = mean(est), mt = mean(truth_vals);
    const cov = est.reduce((s, v, k) => s + (v - me) * (truth_vals[k] - mt), 0);
    const corr = cov / Math.sqrt(est.reduce((s, v) => s + (v - me) ** 2, 0) * truth_vals.reduce((s, v) => s + (v - mt) ** 2, 0));
    expect(corr).toBeGreaterThan(0.98);
    expect((await run(mod, log, [], ['in.nii', '-n4', '-shrink', '0', 'bad.nii'], 'bad.nii')).exitCode).not.toBe(0);
  });
});
