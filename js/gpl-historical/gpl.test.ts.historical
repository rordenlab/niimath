// GPL build tests: the @niivue/niimath/gpl WASM module. Skipped automatically when
// dist/niimath-gpl.js was not built (e.g. a clone without the niimath_gpl submodule).
// Verifies the GPL build adds the SPM coregistration operations on top of everything
// in the BSD build.
import { describe, test, expect, beforeAll } from 'bun:test';
import { GPL_MODULE, gplBuilt, loadModule, makeNifti, run, type EmscriptenModule } from './helpers';

const d = gplBuilt ? describe : describe.skip;

d('GPL build (@niivue/niimath/gpl)', () => {
  let mod: EmscriptenModule;
  let log: string[];

  beforeAll(async () => {
    ({ mod, log } = await loadModule(GPL_MODULE));
  });

  test('basic arithmetic still works (shared core)', async () => {
    const input = { name: 'in.nii', data: makeNifti(8, (i) => i) };
    const r = await run(mod, log, [input], ['in.nii', '-add', '1', 'out.nii', '-odt', 'float'], 'out.nii');
    expect(r.exitCode).toBe(0);
    expect(r.output).not.toBeNull();
  });

  test('-allineate works (BSD allineate present in GPL build too)', async () => {
    const moving = { name: 'a.nii', data: makeNifti(8, (i) => i % 9) };
    const base = { name: 'b.nii', data: makeNifti(8, (i) => i % 9) };
    const r = await run(
      mod, log, [moving, base],
      ['a.nii', '-allineate', 'b.nii', 'out.nii', '-odt', 'float'], 'out.nii',
    );
    expect(r.exitCode).toBe(0);
  });

  test('-spm_coreg is recognized (the GPL delta over BSD)', async () => {
    // Run with the required reference arg present. Self-coregistration should
    // converge (exit 0); regardless, it must NOT be the BSD stub error.
    const moving = { name: 'a.nii', data: makeNifti(8, (i) => i % 9) };
    const ref = { name: 'b.nii', data: makeNifti(8, (i) => i % 9) };
    const r = await run(
      mod, log, [moving, ref],
      ['a.nii', '-spm_coreg', 'b.nii', 'out.nii', '-odt', 'float'], 'out.nii',
    );
    const text = r.log.join('\n');
    // The BSD stub prints a "requires GPL"-style message; the GPL build must not.
    expect(text.toLowerCase()).not.toContain('optional gpl module');
    expect(r.exitCode).toBe(0);
  });
});
