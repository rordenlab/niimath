// BSD build tests: the default @niivue/niimath WASM module. Verifies core math
// works and that the BSD build now ships allineate (-allineate/-deface), which is
// public-domain (AFNI 3dAllineate), not GPL.
import { describe, test, expect, beforeAll } from 'bun:test';
import { BSD_MODULE, loadModule, makeNifti, run, voxel, type EmscriptenModule } from './helpers';

describe('BSD build (@niivue/niimath)', () => {
  let mod: EmscriptenModule;
  let log: string[];

  beforeAll(async () => {
    ({ mod, log } = await loadModule(BSD_MODULE));
  });

  test('basic arithmetic: -add 10 adds to every voxel', async () => {
    const input = { name: 'in.nii', data: makeNifti(8, (i) => i) };
    const r = await run(mod, log, [input], ['in.nii', '-add', '10', 'out.nii', '-odt', 'float'], 'out.nii');
    expect(r.exitCode).toBe(0);
    expect(r.output).not.toBeNull();
    expect(voxel(r.output!, 0)).toBeCloseTo(10, 5);
    expect(voxel(r.output!, 5)).toBeCloseTo(15, 5);
  });

  test('-allineate registers a volume onto a base (BSD, public-domain)', async () => {
    // Register a volume onto itself — must converge with exit 0.
    const moving = { name: 'a.nii', data: makeNifti(8, (i) => i % 9) };
    const base = { name: 'b.nii', data: makeNifti(8, (i) => i % 9) };
    const r = await run(
      mod, log, [moving, base],
      ['a.nii', '-allineate', 'b.nii', 'out.nii', '-odt', 'float'], 'out.nii',
    );
    expect(r.exitCode).toBe(0);
    expect(r.log.join('\n')).toContain('Registration complete');
    expect(r.output).not.toBeNull();
  });

  test('-deface flag is recognized (not an unknown-option error)', async () => {
    // Missing template/mask args, so this is expected to fail — but it must fail
    // because the args are missing, NOT because -deface is an unknown operation.
    // An unknown flag would be the failure mode if allineate were absent.
    const input = { name: 'in.nii', data: makeNifti(8) };
    const r = await run(mod, log, [input], ['in.nii', '-deface', 'out.nii', '-odt', 'float'], 'out.nii');
    const text = r.log.join('\n').toLowerCase();
    expect(text).not.toContain('unknown');
  });

  test('-spm_coreg is NOT functional in the BSD build (GPL-only)', async () => {
    const moving = { name: 'a.nii', data: makeNifti(8) };
    const ref = { name: 'b.nii', data: makeNifti(8) };
    const r = await run(
      mod, log, [moving, ref],
      ['a.nii', '-spm_coreg', 'b.nii', 'out.nii', '-odt', 'float'], 'out.nii',
    );
    // The BSD build carries only a stub for spm_coreg; it must not succeed.
    expect(r.exitCode).not.toBe(0);
  });
});
