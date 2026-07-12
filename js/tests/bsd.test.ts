// BSD build tests: the default @niivue/niimath WASM module. Verifies core math
// works and that the BSD build now ships allineate (-allineate/-deface), which is
// public-domain (AFNI 3dAllineate), not GPL.
import { gzipSync } from 'node:zlib';
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

  test('gzip round-trip: reads .nii.gz input and writes .nii.gz output', async () => {
    // Guards WASM compressed I/O (-DHAVE_ZLIB + emscripten -s USE_ZLIB=1). A build
    // that dropped zlib would fail to READ the gzipped input (nonzero exit) or WRITE
    // a real gzip stream (run() gunzips out.nii.gz, so a plain-bytes file named .gz
    // fails to inflate and leaves output null). Exercising both directions in one
    // run mirrors the user's `in.nii.gz -mul 1 out.nii.gz` case.
    const input = { name: 'in.nii.gz', data: gzipSync(makeNifti(8, (i) => i)) };
    const r = await run(
      mod, log, [input],
      ['in.nii.gz', '-mul', '1', 'out.nii.gz', '-odt', 'float'], 'out.nii.gz',
    );
    expect(r.exitCode).toBe(0);
    expect(r.output).not.toBeNull(); // non-null ⇒ out.nii.gz was a valid gzip stream
    expect(voxel(r.output!, 0)).toBeCloseTo(0, 5);
    expect(voxel(r.output!, 5)).toBeCloseTo(5, 5);
  });

  test('--qc computes a TSV from a T1 and hard segmentation', async () => {
    const side = 8;
    const t1 = { name: 't1.nii', data: makeNifti(side, (i) => 20 + (i % 17)) };
    const seg = {
      name: 'seg.nii',
      data: makeNifti(side, (i) => {
        const x = i % side;
        return x < 2 ? 1 : x < 6 ? 2 : 3;
      }),
    };
    const r = await run(
      mod, log, [t1, seg],
      ['--qc', 't1.nii', '--seg', 'seg.nii', '--csf', '1', '--wm', '3', '--erode', '0', '--out', 'qc.tsv'],
      'qc.tsv',
    );
    expect(r.exitCode).toBe(0);
    expect(r.output).not.toBeNull();
    const rows = new TextDecoder().decode(r.output!).trim().split('\n');
    expect(rows).toHaveLength(2);
    expect(rows[0]).toContain('cjv\tcnr_noair\tsnr_csf');
    expect(rows[0]).toContain('summary_wm_n');
    expect(rows[1].split('\t')).toHaveLength(rows[0].split('\t').length);
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
