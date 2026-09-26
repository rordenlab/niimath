// qc(): argv construction through a fake worker, and the built WASM running --qc --pve.
import { describe, test, expect } from 'bun:test';
import { ImageProcessor } from '../src/core';
import { BSD_MODULE, loadModule, makeNifti } from './helpers';

// Runs qc() against a fake worker that records the posted job and answers with `report`.
async function qcArgv(tissues: Parameters<ImageProcessor['qc']>[0], air?: File, report: object = { cjv: 1 }) {
  let job: { cmd: string[]; outName: string; extraFiles: { name: string }[] } | undefined;
  const worker = {
    onmessage: null as ((e: MessageEvent) => void) | null,
    postMessage(message: typeof job) {
      job = message;
      this.onmessage!({ data: { blob: new Blob([JSON.stringify(report)]), exitCode: 0 } } as MessageEvent);
    },
  };
  const handle = { beginRun: () => worker, settle: () => {}, isCurrent: () => true, fail: () => {} };
  const proc = new ImageProcessor({ handle, file: new File([], 't1.nii'), operators: {} } as never);
  return { report: await proc.qc(tissues, air), job: job! };
}

const file = (name: string) => new File([], name);

describe('qc()', () => {
  test('labels: --qc <in> --seg <s> --csf --wm [--air] --json <out>, report parsed', async () => {
    const { report, job } = await qcArgv({ seg: file('seg.nii'), csf: [3, 4], wm: [1, 5] }, file('tmpl.nii.gz'),
      { cjv: 1, provenance: { air_template: '__nimx1_tmpl.nii.gz' } });
    expect(report).toEqual({ cjv: 1, provenance: { air_template: 'tmpl.nii.gz' } });
    const c = job.cmd;
    expect(c[0]).toBe('--qc');
    expect(c[1]).toContain('t1.nii');
    expect(c[2]).toBe('--seg');
    expect(c.slice(4, 8)).toEqual(['--csf', '3,4', '--wm', '1,5']);
    expect(c[8]).toBe('--air');
    expect(c.slice(-2)).toEqual(['--json', 'qc.json']);
    expect(job.extraFiles.map((f) => f.name)).toEqual([c[3], c[9]]);
  });

  test('pve: --pve <csf> <gm> <wm>, no -odt', async () => {
    const { job } = await qcArgv({ pve: [file('c.nii'), file('g.nii'), file('w.nii')] });
    expect(job.cmd[2]).toBe('--pve');
    expect(job.cmd.slice(3, 6).map((n) => n.replace(/^__nimx\d+_/, ''))).toEqual(['c.nii', 'g.nii', 'w.nii']);
    expect(job.cmd).not.toContain('-odt');
  });
});

test('qc() refuses chain ops before it', async () => {
  const handle = { beginRun: () => { throw new Error('must not run'); }, settle: () => {}, isCurrent: () => true, fail: () => {} };
  const proc = new ImageProcessor({ handle, file: file('t1.nii'), operators: {} } as never);
  (proc as unknown as { commands: string[] }).commands.push('-s', '2');
  await expect(proc.qc({ pve: [file('c'), file('g'), file('w')] })).rejects.toThrow('no chain ops');
});

describe('WASM --qc --pve', () => {
  test('0/1 fractions reproduce --seg --erode 0', async () => {
    const { mod } = await loadModule(BSD_MODULE);
    const side = 12;
    const tissue = (i: number) => Math.floor((i % side) / 4); // three x-slabs: CSF, GM, WM
    mod.FS_createDataFile('.', 't1.nii', makeNifti(side, (i) => [20, 80, 120][tissue(i)] + (i % 7)), true, true);
    mod.FS_createDataFile('.', 'seg.nii', makeNifti(side, (i) => tissue(i) + 1), true, true);
    for (const [t, name] of ['csf', 'gm', 'wm'].entries())
      mod.FS_createDataFile('.', `${name}.nii`, makeNifti(side, (i) => Number(tissue(i) === t)), true, true);
    expect(mod.callMain(['--qc', 't1.nii', '--seg', 'seg.nii', '--csf', '1', '--wm', '3', '--erode', '0', '--json', 'hard.json'])).toBe(0);
    expect(mod.callMain(['--qc', 't1.nii', '--pve', 'csf.nii', 'gm.nii', 'wm.nii', '--json', 'pve.json'])).toBe(0);
    const read = (f: string) => JSON.parse(new TextDecoder().decode(mod.FS_readFile(f)));
    const { provenance: _, ...hard } = read('hard.json');
    const { provenance, ...pve } = read('pve.json');
    expect(provenance.pve).toBe(true);
    expect(pve).toEqual(hard);
  });
});
