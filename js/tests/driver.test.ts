// The FS and ENV exports a pipeline driver needs: nested paths, a per-call
// FSLOUTPUTTYPE, and discovering side outputs the command line wrote.
import { describe, expect, test } from 'bun:test';
import { BSD_MODULE, type EmscriptenModule, makeNifti } from './helpers';

describe('driver exports', () => {
  test('nested paths, ENV output type and side outputs', async () => {
    const factory = (await import(BSD_MODULE)).default as (o: object) => Promise<EmscriptenModule>;
    // libc copies the environment in its startup constructors: ENV must be set in preRun
    const mod = await factory({ preRun: [(m: EmscriptenModule) => { m.ENV.FSLOUTPUTTYPE = 'NIFTI'; }] });
    mod.FS.mkdirTree('/work/sub-01/anat');
    mod.FS_createDataFile('/work/sub-01/anat', 'in.nii', makeNifti(4, () => 1), true, true);
    // extensionless names: niimath picks the extension from FSLOUTPUTTYPE
    expect(mod.callMain(['/work/sub-01/anat/in.nii', '-add', '1', '-save', '/work/sub-01/anat/side', '-mul', '3',
      '/work/sub-01/anat/out'])).toBe(0);
    expect(mod.FS.readdir('/work/sub-01/anat').filter((f) => f !== '.' && f !== '..').sort())
      .toEqual(['in.nii', 'out.nii', 'side.nii']);
    const voxel = (f: string) => new Float32Array(mod.FS_readFile(f).slice(352).buffer)[0];
    expect(voxel('/work/sub-01/anat/side.nii')).toBe(2);
    expect(voxel('/work/sub-01/anat/out.nii')).toBe(6);
  });
});
