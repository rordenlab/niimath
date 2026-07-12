// Fluent -deface API: prove that opts passed to `.deface(tmpl, mask, opts)` are emitted
// AFTER the template/mask operand tokens (so niimath parses them as -deface sub-options,
// not as the template filename), and that `['-cost', 'hel']` is plumbed through to select
// the ordinary engine. The ordinary-vs-fast ENGINE ROUTING of `-cost hel` is covered at the
// C level by gpl-build.yml and the WASI `-cost hel` byte-parity test; here we guard the
// JS-side argv construction that makes that routing reachable.
import { describe, test, expect } from 'bun:test';
import { ImageProcessor } from '../src/core';

function makeProc(): ImageProcessor {
  return new ImageProcessor({
    worker: null as unknown as Worker,
    file: new File([new Uint8Array(0)], 'in.nii'),
    operators: {} as never,
    outputDataType: 'float',
  });
}

const tmpl = new File([new Uint8Array(0)], 'tmpl.nii');
const mask = new File([new Uint8Array(0)], 'mask.nii');

describe('fluent -deface API argv construction', () => {
  test('default (no opts): emits -deface then the two staged operands', () => {
    const p = makeProc().deface(tmpl, mask);
    const cmds = (p as unknown as { commands: string[] }).commands;
    expect(cmds[0]).toBe('-deface');
    expect(cmds.length).toBe(3); // no trailing opts -> niimath defaults to the fast engine
    expect(cmds[1]).toContain('tmpl.nii');
    expect(cmds[2]).toContain('mask.nii');
    // Both operands are staged under generated names that never start with '-'.
    const staged = (p as unknown as { extraFiles: { name: string }[] }).extraFiles.map((f) => f.name);
    expect(staged).toEqual([cmds[1], cmds[2]]);
    expect(cmds[1].startsWith('-')).toBe(false);
    expect(cmds[2].startsWith('-')).toBe(false);
  });

  test("['-cost','hel'] follows the operands (reaches the ordinary engine)", () => {
    const p = makeProc().deface(tmpl, mask, ['-cost', 'hel']);
    const cmds = (p as unknown as { commands: string[] }).commands;
    // The opts MUST come after -deface + the two operands; otherwise niimath would read
    // '-cost' as the template filename and the ordinary engine would never be selected.
    expect(cmds[0]).toBe('-deface');
    expect(cmds.indexOf('-cost')).toBe(3);
    expect(cmds.slice(-2)).toEqual(['-cost', 'hel']);
    expect(cmds.length).toBe(5);
  });

  test('numeric opts are stringified in place', () => {
    const p = makeProc().deface(tmpl, mask, ['-final', 0]);
    const cmds = (p as unknown as { commands: string[] }).commands;
    expect(cmds.slice(-2)).toEqual(['-final', '0']);
  });
});
