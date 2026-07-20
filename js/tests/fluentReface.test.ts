// Fluent -reface API: prove the three REQUIRED file operands are emitted in the exact order
// niimath expects (`-reface <tmpl> <shell> <weight>`), each staged under a unique generated name
// that never starts with '-' (else niimath's parser would consume it as a flag) nor collides, and
// that any opts follow ALL three operands (so they parse as -reface sub-options, not as a filename).
// Unlike deface/allineate this binding stages THREE operands, so operand-order/uniqueness matters.
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
const shell = new File([new Uint8Array(0)], 'shell.nii');
const weight = new File([new Uint8Array(0)], 'weight.nii');

describe('fluent -reface API argv construction', () => {
  test('default (no opts): emits -reface then the three staged operands in order', () => {
    const p = makeProc().reface(tmpl, shell, weight);
    const cmds = (p as unknown as { commands: string[] }).commands;
    expect(cmds[0]).toBe('-reface');
    expect(cmds.length).toBe(4); // flag + tmpl + shell + weight, no trailing opts
    expect(cmds[1]).toContain('tmpl.nii');
    expect(cmds[2]).toContain('shell.nii');
    expect(cmds[3]).toContain('weight.nii');
    // Each operand is staged under a generated name that matches its argv token, never begins with
    // '-', and is unique (three distinct MEMFS entries — no operand shadows another).
    const staged = (p as unknown as { extraFiles: { name: string }[] }).extraFiles.map((f) => f.name);
    expect(staged).toEqual([cmds[1], cmds[2], cmds[3]]);
    expect(new Set(staged).size).toBe(3);
    for (const t of [cmds[1], cmds[2], cmds[3]]) expect(t.startsWith('-')).toBe(false);
  });

  test('opts follow all three operands (parsed as -reface sub-options, not a filename)', () => {
    const p = makeProc().reface(tmpl, shell, weight, ['-cost', 'hel']);
    const cmds = (p as unknown as { commands: string[] }).commands;
    expect(cmds[0]).toBe('-reface');
    expect(cmds.indexOf('-cost')).toBe(4); // after flag + three operands
    expect(cmds.slice(-2)).toEqual(['-cost', 'hel']);
    expect(cmds.length).toBe(6);
  });

  test('numeric opts are stringified in place', () => {
    const p = makeProc().reface(tmpl, shell, weight, ['-p', 1]);
    const cmds = (p as unknown as { commands: string[] }).commands;
    expect(cmds.slice(-2)).toEqual(['-p', '1']);
  });
});
