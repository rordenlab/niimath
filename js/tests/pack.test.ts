// Guards the BSD-only package contract at the tarball level: `npm pack` force-includes
// license-like files regardless of the `files` whitelist, so a stray root LICENSE.GPL-2.0.txt
// would silently ship. This asserts the published tarball contains the BSD LICENSE and NONE of
// the GPL binaries/exports/license material or build/test declarations. One assertion covering
// exports, binaries, and license files (audit follow-up to the GPL de-ship).
import { describe, test, expect } from 'bun:test';
import { execSync } from 'child_process';
import { readFileSync } from 'fs';
import { tmpdir } from 'os';
import { join } from 'path';

const packageRoot = new URL('..', import.meta.url).pathname;
const packageJson = JSON.parse(readFileSync(join(packageRoot, 'package.json'), 'utf8'));

function packedFiles(): string[] {
  // --dry-run computes the file list without writing a tarball; --json emits it structured.
  const out = execSync('npm pack --dry-run --json', {
    cwd: packageRoot,
    encoding: 'utf8',
    stdio: ['ignore', 'pipe', 'ignore'],
    env: { ...process.env, npm_config_cache: join(tmpdir(), 'niimath-npm-pack-cache') },
  });
  const meta = JSON.parse(out);
  return (meta[0]?.files ?? []).map((f: { path: string }) => f.path);
}

describe('npm tarball is BSD-only', () => {
  const files = packedFiles();

  test('ships the BSD LICENSE and package manifest', () => {
    expect(files).toContain('LICENSE');
    expect(files).toContain('package.json');
  });

  test('ships NO GPL binaries, exports, or license material', () => {
    const gpl = files.filter((p) => /gpl|corresponding-source/i.test(p));
    expect(gpl).toEqual([]);
    expect(packageJson.license).toBe('BSD-2-Clause');
    expect(Object.keys(packageJson.exports)).not.toContain('./gpl');
    expect(Object.keys(packageJson.exports)).not.toContain('./niimath-gpl.js');
  });

  test('ships NO test or build-script declarations', () => {
    const decls = files.filter((p) => /\.test\.d\.ts$|check-imports/i.test(p));
    expect(decls).toEqual([]);
  });
});
