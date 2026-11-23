import { expect, test } from 'bun:test';
import { mkdtempSync, mkdirSync, readFileSync, rmSync, writeFileSync } from 'fs';
import { tmpdir } from 'os';
import { join } from 'path';
import { spawnSync } from 'child_process';

const workflow = Bun.YAML.parse(readFileSync(new URL('../../.github/workflows/js.yml', import.meta.url), 'utf8')) as any;
const publish = workflow.jobs.publish;
const script = publish.steps.find((step: any) => step.name === 'Publish to npm if the version is new').run;

test('only the publish job can request OIDC credentials, after a successful master build', () => {
  expect(workflow.permissions?.['id-token']).not.toBe('write');
  for (const [name, job] of Object.entries(workflow.jobs)) {
    expect((job as any).permissions?.['id-token']).toBe(name === 'publish' ? 'write' : undefined);
  }
  expect(publish.needs).toBe('prebuild');
  expect(publish.if).toBe("github.event_name == 'push' && github.ref == 'refs/heads/master'");
});

test.each([
  ['existing version', '["1.0.0", "1.0.1"]', 0, true, false],
  ['new version', '["1.0.0"]', 0, true, true],
  ['registry failure', '', 1, false, false],
  ['registry failure with output', '["1.0.0"]', 1, false, false],
  ['invalid registry JSON', 'unavailable', 0, false, false],
])('npm publishing: %s', (_name, response, status, succeeds, publishes) => {
  const dir = mkdtempSync(join(tmpdir(), 'niimath-publish-'));
  try {
    mkdirSync(join(dir, 'package'));
    writeFileSync(join(dir, 'package/package.json'), JSON.stringify({ version: '1.0.1' }));
    expect(spawnSync('tar', ['-czf', 'package.tgz', 'package'], { cwd: dir }).status).toBe(0);
    rmSync(join(dir, 'package'), { recursive: true });
    // Mock npm only: execute the actual workflow shell, tar extraction and JSON checks.
    const result = spawnSync('bash', ['-e', '-c', `
      npm() {
        printf '%s\n' "$*" >> npm-calls
        if [ "$1" = view ]; then
          printf '%s' "$REGISTRY_RESPONSE"
          return "$REGISTRY_STATUS"
        fi
      }
      ${script}
    `], {
      cwd: dir,
      env: { ...process.env, REGISTRY_RESPONSE: String(response), REGISTRY_STATUS: String(status) },
      encoding: 'utf8',
    });
    expect(result.status === 0).toBe(succeeds);
    const calls = readFileSync(join(dir, 'npm-calls'), 'utf8').trim().split('\n');
    expect(calls).toEqual(publishes ? [
      'view @niivue/niimath versions --json',
      'install -g npm@11',
      'publish ./package.tgz --access public --ignore-scripts',
    ] : ['view @niivue/niimath versions --json']);
  } finally {
    rmSync(dir, { recursive: true, force: true });
  }
});
