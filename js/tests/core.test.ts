// Fake-worker unit tests for the NiimathBase / ImageProcessor worker OWNERSHIP + lifecycle state
// machine in src/core.ts. Unlike the other suites (which drive the built WASM through callMain),
// these inject a fake Worker via NiimathBase's `workerFactory`, so the pure-TypeScript owner —
// init/ready/error, re-init retirement, ready+idle run admission, generation-scoped settlement,
// crash invalidation, dispose settling pending work, and operand staging — is exercised with no
// WASM build. core.ts imports only types, so importing it here has no runtime deps. Every run
// started in a test is settled (a result emitted) or rejected (crash/dispose) — none is left
// pending, so no test leaks owner state into the next.
import { describe, test, expect } from 'bun:test';
import { NiimathBase } from '../src/core';

interface Posted {
  cmd: string[];
  extraFiles: { name: string; data: Blob }[];
}

// Minimal stand-in for a browser Worker: records postMessage payloads and terminations, and lets
// the test drive onmessage/onerror deterministically.
class FakeWorker {
  onmessage: ((e: { data: unknown }) => void) | null = null;
  onerror: ((e: { message: string }) => void) | null = null;
  posted: Posted[] = [];
  terminated = 0;
  postMessage(msg: Posted): void {
    this.posted.push(msg);
  }
  terminate(): void {
    this.terminated++;
  }
  emit(data: unknown): void {
    this.onmessage?.({ data });
  }
  emitError(message: string): void {
    this.onerror?.({ message });
  }
}

// eslint-disable-next-line @typescript-eslint/no-explicit-any
const baseWith = (factory: () => FakeWorker) => new NiimathBase({} as any, () => factory() as unknown as Worker);
const inputNii = () => new File([new Uint8Array([1, 2, 3, 4])], 'in.nii');
const okMsg = () => ({ blob: new Blob([new Uint8Array([0])]), exitCode: 0 });

async function ready(): Promise<{ fake: FakeWorker; niimath: NiimathBase }> {
  const fake = new FakeWorker();
  const niimath = baseWith(() => fake);
  const p = niimath.init();
  fake.emit({ type: 'ready' });
  await p;
  return { fake, niimath };
}

describe('init / dispose lifecycle', () => {
  test('init resolves on ready', async () => {
    const fake = new FakeWorker();
    const p = baseWith(() => fake).init();
    fake.emit({ type: 'ready' });
    expect(await p).toBe(true);
  });

  test('init rejects and terminates on a structured error before ready', async () => {
    const fake = new FakeWorker();
    const p = baseWith(() => fake).init();
    fake.emit({ type: 'error', message: 'wasm instantiate failed' });
    await expect(p).rejects.toThrow('wasm instantiate failed');
    expect(fake.terminated).toBe(1);
  });

  test('init rejects and terminates on a raw worker error', async () => {
    const fake = new FakeWorker();
    const p = baseWith(() => fake).init();
    fake.emitError('failed to load module');
    await expect(p).rejects.toThrow('failed to load');
    expect(fake.terminated).toBe(1);
  });

  test('init rejects (does not throw) when the worker factory throws synchronously', async () => {
    const niimath = baseWith(() => {
      throw new Error('CSP blocked worker construction');
    });
    // Must be a rejected promise, catchable by init().catch(...), not a synchronous throw.
    await expect(niimath.init()).rejects.toThrow('CSP blocked');
  });

  test('re-init retires the prior worker (no leak) and rejects its pending init', async () => {
    const workers = [new FakeWorker(), new FakeWorker()];
    let n = 0;
    const niimath = baseWith(() => workers[n++]);
    const p1 = niimath.init(); // left pending (never readied)
    const p2 = niimath.init(); // supersedes -> must retire+reject worker 0
    await expect(p1).rejects.toThrow('replaced');
    expect(workers[0].terminated).toBe(1);
    workers[1].emit({ type: 'ready' });
    expect(await p2).toBe(true);
    const runP = niimath.image(inputNii()).run('out.nii'); // uses worker 1
    expect(workers[1].posted.length).toBe(1);
    expect(workers[1].terminated).toBe(0);
    workers[1].emit(okMsg());
    await runP;
  });

  test('dispose terminates the worker and is idempotent', async () => {
    const { fake, niimath } = await ready();
    niimath.dispose();
    niimath.dispose();
    expect(fake.terminated).toBe(1);
  });

  test('dispose during a pending init rejects the init', async () => {
    const fake = new FakeWorker();
    const niimath = baseWith(() => fake);
    const p = niimath.init(); // pending
    niimath.dispose();
    await expect(p).rejects.toThrow('disposed');
    expect(fake.terminated).toBe(1);
  });

  test('dispose is safe before init()', () => {
    expect(() => baseWith(() => new FakeWorker()).dispose()).not.toThrow();
  });
});

describe('run() admission and ownership', () => {
  test('run before init rejects clearly', async () => {
    const niimath = baseWith(() => new FakeWorker());
    await expect(niimath.image(inputNii()).run()).rejects.toThrow('not initialized');
  });

  test('run before the worker is READY rejects, and init still resolves', async () => {
    const fake = new FakeWorker();
    const niimath = baseWith(() => fake);
    const initP = niimath.init(); // pending, not ready
    await expect(niimath.image(inputNii()).run('out.nii')).rejects.toThrow('not initialized');
    fake.emit({ type: 'ready' });
    expect(await initP).toBe(true);
  });

  test('a second overlapping run is rejected as busy; the first still completes', async () => {
    const { fake, niimath } = await ready();
    const runA = niimath.image(inputNii()).run('a.nii');
    const runB = niimath.image(inputNii()).run('b.nii'); // A still in flight
    await expect(runB).rejects.toThrow('busy');
    fake.emit(okMsg()); // A's handlers are intact
    expect(await runA).toBeInstanceOf(Blob);
  });

  test('a worker crash during run rejects THAT run and a fresh run then reports not-initialized', async () => {
    const { fake, niimath } = await ready();
    const runP = niimath.image(inputNii()).run('out.nii');
    fake.emitError('abort() called');
    await expect(runP).rejects.toThrow('crashed');
    expect(fake.terminated).toBe(1);
    await expect(niimath.image(inputNii()).run('out2.nii')).rejects.toThrow('not initialized');
  });

  test('a synchronous postMessage failure releases the busy state', async () => {
    const { fake, niimath } = await ready();
    let boom = true;
    const realPost = fake.postMessage.bind(fake);
    fake.postMessage = (msg) => {
      if (boom) {
        boom = false;
        throw new Error('DataCloneError: could not clone payload');
      }
      realPost(msg);
    };
    await expect(niimath.image(inputNii()).run('a.nii')).rejects.toThrow('DataCloneError');
    // The owner must be IDLE again (not stuck "busy") — a fresh run proceeds and completes.
    const runB = niimath.image(inputNii()).run('b.nii');
    expect(fake.posted.length).toBe(1);
    fake.emit(okMsg());
    expect(await runB).toBeInstanceOf(Blob);
  });

  test('dispose during a pending run rejects the run', async () => {
    const { niimath } = await ready();
    const runP = niimath.image(inputNii()).run('out.nii');
    niimath.dispose();
    await expect(runP).rejects.toThrow('disposed');
  });

  test('a late result from a replaced worker does not clear the current run', async () => {
    const workers = [new FakeWorker(), new FakeWorker()];
    let n = 0;
    const niimath = baseWith(() => workers[n++]);
    const [w0, w1] = workers;
    let p = niimath.init();
    w0.emit({ type: 'ready' });
    await p;

    const runA = niimath.image(inputNii()).run('a.nii'); // on w0
    expect(w0.posted.length).toBe(1);
    p = niimath.init(); // dispose rejects runA + retires w0
    await expect(runA).rejects.toThrow('replaced');
    w1.emit({ type: 'ready' });
    await p;

    const runB = niimath.image(inputNii()).run('b.nii'); // on w1
    expect(w1.posted.length).toBe(1);
    w0.emit(okMsg()); // stale result from the retired w0 — must be IGNORED, not settle B
    // Proof B's rejecter survived: dispose must still reject it (would hang if it were cleared).
    niimath.dispose();
    await expect(runB).rejects.toThrow('disposed');
  });

  test('an output name that could alias a staged path is rejected (incl. normalized forms)', async () => {
    const { fake, niimath } = await ready();
    // Bare reserved prefixes AND path forms that MEMFS would normalize onto a staged path.
    for (const bad of ['__nimi_in.nii', '__nimx0_w.nii', './__nimi_in.nii', '/__nimx0_w.nii', 'x/../__nimx0_w.nii']) {
      await expect(niimath.image(inputNii()).run(bad)).rejects.toThrow(/invalid output name/);
    }
    // All rejected up front (no worker acquired), so the instance is still idle — a normal run works.
    const runP = niimath.image(inputNii()).run('out.nii');
    fake.emit(okMsg());
    expect(await runP).toBeInstanceOf(Blob);
  });

  test('run resolves with the output blob on exit code 0', async () => {
    const { fake, niimath } = await ready();
    const runP = niimath.image(inputNii()).run('out.nii');
    fake.emit(okMsg());
    expect(await runP).toBeInstanceOf(Blob);
  });

  test('run rejects on a nonzero exit code', async () => {
    const { fake, niimath } = await ready();
    const runP = niimath.image(inputNii()).run('out.nii');
    fake.emit({ blob: new Blob([]), exitCode: 1 });
    await expect(runP).rejects.toThrow('exit code 1');
  });

  test('a processor created before re-init runs on the NEW worker', async () => {
    const workers = [new FakeWorker(), new FakeWorker()];
    let n = 0;
    const niimath = baseWith(() => workers[n++]);
    let p = niimath.init();
    workers[0].emit({ type: 'ready' });
    await p;
    const proc = niimath.image(inputNii()); // created against worker 0
    p = niimath.init(); // replace with worker 1
    workers[1].emit({ type: 'ready' });
    await p;
    const runP = proc.run('out.nii'); // claims the CURRENT worker, not a cached worker 0
    expect(workers[1].posted.length).toBe(1);
    expect(workers[0].posted.length).toBe(0);
    workers[1].emit(okMsg());
    await runP;
  });
});

describe('allineate -weight binding', () => {
  test('emits -allineate <base> then -weight <img>, staging both operands', async () => {
    const { fake, niimath } = await ready();
    const base = new File([new Uint8Array([5])], 'base.nii');
    const weight = new File([new Uint8Array([6])], 'weight.nii');
    const runP = niimath.image(inputNii()).allineate(base, [], weight).run('out.nii');

    const { cmd, extraFiles } = fake.posted[0];
    const ai = cmd.indexOf('-allineate');
    const wi = cmd.indexOf('-weight');
    expect(ai).toBeGreaterThanOrEqual(0);
    expect(wi).toBeGreaterThan(ai); // -weight comes after -allineate + its base
    expect(cmd[ai + 1]).toContain('base.nii'); // base token immediately follows -allineate
    expect(cmd[wi + 1]).toContain('weight.nii'); // weight token immediately follows -weight
    expect(extraFiles.some((f) => f.name.includes('base.nii'))).toBe(true);
    expect(extraFiles.some((f) => f.name.includes('weight.nii'))).toBe(true);

    fake.emit(okMsg());
    await runP;
  });

  test('omitting the weight emits only -allineate', async () => {
    const { fake, niimath } = await ready();
    const base = new File([new Uint8Array([5])], 'base.nii');
    const runP = niimath.image(inputNii()).allineate(base, []).run('out.nii');

    const { cmd } = fake.posted[0];
    expect(cmd).toContain('-allineate');
    expect(cmd).not.toContain('-weight');

    fake.emit(okMsg());
    await runP;
  });
});
