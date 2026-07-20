import type {
  Operators,
  ImageProcessorMethods,
  MeshOptions,
  BitmapOptions,
  DataType
} from './types';

export type {
  Operators,
  OperatorDefinition,
  ImageProcessorMethods,
  MeshOptions,
  BitmapOptions,
  DataType
} from './types';

export const dataTypes = {
  char: "char" as const,
  short: "short" as const,
  int: "int" as const,
  float: "float" as const,
  double: "double" as const,
  input: "input" as const,
} as const;

interface WorkerReadyMessage {
  type: 'ready';
}

interface WorkerErrorMessage {
  type: 'error';
  message: string;
}

interface WorkerSuccessMessage {
  type?: undefined;
  blob: Blob;
  exitCode: number;
}

type WorkerMessage = WorkerReadyMessage | WorkerErrorMessage | WorkerSuccessMessage;

interface WorkerPostMessage {
  blob: File;
  cmd: string[];
  outName: string;
  extraFiles?: { name: string; data: Blob }[];
}

/**
 * Factory that constructs the WASM-backed Web Worker for a build. The BSD entry
 * point (index.ts) and the GPL entry point (index-gpl.ts) each supply their own
 * factory so that esbuild can statically discover and bundle the correct worker
 * (worker.js vs worker-gpl.js) and its WASM binary. The `new Worker(new URL(...))`
 * literal must live in the entry module for esbuild's worker code-splitting to
 * resolve it, which is why this is injected rather than hard-coded here.
 */
export type WorkerFactory = () => Worker;

/**
 * Niimath runs one WASM worker per instance. **Single-flight contract:** only one `.run()` may
 * be in flight per `Niimath` instance at a time, and the owner ENFORCES this fail-fast rather than
 * cross-wiring handlers: a `.run()` before `init()` has resolved (no ready worker yet) rejects
 * with "Worker not initialized", and a second `.run()` while another is still in flight rejects
 * with "niimath is busy: await the previous run()…". Serialize calls (await the previous `run()`
 * before the next), or create a separate `Niimath` instance per concurrent stream.
 *
 * This base class is shared by the BSD and GPL builds; the only difference is the
 * `WorkerFactory` injected via the constructor. Consumers normally use the concrete
 * `Niimath` exported from the package entry points, not this class directly.
 */
export class NiimathBase {
  // Single owner of the worker and the one in-flight operation. A worker processes one op at a
  // time (the API is awaited sequentially) and the owner ENFORCES that — a second concurrent op
  // is rejected, never silently interleaved. init(), run(), dispose(), and a fatal crash all
  // funnel through this owner, and every state change is scoped to the worker GENERATION so a
  // stale message from a replaced/disposed worker can never settle the current one.
  private worker: Worker | null = null;
  private ready = false; // the current worker has sent 'ready' (init resolved) and is usable
  private pendingReject: ((e: Error) => void) | null = null; // non-null == an op is in flight (busy)
  public readonly operators: Operators;
  private outputDataType: DataType = 'float';
  public readonly dataTypes = dataTypes;
  private readonly workerFactory: WorkerFactory;

  constructor(operators: Operators, workerFactory: WorkerFactory) {
    this.operators = operators;
    this.workerFactory = workerFactory;
  }

  init(): Promise<boolean> {
    // Retire any prior generation FIRST (terminate it and reject its pending op) so a
    // re-init never leaks the previous worker or strands its caller.
    this.dispose('niimath worker replaced by a new init()');
    return new Promise((resolve, reject) => {
      let worker: Worker;
      try {
        // Worker construction can throw synchronously (CSP / bad URL / security). Surface it as
        // a rejection so init().catch() sees EVERY init failure, per the Promise contract.
        worker = this.workerFactory();
      } catch (e) {
        reject(e instanceof Error ? e : new Error(String(e)));
        return;
      }
      this.worker = worker;
      this.ready = false;
      this.pendingReject = reject;
      worker.onmessage = (event: MessageEvent<WorkerMessage>) => {
        if (this.worker !== worker) return; // stale generation (replaced/disposed): ignore
        if (event.data && event.data.type === 'ready') {
          this.ready = true;
          this.pendingReject = null;
          resolve(true);
        } else if (event.data && event.data.type === 'error') {
          // A structured error BEFORE 'ready' (e.g. a WASM fetch/instantiate failure) would
          // otherwise leave init() pending forever. Invalidate and reject via the one path.
          this._fail(worker, new Error(event.data.message || 'niimath worker failed to initialize'));
        }
      };
      // A raw worker-level error during load/init.
      worker.onerror = (error: ErrorEvent) => {
        this._fail(worker, new Error(`Worker failed to load: ${error.message}`));
      };
    });
  }

  // Terminate the worker, release its WASM heap, AND reject any in-flight init()/run() so no
  // caller hangs (Worker.terminate() emits no event). Idempotent: safe before init(), after a
  // failure, or called repeatedly. A processor created earlier becomes non-runnable after this
  // (its next run() sees no ready worker and rejects with "not initialized").
  dispose(reason = 'niimath worker disposed'): void {
    const worker = this.worker;
    const reject = this.pendingReject;
    this.worker = null;
    this.ready = false;
    this.pendingReject = null;
    worker?.terminate();
    reject?.(new Error(reason));
  }

  // Fatal error/crash for `worker`. If it is still the current worker, drop it and reject the
  // in-flight op (visible to EVERY ImageProcessor, since they all read this single owner); a
  // stale/superseded worker is just terminated. The one invalidation path for init and run.
  private _fail(worker: Worker, error: Error): void {
    if (this.worker === worker) {
      const reject = this.pendingReject;
      this.worker = null;
      this.ready = false;
      this.pendingReject = null;
      worker.terminate();
      reject?.(error);
    } else {
      worker.terminate();
    }
  }

  // A capability handle for ImageProcessor: it never holds its own worker reference, so worker
  // ownership stays with this base. All operations are generation-scoped (worker identity) so a
  // stale event cannot settle/clobber a newer generation.
  private _handle(): WorkerHandle {
    return {
      // Fail-fast: a run requires a READY, IDLE worker. Throwing here (with NO state mutation on
      // failure) prevents a pre-ready run, or a second overlapping run, from replacing the
      // in-flight op's handlers/rejecter. On success it registers `reject` and returns the worker.
      beginRun: (reject) => {
        if (this.worker === null || !this.ready) {
          throw new Error('Worker not initialized. Did you await the init() method?');
        }
        if (this.pendingReject !== null) {
          throw new Error('niimath is busy: await the previous run() before starting another');
        }
        this.pendingReject = reject;
        return this.worker;
      },
      // Clear the in-flight op ONLY if `worker` is still current — a late result from a
      // replaced/disposed worker must not clear the new worker's rejecter.
      settle: (worker) => { if (this.worker === worker) this.pendingReject = null; },
      isCurrent: (worker) => this.worker === worker,
      fail: (worker, error) => this._fail(worker, error)
    };
  }

  setOutputDataType(type: DataType): void {
    if (Object.values(this.dataTypes).includes(type)) {
      this.outputDataType = type;
    } else {
      throw new Error(`Invalid data type: ${type}`);
    }
  }

  image(file: File): ImageProcessor {
    return new ImageProcessor({
      handle: this._handle(),
      file,
      operators: this.operators,
      outputDataType: this.outputDataType
    });
  }
}

// The narrow worker-owner capability an ImageProcessor needs. Keeps the worker under a single
// owner (NiimathBase): the processor acquires the live worker per run rather than caching one,
// and every call is generation-scoped by worker identity.
interface WorkerHandle {
  beginRun(reject: (e: Error) => void): Worker; // require a ready+idle worker; register the op; else throw
  settle(worker: Worker): void; // the run finished — clear the op iff `worker` is still current
  isCurrent(worker: Worker): boolean; // is `worker` still the owner's current worker?
  fail(worker: Worker, error: Error): void; // fatal crash: invalidate the worker + reject (if current)
}

interface ImageProcessorConfig {
  handle: WorkerHandle;
  file: File;
  operators: Operators;
  outputDataType?: DataType;
}

class ImageProcessor {
  private handle: WorkerHandle;
  private file: File;
  private operators: Operators;
  private commands: string[] = [];
  private outputDataType: DataType;
  // Files (besides the main input) staged into MEMFS by name for chain ops that
  // take filename argv tokens (e.g. -deface/-spm_deface template + mask).
  private extraFiles: { name: string; data: Blob }[] = [];
  // Monotonic counter for generated staging names (collision-proof argv tokens).
  private stagedCounter = 0;

  // Index signature to allow dynamic method assignment from niimath operators
  [key: string]: unknown;

  constructor({ handle, file, operators, outputDataType }: ImageProcessorConfig) {
    this.handle = handle;
    this.file = file;
    this.operators = operators;
    this.outputDataType = outputDataType ?? 'float'; // default to float
    this._generateMethods();
  }

  private _addCommand(cmd: string, ...args: (string | number)[]): this {
    this.commands.push(cmd, ...args.map(String));
    return this;
  }

  // Chain ops that take input filenames as argv tokens (template/mask/ref). The
  // generated fluent methods only handle scalar args, so these are special-cased.
  // Each File is staged into MEMFS under a GENERATED internal name (a unique
  // prefix + the original name, preserving the extension niimath uses to detect
  // gzip/format) and that name is emitted as the argv token. Generated names keep
  // a template/mask/ref whose File.name collides with the input, output, or
  // another staged file from shadowing or unlinking the wrong MEMFS entry.
  // Extra opts (e.g. '-cost', 'nmi') follow.
  private _addFileCommand(flag: string, files: File[], opts: (string | number)[] = []): this {
    // The leading `__nimx<n>_` prefix makes the token unique and ensures it never
    // begins with '-' (which niimath's option parser would consume) nor with a
    // path separator; sanitizing the original to [A-Za-z0-9._-] strips embedded
    // slashes/spaces while keeping the extension niimath reads for gzip/format.
    const names = files.map(
      (f) => `__nimx${this.stagedCounter++}_${f.name.replace(/[^A-Za-z0-9._-]/g, '_')}`,
    );
    this.commands.push(flag, ...names, ...opts.map(String));
    this.extraFiles.push(...files.map((f, i) => ({ name: names[i], data: f })));
    return this;
  }

  // Affine defacing (BSD allineate): -deface <tmpl> <mask> [opts]
  // Opts follow the template/mask argv tokens, e.g. ['-cost', 'hel'] to select the
  // ordinary AFNI-style engine; omit for the default fast (SPM/FLIRT-inspired) engine.
  deface(tmpl: File, mask: File, opts: (string | number)[] = []): this {
    return this._addFileCommand('-deface', [tmpl, mask], opts);
  }

  // SPM rigid-body defacing (GPL spm_coreg): -spm_deface <tmpl> <mask> [opts]
  spmDeface(tmpl: File, mask: File, opts: (string | number)[] = []): this {
    return this._addFileCommand('-spm_deface', [tmpl, mask], opts);
  }

  // SPM rigid-body coregistration (GPL): -spm_coreg <ref> [opts]
  spmcoreg(ref: File, opts: (string | number)[] = []): this {
    return this._addFileCommand('-spm_coreg', [ref], opts);
  }

  // Affine registration (BSD allineate): -allineate <base> [opts] [-weight <img>]
  // The optional `weight` is a base(fixed)-space GRADED weight image, AFNI 3dAllineate style (its
  // dims + world frame must match `base`): normalized to [0, 1] (divide by max) and used per base
  // voxel — a voxel weighted 0 is excluded, one near 1 dominates. It is NOT an exclusion mask; keep
  // the out-of-ROI head attenuated (nonzero) to anchor global scale (a fully-zeroed exterior lets a
  // cross-modal fit collapse into the scalp). It steers BOTH engines — the ordinary engine
  // (`-cost hel`/`lpc`/`lpa`/`ls`) uses it in place of its manufactured autoweight, the fast engine
  // applies it at the finest 2 mm stage only. It is rejected only with stdin and `-applymat`; when
  // the default fast engine falls back to the ordinary engine, the weight is still honored.
  // Emitted as `-weight <img>` after the base + opts and staged
  // into MEMFS like the other file operands.
  allineate(base: File, opts: (string | number)[] = [], weight?: File): this {
    this._addFileCommand('-allineate', [base], opts);
    if (weight) this._addFileCommand('-weight', [weight]);
    return this;
  }

  // Anonymization by face replacement (BSD allineate/reface): -reface <tmpl> <shell> <weight> [opts].
  // Registers the subject to `tmpl`, back-projects the signed template-space `shell` onto the
  // subject grid, and composites an anonymized image. All three file operands are REQUIRED (the
  // `weight` is reused as the registration weight); opts are the `-cost` tuning as for `deface`.
  // For privacy the coverage diagnostic fails closed (<10% mapped → the run errors, no output).
  reface(tmpl: File, shell: File, weight: File, opts: (string | number)[] = []): this {
    return this._addFileCommand('-reface', [tmpl, shell, weight], opts);
  }

  // Nearest-neighbour reslice of the current image onto another image's grid:
  // -reslice_nn <ref>. (e.g. bring a conformed-space mask back to a native grid.)
  resliceNN(ref: File): this {
    return this._addFileCommand('-reslice_nn', [ref]);
  }

  // Multiply the current image by another image: -mul <img>. The generated `mul`
  // only handles a scalar token; this stages a File operand into MEMFS.
  mulImage(img: File): this {
    return this._addFileCommand('-mul', [img]);
  }

  private _generateMethods(): void {
    Object.keys(this.operators).forEach((methodName) => {
      const definition = this.operators[methodName];

      if (methodName === 'kernel') {
        // Special case for kernels because they have different types with varying arguments
        Object.keys(definition.subOperations!).forEach((subOpName) => {
          const subOpDefinition = definition.subOperations![subOpName];
          const kernelMethodName = `kernel${subOpName.charAt(0).toUpperCase() + subOpName.slice(1)}`;

          this[kernelMethodName] = (...args: (string | number)[]) => {
            if (args.length !== subOpDefinition.args.length) {
              throw new Error(`Expected ${subOpDefinition.args.length} arguments for kernel ${subOpName}, but got ${args.length}`);
            }
            return this._addCommand('-kernel', subOpName, ...args);
          };
        });
      } else if (methodName === 'mesh') {
        // Special case for mesh because it has sub-options that can be passed as an object
        this.mesh = (options: MeshOptions = {}) => {
          const subCommands: (string | number)[] = [];

          Object.keys(options).forEach((subOptionKey) => {
            if (definition.subOperations![subOptionKey]) {
              const subOpDefinition = definition.subOperations![subOptionKey];
              const subOptionValue = options[subOptionKey as keyof MeshOptions];

              if (subOpDefinition.args.length > 0 && subOptionValue === undefined) {
                throw new Error(`Sub-option -${subOptionKey} requires a value.`);
              }

              subCommands.push(`-${subOptionKey}`);

              if (subOpDefinition.args.length > 0) {
                subCommands.push(subOptionValue as string | number);
              }
            } else {
              throw new Error(`Invalid sub-option -${subOptionKey} for mesh.`);
            }
          });

          return this._addCommand('-mesh', ...subCommands);
        };
      } else if (methodName === 'bitmap') {
        // Special case for bitmap because it has sub-options that can be passed as an object
        this.bitmap = (outputPath: string, options: BitmapOptions = {}) => {
          const subCommands: (string | number)[] = [outputPath];

          Object.keys(options).forEach((subOptionKey) => {
            if (definition.subOperations![subOptionKey]) {
              const subOpDefinition = definition.subOperations![subOptionKey];
              const subOptionValue = options[subOptionKey as keyof BitmapOptions];

              if (subOpDefinition.args.length > 0 && subOptionValue === undefined) {
                throw new Error(`Sub-option -${subOptionKey} requires a value.`);
              }

              subCommands.push(`-${subOptionKey}`);

              if (subOpDefinition.args.length > 0) {
                if (Array.isArray(subOptionValue)) {
                  subCommands.push(...subOptionValue);
                } else {
                  subCommands.push(subOptionValue as string | number);
                }
              }
            } else {
              throw new Error(`Invalid sub-option -${subOptionKey} for bitmap.`);
            }
          });

          return this._addCommand('-bitmap', ...subCommands);
        };
      } else {
        // General case for non-kernel, non-mesh, and non-bitmap operations
        this[methodName] = (...args: (string | number)[]) => {
          const expectedArgs = definition.args?.length ?? 0;
          if (args.length < expectedArgs) {
            throw new Error(`Expected ${expectedArgs} arguments for ${methodName}, but got ${args.length}`);
          }
          return this._addCommand(`-${methodName}`, ...args);
        };
      }
    });
  }

  async run(outName: string = 'output.nii'): Promise<Blob> {
    return new Promise((resolve, reject) => {
      // outName is an INTERNAL MEMFS filename (the result is returned as a Blob, not written to
      // the caller's filesystem), so require a plain basename that does not use the reserved
      // staging prefixes. The input is staged as `__nimi_*` and file operands as `__nimx<n>_*` in
      // the same MEMFS the output is written to; a caller-supplied outName that resolves to one of
      // those paths recreates the in-place overwrite/cleanup hazard the prefixes exist to prevent.
      // A bare-prefix check is not enough — MEMFS normalizes paths, so `./__nimi_in.nii`,
      // `/__nimx0_w.nii`, or `x/../__nimx0_w.nii` would slip past it — so also reject any path
      // separator or `..` traversal. Checked up front (no worker acquired yet, so nothing to release).
      if (/[/\\]/.test(outName) || outName.split('/').includes('..') ||
          outName.startsWith('__nimi_') || outName.startsWith('__nimx')) {
        reject(new Error(
          `invalid output name '${outName}': use a plain basename that does not contain a path ` +
          `separator or start with the reserved __nimi_/__nimx prefix`));
        return;
      }
      // Acquire the CURRENT worker from the single owner (never a cached copy) and register this
      // run as the in-flight op in one fail-fast transition: throws "not initialized" if there is
      // no ready worker (before init(), or after a crash/dispose) and "busy" if another op is
      // already in flight — so an overlapping/pre-ready run never clobbers the active op.
      let worker: Worker;
      try {
        worker = this.handle.beginRun(reject);
      } catch (e) {
        reject(e as Error);
        return;
      }

      worker.onmessage = (e: MessageEvent) => {
        if (!this.handle.isCurrent(worker)) return; // stale generation: ignore
        const data = e.data as WorkerMessage;
        if (data.type === 'error') {
          this.handle.settle(worker);
          reject(new Error(data.message));
        } else if ('blob' in data && 'exitCode' in data) {
          // get the output file and the exit code from niimath wasm
          this.handle.settle(worker);
          const { blob, exitCode } = data;
          if (exitCode === 0) {
            // success
            resolve(blob);
          } else {
            // error
            reject(new Error(`niimath processing failed with exit code ${exitCode}`));
          }
        }
      };

      // A raw worker-level crash (WASM abort, OOM) during a run must reject THIS run and
      // invalidate the worker at the OWNER (init()'s onerror only settles init). Routing
      // through handle.fail() clears the base's worker too, so a later image().run() gets the
      // clear "not initialized" error rather than posting to a dead worker and hanging.
      worker.onerror = (error: ErrorEvent) => {
        this.handle.fail(worker, new Error(`niimath worker crashed during run: ${error.message}`));
      };

      // Stage the primary input under a generated internal name (sanitized, extension
      // preserved) rather than the raw file.name — otherwise a caller's file whose name
      // matches the fixed output (e.g. re-running on a prior `defaced.nii.gz`) makes
      // input and output share one MEMFS path (fragile in-place overwrite + cleanup).
      // The `__nimi_` prefix never starts with '-'/'/' and can't collide with outName
      // or the `__nimx<n>_` operand names. Mirrors _addFileCommand's staging.
      try {
        const inName = `__nimi_${this.file.name.replace(/[^A-Za-z0-9._-]/g, '_')}`;
        const inputFile = new File([this.file], inName);
        const args = [inName, ...this.commands, outName, '-odt', this.outputDataType];
        const message: WorkerPostMessage = {
          blob: inputFile,
          cmd: args,
          outName: outName,
          extraFiles: this.extraFiles
        };
        worker.postMessage(message);
      } catch (e) {
        // Acquisition + dispatch is one transaction: if staging/postMessage throws (e.g. a
        // non-cloneable payload), RELEASE the in-flight op we just registered — otherwise the
        // owner stays "busy" forever and every later run() rejects until dispose/re-init.
        this.handle.settle(worker);
        reject(e as Error);
      }
    });
  }
}

// File-operand methods are hand-written on the class (not parsed CLI operators),
// so they are absent from the generated ImageProcessorMethods. Declare them here
// (NOT in the regenerated types.ts) so consumers get a complete typed API.
interface FileOperandMethods {
  deface(tmpl: File, mask: File, opts?: (string | number)[]): this;
  spmDeface(tmpl: File, mask: File, opts?: (string | number)[]): this;
  spmcoreg(ref: File, opts?: (string | number)[]): this;
  allineate(base: File, opts?: (string | number)[], weight?: File): this;
  reface(tmpl: File, shell: File, weight: File, opts?: (string | number)[]): this;
  resliceNN(ref: File): this;
  mulImage(img: File): this;
}

// Use interface merging to add method types to ImageProcessor
interface ImageProcessor extends ImageProcessorMethods, FileOperandMethods {}

export { ImageProcessor };
