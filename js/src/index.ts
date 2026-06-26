import operators from './niimathOperators.json' with { type: 'json' };
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

export class Niimath {
  private worker: Worker | null = null;
  public readonly operators: Operators;
  private outputDataType: DataType = 'float';
  public readonly dataTypes = dataTypes;

  constructor() {
    this.operators = operators as Operators;
  }

  init(): Promise<boolean> {
    this.worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
    return new Promise((resolve, reject) => {
      // Handle worker ready message.
      // This gets reassigned in the run() method,
      // but we need to handle the ready message before that.
      // Maybe there is a less hacky way to do this?
      this.worker!.onmessage = (event: MessageEvent<WorkerMessage>) => {
        if (event.data && event.data.type === 'ready') {
          resolve(true); // Resolve the promise when the worker is ready
        }
      };

      // Handle worker init errors.
      this.worker!.onerror = (error: ErrorEvent) => {
        reject(new Error(`Worker failed to load: ${error.message}`));
      };
    });
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
      worker: this.worker,
      file,
      operators: this.operators,
      outputDataType: this.outputDataType
    });
  }
}

interface ImageProcessorConfig {
  worker: Worker | null;
  file: File;
  operators: Operators;
  outputDataType?: DataType;
}

class ImageProcessor {
  private worker: Worker | null;
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

  constructor({ worker, file, operators, outputDataType }: ImageProcessorConfig) {
    this.worker = worker;
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

  // Affine defacing (BSD allineate): -deface <tmpl> <mask>
  deface(tmpl: File, mask: File): this {
    return this._addFileCommand('-deface', [tmpl, mask]);
  }

  // SPM rigid-body defacing (GPL spm_coreg): -spm_deface <tmpl> <mask> [opts]
  spmDeface(tmpl: File, mask: File, opts: (string | number)[] = []): this {
    return this._addFileCommand('-spm_deface', [tmpl, mask], opts);
  }

  // SPM rigid-body coregistration (GPL): -spmcoreg <ref> [opts]
  spmcoreg(ref: File, opts: (string | number)[] = []): this {
    return this._addFileCommand('-spmcoreg', [ref], opts);
  }

  // Affine registration (BSD allineate): -allineate <base> [opts]
  allineate(base: File, opts: (string | number)[] = []): this {
    return this._addFileCommand('-allineate', [base], opts);
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
      // Check the worker exists BEFORE touching it — otherwise `.run()` before
      // `init()` throws a raw TypeError on `this.worker!.onmessage` instead of
      // this clear message.
      if (this.worker === null) {
        reject(new Error('Worker not initialized. Did you await the init() method?'));
        return;
      }

      this.worker.onmessage = (e: MessageEvent) => {
        const data = e.data as WorkerMessage;
        if (data.type === 'error') {
          reject(new Error(data.message));
        } else if ('blob' in data && 'exitCode' in data) {
          // get the output file and the exit code from niimath wasm
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

      const args = [this.file.name, ...this.commands, outName, '-odt', this.outputDataType];
      const message: WorkerPostMessage = {
        blob: this.file,
        cmd: args,
        outName: outName,
        extraFiles: this.extraFiles
      };
      this.worker.postMessage(message);
    });
  }
}

// Use interface merging to add method types to ImageProcessor
interface ImageProcessor extends ImageProcessorMethods {}