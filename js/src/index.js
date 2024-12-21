import operators from './niimathOperators.json';
export class Niimath {
  constructor() {
    this.worker = null;
    this.operators = operators;
    this.outputDataType = 'float'; // Possible datatypes are: char short int float double input
    this.dataTypes = {char: 'char', short: 'short', int: 'int', float: 'float', double: 'double', input: 'input'};
  }

  init() {
    this.worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
    return new Promise((resolve, reject) => {
      // Handle worker ready message.
      // This gets reassigned in the run() method, 
      // but we need to handle the ready message before that.
      // Maybe there is a less hacky way to do this?
      this.worker.onmessage = (event) => {
        if (event.data && event.data.type === 'ready') {
          resolve(true); // Resolve the promise when the worker is ready
        }
      }

      // Handle worker init errors.
      this.worker.onerror = (error) => {
        reject(new Error(`Worker failed to load: ${error.message}`));
      }
    });
  }

  setOutputDataType(type) {
    if (Object.values(this.dataTypes).includes(type)) {
      this.outputDataType = type;
    } else {
      throw new Error(`Invalid data type: ${type}`);
    }
  }

  image(file) {
    return new ImageProcessor({ worker: this.worker, file, operators: this.operators, outputDataType: this.outputDataType });
  }
}

class ImageProcessor {

  constructor({ worker, file, operators, outputDataType }) {
    this.worker = worker;
    this.file = file;
    this.operators = operators;
    this.commands = [];
    this.outputDataType = outputDataType ? outputDataType : 'float'; // default to float
    this._generateMethods();
  }

  _addCommand(cmd, ...args) {
    this.commands.push(cmd, ...args.map(String));
    return this;
  }

  _generateMethods() {
    Object.keys(this.operators).forEach((methodName) => {
      const definition = this.operators[methodName];

      if (methodName === 'kernel') {
        // Special case for kernels because they have different types with varying arguments
        Object.keys(definition.subOperations).forEach((subOpName) => {
          const subOpDefinition = definition.subOperations[subOpName];
          this[`kernel${subOpName.charAt(0).toUpperCase() + subOpName.slice(1)}`] = (...args) => {
            if (args.length !== subOpDefinition.args.length) {
              throw new Error(`Expected ${subOpDefinition.args.length} arguments for kernel ${subOpName}, but got ${args.length}`);
            }
            return this._addCommand('-kernel', subOpName, ...args);
          };
        });
      } else if (methodName === 'mesh') {
        // Special case for mesh because it has sub-options that can be passed as an object
        this.mesh = (options = {}) => {
          const subCommands = [];

          Object.keys(options).forEach((subOptionKey) => {
            if (definition.subOperations[subOptionKey]) {
              const subOpDefinition = definition.subOperations[subOptionKey];
              const subOptionValue = options[subOptionKey];

              if (subOpDefinition.args.length > 0 && subOptionValue === undefined) {
                throw new Error(`Sub-option -${subOptionKey} requires a value.`);
              }

              subCommands.push(`-${subOptionKey}`);

              if (subOpDefinition.args.length > 0) {
                subCommands.push(subOptionValue);
              }
            } else {
              throw new Error(`Invalid sub-option -${subOptionKey} for mesh.`);
            }
          });

          return this._addCommand('-mesh', ...subCommands);
        };
      } else {
        // General case for non-kernel and non-mesh operations
        this[methodName] = (...args) => {
          if (args.length < definition.args.length || (!definition.optional && args.length > definition.args.length)) {
            throw new Error(`Expected ${definition.args.length} arguments for ${methodName}, but got ${args.length}`);
          }
          return this._addCommand(`-${methodName}`, ...args);
        };
      }
    });
  }

  async run(outName = 'output.nii') {
    return new Promise((resolve, reject) => {
      this.worker.onmessage = (e) => {
        if (e.data.type === 'error') {
          reject(new Error(e.data.message));
        } else {
          // get the output file and the exit code from niimath wasm
          const { blob, exitCode } = e.data;
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
      if (this.worker === null) {
        reject(new Error('Worker not initialized. Did you await the init() method?'));
      }
      this.worker.postMessage({ blob: this.file, cmd: args, outName: outName });
    });
  }
}