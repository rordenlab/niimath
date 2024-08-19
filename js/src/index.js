export class Niimath {
  constructor() {
    this.worker = this.createWorker();
  }

  createWorker() {
    // Create a new worker. The import.meta.url is relevant for module bundlers.
    const worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
    
    // Handle errors from the worker
    worker.onmessage = function (event) {
      if (event.data.type === 'error') {
        console.error('Worker Error:', event.data.message);
        if (event.data.error) {
          console.error('Stack Trace:', event.data.error);
        }
      }
    };

    return worker;
  }

  // Create a new ImageProcessor instance with the given file.
  // This only supports a single file for now. 
  image(file) {
    return new ImageProcessor(this.worker, file);
  }
}

// ImageProcessor class to handle niimath image processing commands.
class ImageProcessor {
  constructor(worker, file) {
    this.worker = worker; // The worker to send commands to. Important to use a separate thread for processing.
    this.file = file;
    this.commands = []; // The list of commands to execute. Will be populated by ImageProcessor niimath operators
  }

  // use _ prefix to indicate a private method (although it's not enforced in JS just like in Python).
  // Could use TypeScript in the future.
  _addCommand(cmd, ...args) {
    // each item must be a string, so we force it here.
    // e.g. niimath.dog(1, 2) will add '-dog', '1', '2' to the commands array.
    this.commands.push(cmd, ...args.map(String));
    return this;
  }

  sobel() {
    return this._addCommand('-sobel');
  }

  dog(sigma1, sigma2) {
    return this._addCommand('-dog', sigma1, sigma2);
  }

  // run command. e.g. niimath.image(file).dog(2, 3.2).run('output.nii')
  async run(outName = 'output.nii') {
    return new Promise((resolve, reject) => {
      this.worker.onmessage = (e) => {
        if (e.data.type === 'error') {
          reject(new Error(e.data.message));
        } else {
          const { blob, exitCode } = e.data;
          if (exitCode === 0) {
            resolve(blob);
          } else {
            reject(new Error(`niimath processing failed with exit code ${exitCode}`));
          }
        }
      };

      const args = [this.file.name, ...this.commands, outName];
      this.worker.postMessage({ blob: this.file, cmd: args, outName: outName });
    });
  }
}