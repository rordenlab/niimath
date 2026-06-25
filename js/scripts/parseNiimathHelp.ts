import { exec } from 'child_process';
import fs from 'fs';

// Define types for the parsed operators
interface SubOperation {
  args: string[];
  help: string;
}

interface OperatorDefinition {
  args: string[];
  help: string;
  subOperations?: Record<string, SubOperation>;
}

interface KernelOperatorDefinition {
  subOperations: Record<string, SubOperation>;
}

type MethodDefinitions = Record<string, OperatorDefinition | KernelOperatorDefinition>;

// Operators whose operands are *filenames* of secondary images. The browser worker
// only stages the primary input file into the Emscripten filesystem, so exposing
// these as high-level methods would guarantee a runtime failure. They remain available
// in the native CLI and via raw `callMain`; exclude them from the generated browser API
// until multi-file staging exists.
const fileOperandOps = new Set<string>([
  'allineate', 'deface', 'skullstrip',
  'reslice', 'reslice_nn', 'reslice_mask',
  'mas', 'restart'
]);

// Function to parse the niimath help text
function parseHelpText(helpText: string): MethodDefinitions {
  const lines = helpText.split('\n');
  const methodDefinitions: MethodDefinitions = {};

  let currentKernel = false;
  let currentMesh = false;
  let currentBitmap = false;

  lines.forEach((line: string) => {
    // Handle kernel operations
    if (line.includes('Kernel operations')) {
      currentKernel = true;
      methodDefinitions.kernel = { subOperations: {} } as KernelOperatorDefinition;
      return;
    }

    // Match lines that start with a dash (indicating a command)
    const match = line.match(/^\s*(-[\w]+)\s*(.*?)\s*:\s*(.*)/);

    if (match) {
      // each -mesh suboption has at least 4 spaces before the dash (DO NOT USE TABS)
      // so we can use this to determine if we are in a mesh suboption after seeing the -mesh string.
      // Store the first 4 characters of the line to determine if we are in a mesh suboption
      const leadingChars = line.substring(0, 4);
      const nSpaces = '    ';
      const command = match[1].trim();
      // Top-level operators are printed with exactly one leading space (" -op ... : ...").
      // Sub-option / continuation lines (e.g. allineate's "-warp XX (...) [default: aff]")
      // are indented far deeper and merely happen to contain a colon; they must NOT become
      // standalone operators. Mesh/bitmap sub-options (4-space indent) are handled above.
      const isTopLevel = line.startsWith(' -');
      // console log below is for debugging
      // console.log(leadingChars, command, leadingChars === nSpaces);
      const argsString = match[2].trim();
      const args = argsString.split(/\s+/).filter(arg => arg.startsWith('<') && arg.endsWith('>'));
      const key = command.replace(/^-+/, ''); // Remove leading dashes

      const helpText = match[3].trim();

      if (currentKernel) {
        // Special handling for kernel operations
        if (key === 'kernel') {
          const subOpMatch = line.match(/-\s*kernel\s*(\w+)\s*(.*?)\s*:/);
          if (subOpMatch) {
            const subOpName = subOpMatch[1].trim();
            const subArgsString = subOpMatch[2].trim();
            const subArgs = subArgsString.split(/\s+/).filter(arg => arg.startsWith('<') && arg.endsWith('>'));
            const kernelDef = methodDefinitions.kernel as KernelOperatorDefinition;
            kernelDef.subOperations[subOpName] = {
              args: subArgs.map(arg => arg.replace(/[<>]/g, '')),
              help: helpText
            };
          }
        }
      } else if (command === '-mesh') {
        // Special handling for the mesh option
        currentMesh = true;
        currentBitmap = false;
        methodDefinitions.mesh = {
          args: args.map(arg => arg.replace(/[<>]/g, '')),
          help: helpText,
          subOperations: {}
        } as OperatorDefinition;
        // check if this is a valid mesh suboption (which is indented in the help text)
      } else if (currentMesh && leadingChars === nSpaces && command.startsWith('-')) {
        // Handling sub-options of the mesh command
        const subKey = command.replace(/^-+/, ''); // Remove leading dashes
        const meshDef = methodDefinitions.mesh as OperatorDefinition;
        meshDef.subOperations![subKey] = {
          args: args.map(arg => arg.replace(/[<>]/g, '')),
          help: helpText
        };
      } else if (command === '-bitmap') {
        // Special handling for the bitmap option
        currentBitmap = true;
        currentMesh = false;
        methodDefinitions.bitmap = {
          args: args.map(arg => arg.replace(/[<>]/g, '')),
          help: helpText,
          subOperations: {}
        } as OperatorDefinition;
        // check if this is a valid bitmap suboption (which is indented in the help text)
      } else if (currentBitmap && leadingChars === nSpaces && command.startsWith('-')) {
        // Handling sub-options of the bitmap command
        const subKey = command.replace(/^-+/, ''); // Remove leading dashes
        const bitmapDef = methodDefinitions.bitmap as OperatorDefinition;
        bitmapDef.subOperations![subKey] = {
          args: args.map(arg => arg.replace(/[<>]/g, '')),
          help: helpText
        };
      } else if (isTopLevel && !fileOperandOps.has(key)) {
        // General case for non-kernel, non-mesh, and non-bitmap top-level operations.
        // File-operand operators are skipped (the browser worker cannot stage their files).
        methodDefinitions[key] = {
          args: args.map(arg => arg.replace(/[<>]/g, '')),
          help: helpText
        } as OperatorDefinition;
        currentMesh = false; // Stop handling mesh sub-options if another main option is encountered
        currentBitmap = false; // Stop handling bitmap sub-options if another main option is encountered
      } else if (isTopLevel) {
        // recognized top-level op, but excluded from the browser API (file operand)
        currentMesh = false;
        currentBitmap = false;
      }
      // deeper-indented lines that are not recognized sub-options are ignored
    }

    // Reset kernel mode when moving past kernel-related lines
    if (line.trim() === '' && currentKernel) {
      currentKernel = false;
    }
  });

  return methodDefinitions;
}

// Function to execute niimath and parse its output
function generateMethodDefinitions(): Promise<MethodDefinitions> {
  return new Promise((resolve, reject) => {
    exec('../src/niimath', (error, stdout, stderr) => {
      if (error) {
        reject(new Error(`Error executing niimath: ${error.message}`));
        return;
      }

      if (stderr) {
        console.error(`stderr: ${stderr}`);
      }

      const methodDefinitions = parseHelpText(stdout);
      resolve(methodDefinitions);
    });
  });
}

// Main function to generate and save the method definitions
async function main(): Promise<void> {
  try {
    const methodDefinitions = await generateMethodDefinitions();
    const jsonString = JSON.stringify(methodDefinitions, null, 2);
    fs.writeFileSync('./src/niimathOperators.json', jsonString);
    console.log('niimathOperators saved to niimathOperators.json');
  } catch (error) {
    console.error(`Failed to generate niimath operators: ${error}`);
  }
}

main();