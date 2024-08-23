const { exec } = require('child_process');
const fs = require('fs');

// Function to parse the niimath help text
function parseHelpText(helpText) {
  const lines = helpText.split('\n');
  const methodDefinitions = {};

  let currentKernel = false;
  let currentMesh = false;

  lines.forEach(line => {
    // Handle kernel operations
    if (line.includes('Kernel operations')) {
      currentKernel = true;
      methodDefinitions.kernel = { subOperations: {} };
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
            methodDefinitions.kernel.subOperations[subOpName] = {
              args: subArgs.map(arg => arg.replace(/[<>]/g, '')),
              help: helpText
            };
          }
        }
      } else if (command === '-mesh') {
        // Special handling for the mesh option
        currentMesh = true;
        methodDefinitions.mesh = {
          args: args.map(arg => arg.replace(/[<>]/g, '')),
          help: helpText,
          subOperations: {}
        };
        // check if this is a valid mesh suboption (which is indented in the help text)
      } else if (currentMesh && leadingChars === nSpaces && command.startsWith('-')) {
        // Handling sub-options of the mesh command
        const subKey = command.replace(/^-+/, ''); // Remove leading dashes
        methodDefinitions.mesh.subOperations[subKey] = {
          args: args.map(arg => arg.replace(/[<>]/g, '')),
          help: helpText
        };
      } else {
        // General case for non-kernel and non-mesh operations
        methodDefinitions[key] = {
          args: args.map(arg => arg.replace(/[<>]/g, '')),
          help: helpText
        };
        currentMesh = false; // Stop handling mesh sub-options if another main option is encountered
      }
    }

    // Reset kernel mode when moving past kernel-related lines
    if (line.trim() === '' && currentKernel) {
      currentKernel = false;
    }
  });

  return methodDefinitions;
}

// Function to execute niimath and parse its output
function generateMethodDefinitions() {
  return new Promise((resolve, reject) => {
    exec('../src/niimath', (error, stdout, stderr) => {
      if (error) {
        reject(`Error executing niimath: ${error.message}`);
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
async function main() {
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