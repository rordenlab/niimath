const { exec } = require('child_process');
const fs = require('fs');

// Function to parse the niimath help text
function parseHelpText(helpText) {
  const lines = helpText.split('\n');
  const methodDefinitions = {};

  let currentKernel = false;

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
      const command = match[1].trim();
      const argsString = match[2].trim();
      const args = argsString.split(/\s+/).filter(arg => arg.startsWith('<') && arg.endsWith('>'));
      const key = command.replace(/^-+/, ''); // Remove leading dashes

      const helpText = line.trim();

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
      } else {
        // General case for non-kernel operations
        methodDefinitions[key] = {
          args: args.map(arg => arg.replace(/[<>]/g, '')),
          help: helpText
        };
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