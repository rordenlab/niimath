/* Build the Node smoke bundle, run it under Node, and remove the temporary artifact.
 * Kept as a Bun orchestration script so package.json needs no POSIX shell syntax or /tmp path. */
import { spawnSync } from "node:child_process";
import { rmSync } from "node:fs";
import { tmpdir } from "node:os";
import { join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const jsRoot = resolve(fileURLToPath(new URL("..", import.meta.url)));
const output = join(tmpdir(), `niimath-wasi-node-smoke-${process.pid}.js`);
let exitCode = 1;

try {
  const build = spawnSync(
    process.execPath,
    ["build", "scripts/wasi-node-smoke.mjs", "--target=node", `--outfile=${output}`],
    { cwd: jsRoot, stdio: "inherit" },
  );
  if (build.status !== 0) {
    exitCode = build.status ?? 1;
  } else {
    const run = spawnSync("node", [output], {
      cwd: jsRoot,
      stdio: "inherit",
      env: { ...process.env, NIIMATH_ROOT: resolve(jsRoot, "..") },
    });
    if (run.error) console.error(`Node smoke failed to start: ${run.error.message}`);
    else exitCode = run.status ?? 1;
  }
} finally {
  rmSync(output, { force: true });
}
process.exit(exitCode);
