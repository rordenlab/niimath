// emcc -O2 -s ALLOW_MEMORY_GROWTH -s MAXIMUM_MEMORY=4GB -s TOTAL_MEMORY=268435456 -s WASM=1 -DUSING_WASM -I. core32.c nifti2_wasm.c core.c walloc.c -o funcx.js; node test_NEW.js

const fs = require('fs')
let instance = null;

let granule_size = 8;
let bits_per_byte = 8;
let bits_per_byte_log2 = 3;

function assert(c, msg) { if (!c) throw new Error(msg); }
function power_of_two(x) { return x && (x & (x - 1)) == 0; }
function assert_power_of_two(x) {
    assert(power_of_two(x), `not power of two: ${x}`);
}
function aligned(x, y) {
    assert_power_of_two(y);
    return (x & (y - 1)) == 0;
}

function assert_aligned(x, y) {
    assert(aligned(x, y), `bad alignment: ${x} % ${y}`);
}

class HeapVerifier {
    constructor(maxbytes) {
        this.maxwords = maxbytes / granule_size;
        this.state = new Uint8Array(this.maxwords / bits_per_byte);
        this.allocations = new Map;
    }
    acquire(offset, len) {
        assert_aligned(offset, granule_size);
        for (let i = 0; i < len; i += granule_size) {
            let bit = (offset + i) / granule_size;
            let byte = bit >> bits_per_byte_log2;
            let mask = 1 << (bit & (bits_per_byte - 1));
            assert((this.state[byte] & mask) == 0, "word in use");
            this.state[byte] |= mask;
        }
        this.allocations.set(offset, len);
    }
    release(offset) {
        assert(this.allocations.has(offset))
        let len = this.allocations.get(offset);
        this.allocations.delete(offset);
        for (let i = 0; i < len; i += granule_size) {
            let bit = (offset + i) / granule_size;
            let byte = bit >> bits_per_byte_log2;
            let mask = 1 << (bit & (bits_per_byte - 1));
            this.state[byte] &= ~mask;
        }
    }
}

class LinearMemory {
    constructor({initial = 2048, maximum = 2048}) {
    //https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/WebAssembly/Memory
    //shared is problematic if pointer moved (e.g. pointer not handle)
        //this.memory = new WebAssembly.Memory({ initial, maximum, shared: true });
        this.memory = new WebAssembly.Memory({ initial, maximum});
        this.verifier = new HeapVerifier(maximum * 65536);
    }
    record_malloc(ptr, len) { this.verifier.acquire(ptr, len); }
    record_free(ptr) { this.verifier.release(ptr); }
    read_string(offset) {
        let view = new Uint8Array(this.memory.buffer);
        let bytes = []
        for (let byte = view[offset]; byte; byte = view[++offset])
            bytes.push(byte);
        return String.fromCharCode(...bytes);
    }
    log(str)      { console.log(`wasm log: ${str}`) }
    log_i(str, i) { console.log(`wasm log: ${str}: ${i}`) }
    env() {
        return {
            memory: this.memory,
            wasm_log: (off) => this.log(this.read_string(off)),
            wasm_log_i: (off, i) => this.log_i(this.read_string(off), i)
        }
    }
}

function niimath_shim_u8(instance, memory, cmd) {
  let {niimath, walloc, wfree} = instance.exports;
  let datatype = 2; //2=DT_UNSIGNED_CHAR, 4=DT_SIGNED_SHORT, 16=DT_FLOAT
  const nx = 256 //number of columns
  const ny = 256 //number of rows
  const nz = 256 //number of slices
  const nt = 1 //number of volumes
  const dx = 2 //space between columns
  const dy = 2 //space between rows
  const dz = 2 //space between slices
  const dt = 3 //time between volumes
  let nvox = nx * ny * nz * nt; 
  let ptr = walloc(nvox);
  memory.record_malloc(ptr, nvox);
  let jsimg = new Uint8Array(nvox);
  for (let i = 0; i < nvox; i++) 
      jsimg[i] = i % 4;
  console.log("");
  console.log(`dims ${nx}x${ny}x${nz}x${nt}`);
  console.log("input ", jsimg[0],",", jsimg[1],",", jsimg[2],",", jsimg[3],",", jsimg[4],",", jsimg[5],"...")
  let cimg = new Uint8Array(instance.exports.memory.buffer, ptr, nvox);
  // Copy JavaScript data in to be used by WebAssembly.
  cimg.set(jsimg);
  //copy command string
  let cptr = walloc(cmd.length + 1);
  memory.record_malloc(cptr, cmd.length + 1);
  let cmdstr = new Uint8Array(cmd.length + 1);
  for (let i = 0; i < cmd.length; i++)
    cmdstr[i] = cmd.charCodeAt(i);
  let cstr = new Uint8Array(instance.exports.memory.buffer, cptr, cmd.length + 1);
  cstr.set(cmdstr);
  //run WASM
  startTime = new Date();
  let ok = niimath(ptr, datatype, nx, ny, nz, nt, dx, dy, dz, dt, cptr)
  if (ok != 0) {
    console.log(" -> '", cmd, " generated a fatal error: ", ok);
    return
  }
  console.log("  ", Math.round((new Date()) - startTime) + "ms -> '", cmd, "' -> output ", cimg[0],",", cimg[1],",", cimg[2],",", cimg[3],",", cimg[4],",", cimg[5],"...")
  // Copy data out to JavaScript.
  // jsimg.set(cimg);
  // console.log("out ", jsimg[0],",", jsimg[1],",", jsimg[2],",", jsimg[3],",", jsimg[4],",", jsimg[5],"...")
  //Free
  memory.record_free(cptr)
  wfree(cptr)
  memory.record_free(ptr)
  wfree(ptr)
}

async function main() {
  // Load the wasm into a buffer.
  const buf = fs.readFileSync('./funcx.wasm')
  let mod = new WebAssembly.Module(buf);
  let memory = new LinearMemory({ initial: 2048, maximum: 2048 });
  let instance = new WebAssembly.Instance(mod, { env: memory.env() });
  niimath_shim_u8(instance, memory, "-mul 7")
  niimath_shim_u8(instance, memory, "-add 3 -mul 2")
  niimath_shim_u8(instance, memory, "-s 2")
  niimath_shim_u8(instance, memory, "-dehaze 5")
  niimath_shim_u8(instance, memory, "-otsu 5")
}
main().then(() => console.log('Done'))
