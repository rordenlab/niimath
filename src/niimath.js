//Install dependencies
// npm install nifti-reader-js
//Compile funcx.wasm
// emcc -O2 -s ALLOW_MEMORY_GROWTH -s MAXIMUM_MEMORY=4GB -s TOTAL_MEMORY=268435456 -s WASM=1 -DUSING_WASM -I. core32.c nifti2_wasm.c core.c walloc.c -o funcx.js
// emcc -O2 -s ALLOW_MEMORY_GROWTH -s MAXIMUM_MEMORY=4GB -s WASM=1 -DUSING_WASM -I. core32.c nifti2_wasm.c core.c walloc.c -o funcx.js

//Test on image
// node niimath.js T1.nii -sqr sT1.nii

const fs = require('fs')
const nifti = require('nifti-reader-js')
let instance = null

let granule_size = 8
let bits_per_byte = 8
let bits_per_byte_log2 = 3

function assert(c, msg) { if (!c) throw new Error(msg); }
function power_of_two(x) { return x && (x & (x - 1)) == 0; }
function assert_power_of_two(x) {
    assert(power_of_two(x), `not power of two: ${x}`)
}
function aligned(x, y) {
    assert_power_of_two(y)
    return (x & (y - 1)) == 0
}

function assert_aligned(x, y) {
    assert(aligned(x, y), `bad alignment: ${x} % ${y}`)
}

class HeapVerifier {
    constructor(maxbytes) {
        this.maxwords = maxbytes / granule_size
        this.state = new Uint8Array(this.maxwords / bits_per_byte)
        this.allocations = new Map
    }
    acquire(offset, len) {
        assert_aligned(offset, granule_size)
        for (let i = 0; i < len; i += granule_size) {
            let bit = (offset + i) / granule_size
            let byte = bit >> bits_per_byte_log2
            let mask = 1 << (bit & (bits_per_byte - 1))
            assert((this.state[byte] & mask) == 0, "word in use")
            this.state[byte] |= mask
        }
        this.allocations.set(offset, len)
    }
    release(offset) {
        assert(this.allocations.has(offset))
        let len = this.allocations.get(offset)
        this.allocations.delete(offset)
        for (let i = 0; i < len; i += granule_size) {
            let bit = (offset + i) / granule_size
            let byte = bit >> bits_per_byte_log2
            let mask = 1 << (bit & (bits_per_byte - 1))
            this.state[byte] &= ~mask
        }
    }
}

class LinearMemory {
    constructor({initial = 2048, maximum = 2048}) {
    //https://developer.mozilla.org/en-US/docs/Web/JavaScript/Reference/Global_Objects/WebAssembly/Memory
    //shared is problematic if pointer moved (e.g. pointer not handle)
        //this.memory = new WebAssembly.Memory({ initial, maximum, shared: true })
        this.memory = new WebAssembly.Memory({ initial, maximum})
        this.verifier = new HeapVerifier(maximum * 65536)
    }
    record_malloc(ptr, len) { this.verifier.acquire(ptr, len); }
    record_free(ptr) { this.verifier.release(ptr); }
    read_string(offset) {
        let view = new Uint8Array(this.memory.buffer)
        let bytes = []
        for (let byte = view[offset]; byte; byte = view[++offset])
            bytes.push(byte)
        return String.fromCharCode(...bytes)
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

function niimath_shim(instance, memory, hdr, img8, cmd) {
    datatype = hdr.datatypeCode
    if ((datatype != 2) && (datatype != 4) && (datatype != 16)) {
      console.log('Only dataypes 2,4,16 supported')
      return
    }
    let {niimath, walloc, wfree} = instance.exports
    let nx = hdr.dims[1] //number of columns
    let ny = hdr.dims[2] //number of rows
    let nz = hdr.dims[3] //number of slices
    let nt = hdr.dims[4] //number of volumes
    let dx = hdr.pixDims[1] //space between columns
    let dy = hdr.pixDims[2] //space between rows
    let dz = hdr.pixDims[3] //space between slices
    let dt = hdr.pixDims[4] //time between volumes
    let bpv = Math.floor(hdr.numBitsPerVoxel/8) 
    console.log(`dims: ${nx}x${ny}x${nz}x${nt} nbyper: ${bpv}`)
    //allocate WASM null-terminated command string
    let cptr = walloc(cmd.length + 1)
    memory.record_malloc(cptr, cmd.length + 1)
    let cmdstr = new Uint8Array(cmd.length + 1)
    for (let i = 0; i < cmd.length; i++)
        cmdstr[i] = cmd.charCodeAt(i)
    let cstr = new Uint8Array(instance.exports.memory.buffer, cptr, cmd.length + 1)
    cstr.set(cmdstr)
    //allocate WASM image data
    let nvox = nx * ny * nz * nt; 
    let ptr = walloc(nvox * bpv)
    memory.record_malloc(ptr, nvox * bpv)
    const isSharedMem = true
    let cimg
    if (isSharedMem) {
        //https://stackoverflow.com/questions/54112373/how-to-cast-an-arraybuffer-to-a-sharedarraybuffer-in-javascript
        cimg = new Uint8Array(new SharedArrayBuffer(instance.exports.memory.buffer, ptr, nvox * bpv))
        //cimg.set(img8) // <- does not work???
        for (let i = 0; i < nvox * bpv; ++i)
            cimg[i] = img8[i]
    } else {
        cimg = new Uint8Array(instance.exports.memory.buffer, ptr, nvox * bpv)
        cimg.set(img8)
    }
    //run WASM
    startTime = new Date()
    let ok = niimath(ptr, datatype, nx, ny, nz, nt, dx, dy, dz, dt, cptr)
    if (ok != 0) {
        console.log(" -> '", cmd, " generated a fatal error: ", ok)
        return
    }
    console.log("'", cmd, "' required ", Math.round((new Date()) - startTime), "ms")
    //copy WASM image data to JavaScript image data
    img8.set(cimg)
    //free WASM memory
    memory.record_free(cptr)
    wfree(cptr)
    memory.record_free(ptr)
    wfree(ptr)
    return
}

async function main() {
    // Load the wasm into a buffer.
    const buf = fs.readFileSync('./funcx.wasm')
    let mod = new WebAssembly.Module(buf)
    let memory = new LinearMemory({ initial: 2048, maximum: 2048 })
    let instance = new WebAssembly.Instance(mod, { env: memory.env() })
    //check inputs
    var argv = process.argv.slice(2)
    argc = argv.length
    if (argc < 3) {
        console.log("At least three arguments required: 'node niimath.js T1 -sqr sT1.nii'")
        return
    }
    //check input filename
    fnm = argv[0]
    if (!fs.existsSync(fnm)) {
        ifnm = fnm
        fnm += '.nii'
        if (!fs.existsSync(fnm)) {
            fnm += '.gz'
            if (!fs.existsSync(fnm)) {
                console.log("Unable to find NIfTI: "+ifnm)
                return
            }
        }
    }
    //parse command string, e.g. '-dog 2 3.2'
    cmd = ''
    for (let i = 1; i < (argc -1); i++)
        cmd += argv[i]+' '
    //load image data
    var data = nifti.Utils.toArrayBuffer(fs.readFileSync(fnm))
    if (nifti.isCompressed(data)) {
        data = nifti.decompress(data)
    }
    if (!nifti.isNIFTI(data)) {
        console.log("not NIfTI ", fnm)
        return
    }
    //load header data
    hdr = nifti.readHeader(data)
    if (hdr.vox_offset < 352) { //minimal size for attached NIfTI1 header
        console.log("Only supports attached header (.nii), not detached headers (.hdr;.img)")
        return
    }
    let img8 = new Uint8Array(nifti.readImage(hdr, data))
    //process data
    niimath_shim(instance, memory, hdr, img8, cmd)
    //write output
    ofnm = argv[argc-1]
    //get file extension
    var re = /(?:\.([^.]+))?$/
    let ext = re.exec(ofnm)[1]
    if  ((ext !== undefined) && (ext.toLowerCase() === "gz")) {
        console.log('gz output not (yet) supported. See https://github.com/rii-mango/NIFTI-Reader-JS/issues/4')
        return
    }
    if ((ext === undefined) || (ext.toLowerCase() !== "nii"))
        ofnm += '.nii'
    var ohdr = new Uint8Array(data,0,hdr.vox_offset)
    //var odata = new Uint8Array(ohdr.length + img8.length)
    var odata = new Uint8Array(ohdr.length + img8.length)
    odata.set(ohdr)
    odata.set(img8, ohdr.length)
    fs.writeFile(ofnm, Buffer.from(odata), "binary", function(err) {
        if(err) {
            console.log(err)
        }
    })
}
main().then(() => console.log('Done'))
