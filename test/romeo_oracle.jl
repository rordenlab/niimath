#!/usr/bin/env julia
#
# ROMEO parity oracle for niimath's `-romeo` port (plan M0).
#
# Runs the PINNED upstream Julia implementation (ROMEO.jl + MriResearchTools) and dumps every
# intermediate the C port has to reproduce, as raw little-endian binaries plus a text manifest of
# exact scalars. The C side emits the same raw binaries via the hidden `-romeo-dump <dir>` option,
# so parity is checked by byte/ULP comparison of arrays rather than by NIfTI container bytes.
#
# Usage:  julia --startup-file=no romeo_oracle.jl <ROMEO.jl-checkout> <outdir> [--skip-env-check]
#
# Always invoke through test/romeo_oracle.sh, which pins JULIA_NUM_THREADS=1 and verifies the
# environment described in romeo_plan.md §2.1.  Nothing here mutates the upstream checkout.

using Printf
using SHA

const ARGS_CHECKOUT = length(ARGS) >= 1 ? abspath(ARGS[1]) : error("usage: romeo_oracle.jl <ROMEO.jl-checkout> <outdir>")
const OUTDIR        = length(ARGS) >= 2 ? abspath(ARGS[2]) : error("usage: romeo_oracle.jl <ROMEO.jl-checkout> <outdir>")
const SKIP_ENV      = "--skip-env-check" in ARGS

import Pkg
Pkg.activate(ARGS_CHECKOUT; io = devnull)
using ROMEO, MriResearchTools, Statistics, StatsBase

const R    = ROMEO
const MRT  = MriResearchTools
const DATA = joinpath(ARGS_CHECKOUT, "romeo")

# ---------------------------------------------------------------------------------------------
# Environment pin (romeo_plan.md §2.1)
# ---------------------------------------------------------------------------------------------

const PIN_JULIA         = "1.12.3"
const PIN_ROMEO_COMMIT  = "60d83fbb69669560d227c252dcd844afbb6648e1"
const PIN_MANIFEST_SHA  = "f6718247536608251ddef803f9674656cdb20b3141e204f250b90e6a6c368065"

filesha(p) = bytes2hex(open(sha256, p))

function git_head(dir)
    try
        return strip(read(`git -C $dir rev-parse HEAD`, String))
    catch
        return "unknown"
    end
end

function check_environment()
    manifest = joinpath(ARGS_CHECKOUT, "Manifest.toml")
    head     = git_head(ARGS_CHECKOUT)
    msha     = isfile(manifest) ? filesha(manifest) : "missing"
    bad = String[]
    string(VERSION) == PIN_JULIA          || push!(bad, "julia $(VERSION) != $PIN_JULIA")
    head == PIN_ROMEO_COMMIT              || push!(bad, "ROMEO.jl HEAD $head != $PIN_ROMEO_COMMIT")
    msha == PIN_MANIFEST_SHA              || push!(bad, "Manifest.toml sha256 $msha != $PIN_MANIFEST_SHA")
    if !isempty(bad)
        msg = "oracle environment mismatch (romeo_plan.md §2.1):\n  " * join(bad, "\n  ")
        SKIP_ENV ? @warn(msg) : error(msg)
    end
    return (julia = string(VERSION), commit = head, manifest = msha)
end

# ---------------------------------------------------------------------------------------------
# Raw dump helpers.  Every array is written in Julia's column-major (== NIfTI x-fastest) order,
# which is exactly the C linear order, so the C port dumps the same bytes with a plain fwrite.
# ---------------------------------------------------------------------------------------------

const MANIFEST = String[]

note(line) = push!(MANIFEST, line)

function dumpraw(name, A::AbstractArray{T}) where {T}
    path = joinpath(OUTDIR, name)
    open(path, "w") do io
        write(io, convert(Array{T}, collect(A)))
    end
    note(@sprintf("file %-34s %-8s n=%-10d sha256=%s", name, string(T), length(A), filesha(path)))
    return path
end

dump_u8(name, A)  = dumpraw(name, UInt8.(A))
dump_f32(name, A) = dumpraw(name, Float32.(A))

# exact scalar record: decimal repr plus the raw bit pattern, so the C comparison is unambiguous
bits(x::Float64) = @sprintf("0x%016x", reinterpret(UInt64, x))
bits(x::Float32) = @sprintf("0x%08x", reinterpret(UInt32, x))

function scalar(name, x)
    if x isa AbstractFloat
        note(@sprintf("scalar %-40s %-24s %s %s", name, repr(x), string(typeof(x)), bits(x)))
    else
        note(@sprintf("scalar %-40s %-24s %s", name, repr(x), string(typeof(x))))
    end
    return x
end

# ---------------------------------------------------------------------------------------------
# robustmask, decomposed into the four observable stages (MriResearchTools masking.jl).
# The final stage is asserted equal to the upstream one-shot robustmask().
# ---------------------------------------------------------------------------------------------

function robustmask_stages(weight; factor = 1, threshold = nothing, tag = "")
    thr = threshold
    if thr === nothing
        w = MRT.sample(weight)
        q05, q15, q8, q99 = quantile.(Ref(w), (0.05, 0.15, 0.8, 0.99))
        high_intensity = mean(filter(isfinite, w[q8 .<= w .<= q99]))
        noise = mean(filter(isfinite, w[w .<= q15]))
        noise_stage = 1
        if noise > high_intensity / 10
            noise = mean(filter(isfinite, w[w .<= q05]))
            noise_stage = 2
            if noise > high_intensity / 10
                noise = 0
                noise_stage = 3
            end
        end
        thr = max(5noise, high_intensity / 5)
        scalar("$(tag)sample_len",     length(w))
        scalar("$(tag)q05",            q05)
        scalar("$(tag)q15",            q15)
        scalar("$(tag)q8",             q8)
        scalar("$(tag)q99",            q99)
        scalar("$(tag)high_intensity", high_intensity)
        scalar("$(tag)noise",          noise)
        scalar("$(tag)noise_stage",    noise_stage)
        scalar("$(tag)threshold",      thr)
    else
        scalar("$(tag)threshold_given", thr)
    end

    s1 = weight .> (thr * factor)
    boxsizes = [[5] for _ in 1:ndims(weight)]
    sm1 = gaussiansmooth3d(s1; nbox = 1, boxsizes)
    s2 = sm1 .> 0.4
    s3 = MRT.fill_holes(s2)
    boxsizes = [[3, 3] for _ in 1:ndims(weight)]
    sm2 = gaussiansmooth3d(s3; nbox = 2, boxsizes)
    s4 = sm2 .> 0.6

    if threshold === nothing
        @assert s4 == robustmask(weight; factor) "stage decomposition != robustmask()"
    else
        @assert s4 == robustmask(weight; factor, threshold) "stage decomposition != robustmask()"
    end
    return (s1 = s1, sm1 = sm1, s2 = s2, s3 = s3, sm2 = sm2, s4 = s4, threshold = thr)
end

# ---------------------------------------------------------------------------------------------
# Seed scalars: replicate the first two statements of ROMEO.getseedfunction's addseed! on a fresh
# `visited`, using ROMEO's own functions.  Nothing upstream is modified.
# ---------------------------------------------------------------------------------------------

function seed_scalars(weights, sz, tag)
    visited = zeros(UInt8, sz)
    sq = R.getseedqueue(weights)
    seed = R.findseed!(sq, weights, visited)
    scalar("$(tag)seed_index", seed)
    if seed != 0
        sw = weights[R.getedgeindex.(seed, 1:3)]
        nst = R.NBINS - div(R.NBINS - sum(sw) / 3, 2)
        scalar("$(tag)seed_w1", Int(sw[1]))
        scalar("$(tag)seed_w2", Int(sw[2]))
        scalar("$(tag)seed_w3", Int(sw[3]))
        scalar("$(tag)new_seed_thresh", nst)
    end
    return seed
end

# ---------------------------------------------------------------------------------------------
# Weight-flag selections dumped for every case (plan M2 exit criterion).
# ---------------------------------------------------------------------------------------------

flagvec(bits::String) = Bool[c == '1' for c in bits]

const WEIGHT_SELECTIONS = [
    ("romeo3", :romeo3),
    ("romeo2", :romeo2),
    ("romeo4", :romeo4),
    ("romeo6", :romeo6),
    ("f100000", flagvec("100000")),
    ("f110000", flagvec("110000")),
    ("f101000", flagvec("101000")),
    ("f000100", flagvec("000100")),
    ("f111111", flagvec("111111")),
    ("f010000", flagvec("010000")),
]

# ---------------------------------------------------------------------------------------------
# One validation case.
# ---------------------------------------------------------------------------------------------

function run_case(tag; phasefile, magfile, TEs, template = 1)
    note("")
    note("# ===== case $tag  phase=$(basename(phasefile)) mag=$(magfile === nothing ? "none" : basename(magfile)) TEs=$TEs =====")
    note(@sprintf("input %-20s sha256=%s", basename(phasefile), filesha(phasefile)))
    magfile !== nothing && note(@sprintf("input %-20s sha256=%s", basename(magfile), filesha(magfile)))

    phasevol = readphase(phasefile; mmap = false)
    phase = Float32.(phasevol)             # scaled, exactly as the app sees it
    scalar("$(tag)_phase_slope", phasevol.header.scl_slope)
    scalar("$(tag)_phase_inter", phasevol.header.scl_inter)
    dump_f32("$(tag)_phase_rescaled.f32", phase)

    mag = magfile === nothing ? nothing : Float32.(readmag(magfile; mmap = false))
    neco = size(phase, 4)
    sz3 = size(phase)[1:3]
    scalar("$(tag)_dims", collect(sz3))
    scalar("$(tag)_neco", neco)

    # ---- mask (default: robustmask on the template echo of the magnitude) ------------------
    mask = nothing
    if mag !== nothing
        template_echo = min(template, size(mag, 4))
        magt = ndims(mag) == 4 ? mag[:, :, :, template_echo] : mag
        st = robustmask_stages(magt; tag = "$(tag)_rm_")
        dump_u8("$(tag)_mask_s1_thresh.u8",   st.s1)
        dump_f32("$(tag)_mask_sm1.f32",       st.sm1)
        dump_u8("$(tag)_mask_s2_smooth1.u8",  st.s2)
        dump_u8("$(tag)_mask_s3_fill.u8",     st.s3)
        dump_f32("$(tag)_mask_sm2.f32",       st.sm2)
        dump_u8("$(tag)_mask_s4_final.u8",    st.s4)
        mask = st.s4
    end

    # ---- keyargs, mirroring RomeoApp.get_keyargs at CLI defaults ----------------------------
    key = Dict{Symbol,Any}()
    mag  !== nothing && (key[:mag]  = mag)
    mask !== nothing && (key[:mask] = mask)
    key[:TEs] = TEs
    key[:correctglobal] = false
    key[:maxseeds] = 1
    key[:merge_regions] = false
    key[:correct_regions] = false
    key[:wrap_addition] = 0.0
    key[:temporal_uncertain_unwrapping] = 0.0
    key[:individual] = false
    key[:template] = template
    # MUST be set explicitly: the CLI resolves "romeo" to romeo3 (with a magnitude) or romeo4
    # (without) in load_data_and_resolve_args!, whereas the LIBRARY default for calculateweights
    # is :romeo (flags 1-4).  Omitting it silently unwraps with a different weight set.
    key[:weights] = mag === nothing ? :romeo4 : :romeo3

    # For 4D the library reduces to the template echo before computing weights; replicate so the
    # dumped weights correspond to exactly the array the C port builds.
    wkey = Dict{Symbol,Any}(key)
    if neco > 1
        p2ref = template == 1 ? 2 : template - 1
        wkey[:phase2] = phase[:, :, :, p2ref]
        wkey[:TEs] = TEs[[template, p2ref]]
        mag !== nothing && (wkey[:mag] = mag[:, :, :, template])
        scalar("$(tag)_p2ref", p2ref)
        scalar("$(tag)_weight_TEs", wkey[:TEs])
    end
    ptemplate = neco > 1 ? phase[:, :, :, template] : phase

    # maxmag (Float64 quantile of mag .* mask) — the single most fragile scalar in the weights
    if mag !== nothing
        mm = wkey[:mag]
        mm = mask === nothing ? mm : mm .* mask
        scalar("$(tag)_maxmag", quantile(mm[isfinite.(mm)], 0.95))
    end

    # ---- weights, every selection -----------------------------------------------------------
    wref = nothing
    # NOTE: `weights` MUST come AFTER the splat. Julia's keyword splatting is rightmost-wins
    # (ROMEO relies on this itself: `unwrap!(view(...); args..., weights)`), so putting the
    # per-selection value first lets wkey[:weights] silently override every selection.
    for (wname, wsel) in WEIGHT_SELECTIONS
        w = R.calculateweights(ptemplate; wkey..., weights = wsel)
        dump_u8("$(tag)_weights_$(wname).u8", w)
        # capture the DEFAULT selection's weights for the seed scalars below -- hardcoding
        # "romeo3" here silently used the wrong weight set for every magnitude-free case.
        wname == string(mag === nothing ? :romeo4 : :romeo3) && (wref = w)
    end
    # the resolved default: `romeo` -> romeo3 with magnitude, romeo4 without
    wdefault = mag === nothing ? :romeo4 : :romeo3
    scalar("$(tag)_weights_default", wdefault)

    # pre-rescale Float64 weights (sampled edges) for the <=4 ULP check
    let flags = R.updateflags(wdefault == :romeo3 ? Bool[1, 1, 0, 1, 0, 0] : Bool[1, 1, 1, 1, 0, 0],
                              get(wkey, :phase2, nothing), get(wkey, :TEs, nothing), get(wkey, :mag, nothing))
        note("flags_active " * join(Int.(flags), ""))
        wf = R.calculateweights(ptemplate; wkey..., weights = wdefault, type = Float64, rescale = identity)
        dumpraw("$(tag)_weights_prerescale.f64", wf)
    end

    # ---- quality map on the still-WRAPPED phase, WITHOUT a mask.  This is what set_mask!'s
    #      `-k qualitymask` branch sees (data["mask"] does not exist yet at that point).
    let k0 = Dict{Symbol,Any}(key)
        delete!(k0, :mask)
        dump_f32("$(tag)_qmap_wrapped.f32", voxelquality(phase; k0...))
    end

    # ---- seed --------------------------------------------------------------------------------
    seed_scalars(wref === nothing ? R.calculateweights(ptemplate; weights = wdefault, wkey...) : wref,
                 sz3, "$(tag)_")

    # ---- unwrapping --------------------------------------------------------------------------
    regions = zeros(UInt8, sz3)
    unwrapped = copy(phase)
    unwrap!(unwrapped; key..., regions)
    dump_u8("$(tag)_visited.u8", regions)
    dump_f32("$(tag)_unwrapped.f32", unwrapped)
    scalar("$(tag)_unwrapped_extrema", collect(extrema(unwrapped)))

    # ---- quality maps, exactly as write_qualitymap does: AFTER unwrapping, on data["phase"],
    #      with the full keyargs (mask included).
    qmap = voxelquality(unwrapped; key...)
    dump_f32("$(tag)_qmap.f32", qmap)
    for i in 1:6
        flags = falses(6); flags[i] = true
        qm = voxelquality(unwrapped; key..., weights = flags)
        dump_f32("$(tag)_qmap_$(i).f32", qm)
        note("qmap_$(i)_all_one " * string(all(qm[1:end-1, 1:end-1, 1:end-1] .== 1.0f0)))
    end

    # ---- B0 (plan M8). The app's computeB0 runs on the UNWRAPPED phase, and substitutes a
    #      Float64 exp(-TE/20) T2* decay for the magnitude when none is supplied. Equivalent to
    #      `romeo --compute-B0 --phase-offset-correction off` (MCPC-3D-S is out of scope).
    for wmode in (:phase_snr, :phase_var, :average, :TEs, :mag, :simulated_mag)
        b0mag = mag === nothing ? MRT.to_dim(exp.(-TEs / 20), 4) : mag
        b0 = calculateB0_unwrapped(unwrapped, b0mag, TEs, wmode)
        snr = get_B0_snr(b0mag, TEs, wmode)
        dump_f32("$(tag)_b0_$(wmode).f32", b0)
        dump_f32("$(tag)_b0_snr_$(wmode).f32", ndims(snr) >= 4 ? dropdims(snr; dims=4) : snr)
    end

    # ---- variants required by later milestones ----------------------------------------------
    if neco > 1
        let u = copy(phase)
            unwrap!(u; key..., temporal_uncertain_unwrapping = 0.5)
            dump_f32("$(tag)_unwrapped_tuu05.f32", u)
        end
        let u = copy(phase)
            unwrap!(u; key..., individual = true)
            dump_f32("$(tag)_unwrapped_individual.f32", u)
        end
        if neco >= 2
            # The app rebuilds the mask for the SELECTED template: set_mask! uses
            # `robustmask(mag[:,:,:,min(settings["template"], size(mag,4))])`.  Reusing the
            # template-1 mask here would make this variant a run the CLI can never produce.
            let u = copy(phase), k2 = Dict{Symbol,Any}(key)
                k2[:template] = 2
                if mag !== nothing
                    te2 = min(2, size(mag, 4))
                    k2[:mask] = robustmask(ndims(mag) == 4 ? mag[:, :, :, te2] : mag)
                end
                unwrap!(u; k2...)
                dump_f32("$(tag)_unwrapped_template2.f32", u)
            end
        end
    end
    let u = copy(phase)
        unwrap!(u; key..., correctglobal = true)
        dump_f32("$(tag)_unwrapped_correctglobal.f32", u)
    end
    return nothing
end

# ---------------------------------------------------------------------------------------------
# Synthetic fixtures: tiny, non-cubic, boundary-heavy volumes that exercise the linear-index
# semantics, the half-even ties, NaN/Inf handling and the documented upstream weight cases.
# Written as .nii by the oracle so both sides read identical inputs.
# ---------------------------------------------------------------------------------------------

"deterministic, portable pseudo-random sequence (no RNG dependency across Julia versions)"
function detseq(n, seed)
    x = UInt64(seed)
    out = Vector{Float64}(undef, n)
    for i in 1:n
        x = (x * 6364136223846793005 + 1442695040888963407) % UInt64
        out[i] = ((x >> 11) % UInt64) / Float64(1 << 53)
    end
    return out
end

function synth_phase(sz, seed)
    n = prod(sz)
    v = detseq(n, seed)
    # smooth ramp + wrapped detail so unwrapping is non-trivial but deterministic
    p = Array{Float32}(undef, sz)
    idx = 0
    for k in 1:sz[3], j in 1:sz[2], i in 1:sz[1]
        idx += 1
        ramp = 0.7 * i + 0.35 * j + 0.21 * k
        p[i, j, k] = Float32(rem2pi(ramp + 0.3 * (v[idx] - 0.5), RoundNearest))
    end
    return p
end

function synth_mag(sz, seed)
    n = prod(sz)
    v = detseq(n, seed)
    m = Array{Float32}(undef, sz)
    idx = 0
    c = (sz .+ 1) ./ 2
    for k in 1:sz[3], j in 1:sz[2], i in 1:sz[1]
        idx += 1
        r = sqrt(((i - c[1]) / c[1])^2 + ((j - c[2]) / c[2])^2 + ((k - c[3]) / c[3])^2)
        base = r < 0.75 ? 1000.0 : 12.0
        m[i, j, k] = Float32(base * (0.85 + 0.3 * v[idx]))
    end
    return m
end

function write_fixtures(dir)
    mkpath(dir)
    made = String[]
    # NOTE: every fixture must have >= 20 voxels.  MriResearchTools' fill_holes passes
    # (1, length(mask)/20) to ImageMorphology.imfill, which throws a DomainError when the upper
    # bound falls below the lower one.  That is an upstream limitation, mirrored by the C port.
    for (name, sz, seed) in (("small_a", (7, 5, 3), 12345), ("small_b", (5, 9, 4), 999),
                             ("line_x", (23, 1, 1), 7), ("plane_xy", (6, 4, 1), 8))
        p = synth_phase(sz, seed)
        m = synth_mag(sz, seed + 1)
        savenii(p, joinpath(dir, "$(name)_phase.nii"))
        savenii(m, joinpath(dir, "$(name)_mag.nii"))
        push!(made, name)
    end
    # two-echo synthetic
    sz = (7, 5, 3)
    p1 = synth_phase(sz, 4242)
    p2 = Float32.(rem2pi.(Float64.(p1) .* (38.56 / 16.8), RoundNearest))
    p4 = cat(p1, p2; dims = 4)
    m1 = synth_mag(sz, 4243)
    m4 = cat(m1, m1 .* 0.8f0; dims = 4)
    savenii(p4, joinpath(dir, "me_a_phase.nii"))
    savenii(m4, joinpath(dir, "me_a_mag.nii"))
    push!(made, "me_a")
    # degenerate magnitudes: zeros, equal values, one huge outlier
    md = fill(5.0f0, sz); md[1] = 0f0; md[end] = 1f6
    savenii(synth_phase(sz, 31337), joinpath(dir, "degen_phase.nii"))
    savenii(md, joinpath(dir, "degen_mag.nii"))
    push!(made, "degen")
    return made
end

# ---------------------------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------------------------

function main()
    mkpath(OUTDIR)
    env = check_environment()
    note("# ROMEO parity oracle — regenerate with test/romeo_oracle.sh")
    note("julia_version   $(env.julia)")
    note("romeo_commit    $(env.commit)")
    note("manifest_sha256 $(env.manifest)")
    note("project_toml_sha256 $(filesha(joinpath(ARGS_CHECKOUT, "Project.toml")))")
    note("nthreads        $(Threads.nthreads())")
    note("NBINS           $(R.NBINS)")

    fixdir = joinpath(OUTDIR, "fixtures")
    made = write_fixtures(fixdir)
    note("fixtures        " * join(made, " "))

    # --- real validation cases (romeo_plan.md §6) ---
    if isdir(DATA)
        run_case("e0"; phasefile = joinpath(DATA, "phase0.nii.gz"), magfile = joinpath(DATA, "mag0.nii.gz"), TEs = [16.8])
        run_case("e1"; phasefile = joinpath(DATA, "phase1.nii.gz"), magfile = joinpath(DATA, "mag1.nii.gz"), TEs = [38.56])
        run_case("me"; phasefile = joinpath(DATA, "phase.nii"),     magfile = joinpath(DATA, "mag.nii"),     TEs = [16.8, 38.56])
        run_case("e0n"; phasefile = joinpath(DATA, "phase0.nii.gz"), magfile = nothing, TEs = [16.8])
    else
        @warn "validation data not found at $DATA — synthetic cases only"
    end

    # --- synthetic cases ---
    for name in ("small_a", "small_b", "line_x", "plane_xy", "degen")
        run_case(name; phasefile = joinpath(fixdir, "$(name)_phase.nii"),
                       magfile   = joinpath(fixdir, "$(name)_mag.nii"), TEs = [16.8])
        run_case("$(name)_nomag"; phasefile = joinpath(fixdir, "$(name)_phase.nii"),
                                  magfile = nothing, TEs = [16.8])
    end
    run_case("me_a"; phasefile = joinpath(fixdir, "me_a_phase.nii"),
                     magfile   = joinpath(fixdir, "me_a_mag.nii"), TEs = [16.8, 38.56])

    # --- scalar unit oracles: gamma / rem2pi / rescale / unwrapvoxel / quantile -------------
    note("")
    note("# ===== scalar unit oracles =====")
    gam_in = Float32[0, 1, -1, 3.1415925f0, 3.1415927f0, -3.1415925f0, -3.1415927f0,
                     3.2f0, -3.2f0, 6.2831855f0, -6.2831855f0, 1f-30, -1f-30, 1f10, NaN32, Inf32, -Inf32]
    for x in gam_in
        note(@sprintf("gamma %-16s -> %-16s", repr(x), repr(R.γ(x))))
    end
    r2p_in = Float32[0, 1, 3.1415925f0, 3.1415927f0, -3.1415927f0, 6.2831855f0, -6.2831855f0,
                     9.42477f0, 12.566371f0, 100f0, 1f5, 1f7, 1f10, -1f10, 1.5f0, -1.5f0,
                     prevfloat(Float32(pi)), nextfloat(Float32(pi)), 2f0 * Float32(pi), 4f0 * Float32(pi)]
    for x in r2p_in
        note(@sprintf("rem2pi32 %-18s -> %-18s %s", repr(x), repr(rem2pi(x, RoundNearest)), bits(rem2pi(x, RoundNearest))))
    end
    r2p64 = Float64[0, Float64(pi), -Float64(pi), nextfloat(Float64(pi)), prevfloat(Float64(pi)), 2pi, -2pi, 3pi/2, 7pi/4, 1e6, 1e7,
                    1.6e6, 1.7e6, 1e10, -1e10, 1e15, 1e20, 1.5, -1.5]
    for x in r2p64
        note(@sprintf("rem2pi64 %-24s -> %-24s %s", repr(x), repr(rem2pi(x, RoundNearest)), bits(rem2pi(x, RoundNearest))))
    end
    for w in (0.0, 1.0, 0.5, 0.5 / 255, 1.5 / 255, 2.5 / 255, 1 - 0.5/255, 1 - 1.5/255,
              1 + eps(), -eps(), NaN, Inf, -Inf, 0.9980392156862745, 0.99607843137254903)
        note(@sprintf("rescale %-24s -> %d", repr(w), R.rescale(w)))
    end
    for (nw, od) in ((1.0f0, 1.0f0), (3.0f0, -3.0f0), (0.0f0, 3.1415927f0), (0.0f0, -3.1415927f0),
                     (1f0, 1f0 + 6.2831855f0), (-1f0, 100f0), (0f0, 1f10))
        note(@sprintf("unwrapvoxel %-14s %-14s -> %-24s", repr(nw), repr(od), repr(R.unwrapvoxel(nw, od))))
    end
    # type-7 quantile oracle on a small exact vector
    let v = Float32.([3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5])
        for p in (0.0, 0.05, 0.15, 0.5, 0.8, 0.95, 0.99, 1.0)
            note(@sprintf("quantile11 p=%-6s -> %-24s %s", repr(p), repr(quantile(v, p)), bits(quantile(v, p))))
        end
    end
    # sample() index geometry for the validation volume size
    let n = 76 * 76 * 46
        len = ceil(Int, sqrt(min(1e5, n)))
        si = round.(Int, range(0, n - len; length = len))
        note(@sprintf("sample_geometry n=%d len=%d first=%d second=%d last=%d total=%d",
                      n, len, si[1], si[2], si[end], len * len))
    end

    # --- primitive parity tables. MUST stay in sync with RM_PRIM_* in src/romeo.c: the C side
    #     evaluates the identical inputs under -romeo-dump and the two raw outputs are compared
    #     byte-for-byte. This is the only coverage of the Payne-Hanek branch of rem2pi.
    prim_d = Float64[
        0.0, 1.0, -1.0, 3.141592653589793, -3.141592653589793,
        1.5707963267948966, -1.5707963267948966, 4.71238898038469, 6.283185307179586, -6.283185307179586,
        2.5, -2.5, 3.5, -3.5, 5.0, -5.0, 7.0, -7.0,
        100.0, -100.0, 1000.0, 100000.0, 1000000.0, 1600000.0,
        1650000.0, -1650000.0, 1.0e7, -1.0e7, 1.0e10, -1.0e10, 1.0e15, 1.0e20, 1.0e30, 1.0e100,
        0.5, -0.5, 1.5, 2.5000000000000004,
        5e-324, -5e-324, 1e-310, 2.2250738585072011e-308]
    prim_f = Float32[
        0.0, 1.0, -1.0, 3.1415925f0, 3.1415927f0, -3.1415925f0, -3.1415927f0,
        3.2f0, -3.2f0, 6.2831855f0, -6.2831855f0, 9.42477f0, 12.566371f0,
        1.0f-30, -1.0f-30, 100.0f0, 1.0f5, 1.0f7, 1.0f10, -1.0f10, 1.0f20, 4.0f0, -4.0f0,
        1.4012984643f-45, -1.4012984643f-45, 1.1754942107f-38]
    prim_w = Float64[
        0.0, 1.0, 0.5, 0.5 / 255, 1.5 / 255, 2.5 / 255,
        1 - 0.5 / 255, 1 - 1.5 / 255, 1.0000000000000002, -1.0e-17,
        0.9980392156862745, 0.99607843137254903, 0.25, 0.75]
    prim_uv = Tuple{Float32,Float32}[
        (1.0f0, 1.0f0), (3.0f0, -3.0f0), (0.0f0, 3.1415927f0), (0.0f0, -3.1415927f0),
        (1.0f0, 7.2831855f0), (-1.0f0, 100.0f0), (0.0f0, 1.0f10), (3.1415927f0, -3.1415927f0),
        (2.0f0, 2.0f0 + 3.1415927f0), (-2.0f0, -2.0f0 - 3.1415927f0)]

    dumpraw("prim_rem2pi64.f64", [rem2pi(x, RoundNearest) for x in prim_d])
    dumpraw("prim_rem2pi32_gamma.f32",
            vcat(Float32[rem2pi(x, RoundNearest) for x in prim_f], Float32[R.γ(x) for x in prim_f]))
    dumpraw("prim_rescale.u8", UInt8[R.rescale(w) for w in prim_w])
    # unwrapvoxel with a Float32 `old` (the pass-through branch of unwrapedge!) and with a
    # Float64 `old` (the threshold branches / the temporal path) -- DIFFERENT subtraction widths.
    dumpraw("prim_unwrapvoxel.f32",
            vcat(Float32[R.unwrapvoxel(n, o) for (n, o) in prim_uv],
                 Float32[R.unwrapvoxel(n, Float64(o)) for (n, o) in prim_uv]))

    open(joinpath(OUTDIR, "manifest.txt"), "w") do io
        for l in MANIFEST
            println(io, l)
        end
    end
    println("wrote $(length(MANIFEST)) manifest lines to $(joinpath(OUTDIR, "manifest.txt"))")
    return 0
end

main()
