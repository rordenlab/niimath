# Prior art and patent notes: MCPC-3D-S in niimath's `--medic`

Raised during development: [korbinian90/ASPIRE](https://github.com/korbinian90/ASPIRE) carries a patent notice and links to [US10605885B2](https://patents.google.com/patent/US10605885B2/en). Does that patent read on the MCPC-3D-S phase-offset step in `src/medic.c`?

**This is a factual reading of the published claim, not legal advice.** Claim scope beyond claim 1, continuations, and jurisdiction are questions for counsel.

## The patent

| field | value |
| --- | --- |
| number | US10605885B2 |
| title | Phase offset determination in magnetic resonance imaging |
| assignee | Medizinische Universität Wien (Medical University of Vienna) |
| priority | 2016-08-08 |
| filed | 2017-06-09 |
| published | 2020-03-31 |
| status | active, anticipated expiry 2037-06-09 |

Neither "MCPC-3D-S" nor "ASPIRE" appears in the patent text. The ASPIRE README's notice reads: *"There is a patent on ASPIRE and a license is required for commercial use."* — it names **ASPIRE**, not MCPC-3D-S.

### Independent claim 1 (verbatim)

> A method for determining phase offsets in a complex-valued image acquired with a receiver coil at an echo time following an excitation by a transmitter coil in Magnetic Resonance Imaging, each pixel of said image representing a volume element of a 3-dimensional object, comprising: immobilising the object and acquiring a first image of the object at a predetermined first echo time, the first image being separated into a first magnitude image and a first phase image, and a second image of the object at a predetermined second echo time, the second image being separated into a second magnitude image and a second phase image, **wherein a ratio between said first echo time and said second echo time is chosen to be n:(n+1), n being a positive integer**; generating, pixel by pixel, a phase evolution image representing phase changes from the first phase image to the second phase image; and subtracting, pixel by pixel, **an n-fold** of the phase evolution image from the first phase image to obtain a phase offset image containing said phase offsets.

(emphasis added)

The `n:(n+1)` echo-time constraint **is** the ASPIRE invention: it forces the extrapolation factor to an integer `n`, which is what lets the phase offset be recovered **without any unwrapping**. That is the whole point of the method and the thing being protected.

## What niimath implements

`md_mcpc3ds()` in `src/medic.c`, on **coil-combined** input:

```text
hip    = m1*m2 * exp(i*(phi2 - phi1))
d_uw   = ROMEO_unwrap(angle(hip), mag=|hip|)          <-- an unwrapping step
offset = wrap( phi1 - TE1/(TE2-TE1) * d_uw )
```

## Claim-by-limitation comparison

| claim 1 limitation | niimath `--medic` |
| --- | --- |
| complex image **per receiver coil**; multi-channel combination | **coil-combined magnitude and phase only.** Coil combination and 5-D coil-channel input are explicit non-goals (plan §1 scope, §11 non-goals). No per-coil data is ever read. |
| echo times chosen so **TE₁:TE₂ = n:(n+1)**, n a positive integer | **Not required and not the case.** The demo data is 16.8 : 38.56 ms = 1 : 2.295. `n:(n+1)` ratios lie strictly between 1 and 2; 2.295 is outside that range for every integer `n`. Arbitrary strictly-increasing echo times are accepted. |
| subtract an **n-fold** (integer multiple) of the phase-evolution image | **Subtracts a non-integer factor** `TE₁/(TE₂−TE₁)` = 16.8/21.76 = **0.7721** on the demo data; in general a real number, not an integer. |
| (the purpose: obtain the offset **without unwrapping**) | **Explicitly unwraps** the phase-evolution image with ROMEO before extrapolating — the step the patented method exists to avoid. |
| "immobilising the object" | Not applicable; `--medic` re-estimates the field per frame precisely *because* the subject moves. |

Every limitation of a claim must be present for infringement. Four are absent here, and the two that define the invention — the `n:(n+1)` ratio and the integer `n`-fold subtraction — are absent structurally, not incidentally.

## Prior art

MCPC-3D-S predates the patent's 2016-08-08 priority date:

- Robinson SD, Grabner G, Witoszynskyj S, Trattnig S. **Combining phase images from multi-channel RF coils using 3D phase offset maps derived from a dual-echo scan.** *Magn Reson Med* 2011;65(6):1638–1648. doi:10.1002/mrm.22753 — the original MCPC-3D.
- Robinson SD, Bredies K, Khabipova D, Dymerska B, Marques JP, Schweser F. **An illustrated comparison of processing methods for MR phase imaging and QSM: combining array coil signals and phase unwrapping.** *NMR Biomed* 2017;30(4):e3601. doi:10.1002/nbm.3601 — describes MCPC-3D-S by that name.
- Eckstein K, Dymerska B, Bachrata B, et al. **Computationally efficient combination of multi-channel phase data from multi-echo acquisitions (ASPIRE).** *Magn Reson Med* 2018;79(6):2996–3006. doi:10.1002/mrm.26963 — ASPIRE, the patented method, presented as an *improvement over* MCPC-3D-S.

ASPIRE is positioned in its own paper as faster than MCPC-3D-S. The predecessor is not covered by a patent claiming the successor's distinguishing feature.

## Residual caveats

1. **MIT has no express patent grant.** ROMEO.jl and MriResearchTools.jl are MIT-licensed, which conveys copyright permissions only. Copying MIT code does not, by itself, convey any patent licence. (This is a general property of MIT, not specific to these projects.)
2. **Only claim 1 was analysed.** Dependent claims, other independent claims, continuations, and non-US family members were not reviewed.
3. **Research use is unencumbered either way.** The ASPIRE README states scientific use is free and requires no licence; the restriction it names is commercial use of ASPIRE.

If niimath is ever distributed into a context where this matters commercially, the check to run is a proper freedom-to-operate opinion on the full family — not this note.

## Mitigation already in the code

`--medic --phase-offset none` disables the MCPC-3D-S stage entirely; the rest of the pipeline (ROMEO unwrapping, weighted regression, temporal correction, low-rank filter, inversion, displacement) runs unchanged. Removing the stage would degrade multi-echo unwrapping quality on data with a large receive-phase offset, but it is a single flag, not a rework.

Related: [test/medic_reference_manifest.md](test/medic_reference_manifest.md) §4 records the measurement establishing the MCPC formula, and §2 the clean-room boundary observed throughout.
