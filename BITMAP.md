# niimath bitmap functions

This document describes the `-bitmap` subcommand and its arguments for the `niimath` tool (nim → PNG bitmap exporter).  It documents the inputs, options, and examples.  This function takes inspiration from FSL's slicer and AFNI's snapshot_volreg. Whereas slicer does not reorient images, niimath's bitmap function always rotates images a common orientation. 

---

## Basic usage

```
niimath <input_nifti> [optional niimath commands] -bitmap [bitmap options] <output.png>
```

Examples:

```
niimath PSR.nii -bitmap -s 1.5 -n -f -u -m out.png
niimath RAS.nii.gz -dehaze -5 -dog 2 3.2 -bitmap -x 0.5 -y 0.5 -z 0.5 dog.png
niimath AIR.nii.gz -bitmap -m -r -a outer.png
niimath ~/Desktop/fslmean -bitmap ~/Desktop/fslt -m -T 4 7  t4_7.png
niimath ~/Desktop/fslmean -bitmap ~/Desktop/fslt -m -e  edge.png 
niimath ~/Desktop/fslmean -bitmap ~/Desktop/fslt -x 0.25 -x 0.5 -r -y 0.5 -z 0.5  bmp.png
niimath fslmean -bitmap fslt  -c gray -T 2 7 -s 3 -C redyellow 0.5 -N -a redyell.png
```

- The last argument following `-bitmap` is the output PNG filename.
- `-bitmap` may accept an optional second input filename (overlay). If provided, that image is read and used as an overlay.

---

## Tile / layout options (specify slices to include)

These options build a tiled layout of slices (sagittal/coronal/axial).

- `-x <val> [<val> ...]` : sagittal slices. `<val>` may be:
  - Fraction `0..1` (e.g. `0.5` = middle) — computed as `round(frac * (nx-1))`
  - Negative integer `-N` to request absolute slice number `N` (FSL style)
  - Multiple values allowed after a single `-x`
- `-y <val> [<val> ...]` : coronal slices (same semantics as `-x`)
- `-z <val> [<val> ...]` : axial slices (same semantics as `-x`)
- `X`, `Y`, `Z` act like their lowercase counterparts but show white cross lines marking the position of other tiles.
- `-o <val> [<val> ...]` : optimal (largest) slices (same semantics as `-x`)
- `-r` : row separator — forces following tiles onto a new row
- `-a` : alias — adds three slices (middle of x,y,z): equivalent to `-x 0.5 -y 0.5 -z 0.5`
- `-m` : mosaic alias — builds a 3×3 mosaic:
  equivalent to `-x 0.25 0.5 0.75 -n -y 0.25 0.5 0.75 -n -z 0.25 0.5 0.75`
  (where `-n` here is a separator within the alias expansion; note `-n` may be used elsewhere)

---

## Global image options (apply to entire output)

These are image-wide and may appear anywhere in the `-bitmap` argument list.

- `-f [0|1]` : radiological orientation (default `1`). If omitted the shorthand `-f` sets to `0` (neurological).  
  `isRadiological = 1` means radiological (left appears on image right).
- `-u [0|1]` : draw left/right labels (default `1`). Omit argument to set to `0`.
- `-n [0|1]` : interpolation order: `0 = nearest`, `1 = linear` (default `1`). If alone, `-n` sets `0`.
- `-N [0|1]` : turn negative colormap on or off.
- `-s <scale>` : output scale factor (float > 0). Example `-s 1.5`.
- `-t <min> <max>` : set base image calibration range (cal_min, cal_max).
- `-T <min> <max>` : set overlay image calibration range (cal_min2, cal_max2).
- `-e` : enable edge rendering (you can supply a threshold elsewhere); presence of `-e` in later args sets `isEdge = 1`.

---

## Color / LUT options

You can either pass a numeric RGBA quadruple (R G B A, values `0..1`) or a named LUT with optional alpha.

Available named LUTs:
- `none`, `gray`, `red`, `green`, `blue`, `cyan`, `yellow`, `bluegreen`, `redyellow`, `viridis`, `magma`, `inferno`, `plasma`

**Base (tint) color / LUT**
- `-c R G B A` : numeric base tint (applies multiplicatively to grayscale intensities)
- `-c <name> [alpha]` : use named `LUT` for base; `alpha` optional (0..1). If `alpha` omitted, `1.0` is assumed.

Values are stored in `opts->RGBA` and `opts->LUT`.

**Overlay (second) color / LUT**
- `-C R G B A` or `-C <name> [alpha]` : sets overlay LUT / RGBA (stored in `opts->RGBA2` and `opts->LUT2`)

**Background color**
- `-b R G B A` : set the background fill color for empty pixels (values `0..1`) — stored in `opts->RGBA2` if needed.

Notes:
- If you choose named LUTs, the implementation builds a 256-entry LUT by linearly interpolating RGB from a low endpoint to a high endpoint (alpha is constant per LUT as specified).
- Example: `-c red 0.5` will set the base LUT to red with alpha 0.5; `-C blue` sets overlay LUT2 to blue.

---

## Calibration and tint behavior

- `cal_min`/`cal_max` (`-t`) map voxel intensity linearly to `[0..255]` for LUT lookup or grayscale conversion.
- Overlay `cal_min2`/`cal_max2` (`-T`) is used for the overlay image if one is provided.
- Base color tint: numeric base RGBA multiplies the grayscale intensity (0..255) — the implementation builds a base LUT from `RGBA`.
- Overlay color: numeric overlay RGBA multiplies the overlay intensity and the LUT stores the alpha component for blending.

---

## Labels

- Left/right labels drawn as 5×4 pixel white marks. Labels are drawn for non-`x` axes only (i.e., coronal and axial tiles).
- Radiological option (`-f`) inverts which side is drawn as `L`/`R`.

---

## Output and scaling

- By default, the renderer writes the generated image directly to PNG.
- If `-s` (scale) is specified and not 1.0, the implementation uses a two-pass resampler (from Graphics Gems) to rescale the image:
  - `isLinear = 1` → triangle (tent) filter (bilinear-like)
  - `isLinear = 0` → box filter (nearest)

---

## Advanced Examples

```
# simple middle slices into out.png
niimath in.nii -bitmap -a out.png

# 3x3 mosaic, tinted red overlay at 50% opacity, scaled 2x
niimath in.nii -bitmap -m -C red 0.5 -s 2.0 out.png

# produce axial slices 10 and 20 using absolute indices, set base range
niimath in.nii -bitmap -z -10 -z -20 -t 0 100 out.png

# overlay image with separate calibration
niimath main.nii -bitmap overlay.nii -T 0 200 -C cyan 0.4 out.png
```

## License / Attribution

This README documents behavior of the `niimath` `-bitmap` exporter implemented by you. Portions of the scaling code are derived from Graphics Gems (resampling algorithm); `stb_image_write` is used to write PNG files.

## Alternatives

- Python and Matlab NIfTI-Image-Converter
 [nii2png](https://github.com/alexlaurence/NIfTI-Image-Converter).
- [FSL](https://fsl.fmrib.ox.ac.uk/fsl/docs/)  includes a utility `slicer` for generating images.
- [med2image](https://github.com/FNNDSC/med2image) is a Python script for generating NG/JPG bitmaps from DICOM/NIfTI volumes.