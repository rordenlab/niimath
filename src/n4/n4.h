/* niimath -n4: C port of ITK/ANTs N4 bias-field correction.
 *
 *  Derived from the Insight Toolkit (Copyright 1999-2019 Insight Software Consortium,
 *  Copyright 2020-present NumFOCUS) and ANTs (Copyright the ANTs contributors).
 *
 *  Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file
 *  except in compliance with the License. You may obtain a copy of the License at
 *      https://www.apache.org/licenses/LICENSE-2.0.txt  (also LICENSE in this directory)
 *  Unless required by applicable law or agreed to in writing, software distributed under the
 *  License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND,
 *  either express or implied. See the License for the specific language governing
 *  permissions and limitations under the License. See NOTICE in this directory.
 */
#ifndef NIIMATH_N4_H
#define NIIMATH_N4_H

#include <stdint.h>

#define N4_MAX_LEVELS 8
#define N4_MAX_LATTICE (1 << 24) // control points after the last refinement (64 MB of floats)
#define N4_MAX_PARTIAL (1 << 26) // floats in the fit's per-slice partial lattices (256 MB each for omega, delta)
#define N4_SPLINE_ORDER 3        // cubic, as ANTs; the B-spline code is written for this order

// Defaults are the ANTs N4BiasFieldCorrection CLI's (see plan_n4.md).
typedef struct {
	int shrink;                    // subsampling factor per axis
	int levels;                    // fitting levels
	int iterations[N4_MAX_LEVELS]; // maximum iterations per level
	double convergence;            // coefficient-of-variation threshold (0: run every iteration)
	int control_points[3];         // initial lattice size per axis (spans + N4_SPLINE_ORDER)
	int bins;                      // histogram bins for sharpening
	double fwhm, wiener;           // bias-field FWHM (histogram units) and Wiener filter noise
	double spline_distance;        // >0: mm between knots; pads the image to whole spans (ANTs -b d)
	int rescale;                   // with a mask: rescale inside it to the input's range (ANTs -r)
} N4Params;

void n4_default_params(N4Params *p);

// Corrects data (nx*ny*nz floats) in place, as ANTs N4BiasFieldCorrection does. mask and confidence
// may be NULL (all ones); a voxel takes part in the estimate iff mask != 0 and confidence > 0. With
// a mask (ANTs -x), voxels outside it keep the input and, if p->rescale, voxels inside are rescaled
// to the input's range there. field_out (optional, nx*ny*nz) receives the full-resolution log bias
// field. Returns 0 on success.
int n4_correct(float *data, const uint8_t *mask, const float *confidence, int nx, int ny, int nz,
               float dx, float dy, float dz, const N4Params *p, float *field_out);

#endif // NIIMATH_N4_H
