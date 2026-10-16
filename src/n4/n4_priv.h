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
#ifndef NIIMATH_N4_PRIV_H
#define NIIMATH_N4_PRIV_H

#include <stddef.h>
#include <stdint.h>
#include <math.h>
#ifndef M_PI // MSVC gates M_PI behind _USE_MATH_DEFINES
#define M_PI 3.14159265358979323846
#endif

#ifdef N4_DUMP // parity-bench builds only (-DN4_DUMP): at run time N4_DUMP=<dir> receives <name>.f32 stage arrays
void n4_dump(const char *name, const float *v, size_t n);
#else
#define n4_dump(name, v, n) ((void)0)
#endif

// In-place complex FFT of power-of-two length n (both directions unscaled, as vnl_fft_1d).
void n4_fft(double *re, double *im, int n, int inverse);

// ITK's float quotient a/b, correctly rounded even under -freciprocal-math (niimath is built
// -ffast-math): a double quotient rounded to float. The volatile stops clang narrowing it back to
// a float divide (which it may then turn into a multiply by a reciprocal); a reciprocal applied in
// double is far below float rounding.
static inline float n4_fdiv(float a, float b) {
	volatile double q = (double)a / (double)b;
	return (float)q;
}

// itk::BSplineScatteredDataPointSetToImageFilter, one level: fit a cubic B-spline lattice
// (nc[0]*nc[1]*nc[2] floats, x fastest) to the included voxels of `resid` on an m[0]*m[1]*m[2] grid
// of spacing sp, each weighted by conf (NULL: 1). Returns 0 on success.
int n4_bspline_fit(const float *resid, const uint8_t *incl, const float *conf, const int m[3],
                   const double sp[3], const int nc[3], float *phi);

// itk::BSplineControlPointImageFilter::RefineControlPointLattice for one doubling of a cubic
// lattice: nc[d] -> 2*nc[d]-3 control points per axis.
void n4_bspline_refine(const float *phi, const int nc[3], float *out);

// itk::BSplineControlPointImageFilter: evaluate the cubic lattice on an m[0]*m[1]*m[2] grid
// spanning the lattice's parametric domain (x fastest). Returns 0 on success.
int n4_bspline_eval(const float *phi, const int nc[3], const int m[3], const double sp[3], float *out);

#endif // NIIMATH_N4_PRIV_H
