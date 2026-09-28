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
#include "n4_priv.h"

#include <math.h>

// Iterative radix-2 transform. Forward uses exp(+2*pi*i*j*k/n) and neither direction is scaled,
// matching vnl_fft_1d's fwd_transform/bwd_transform as ITK's N4 calls them.
void n4_fft(double *re, double *im, int n, int inverse) {
	for (int i = 1, j = 0; i < n; i++) { // bit-reversal permutation
		int bit = n >> 1;
		for (; j & bit; bit >>= 1) j ^= bit;
		j ^= bit;
		if (i < j) {
			double t = re[i]; re[i] = re[j]; re[j] = t;
			t = im[i]; im[i] = im[j]; im[j] = t;
		}
	}
	const double sign = inverse ? -1.0 : 1.0;
	for (int len = 2; len <= n; len <<= 1) {
		int half = len >> 1;
		for (int k = 0; k < half; k++) {
			double ang = sign * 2.0 * M_PI * k / len, wr = cos(ang), wi = sin(ang);
			for (int i = k; i < n; i += len) {
				int j = i + half;
				double tr = re[j] * wr - im[j] * wi, ti = re[j] * wi + im[j] * wr;
				re[j] = re[i] - tr; im[j] = im[i] - ti;
				re[i] += tr; im[i] += ti;
			}
		}
	}
}
