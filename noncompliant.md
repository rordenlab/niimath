# Archived non-compliant separable box filters

## Status

This implementation was removed from `src/coreFLT.c`. It is retained here for possible future work and must not be copied back without resolving the non-finite semantics described below.

The code accelerated uniform rectangular box kernels:

- `-fmean` and `-fmeanu`: three separable running-sum passes.
- `-dilF` and `-eroF`: separable monotonic-deque min/max passes.
- `-ero`: separable maximum over a zero-indicator mask.

These operations are uncommon. Keeping both a fast finite-only algorithm and a reference fallback would add detection, dispatch, testing, and maintenance complexity that is not justified by their usage.

## Compliance failure

FSL 6.0.7.21 preserves NaN, positive infinity, and negative infinity on input. For a central non-finite voxel and a `boxv 3` kernel:

- `-fmean` and `-fmeanu` propagate the value only to the local 3x3x3 neighborhood.
- `-dilF` ignores NaN and negative infinity, but propagates positive infinity locally.
- `-eroF` ignores NaN and positive infinity, but propagates negative infinity locally.
- `-nan` converts all three non-finite forms to zero; `-nanm` marks all three.

The running sums below cannot recover after NaN enters an accumulator. Removing infinity from a sliding sum computes `Inf - Inf`, producing NaN. One central NaN therefore contaminated 125 voxels instead of 27; one infinity produced 27 signed infinities plus 98 NaNs. The monotonic deques also changed NaN handling and location.

The implementation was correct for finite images and delivered large speedups, but it was not FSL-compatible for valid floating-point NIfTI data. The original local gather implementation has been restored. This removes the widened/non-local contamination and restores historical niimath behavior. The historical `-dilF`/`-eroF` gather still retains one NaN where FSL ignores it; that pre-existing rare-case difference was not expanded as part of this rollback.

The current compliant optimization does not reuse sums, extrema, or erosion state across voxels. It removes redundant bounds checks only for proven-interior voxels and evaluates eight adjacent outputs independently for SIMD. Each output still visits its own kernel taps in the original order, so NaN and infinity retain the historical local behavior for `-fmean`, `-fmeanu`, `-dilF`, `-eroF`, and `-ero`.

## Archived functions

### `kernel_uniform_box`

```c
staticx int kernel_uniform_box(const int *kernel, int nkernel, int nx, int nxy,
								int *pxlo, int *pxhi, int *pylo, int *pyhi, int *pzlo, int *pzhi) {
	if (nkernel < 1)
		return 0;
	int w0 = kernel[nkernel * 3];
	int xlo = INT_MAX, xhi = INT_MIN, ylo = INT_MAX, yhi = INT_MIN, zlo = INT_MAX, zhi = INT_MIN;
	for (int k = 0; k < nkernel; k++) {
		if (kernel[k + nkernel * 3] != w0)
			return 0;
		int dx = kernel[k + nkernel];
		int dy = kernel[k + nkernel + nkernel];
		int dz = (kernel[k] - dy * nx - dx) / nxy;
		if (dx < xlo) xlo = dx;
		if (dx > xhi) xhi = dx;
		if (dy < ylo) ylo = dy;
		if (dy > yhi) yhi = dy;
		if (dz < zlo) zlo = dz;
		if (dz > zhi) zhi = dz;
	}
	long long prod = (long long)(xhi - xlo + 1) * (yhi - ylo + 1) * (zhi - zlo + 1);
	if (prod != nkernel)
		return 0;
	int gx = xhi - xlo + 1, gy = yhi - ylo + 1;
	unsigned char *seen = (unsigned char *)calloc((size_t)prod, 1);
	if (!seen)
		return 0;
	int ok = 1;
	for (int k = 0; k < nkernel; k++) {
		int dx = kernel[k + nkernel];
		int dy = kernel[k + nkernel + nkernel];
		int dz = (kernel[k] - dy * nx - dx) / nxy;
		size_t idx = ((size_t)(dz - zlo) * gy + (dy - ylo)) * gx + (dx - xlo);
		if (seen[idx]) { ok = 0; break; }
		seen[idx] = 1;
	}
	free(seen);
	if (!ok)
		return 0;
	*pxlo = xlo; *pxhi = xhi; *pylo = ylo; *pyhi = yhi; *pzlo = zlo; *pzhi = zhi;
	return 1;
}
```

### `boxsum_x_fd`

```c
staticx void boxsum_x_fd(const flt *src, double *dst, int nx, int ny, int nz, int lo, int hi) {
	size_t nxy = (size_t)nx * ny;
	for (int z = 0; z < nz; z++)
		for (int y = 0; y < ny; y++) {
			const flt *s = src + (size_t)z * nxy + (size_t)y * nx;
			double *d = dst + (size_t)z * nxy + (size_t)y * nx;
			double S = 0.0;
			int a = MAX(0, lo), b = MIN(nx - 1, hi);
			for (int j = a; j <= b; j++)
				S += s[j];
			d[0] = S;
			for (int x = 1; x < nx; x++) {
				int rem = (x - 1) + lo;
				if (rem >= 0)
					S -= s[rem];
				int add = x + hi;
				if (add < nx)
					S += s[add];
				d[x] = S;
			}
		}
}
```

### `boxsum_y_dd`

```c
staticx int boxsum_y_dd(const double *src, double *dst, int nx, int ny, int nz, int lo, int hi) {
	size_t nxy = (size_t)nx * ny;
	double *R = (double *)malloc((size_t)nx * sizeof(double));
	if (!R) return 1;
	for (int z = 0; z < nz; z++) {
		const double *base = src + (size_t)z * nxy;
		double *dbase = dst + (size_t)z * nxy;
		for (int x = 0; x < nx; x++)
			R[x] = 0.0;
		int a = MAX(0, lo), b = MIN(ny - 1, hi);
		for (int yy = a; yy <= b; yy++) {
			const double *row = base + (size_t)yy * nx;
			for (int x = 0; x < nx; x++)
				R[x] += row[x];
		}
		for (int x = 0; x < nx; x++)
			dbase[x] = R[x];
		for (int y = 1; y < ny; y++) {
			int rem = (y - 1) + lo;
			if (rem >= 0) {
				const double *row = base + (size_t)rem * nx;
				for (int x = 0; x < nx; x++)
					R[x] -= row[x];
			}
			int add = y + hi;
			if (add < ny) {
				const double *row = base + (size_t)add * nx;
				for (int x = 0; x < nx; x++)
					R[x] += row[x];
			}
			double *drow = dbase + (size_t)y * nx;
			for (int x = 0; x < nx; x++)
				drow[x] = R[x];
		}
	}
	free(R);
	return 0;
}
```

### `boxsum_z_dd`

```c
staticx int boxsum_z_dd(const double *src, double *dst, int nx, int ny, int nz, int lo, int hi) {
	size_t nxy = (size_t)nx * ny;
	double *R = (double *)malloc(nxy * sizeof(double));
	if (!R) return 1;
	for (size_t p = 0; p < nxy; p++)
		R[p] = 0.0;
	int a = MAX(0, lo), b = MIN(nz - 1, hi);
	for (int zz = a; zz <= b; zz++) {
		const double *sl = src + (size_t)zz * nxy;
		for (size_t p = 0; p < nxy; p++)
			R[p] += sl[p];
	}
	for (size_t p = 0; p < nxy; p++)
		dst[p] = R[p];
	for (int z = 1; z < nz; z++) {
		int rem = (z - 1) + lo;
		if (rem >= 0) {
			const double *sl = src + (size_t)rem * nxy;
			for (size_t p = 0; p < nxy; p++)
				R[p] -= sl[p];
		}
		int add = z + hi;
		if (add < nz) {
			const double *sl = src + (size_t)add * nxy;
			for (size_t p = 0; p < nxy; p++)
				R[p] += sl[p];
		}
		double *d = dst + (size_t)z * nxy;
		for (size_t p = 0; p < nxy; p++)
			d[p] = R[p];
	}
	free(R);
	return 0;
}
```

### `box_sum_separable`

```c
staticx double *box_sum_separable(int nx, int ny, int nz, const flt *inf32,
								  int xlo, int xhi, int ylo, int yhi, int zlo, int zhi) {
	size_t nvox = (size_t)nx * ny * (size_t)nz;
	double *dA = (double *)malloc(nvox * sizeof(double));
	double *dB = (double *)malloc(nvox * sizeof(double));
	if (!dA || !dB) { free(dA); free(dB); return NULL; }
	boxsum_x_fd(inf32, dA, nx, ny, nz, xlo, xhi);
	if (boxsum_y_dd(dA, dB, nx, ny, nz, ylo, yhi) ||
		boxsum_z_dd(dB, dA, nx, ny, nz, zlo, zhi)) {
		free(dA); free(dB); return NULL;
	}
	free(dB);
	return dA;
}
```

### `box_mean_separable`

```c
staticx int box_mean_separable(int nx, int ny, int nz, flt *f32, const flt *inf32,
							 int xlo, int xhi, int ylo, int yhi, int zlo, int zhi) {
	size_t nxy = (size_t)nx * ny;
	double *sum = box_sum_separable(nx, ny, nz, inf32, xlo, xhi, ylo, yhi, zlo, zhi);
	if (!sum) return 1;
	int *cx = (int *)malloc((size_t)nx * sizeof(int));
	int *cy = (int *)malloc((size_t)ny * sizeof(int));
	int *cz = (int *)malloc((size_t)nz * sizeof(int));
	if (!cx || !cy || !cz) { free(cx); free(cy); free(cz); free(sum); return 1; }
	for (int x = 0; x < nx; x++)
		cx[x] = MIN(nx - 1, x + xhi) - MAX(0, x + xlo) + 1;
	for (int y = 0; y < ny; y++)
		cy[y] = MIN(ny - 1, y + yhi) - MAX(0, y + ylo) + 1;
	for (int z = 0; z < nz; z++)
		cz[z] = MIN(nz - 1, z + zhi) - MAX(0, z + zlo) + 1;
	for (int z = 0; z < nz; z++)
		for (int y = 0; y < ny; y++) {
			const double *s = sum + (size_t)z * nxy + (size_t)y * nx;
			flt *d = f32 + (size_t)z * nxy + (size_t)y * nx;
			double cyz = (double)cy[y] * (double)cz[z];
			for (int x = 0; x < nx; x++)
				d[x] = (flt)(s[x] / (cyz * (double)cx[x]));
		}
	free(cx); free(cy); free(cz); free(sum);
	return 0;
}
```

### `box_meanu_separable`

```c
staticx int box_meanu_separable(int nx, int ny, int nz, flt *f32, const flt *inf32,
								int xlo, int xhi, int ylo, int yhi, int zlo, int zhi, double scale) {
	size_t nvox = (size_t)nx * ny * (size_t)nz;
	double *sum = box_sum_separable(nx, ny, nz, inf32, xlo, xhi, ylo, yhi, zlo, zhi);
	if (!sum) return 1;
	for (size_t p = 0; p < nvox; p++)
		f32[p] = (flt)(sum[p] * scale);
	free(sum);
	return 0;
}
```

### `minmax_axis`

```c
staticx void minmax_axis(const flt *src, flt *dst, int nx, int ny, int nz,
						  int axis, int lo, int hi, int isMax, int *dq) {
	size_t nxy = (size_t)nx * ny;
	int n; size_t stride; int oa, ob; size_t sa, sb;
	if (axis == 0) { n = nx; stride = 1;   oa = ny; sa = nx;  ob = nz; sb = nxy; }
	else if (axis == 1) { n = ny; stride = nx; oa = nx; sa = 1; ob = nz; sb = nxy; }
	else { n = nz; stride = nxy; oa = nx; sa = 1; ob = ny; sb = nx; }
	for (int b = 0; b < ob; b++)
		for (int a = 0; a < oa; a++) {
			size_t base = (size_t)b * sb + (size_t)a * sa;
			int head = 0, tail = 0, added = 0;
			for (int x = 0; x < n; x++) {
				int rlim = x + hi; if (rlim > n - 1) rlim = n - 1;
				while (added <= rlim) {
					flt v = src[base + (size_t)added * stride];
					if (isMax)
						while (tail > head && src[base + (size_t)dq[tail - 1] * stride] <= v) tail--;
					else
						while (tail > head && src[base + (size_t)dq[tail - 1] * stride] >= v) tail--;
					dq[tail++] = added; added++;
				}
				int llim = x + lo; if (llim < 0) llim = 0;
				while (head < tail && dq[head] < llim) head++;
				dst[base + (size_t)x * stride] = src[base + (size_t)dq[head] * stride];
			}
		}
}
```

### `box_minmax_separable`

```c
staticx int box_minmax_separable(int nx, int ny, int nz, flt *f32, flt *inf32,
								  int xlo, int xhi, int ylo, int yhi, int zlo, int zhi, int isMax) {
	int maxn = MAX(nx, MAX(ny, nz));
	int *dq = (int *)malloc((size_t)maxn * sizeof(int));
	if (!dq) return 1;
	minmax_axis(inf32, f32, nx, ny, nz, 0, xlo, xhi, isMax, dq);
	minmax_axis(f32, inf32, nx, ny, nz, 1, ylo, yhi, isMax, dq);
	minmax_axis(inf32, f32, nx, ny, nz, 2, zlo, zhi, isMax, dq);
	free(dq);
	return 0;
}
```

### `box_ero_separable`

```c
staticx int box_ero_separable(int nx, int ny, int nz, flt *f32, const flt *inf32,
							  int xlo, int xhi, int ylo, int yhi, int zlo, int zhi) {
	size_t nvox = (size_t)nx * ny * nz;
	int maxn = MAX(nx, MAX(ny, nz));
	int *dq = (int *)malloc((size_t)maxn * sizeof(int));
	flt *ind = (flt *)malloc(nvox * sizeof(flt));
	flt *hz = (flt *)malloc(nvox * sizeof(flt));
	if (!dq || !ind || !hz) { free(dq); free(ind); free(hz); return 1; }
	for (size_t p = 0; p < nvox; p++)
		ind[p] = (inf32[p] == 0.0f) ? 1.0f : 0.0f;
	minmax_axis(ind, hz, nx, ny, nz, 0, xlo, xhi, 1, dq);
	minmax_axis(hz, ind, nx, ny, nz, 1, ylo, yhi, 1, dq);
	minmax_axis(ind, hz, nx, ny, nz, 2, zlo, zhi, 1, dq);
	for (size_t p = 0; p < nvox; p++)
		f32[p] = (inf32[p] != 0.0f && hz[p] > 0.5f) ? 0.0f : inf32[p];
	free(dq); free(ind); free(hz);
	return 0;
}
```
