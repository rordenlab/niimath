/* filter.c - multi-component aware resampling using Graphics Gems single-channel code.
   Public-domain algorithm adapted from Dale Schumacher, "General Filtered Image Rescaling",
   Graphics Gems III (1992). This implementation wraps a robust single-channel resampler
   and adds a multi-component loop to handle interleaved RGB(A)/gray images.

   Build: compile and link this file together with your project. Keep one copy only.
*/

#include "filter.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* small macros */
#ifndef CLAMP
#define CLAMP(v, lo, hi) (((v) < (lo)) ? (lo) : (((v) > (hi)) ? (hi) : (v)))
#endif

/* create/destroy Image: allocate interleaved buffer (xsize * ysize * components) */
Image *new_image(int xsize, int ysize, int components) {
	if (xsize <= 0 || ysize <= 0 || components <= 0)
		return NULL;
	Image *image = (Image *)malloc(sizeof(Image));
	if (!image)
		return NULL;
	size_t count = (size_t)xsize * (size_t)ysize * (size_t)components;
	image->data = (Pixel *)calloc(count, sizeof(Pixel));
	if (!image->data) {
		free(image);
		return NULL;
	}
	image->xsize = xsize;
	image->ysize = ysize;
	image->components = components;
	image->span = xsize * components;
	return image;
}

void free_image(Image *image) {
	if (!image)
		return;
	free(image->data);
	free(image);
}

/* Basic access helpers operate on interleaved storage.
   They work per-component (c): 0 <= c < components.
*/
Pixel get_pixel(const Image *image, int x, int y, int c) {
	if (!image)
		return 0;
	if (x < 0 || x >= image->xsize || y < 0 || y >= image->ysize || c < 0 || c >= image->components)
		return 0;
	size_t idx = (size_t)y * (size_t)image->span + (size_t)x * (size_t)image->components + (size_t)c;
	return image->data[idx];
}

void put_pixel(Image *image, int x, int y, int c, Pixel val) {
	if (!image)
		return;
	if (x < 0 || x >= image->xsize || y < 0 || y >= image->ysize || c < 0 || c >= image->components)
		return;
	size_t idx = (size_t)y * (size_t)image->span + (size_t)x * (size_t)image->components + (size_t)c;
	image->data[idx] = val;
}

/* copy one channel out of interleaved row (row = output buffer) */
void get_row(Pixel *row, const Image *image, int y, int c) {
	if (!image || !row)
		return;
	if (y < 0 || y >= image->ysize)
		return;
	const Pixel *src = image->data + (size_t)y * (size_t)image->span + (size_t)c;
	for (int x = 0; x < image->xsize; ++x) {
		row[x] = src[x * image->components];
	}
}

/* copy one column out of interleaved image into col[] for channel c */
void get_column(Pixel *col, const Image *image, int x, int c) {
	if (!image || !col)
		return;
	if (x < 0 || x >= image->xsize)
		return;
	const Pixel *p = image->data + (size_t)x * (size_t)image->components + (size_t)c;
	int span = image->span;
	for (int y = 0; y < image->ysize; ++y, p += span) {
		col[y] = *p;
	}
}

/* Filter support functions: box and triangle (tent) */
double box_filter(double t) {
	if ((t > -0.5) && (t <= 0.5))
		return 1.0;
	return 0.0;
}
double triangle_filter(double t) {
	if (t < 0.0)
		t = -t;
	if (t < 1.0)
		return 1.0 - t;
	return 0.0;
}

/* ----------------- Single-channel zoom implementation (static) -----------------
   This is essentially the Graphics Gems two-pass implementation adapted to use
   helper functions that read/write rows/columns (get_row/put_pixel) for single
   channel buffers. We'll implement zoom_single that operates on "ImageSC",
   a tiny single-channel wrapper.
*/

typedef struct {
	int xsize;
	int ysize;
	Pixel *data;
	int span; /* bytes per row (== xsize) for single-channel buffer */
} ImageSC;

/* allocate single-channel image (xsize * ysize bytes) */
static ImageSC *new_image_sc(int xsize, int ysize) {
	ImageSC *img = (ImageSC *)malloc(sizeof(ImageSC));
	if (!img)
		return NULL;
	img->data = (Pixel *)calloc((size_t)xsize * (size_t)ysize, sizeof(Pixel));
	if (!img->data) {
		free(img);
		return NULL;
	}
	img->xsize = xsize;
	img->ysize = ysize;
	img->span = xsize;
	return img;
}
static void free_image_sc(ImageSC *img) {
	if (!img)
		return;
	free(img->data);
	free(img);
}

/* Access helpers for ImageSC (simple contiguous single-channel layout) */
static inline Pixel sc_get_pixel(const ImageSC *img, int x, int y) {
	if (!img)
		return 0;
	if (x < 0 || x >= img->xsize || y < 0 || y >= img->ysize)
		return 0;
	return img->data[y * img->span + x];
}
static void sc_get_row(Pixel *row, const ImageSC *img, int y) {
	if (!img || !row)
		return;
	memcpy(row, img->data + (size_t)y * (size_t)img->span, (size_t)img->xsize);
}
static void sc_get_column(Pixel *col, const ImageSC *img, int x) {
	if (!img || !col)
		return;
	const Pixel *p = img->data + x;
	for (int i = 0; i < img->ysize; ++i, p += img->span)
		col[i] = *p;
}
static void sc_put_pixel(ImageSC *img, int x, int y, Pixel v) {
	if (!img)
		return;
	if (x < 0 || x >= img->xsize || y < 0 || y >= img->ysize)
		return;
	img->data[y * img->span + x] = v;
}

/* Contribution lists used by single-channel zoom */
typedef struct {
	int pixel;
	double weight;
} CONTRIB;
typedef struct {
	int n;
	CONTRIB *p;
} CLIST;
static CLIST *contrib = NULL; /* allocated by zoom_single */

/* zoom_single: two-pass separable resampler for ImageSC */
static void zoom_single(ImageSC *dst, ImageSC *src, double (*filterf)(double), double fwidth) {
	if (!dst || !src || !filterf)
		return;
	ImageSC *tmp;
	double xscale = (double)dst->xsize / (double)src->xsize;
	double yscale = (double)dst->ysize / (double)src->ysize;

	tmp = new_image_sc(dst->xsize, src->ysize);
	if (!tmp)
		return;

	contrib = (CLIST *)calloc((size_t)dst->xsize, sizeof(CLIST));
	if (!contrib) {
		free_image_sc(tmp);
		return;
	}

	if (xscale < 1.0) {
		double width = fwidth / xscale;
		double fscale = 1.0 / xscale;
		for (int i = 0; i < dst->xsize; ++i) {
			contrib[i].n = 0;
			int pcount = (int)(width * 2 + 1);
			contrib[i].p = (CONTRIB *)calloc((size_t)pcount, sizeof(CONTRIB));
			double center = (double)i / xscale;
			int left = (int)ceil(center - width);
			int right = (int)floor(center + width);
			for (int j = left; j <= right; ++j) {
				double weight = center - (double)j;
				weight = (*filterf)(weight / fscale) / fscale;
				int n;
				if (j < 0)
					n = -j;
				else if (j >= src->xsize)
					n = (src->xsize - j) + src->xsize - 1;
				else
					n = j;
				int k = contrib[i].n++;
				contrib[i].p[k].pixel = n;
				contrib[i].p[k].weight = weight;
			}
		}
	} else {
		for (int i = 0; i < dst->xsize; ++i) {
			contrib[i].n = 0;
			int pcount = (int)(fwidth * 2 + 1);
			contrib[i].p = (CONTRIB *)calloc((size_t)pcount, sizeof(CONTRIB));
			double center = (double)i / xscale;
			int left = (int)ceil(center - fwidth);
			int right = (int)floor(center + fwidth);
			for (int j = left; j <= right; ++j) {
				double weight = center - (double)j;
				weight = (*filterf)(weight);
				int n;
				if (j < 0)
					n = -j;
				else if (j >= src->xsize)
					n = (src->xsize - j) + src->xsize - 1;
				else
					n = j;
				int k = contrib[i].n++;
				contrib[i].p[k].pixel = n;
				contrib[i].p[k].weight = weight;
			}
		}
	}

	/* horizontal pass: src -> tmp */
	Pixel *raster = (Pixel *)calloc((size_t)src->xsize, sizeof(Pixel));
	if (!raster) {
		for (int i = 0; i < dst->xsize; ++i)
			free(contrib[i].p);
		free(contrib);
		free_image_sc(tmp);
		return;
	}
	for (int k = 0; k < tmp->ysize; ++k) {
		sc_get_row(raster, src, k);
		for (int i = 0; i < tmp->xsize; ++i) {
			double weight = 0.0;
			for (int j = 0; j < contrib[i].n; ++j) {
				weight += raster[contrib[i].p[j].pixel] * contrib[i].p[j].weight;
			}
			sc_put_pixel(tmp, i, k, (Pixel)CLAMP(weight, 0, 255));
		}
	}
	free(raster);

	for (int i = 0; i < tmp->xsize; ++i)
		free(contrib[i].p);
	free(contrib);

	/* vertical contribution lists */
	contrib = (CLIST *)calloc((size_t)dst->ysize, sizeof(CLIST));
	if (!contrib) {
		free_image_sc(tmp);
		return;
	}

	if (yscale < 1.0) {
		double width = fwidth / yscale;
		double fscale = 1.0 / yscale;
		for (int i = 0; i < dst->ysize; ++i) {
			contrib[i].n = 0;
			int pcount = (int)(width * 2 + 1);
			contrib[i].p = (CONTRIB *)calloc((size_t)pcount, sizeof(CONTRIB));
			double center = (double)i / yscale;
			int left = (int)ceil(center - width);
			int right = (int)floor(center + width);
			for (int j = left; j <= right; ++j) {
				double weight = center - (double)j;
				weight = (*filterf)(weight / fscale) / fscale;
				int n;
				if (j < 0)
					n = -j;
				else if (j >= tmp->ysize)
					n = (tmp->ysize - j) + tmp->ysize - 1;
				else
					n = j;
				int k = contrib[i].n++;
				contrib[i].p[k].pixel = n;
				contrib[i].p[k].weight = weight;
			}
		}
	} else {
		for (int i = 0; i < dst->ysize; ++i) {
			contrib[i].n = 0;
			int pcount = (int)(fwidth * 2 + 1);
			contrib[i].p = (CONTRIB *)calloc((size_t)pcount, sizeof(CONTRIB));
			double center = (double)i / yscale;
			int left = (int)ceil(center - fwidth);
			int right = (int)floor(center + fwidth);
			for (int j = left; j <= right; ++j) {
				double weight = center - (double)j;
				weight = (*filterf)(weight);
				int n;
				if (j < 0)
					n = -j;
				else if (j >= tmp->ysize)
					n = (tmp->ysize - j) + tmp->ysize - 1;
				else
					n = j;
				int k = contrib[i].n++;
				contrib[i].p[k].pixel = n;
				contrib[i].p[k].weight = weight;
			}
		}
	}

	/* vertical pass: tmp -> dst */
	raster = (Pixel *)calloc((size_t)tmp->ysize, sizeof(Pixel));
	if (!raster) {
		for (int i = 0; i < dst->ysize; ++i)
			free(contrib[i].p);
		free(contrib);
		free_image_sc(tmp);
		return;
	}
	for (int k = 0; k < dst->xsize; ++k) {
		sc_get_column(raster, tmp, k);
		for (int i = 0; i < dst->ysize; ++i) {
			double weight = 0.0;
			for (int j = 0; j < contrib[i].n; ++j) {
				weight += raster[contrib[i].p[j].pixel] * contrib[i].p[j].weight;
			}
			sc_put_pixel(dst, k, i, (Pixel)CLAMP(weight, 0, 255));
		}
	}
	free(raster);

	for (int i = 0; i < dst->ysize; ++i)
		free(contrib[i].p);
	free(contrib);

	free_image_sc(tmp);
}

/* ----------------- Multi-component wrapper (public zoom) -----------------
   Extract each channel from interleaved src into a single-channel buffer,
   call zoom_single for that channel, and write result back into interleaved dst.
*/
static unsigned char *extract_channel_interleaved(const Image *src, int c) {
	size_t n = (size_t)src->xsize * (size_t)src->ysize;
	unsigned char *buf = (unsigned char *)malloc(n);
	if (!buf)
		return NULL;
	const Pixel *p = src->data + (size_t)c;
	for (int y = 0; y < src->ysize; ++y) {
		const Pixel *prow = p + (size_t)y * (size_t)src->span;
		unsigned char *out = buf + (size_t)y * (size_t)src->xsize;
		for (int x = 0; x < src->xsize; ++x) {
			out[x] = prow[x * src->components];
		}
	}
	return buf;
}

static int write_channel_interleaved(Image *dst, int c, const ImageSC *img) {
	if (!dst || !img)
		return -1;
	for (int y = 0; y < dst->ysize; ++y) {
		Pixel *dprow = dst->data + (size_t)y * (size_t)dst->span;
		const Pixel *srcrow = img->data + (size_t)y * (size_t)img->span;
		for (int x = 0; x < dst->xsize; ++x) {
			dprow[x * dst->components + c] = srcrow[x];
		}
	}
	return 0;
}

int zoom(Image *dst, const Image *src, double (*filterf)(double), double fwidth) {
	if (!dst || !src || !filterf)
		return -1;
	if (dst->xsize <= 0 || dst->ysize <= 0 || src->xsize <= 0 || src->ysize <= 0)
		return -1;
	if (dst->components <= 0 || src->components <= 0)
		return -1;
	if (dst->components != src->components) {
		fprintf(stderr, "zoom: component mismatch dst=%d src=%d\n", dst->components, src->components);
		return -1;
	}

	int comp = dst->components;
	if (comp == 1) {
		/* Fast path: single-channel -- create ImageSC wrappers pointing to contiguous single-channel memory
		   If underlying storage already single-channel contiguous, we can wrap without copy.
		   However our Image is interleaved, so for safety we still extract into an ImageSC, except when span==xsize and components==1.
		*/
		if (src->components == 1 && src->span == src->xsize && dst->components == 1 && dst->span == dst->xsize) {
			ImageSC srcsc = {src->xsize, src->ysize, src->data, src->span};
			ImageSC dstsc = {dst->xsize, dst->ysize, dst->data, dst->span};
			zoom_single(&dstsc, &srcsc, filterf, fwidth);
			return 0;
		}
	}

	/* For multi-component images, process each channel separately */
	for (int c = 0; c < comp; ++c) {
		unsigned char *src_plane = extract_channel_interleaved(src, c);
		if (!src_plane) {
			fprintf(stderr, "zoom: OOM extracting channel %d\n", c);
			return -1;
		}
		ImageSC *srcsc = new_image_sc(src->xsize, src->ysize);
		if (!srcsc) {
			free(src_plane);
			return -1;
		}
		/* copy extracted data into srcsc->data */
		memcpy(srcsc->data, src_plane, (size_t)src->xsize * (size_t)src->ysize);
		free(src_plane);

		ImageSC *dstsc = new_image_sc(dst->xsize, dst->ysize);
		if (!dstsc) {
			free_image_sc(srcsc);
			return -1;
		}

		/* run single-channel resampler */
		zoom_single(dstsc, srcsc, filterf, fwidth);

		/* write back into interleaved dst */
		write_channel_interleaved(dst, c, dstsc);

		free_image_sc(srcsc);
		free_image_sc(dstsc);
	}

	return 0;
}
