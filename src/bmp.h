#ifndef BMPGEN_H
#define BMPGEN_H

#include <nifti2_io.h>
#include <stdint.h>

/**
 * Create a test PNG (and JPG) image using stb_image_write.
 *
 * @param filename_png  Path to output PNG file
 * @param filename_jpg  Path to output JPG file
 * @param width         Image width in pixels
 * @param height        Image height in pixels
 * @return 0 on success, nonzero on failure
 */

int nim2png(nifti_image *nim, nifti_image *nim2, size_t argc, const char *argv[], const char *out_path);

#endif // BMPGEN_H
