#ifndef MINDGRAB_H
#define MINDGRAB_H

#ifdef __cplusplus
extern "C" {
#endif

// MindGrab skull stripping (-mindgrab). Clean-room C implementation of the MeshNet
// inference the brainchop-cli reference performs with tinygrad; the weights themselves
// are the upstream MIT-licensed MindGrab model (src/mindgrab.LICENSE).
//
// mindgrab_segment() is the whole conformed-space job: quantile normalisation, the
// 25-block dilated MeshNet, the 2-class argmax and the largest 26-connected component.
//
//   conformed  256^3 uint8, x fastest, as produced by `-conform ... -odt char`
//   mask       256^3 float, caller-allocated, written 0.0f / 1.0f
//
// Returns 0 on success; nonzero (with mask untouched) if one of ITS OWN allocations fails. The
// final connected-component pass calls bwlabel(), whose allocator is fail-closed and calls
// exit(EXIT_FAILURE) rather than returning -- matching niimath's nii_calloc policy, but it does
// mean an OOM there terminates the process instead of unwinding. Allocates 2.15 GB
// of activation buffers (measured whole-process peak for -mindgrab is about 2.5 GB).
// Deterministic: no reduction crosses voxels and the two
// per-channel moments are reduced in a fixed slice order, so the result is identical for
// any thread count.
int mindgrab_segment(const unsigned char *conformed, float *mask);

#ifdef __cplusplus
}
#endif

#endif // MINDGRAB_H
