#ifndef EX_BWLABEL_H
#define EX_BWLABEL_H

#ifdef  __cplusplus
extern "C" {
#endif

#include <nifti2_io.h>

int bwlabel(nifti_image *nim, int conn);

//int butter_design(int order, double fl, double fh, double ** a, double ** b, double ** IC);
//void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])


#endif // EX_BUTTER_H
