#ifndef NII_DTIFIT_H
#define NII_DTIFIT_H

#ifdef __cplusplus
extern "C" {
#endif

// Entry point for "niimath --dtifit ...". Parses its own argv (FSL dtifit-style
// flags) and writes FA/MD/L*/V*/S0/MO/tensor maps. Returns EXIT_SUCCESS/FAILURE.
int nii_dtifit(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif // NII_DTIFIT_H
