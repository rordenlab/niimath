#ifndef NII_QC_H
#define NII_QC_H

#ifdef __cplusplus
extern "C" {
#endif

// Entry point for "niimath --qc ...". Given a T1-weighted image, a matching
// integer tissue segmentation, and the label values that denote CSF and white
// matter, computes MRIQC-style hard-mask IQMs (CJV, no-air CNR, within-tissue
// SNR, WM2MAX, brain-only EFC, ICV fractions, per-tissue summary stats) and writes
// them to a wide TSV. Parses its own argv. Returns EXIT_SUCCESS/FAILURE.
int nii_qc(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif // NII_QC_H
