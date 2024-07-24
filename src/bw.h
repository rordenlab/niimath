#ifndef EX_BUTTER_H
#define EX_BUTTER_H

#ifdef  __cplusplus
extern "C" {
#endif

double *dcof_bwlp( int n, double fcf );
double *dcof_bwhp( int n, double fcf );
double *dcof_bwbp( int n, double f1f, double f2f );
double *dcof_bwbs( int n, double f1f, double f2f );
int *ccof_bwlp( int n );
int *ccof_bwhp( int n );
int *ccof_bwbp( int n );
double *ccof_bwbs( int n, double f1f, double f2f );
double sf_bwlp( int n, double fcf );
double sf_bwhp( int n, double fcf );
double sf_bwbp( int n, double f1f, double f2f );
double sf_bwbs( int n, double f1f, double f2f );
void Filt(double *X, int nX, double *a, double *b, int order, double *Z);
void FiltRev(double *X, int nX, double *a, double *b, int order, double *Z);
int butter_design(int order, double fl, double fh, double ** a, double ** b, double ** IC);

#ifdef  __cplusplus
}
#endif

#endif // EX_BUTTER_H
