/* Statistics for Python
 * Copyright (C) 2006 Michiel Jan Laurens de Hoon.
 *
 * This library was written at the Center for Computational Biology and
 * Bioinformatics, Columbia University, 1130 Saint Nicholas Avenue,
 * New York, NY 10025, United States.
 * Contact: mdehoon@c2b2.columbia.edu
 *
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 * 
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 * 
 */



/* Kernel-based probability density estimation */
void pdf(int ndata, double data[], int dstride, double weight[], int wstride,
  double h, char kernel, int n, double x[], int xstride, double y[]);
void cpdf(int ndata, double data[], int dstride, double h, char kernel,
  int n, double x[], int xstride, double y[]);
void cpdfc(int ndata, double data[], int dstride, double h, char kernel,
  int n, double x[], int xstride, double y[]);
double bandwidth (int ndata, double data[], int stride,
  double weight[], int wstride, char kernel);

/* Kernel-based local regression */
void localfit(int ndata, double* data, int strides[2], double h, char kernel,
  int n, double x[], int xstride, double y[]);

/* Utility routines, currently undocumented */
double mean(int n, double x[]);
double variance(int n, double data[], int stride, char mode);
double median(int n, double x[]);
int covariance(int n, int m, double data[], int strides[2], char mode, double matrix[]);
int pearson(int n, int m, double data[], int strides[2], double matrix[]);
int spearman(int n, int m, double data[], double matrix[]);
int intraclass(int n, int m, double data[], int strides[2], double matrix[]);
int regression(int n, double x[], double y[], double* a, double* b);
