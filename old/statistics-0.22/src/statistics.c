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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257390
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

#ifndef M_PI_4
#define M_PI_4 0.78539816339744830962
#endif

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#ifndef min
#define min(x, y)       ((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define max(x, y)       ((x) > (y) ? (x) : (y))
#endif

#include <statistics.h>

static const double* sortdata = NULL; /* used in the quicksort algorithm */

/* ************************************************************************ */

double mean(int n, double x[])
{ double result = 0.;
  int i;
  for (i = 0; i < n; i++) result += x[i];
  result /= n;
  return result;
}

/* ************************************************************************ */

double median (int n, double x[])
/*
Find the median of X(1), ... , X(N), using as much of the quicksort
algorithm as is needed to isolate it.
N.B. On exit, the array X is partially ordered.
Based on Alan J. Miller's median.f90 routine.
*/

{ int i, j;
  int nr = n / 2;
  int nl = nr - 1;
  int even = 0;
  /* hi & lo are position limits encompassing the median. */
  int lo = 0;
  int hi = n-1;

  if (n==2*nr) even = 1;
  if (n<3)
  { if (n<1) return 0.;
    if (n == 1) return x[0];
    return 0.5*(x[0]+x[1]);
  }

  /* Find median of 1st, middle & last values. */
  do
  { int loop;
    int mid = (lo + hi)/2;
    double result = x[mid];
    double xlo = x[lo];
    double xhi = x[hi];
    if (xhi<xlo)
    { double temp = xlo;
      xlo = xhi;
      xhi = temp;
    }
    if (result>xhi) result = xhi;
    else if (result<xlo) result = xlo;
    /* The basic quicksort algorithm to move all values <= the sort key (XMED)
     * to the left-hand end, and all higher values to the other end.
     */
    i = lo;
    j = hi;
    do
    { while (x[i]<result) i++;
      while (x[j]>result) j--;
      loop = 0;
      if (i<j)
      { double temp = x[i];
        x[i] = x[j];
        x[j] = temp;
        i++;
        j--;
        if (i<=j) loop = 1;
      }
    } while (loop); /* Decide which half the median is in. */

    if (even)
    { if (j==nl && i==nr)
        /* Special case, n even, j = n/2 & i = j + 1, so the median is
         * between the two halves of the series.   Find max. of the first
         * half & min. of the second half, then average.
         */
        { int k;
          double xmax = x[0];
          double xmin = x[n-1];
          for (k = lo; k <= j; k++) xmax = max(xmax,x[k]);
          for (k = i; k <= hi; k++) xmin = min(xmin,x[k]);
          return 0.5*(xmin + xmax);
        }
      if (j<nl) lo = i;
      if (i>nr) hi = j;
      if (i==j)
      { if (i==nl) lo = nl;
        if (j==nr) hi = nr;
      }
    }
    else
    { if (j<nr) lo = i;
      if (i>nr) hi = j;
      /* Test whether median has been isolated. */
      if (i==j && i==nr) return result;
    }
  }
  while (lo<hi-1);

  if (even) return (0.5*(x[nl]+x[nr]));
  if (x[lo]>x[hi])
  { double temp = x[lo];
    x[lo] = x[hi];
    x[hi] = temp;
  }
  return x[nr];
}

/* ************************************************************************ */

double variance(int n, double data[], int stride, char mode)
/* This algorithm is described in:
 * B.P. Welford:
 * "Note on a method for calculating corrected sums of squares and products."
 * Technometrics 4(3): 419-420 (1962).
 * Also see:
 * Peter M. Neely:
 * "Comparison of several algorithms for computation of means, standard
 * deviations and correlation coefficients."
 * Communications of the ACM 9(7): 496-499 (1966).
 */
{ int i;
  double mean = 0.0;
  double result = 0.0;
  double value;
  double delta;
  double* p = data;
  for (i = 1; i <= n; i++, p+=stride)
  { value = *p;
    delta = value - mean;
    mean += delta/i;
    result += delta * (value-mean);
  }
  if (mode=='u') result /= (n-1); /* Unbiased estimate */
  else result /= n;               /* Maximum-likelihood estimate */
  return result;
}

/* ************************************************************************ */

int covariance(int n, int m, double data[], int strides[2], char mode, double matrix[])
/* This algorithm is described in:
 * B.P. Welford:
 * "Note on a method for calculating corrected sums of squares and products."
 * Technometrics 4(3): 419-420 (1962).
 * Also see:
 * Peter M. Neely:
 * "Comparison of several algorithms for computation of means, standard
 * deviations and correlation coefficients."
 * Communications of the ACM 9(7): 496-499 (1966).
 */
{ int i, j, k;
  double* p;
  double* q;
  double* q1;
  double* q2;
  double x1;
  double x2;
  const int stride1 = strides[0];
  const int stride2 = strides[1];
  const int denominator = (mode=='u') ? n-1 : n;
  /* 'u': Unbiased estimate; otherwise maximum-likelihood estimate */
  double* average = malloc(m*sizeof(double));
  if(!average) return 0;

  /* Initialize to zero */
  p = matrix;
  for (i = 0; i < m; i++)
  { average[i] = 0.0;
    for (j = 0; j < m; j++, p++) *p = 0.0;
  }

  /* Calculate the sums of squares */
  for (k = 0; k < n; k++)
  { const double scale = k+1.0;
    const double factor = k/scale;
    q = data + k*stride1;
    q1 = q;
    for (i = 0; i < m; i++, q1+=stride2)
    { p = matrix + i*m;
      x1 = (*q1) - average[i];
      q2 = q;
      for (j = 0; j <= i; j++, p++, q2+=stride2)
      { x2 = (*q2) - average[j];
        *p += factor*x1*x2;
      }
    }
    for (i = 0; i < m; i++, q+=stride2)
      average[i] = factor * average[i] + (*q) / scale;
  }
  free(average);
 
  /* Scale, and copy to the upper half of the matrix */
  for (i = 0; i < m; i++)
  { p = matrix + i*m;
    q = matrix + i + (i+1)*m;
    p[i] /= denominator;
    for (j = i+1; j < m; j++, q+=m)
    { (*q) /= denominator;
      p[j] = (*q);
    }
  }
  return 1;
}

/* ************************************************************************ */

int pearson(int n, int m, double data[], int strides[2], double matrix[])
/* This algorithm is described in:
 * B.P. Welford:
 * "Note on a method for calculating corrected sums of squares and products."
 * Technometrics 4(3): 419-420 (1962).
 * Also see:
 * Peter M. Neely:
 * "Comparison of several algorithms for computation of means, standard
 * deviations and correlation coefficients."
 * Communications of the ACM 9(7): 496-499 (1966).
 */
{ int i, j, k;
  double* p;
  double* q;
  double* q1;
  double* q2;
  double x1;
  double x2;
  const int stride1 = strides[0];
  const int stride2 = strides[1];
  double scale, factor;
  double* average = malloc(m*sizeof(double));
  if(!average) return 0;

  /* Initialize to zero */
  p = matrix;
  for (i = 0; i < m; i++)
  { average[i] = 0.0;
    for (j = 0; j < m; j++, p++) *p = 0.0;
  }

  /* Calculate the sums of squares */
  for (k = 0; k < n; k++)
  { scale = k+1.0;
    factor = k/scale;
    q = data + k*stride1;
    q1 = q;
    for (i = 0; i < m; i++, q1+=stride2)
    { p = matrix + i*m;
      x1 = (*q1) - average[i];
      q2 = q;
      for (j = 0; j <= i; j++, p++, q2+=stride2)
      { x2 = (*q2) - average[j];
        *p += factor*x1*x2;
      }
    }
    for (i = 0; i < m; i++, q+=stride2)
      average[i] = factor * average[i] + (*q) / scale;
  }
  free(average);
 
  /* Divide by the standard deviation, and
   * copy to the upper half of the matrix */
  p = matrix;
  k = m+1;
  for (i = 0; i < m; i++, p+=k) *p = sqrt(*p);

  for (i = 0; i < m; i++)
  { p = matrix + i*m;
    scale = p[i];
    q1 = matrix + (i+1)*m + i;
    q2 = q1 + 1;
    for (j = i+1; j < m; j++, q1+=m, q2+=k)
    { (*q1) /= ((*q2)*scale);
      p[j] = (*q1);
    }
  }

  p = matrix;
  for (i = 0; i < m; i++, p+=k) *p = 1.0;

  return 1;
}

/* ********************************************************************** */

static
int compare(const void* a, const void* b)
/* Helper function for sort. */
{ const int i1 = *(const int*)a;
  const int i2 = *(const int*)b;
  const double term1 = sortdata[i1];
  const double term2 = sortdata[i2];
  if (term1 < term2) return -1;
  if (term1 > term2) return +1;
  return 0;
}

/* ---------------------------------------------------------------------- */

void sort(int n, const double data[], int index[])
/* Sets up an index table given the data, such that data[index[]] is in
 * increasing order. Sorting is done on the indices; the array data
 * is unchanged.
 */
{ int i;
  sortdata = data;
  for (i = 0; i < n; i++) index[i] = i;
  qsort(index, n, sizeof(int), compare);
}

/* ********************************************************************** */

int spearman(int n, int m, double data[], double matrix[])
/* Find the Spearman rank correlation between the n rows in data.
 * Each row has m columns.
 */
{ int i, j, k;
  double* p;
  double* q1;
  double* q2;
  double x1;
  double x2;
  double scale;
  const double average = 0.5*(m-1); /* Average rank */
  /* If two elements have the same rank, the squared sum of
   * their ranks will change. The squared sum of the ranks
   * can therefore not be precalculated.
   */

  /* Convert all values to ranks */
  int* index = malloc(m*sizeof(int));
  if (!index) return 0;
  p = data;
  for (i = 0; i < n; i++, p+=m)
  { double value, rank;
    j = 0;
    k = 0;
    /* Call sort to get an index table */
    sort(m, p, index);
    /* Build a rank table */
    while (j < m)
    {
      value = p[index[j]];
      do j++; while(j < m && p[index[j]]==value);
      rank = 0.5*(j+k)-0.5;
      for ( ; k < j; k++) p[index[k]] = rank;
    }
  }
  free(index);

  /* Initialize to zero */
  for (i = 0; i < n*n; i++) matrix[i] = 0.0;

  /* Calculate the sums of squares */
  for (k = 0; k < m; k++)
  { q1 = data + k;
    for (i = 0; i < n; i++, q1+=m)
    { p = matrix + i*n;
      x1 = (*q1) - average;
      q2 = data + k;
      for (j = 0; j <= i; j++, p++, q2+=m)
      { x2 = (*q2) - average;
        *p += x1*x2;
      }
    }
  }
 
  /* Divide by the standard deviation, and
   * copy to the upper half of the matrix */
  p = matrix;
  k = n+1;
  for (i = 0; i < n; i++, p+=k) *p = sqrt(*p);

  for (i = 0; i < n; i++)
  { p = matrix + i*n;
    scale = p[i];
    q1 = matrix + (i+1)*n + i;
    q2 = q1 + 1;
    for (j = i+1; j < n; j++, q1+=n, q2+=k)
    { (*q1) /= ((*q2)*scale);
      p[j] = (*q1);
    }
  }

  p = matrix;
  for (i = 0; i < n; i++, p+=k) *p = 1.0;

  return 1;
}

/* ************************************************************************ */

int intraclass(int n, int m, double data[], int strides[2], double matrix[])
/* The intraclass correlation is described in:
 * Ronald A. Fisher: Statistical Methods for Research Workers, chapter VII.
 * Oliver and Boyd, Edinburgh/London (1925).
 * The recursion algorithm used here is related to the Welford algorithm:
 * B.P. Welford:
 * "Note on a method for calculating corrected sums of squares and products."
 * Technometrics 4(3): 419-420 (1962).
 */
{ int i, j, k;
  double* p;
  double* q;
  double* q1;
  double* q2;
  double x1;
  double x2;
  double x3;
  const int stride1 = strides[0];
  const int stride2 = strides[1];
  double scale, factor;
  double* average = malloc(m*sizeof(double));
  if(!average) return 0;

  /* Initialize to zero */
  p = matrix;
  for (i = 0; i < m; i++)
  { average[i] = 0.0;
    for (j = 0; j < m; j++, p++) *p = 0.0;
  }

  /* Calculate the sums of squares */
  for (k = 0; k < n; k++)
  { scale = k+1.0;
    factor = k/scale;
    q = data + k*stride1;
    q1 = q;
    for (i = 0; i < m; i++, q1+=stride2)
    { p = matrix + i*m;
      x1 = (*q1) - average[i];
      q2 = q;
      for (j = 0; j < i; j++, p++, q2+=stride2)
      { x2 = x1 + (*q2) - average[j];
        x3 = (*q1) - (*q2);
        *p += factor*x2*x2 - x3*x3;
      }
      *p += factor*x1*x1;
    }
    for (i = 0; i < m; i++, q+=stride2)
      average[i] = factor * average[i] + (*q) / scale;
  }
 
  /* Divide by the standard deviation, and
   * copy to the upper half of the matrix */
  p = matrix;
  k = m+1;

  for (i = 0; i < m; i++)
  { p = matrix + i*m;
    scale = p[i];
    q1 = matrix + (i+1)*m + i;
    q2 = q1 + 1;
    for (j = i+1; j < m; j++, q1+=m, q2+=k)
    { x3 = average[i] - average[j];
      (*q1) /= (2*(p[i]+(*q2)) + n*x3*x3);
      p[j] = (*q1);
    }
  }
  free(average);

  p = matrix;
  for (i = 0; i < m; i++, p+=k) *p = 1.0;

  return 1;
}

/* ********************************************************************** */

int regression(int n, double x[], double y[], double* a, double* b)
/* Fit a linear regression line to the vectors x and y, each consisting
 * of n elements. The output variable a is the intercept; b is the slope.
 * The function returns 0 if an error occurs (in particular, if the variance
 * in x is zero), and returns 1 otherwise.
 */
{ int i;
  double sx = 0.0;
  double sy = 0.0;
  double sxx = 0.0;
  double sxy = 0.0;
  double numerator, denominator;
  for (i = 0; i < n; i++)
  { double xv = x[i];
    double yv = y[i];
    sx += xv;
    sy += yv;
    sxx += xv*xv;
    sxy += xv*yv;
  }
  numerator = n*sxy - sx*sy;
  denominator = n*sxx - sx*sx;
  if (denominator==0.0) return 0; /* Error */
  *b = numerator / denominator;
  *a = (sy - (*b)*sx)/n;
  return 1;
}

/* ********************************************************************* */

void pdf(int ndata, double data[], int dstride, double weight[], int wstride,
         double h, char kernel, int n, double x[], int xstride, double y[])
{ int i, j;
  double factor;
  if (weight)
  { factor = 0.0;
    for (j = 0; j < ndata; j++) factor += weight[j];
    factor = 1.0 / (factor*h);
    switch (kernel)
    { case 'e': /* Epanechnikov */
      { factor *= 0.75;
        for (i = 0; i < n; i++)
        { double temp = x[i*xstride];
          y[i] = 0.0;
          for (j = 0; j < ndata; j++)
          { double u = (data[j*dstride] - temp) / h;
            if (u < -1.0 || u > 1.0) continue;
            y[i] += weight[j*wstride] * (1-u*u);
          }
        }
        break;
      }
      case 'u': /* Uniform */
      { for (i = 0; i < n; i++)
        { double temp = x[i*xstride];
          y[i] = 0.0;
          for (j = 0; j < ndata; j++)
          { double u = data[j*dstride] - temp;
            if (u < -h || u > h) continue;
            y[i] += 0.5 * weight[j*wstride];
          }
        }
        break;
      }
      case 't': /* Triangle */
      { for (i = 0; i < n; i++)
        { double temp = x[i*xstride];
          y[i] = 0.0;
          for (j = 0; j < ndata; j++)
          { double u = fabs(data[j*dstride] - temp) / h;
            if (u > 1) continue;
            y[i] += (1-u) * weight[j*wstride];
          }
        }
        break;
      }
      case 'g': /* Gaussian */
      { factor *= 0.5 * M_SQRT1_2 * M_2_SQRTPI;
        for (i = 0; i < n; i++)
        { double temp = x[i*xstride];
          y[i] = 0.0;
          for (j = 0; j < ndata; j++)
          { double u = (data[j*dstride] - temp) / h;
            y[i] += exp(-0.5*u*u) * weight[j*wstride];
          }
        }
        break;
      }
      case 'b': /* Biweight or quartic */
      { factor *= 0.9375;
        for (i = 0; i < n; i++)
        { double temp = x[i*xstride];
          y[i] = 0.0;
          for (j = 0; j < ndata; j++)
          { double u = (data[j*dstride] - temp) / h;
            if (u < -1.0 || u > 1.0) continue;
            y[i] += (1-u*u)*(1-u*u) * weight[j*wstride];
          }
        }
        break;
      }
      case '3': /* Triweight */
      { factor *= 1.09375;
        for (i = 0; i < n; i++)
        { double temp = x[i*xstride];
          y[i] = 0.0;
          for (j = 0; j < ndata; j++)
          { double u = (data[j*dstride] - temp) / h;
            if (u < -1.0 || u > 1.0) continue;
            y[i] += (1-u*u)*(1-u*u)*(1-u*u) * weight[j*wstride];
          }
        }
        break;
      }
      case 'c': /* Cosine */
      { factor *= M_PI_4;
        for (i = 0; i < n; i++)
        { double temp = x[i*xstride];
          y[i] = 0.0;
          for (j = 0; j < ndata; j++)
          { double u = (data[j*dstride] - temp) / h;
            if (u < -1.0 || u > 1.0) continue;
            y[i] += cos(M_PI_2*u) * weight[j*wstride];
          }
        }
        break;
      }
    }
  }
  else
  { factor = 1.0/(ndata*h);
    switch (kernel)
    { case 'e': /* Epanechnikov */
      { factor *= 0.75;
        for (i = 0; i < n; i++)
        { double temp = x[i*xstride];
          y[i] = 0.0;
          for (j = 0; j < ndata; j++)
          { double u = (data[j*dstride] - temp) / h;
            if (u < -1.0 || u > 1.0) continue;
            y[i] += 1-u*u;
          }
        }
        break;
      }
      case 'u': /* Uniform */
      { for (i = 0; i < n; i++)
        { double temp = x[i*xstride];
          y[i] = 0.0;
          for (j = 0; j < ndata; j++)
          { double u = data[j*dstride] - temp;
            if (u < -h || u > h) continue;
            y[i] += 0.5;
          }
        }
        break;
      }
      case 't': /* Triangle */
      { for (i = 0; i < n; i++)
        { double temp = x[i*xstride];
          y[i] = 0.0;
          for (j = 0; j < ndata; j++)
          { double u = fabs(data[j*dstride] - temp) / h;
            if (u > 1) continue;
            y[i] += 1-u;
          }
        }
        break;
      }
      case 'g': /* Gaussian */
      { factor *= 0.5 * M_SQRT1_2 * M_2_SQRTPI;
        for (i = 0; i < n; i++)
        { double temp = x[i*xstride];
          y[i] = 0.0;
          for (j = 0; j < ndata; j++)
          { double u = (data[j*dstride] - temp) / h;
            y[i] += exp(-0.5*u*u);
          }
        }
        break;
      }
      case 'b': /* Biweight or quartic */
      { factor *= 0.9375;
        for (i = 0; i < n; i++)
        { double temp = x[i*xstride];
          y[i] = 0.0;
          for (j = 0; j < ndata; j++)
          { double u = (data[j*dstride] - temp) / h;
            if (u < -1.0 || u > 1.0) continue;
            y[i] += (1-u*u)*(1-u*u);
          }
        }
        break;
      }
      case '3': /* Triweight */
      { factor *= 1.09375;
        for (i = 0; i < n; i++)
        { double temp = x[i*xstride];
          y[i] = 0.0;
          for (j = 0; j < ndata; j++)
          { double u = (data[j*dstride] - temp) / h;
            if (u < -1.0 || u > 1.0) continue;
            y[i] += (1-u*u)*(1-u*u)*(1-u*u);
          }
        }
        break;
      }
      case 'c': /* Cosine */
      { factor *= M_PI_4;
        for (i = 0; i < n; i++)
        { double temp = x[i*xstride];
          y[i] = 0.0;
          for (j = 0; j < ndata; j++)
          { double u = (data[j*dstride] - temp) / h;
            if (u < -1.0 || u > 1.0) continue;
            y[i] += cos(M_PI_2*u);
          }
        }
        break;
      }
    }
  }
  for (i = 0; i < n; i++) y[i] *= factor;
}

/* ********************************************************************* */

void cpdf(int ndata, double data[], int dstride, double h, char kernel,
          int n, double x[], int xstride, double y[])
{ int i, j;
  switch (kernel)
  { case 'e': /* Epanechnikov */
    { const double factor = 0.25/ndata;
      for (i = 0; i < n; i++)
      { double ordinate = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (ordinate - data[j*dstride]) / h;
          if (u < -1.0) continue;
          if (u > +1.0) y[i] += 4.0;
          else y[i] += (-u*u*u+3*u+2);
        }
      }
      for (i = 0; i < n; i++) y[i] *= factor;
      return;
    }
    case 'u': /* Uniform */
    { const double factor = 0.5/ndata;
      for (i = 0; i < n; i++)
      { double ordinate = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (ordinate - data[j*dstride]) / h;
          if (u < -1.0) continue;
          if (u > +1.0) y[i] += 2.0;
          else y[i] += u+1.0;
        }
      }
      for (i = 0; i < n; i++) y[i] *= factor;
      return;
    }
    case 't': /* Triangle */
    { for (i = 0; i < n; i++)
      { double ordinate = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (ordinate - data[j*dstride]) / h;
          if (u < -1.0) continue;
          else if (u < 0.0) y[i] +=  0.5*u*u + u + 0.5;
          else if (u < 1.0) y[i] += -0.5*u*u + u + 0.5;
          else y[i] += 1.0;
        }
      }
      for (i = 0; i < n; i++) y[i] /= ndata;
      return;
    }
#ifdef HAVE_ERF
    case 'g': /* Gaussian */
    { const double factor = 0.5/ndata;
      for (i = 0; i < n; i++)
      { double ordinate = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (ordinate - data[j*dstride]) / h;
          y[i] += erf(u*M_SQRT1_2);
        }
      }
      for (i = 0; i < n; i++) y[i] = y[i] * factor + 0.5;
      return;
    }
#endif
    case 'b': /* Biweight or quartic */
    { const double factor = 0.0625/ndata;
      for (i = 0; i < n; i++)
      { double ordinate = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (ordinate - data[j*dstride]) / h;
          if (u < -1.0) continue;
          if (u > +1.0) y[i] += 16.0;
          else y[i] += (8.0+15.0*u-10.0*u*u*u+3.0*u*u*u*u*u);
        }
      }
      for (i = 0; i < n; i++) y[i] *= factor;
      return;
    }
    case '3': /* Triweight */
    { const double factor = 0.03125/ndata;
      for (i = 0; i < n; i++)
      { double ordinate = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (ordinate - data[j*dstride]) / h;
          if (u < -1.0) continue;
          if (u > +1.0) y[i] += 32.0;
          else y[i] += (16.0+35.0*u-35.0*u*u*u+21.0*u*u*u*u*u-5.0*u*u*u*u*u*u*u);
        }
      }
      for (i = 0; i < n; i++) y[i] *= factor;
      return;
    }
    case 'c': /* Cosine */
    { const double factor = 0.5/ndata;
      for (i = 0; i < n; i++)
      { double ordinate = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (ordinate - data[j*dstride]) / h;
          if (u < -1.0) continue;
          if (u > +1.0) y[i] += 2.0;
          else y[i] += sin(M_PI_2*u) + 1.0;
        }
      }
      for (i = 0; i < n; i++) y[i] *= factor;
      return;
    }
  }
}

/* ********************************************************************* */

void cpdfc(int ndata, double data[], int dstride, double h, char kernel,
           int n, double x[], int xstride, double y[])
{ int i, j;
  switch (kernel)
  { case 'e': /* Epanechnikov */
    { const double factor = 0.25/ndata;
      for (i = 0; i < n; i++)
      { double ordinate = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (ordinate - data[j*dstride]) / h;
          if (u > +1.0) continue;
          if (u < -1.0) y[i] += 4.0;
          else y[i] += (2+u*u*u-3*u);
        }
      }
      for (i = 0; i < n; i++) y[i] *= factor;
      return;
    }
    case 'u': /* Uniform */
    { const double factor = 0.5/ndata;
      for (i = 0; i < n; i++)
      { double ordinate = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (ordinate - data[j*dstride]) / h;
          if (u > +1.0) continue;
          if (u < -1.0) y[i] += 2.0; 
          else y[i] += 1.0-u;
        }
      }
      for (i = 0; i < n; i++) y[i] *= factor;
      return;
    }
    case 't': /* Triangle */
    { for (i = 0; i < n; i++)
      { double ordinate = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (ordinate - data[j*dstride]) / h;
          if (u > 1.0) continue;
          else if (u >  0.0) y[i] += +0.5*u*u - u + 0.5;
          else if (u > -1.0) y[i] += -0.5*u*u - u + 0.5;
          else y[i] += 1.0;
        }
      }
      for (i = 0; i < n; i++) y[i] /= ndata;
      return;
    }
#ifdef HAVE_ERFC
    case 'g': /* Gaussian */
    { const double factor = 0.5/ndata;
      for (i = 0; i < n; i++)
      { double ordinate = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (ordinate - data[j*dstride]) / h;
          y[i] += erfc(u*M_SQRT1_2);
        }
      }
      for (i = 0; i < n; i++) y[i] = y[i] * factor;
      return;
    }
#endif
    case 'b': /* Biweight or quartic */
    { const double factor = 0.0625/ndata;
      for (i = 0; i < n; i++)
      { double ordinate = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (ordinate - data[j*dstride]) / h;
          if (u > +1.0) continue;
          if (u < -1.0) y[i] += 16.0;
          else y[i] += (8.0-15.0*u+10.0*u*u*u-3.0*u*u*u*u*u);
        }
      }
      for (i = 0; i < n; i++) y[i] *= factor;
      return;
    }
    case '3': /* Triweight */
    { const double factor = 0.03125/ndata;
      for (i = 0; i < n; i++)
      { double ordinate = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (ordinate - data[j*dstride]) / h;
          if (u > +1.0) continue;
          if (u < -1.0) y[i] += 32.0;
          else y[i] += 16.0-35.0*u+35.0*u*u*u-21.0*u*u*u*u*u+5.0*u*u*u*u*u*u*u;
        }
      }
      for (i = 0; i < n; i++) y[i] *= factor;
      return;
    }
    case 'c': /* Cosine */
    { const double factor = 0.5/ndata;
      for (i = 0; i < n; i++)
      { double ordinate = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (ordinate - data[j*dstride]) / h;
          if (u > +1.0) continue;
          if (u < -1.0) y[i] += 2.0;
          else y[i] += 1.0 - sin(M_PI_2*u);
        }
      }
      for (i = 0; i < n; i++) y[i] *= factor;
      return;
    }
  }
}

/* ********************************************************************* */

void localfit(int ndata, double* data, int strides[2], double h, char kernel,
              int n, double x[], int xstride, double y[])
{ int i, j;
  int dstride = strides[0];
  switch (kernel)
  { case 'e': /* Epanechnikov */
    { for (i = 0; i < n; i++)
      { double temp = x[i*xstride];
        double numerator = 0.0;
        double denominator = 0.0;
        for (j = 0; j < ndata; j++)
        { int index = j*strides[0];
          double u = (data[index] - temp) / h;
          double value = data[index+strides[1]];
          if (u < -1.0 || u > 1.0) continue;
          numerator += (1-u*u) * value;
          denominator += 1-u*u;
        }
        if (denominator > 0.0) y[i] = numerator / denominator;
        else y[i] = 0.0;
      }
      break;
    }
    case 'u': /* Uniform */
    { for (i = 0; i < n; i++)
      { double temp = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = data[j*dstride] - temp;
          if (u < -h || u > h) continue;
          y[i] += 0.5;
        }
      }
      break;
    }
    case 't': /* Triangle */
    { for (i = 0; i < n; i++)
      { double temp = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = fabs(data[j*dstride] - temp) / h;
          if (u > 1) continue;
          y[i] += u;
        }
      }
      break;
    }
    case 'g': /* Gaussian */
    { for (i = 0; i < n; i++)
      { double temp = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (data[j*dstride] - temp) / h;
          y[i] += exp(-0.5*u*u);
        }
      }
      break;
    }
    case 'b': /* Biweight or quartic */
    { for (i = 0; i < n; i++)
      { double temp = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (data[j*dstride] - temp) / h;
          if (u < -1.0 || u > 1.0) continue;
          y[i] += (1-u*u)*(1-u*u);
        }
      }
      break;
    }
    case '3': /* Triweight */
    { for (i = 0; i < n; i++)
      { double temp = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (data[j*dstride] - temp) / h;
          if (u < -1.0 || u > 1.0) continue;
          y[i] += (1-u*u)*(1-u*u)*(1-u*u);
        }
      }
      break;
    }
    case 'c': /* Cosine */
    { for (i = 0; i < n; i++)
      { double temp = x[i*xstride];
        y[i] = 0.0;
        for (j = 0; j < ndata; j++)
        { double u = (data[j*dstride] - temp) / h;
          if (u < -1.0 || u > 1.0) continue;
          y[i] += cos(M_PI_2*u);
        }
      }
      break;
    }
  }
}

double bandwidth(int n, double data[], int dstride,
  double weight[], int wstride, char kernel)
{ int i;
  double mean = 0.0;
  double std = 0.0;
  double nd = 0.0;
  if (weight)
  { for (i=0; i < n; i++)
    { const double x = data[i*dstride];
      const double w = weight[i*wstride];
      mean += x * w;
      std += x*x * w;
      nd += w;
    }
  }
  else
  { for (i=0; i < n; i++)
    { const double x = data[i*dstride];
      mean += x;
      std += x*x;
    }
    nd = n;
  }
  if (nd <= 1.0)
      return -1.0; /* Not enough data */
  std -= mean*mean/nd;
  std = sqrt(std/(nd-1)); /* standard deviation */
  switch (kernel)
  { case 'e': return std * pow( 80.0/(M_2_SQRTPI*nd), 0.2);
    case 'u': return std * pow( 24.0/(M_2_SQRTPI*nd), 0.2);
    case 't': return std * pow(128.0/(M_2_SQRTPI*nd), 0.2);
    case 'g': return std * pow(4.0/(3.0*nd), 0.2);
    case 'b': return std * pow(560.0/(3.0*M_2_SQRTPI*nd), 0.2);
    case '3': return std * pow(50400.0/(143.0*M_2_SQRTPI*nd), 0.2);
    case 'c':
    { const double k = M_PI*M_PI-8;
      return std * pow(M_PI,0.65)/pow(6.0*k*k*nd,0.2);
    }
  }
  return -1.0; /* Unknown kernel type */
}
