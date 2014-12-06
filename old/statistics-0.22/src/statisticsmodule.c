#include "Python.h"
#include "numpy/arrayobject.h"
#include "statistics.h"
#include <float.h>

static char buffer[512];

/* ========================================================================== */
/* -- Helper routines ------------------------------------------------------- */
/* ========================================================================== */

/* -- data ------------------------------------------------------------------ */

static PyArrayObject* parse_data(PyObject* object)
/* Takes the Python object from the argument list, and finds the 1D array
 * containing the data set. In case of an error, the function returns NULL.
 * The returned PyArrayObject* is an owned reference that should be
 * PyDECREF'ed */
{ PyArrayObject* array = NULL;
#if PY_MAJOR_VERSION < 3
  if(PyFloat_Check(object) || PyInt_Check(object) || PyLong_Check(object))
#else
  if(PyFloat_Check(object) || PyLong_Check(object))
#endif
  /* User passed a number instead of an array */
  { double value = PyFloat_AS_DOUBLE(object);
    npy_intp dims = 1;
    object = PyArray_SimpleNew(1, &dims, NPY_DOUBLE);
    if (object==NULL) return NULL;
    array = (PyArrayObject*) object;
    *((double*)PyArray_DATA(array)) = value;
    return array;
  }
  if(!PyArray_Check(object))
  /* Try to convert object to a 1D double array */
  { object = PyArray_FromObject(object, NPY_DOUBLE, 1, 1);
    if (object==NULL) return NULL;
    array = (PyArrayObject*) object;
  }
  else /* User passed an array */
  { array = (PyArrayObject*) object;
    if (PyArray_TYPE(array) != NPY_DOUBLE) /* Cast to type double */
    { object = PyArray_Cast(array, NPY_DOUBLE);
      if (!object) return NULL;
      array = (PyArrayObject*) object;
    } 
    else Py_INCREF(object);
    if (PyArray_NDIM(array) != 1) /* Checking number of dimensions */
    { char* message = strchr(buffer, '\0');
      sprintf(message,
              "array has incorrect rank (%d expected 1)",
              PyArray_NDIM(array));
      PyErr_SetString(PyExc_ValueError, buffer);
      Py_DECREF(object);
      return NULL;
    }
  }
  if (PyArray_DIM(array, 0) < 1)
  { PyErr_SetString(PyExc_ValueError, "array is empty");
    Py_DECREF(object);
    return NULL;
  }
  if (PyArray_DIM(array, 0) > INT_MAX)
  { char* message = strchr(buffer, '\0');
    sprintf(message,
            "received too many data (%" NPY_INTP_FMT " data points received)",
            PyArray_DIM(array, 0));
    PyErr_SetString (PyExc_ValueError, buffer);
    Py_DECREF(object);
    return NULL;
  }
  return array;
}

static PyArrayObject* parse_dataxy(PyObject* object)
/* Takes the Python object from the argument list, and finds the 2D array
 * containing the data set. In case of an error, the function returns NULL.
 * The returned PyArrayObject* is an owned reference that should be
 * PyDECREF'ed */
{ PyArrayObject* array = NULL;
  if(!PyArray_Check (object)) /* Try to convert object to a 2D double array */
  { object = PyArray_FromObject(object, NPY_DOUBLE, 2, 2);
    if (object==NULL) return NULL;
    array = (PyArrayObject*) object;
  }
  else /* User passed an array */
  { array = (PyArrayObject*) object;
    if (PyArray_TYPE(array) != NPY_DOUBLE) /* Cast to type double */
    { object = PyArray_Cast(array, NPY_DOUBLE);
      if (!object) return NULL;
      array = (PyArrayObject*) object;
    } 
    else Py_INCREF(object);
    if (PyArray_NDIM(array) != 2) /* Checking number of dimensions */
    { char* message = strchr(buffer, '\0');
      sprintf(message,
              "array has incorrect rank (%d expected 2)",
              PyArray_NDIM( array));
      PyErr_SetString(PyExc_ValueError, buffer);
      Py_DECREF(object);
      return NULL;
    }
  }
  if (PyArray_DIM(array, 0) < 1)
  { PyErr_SetString(PyExc_ValueError, "array is empty");
    Py_DECREF(object);
    return NULL;
  }
  if (PyArray_DIM(array, 1)!=2)
  { PyErr_SetString(PyExc_ValueError, "array should have exactly two columns");
    Py_DECREF(object);
    return NULL;
  }
  return array;
}

static char
find_mnemonic(const char arg[], const char* names[], const char mnemonics[])
{ int i;
  const char* name;
  const int n = strlen(arg);
  char result = '\0';
  char* temp = malloc((n+1)*sizeof(char));
  strcpy(temp, arg);
  for (i = 0; i < n; i++) temp[i] = tolower((int)arg[i]);
  for (i = 0; (name = names[i]); i++)
  { if (!strcmp(temp, name))
    { result = mnemonics[i];
      break;
    }
  }
  if (!result && n==1)
  { const char* p;
    const char first = temp[0];
    for (p = mnemonics ; (result = *p); p++)
      if (result==first) break;
  }
  free(temp);
  return result;
}

static char
parse_kernel(const char arg[])
{ const char* names[] = {"epanechnikov",
                         "uniform",
                         "triangle",
                         "gaussian",
                         "biweight",
                         "triweight",
                         "cosine",
                         NULL};
  const char mnemonics[] = {'e', 'u', 't', 'g', 'b', '3', 'c', '\0'};
  return find_mnemonic(arg, names, mnemonics);
}

static npy_intp parse_integer(PyObject* object)
{ long n;
#if PY_MAJOR_VERSION < 3
  if(PyInt_Check(object)) n = PyInt_AS_LONG(object);
  else if(PyLong_Check(object)) n = PyLong_AsLong(object);
#else
  if(PyLong_Check(object)) n = PyLong_AsLong(object);
  /* returns -1 if an exception is raised */
#endif
  else
  { strcat(buffer, "n should be an integer");
    PyErr_SetString (PyExc_TypeError, buffer);
    return -1;
  }
  if (n > NPY_MAX_INTP)
  { char* message = strchr(buffer, '\0');
    sprintf(message,
            "n was chosen too large (maximum is %" NPY_INTP_FMT ")",
            (npy_intp) NPY_MAX_INTP);
    PyErr_SetString (PyExc_TypeError, buffer);
    return -1;
  }
  return (npy_intp)n;
}

/* ========================================================================== */
/* -- Methods --------------------------------------------------------------- */
/* ========================================================================== */

/* pdf */
static char pdf__doc__[] =
"y, x = pdf(data, weight = None, h = None, kernel = 'Epanechnikov', n = 100)\n"
"or\n"
"y = pdf(data, x, weight = None, h = None, kernel = 'Epanechnikov')\n"
"\n"
"This function estimates the probability density function from the random\n"
"numbers in the array data, using the bandwidth h and the specified kernel\n"
"function.\n"
"\n"
"You can either explicitly specify the points x at which the probability\n"
"density function is to be calculated, or let the function do it for you.\n"
"In this case, the function will estimate the probability density function\n"
"at n equally spaced points (defaulting to 100).\n"
"\n"
"If you wish to estimate the probability density function at specific\n"
"values of x, you can pass an array of numbers via the argument x.\n"
"\n"
"Arguments:\n"
"o) The one-dimensional array data contains the sample data from which the\n"
"   probability density function is to be estimated.\n"
"o) The one-dimensional array weight contains weights for the sample data.\n"
"   If weight==None, then each data point receives an equal weight of 1.\n"
"o) Use the keyword argument 'x' to specify for which value of x you wish\n"
"   to estimate the probability density function. You can either pass a\n"
"   single value for x, or a 1D array of values. If you don't specify x,\n"
"   the function will create x as a 1D array of n values for you and return\n"
"   it together with the estimated probability density function.\n"
"o) The keyword argument 'h' specifies the bandwidth to be used for the\n"
"   estimation. If h is None (and also if the user specifies a zero or\n"
"   negative h), the optimal bandwidth is used (which can be calculated\n"
"   explicitly by the function 'bandwidth').\n"
"o) The keyword argument 'kernel' specifies the kernel function:\n"
"   -'E' or 'Epanechnikov' : Epanechnikov kernel (default)\n"
"   -'U' or 'Uniform'      : Uniform kernel\n"
"   -'T' or 'Triangle'     : Triangle kernel\n"
"   -'G' or 'Gaussian'     : Gaussian kernel\n"
"   -'B' or 'Biweight'     : Quartic/biweight kernel\n"
"   -'3' or 'Triweight'    : Triweight kernel\n"
"   -'C' or 'Cosine'       : Cosine kernel\n"
"   (case is ignored)\n"
"o) The keyword argument 'n' specifies the number of points for which the\n"
"   probability density function is to be estimated. This argument is\n"
"   meaningful only if you don't specify x explicitly; passing both x and\n"
"   n raises an error. Default value of n is 100.\n"
"\n"
"Return values:\n"
"o) If you specified x explicitly: the probability density, estimated at\n"
"   the values in x.\n"
"o) If you did not specify x explicitly: the estimated probability density,\n"
"   as well as the corresponding values of x.\n";

static PyObject*
py_pdf(PyObject* self, PyObject* args, PyObject* keywords)
{ PyObject *DATA = NULL;
  PyArrayObject* data = NULL;
  npy_intp ndata = 0;
  PyObject *WEIGHT = NULL;
  PyArrayObject* weight = NULL;
  PyObject *X = NULL;
  PyArrayObject* x;
  PyArrayObject* y;
  double h = -1.0;
  char* KERNEL = NULL;
  char kernel = 'e';
  npy_intp n = 100;
  PyObject* N = NULL;

  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = {"data", "x", "weight", "h", "kernel", "n", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|OOdsO", kwlist,
                                  &DATA, &X, &WEIGHT, &h, &KERNEL, &N))
      return NULL;

  /* Set the function name for error messages */
  strcpy (buffer, "pdf: ");

  /* -- Check the kernel variable ---------------------------------------- */
  if (KERNEL) 
  { kernel = parse_kernel(KERNEL);
    if(!kernel)
    { strcat(buffer, "unknown kernel specified");
      PyErr_SetString(PyExc_ValueError, buffer);
      return NULL;
    }
  }

  /* -- Check the data input array --------------------------------------- */
  data = parse_data(DATA);
  if (!data) return NULL;
  ndata = PyArray_DIM(data, 0);
  if (ndata < 2)
  { char* message = strchr(buffer, '\0');
    sprintf(message,
            "insufficient data (need at least two data points, received %" NPY_INTP_FMT ")",
            PyArray_DIM(data, 0));
    PyErr_SetString (PyExc_ValueError, buffer);
    Py_DECREF((PyObject*) data);
    return NULL;
  }

  /* -- Check the weight input array ------------------------------------- */
  if (WEIGHT && WEIGHT!=Py_None)
  { weight = parse_data(WEIGHT);
    if (!weight)
    { Py_DECREF((PyObject*) data);
      return NULL;
    }
    if (PyArray_DIM(weight, 0) != ndata)
    { sprintf(buffer,
             "pdf: the array weight has an incorrect size (expected %" NPY_INTP_FMT ", received %" NPY_INTP_FMT ")", ndata, PyArray_DIM(weight, 0));
      PyErr_SetString (PyExc_ValueError, buffer);
      Py_DECREF((PyObject*) data);
      return NULL;
    }
  }

  /* -- Check the bandwidth argument ------------------------------------- */
  if (h <= 0.0)
  { /* Need to choose the bandwidth ourselves */
    h = bandwidth(PyArray_DIM(data, 0),
                  (double*) PyArray_DATA(data),
                  PyArray_STRIDE(data, 0) / sizeof(double),
                  weight ? (double*) PyArray_DATA(weight) : NULL,
                  weight ? PyArray_STRIDE(weight, 0) / sizeof(double) : 0,
                  kernel);
    if (h < 0.0)
    { PyErr_SetString(PyExc_ValueError, "Not enough data to calculate the bandwidth");
      Py_DECREF((PyObject*) data);
      return NULL;
    }
  }
 
  /* -- Check if the user passed x ------------------------------------------ */
  if (X)
  { if (N)
    { PyErr_SetString(PyExc_TypeError, "Don't specify n if you specify x");
      Py_DECREF((PyObject*) data);
      return NULL;
    }
    x = parse_data(X);
    if (!x)
    { Py_DECREF((PyObject*) data);
      return NULL;
    }
    n = PyArray_DIM(x, 0);
  }
  else /* Create the x array */
  { npy_intp i;
    double xmin = DBL_MAX;
    double xmax = DBL_MIN;
    double dx;
    /* -- Create the x array ------------------------------------------------ */
    double* p = (double*) PyArray_DATA(data);
    npy_intp stride = PyArray_STRIDE(data, 0)/sizeof(double);
    for (i = 0; i < ndata; i++, p+=stride)
    { const double value = *p;
      if (value < xmin) xmin = value; 
      if (value > xmax) xmax = value; 
    }
    if (kernel=='g')
    { xmax += 3*h;
      xmin -= 3*h;
    }
    else
    { xmax += h;
      xmin -= h;
    }

    if (N)
    { n = parse_integer(N);
      if (n < 0)
      { Py_DECREF((PyObject*) data);
        return NULL;
      }
    }

    x = (PyArrayObject*) PyArray_SimpleNew(1, &n, NPY_DOUBLE);
    if (!x)
    { Py_DECREF((PyObject*) data);
      return NULL;
    }

    p = (double*) PyArray_DATA(x);
    dx = (xmax-xmin)/(n-1);
    for (i = 0; i < n; i++) p[i] = xmin + i*dx;
  }
  
  y = (PyArrayObject*) PyArray_SimpleNew(1, &n, NPY_DOUBLE);
  if (!y)
  { Py_DECREF((PyObject*) data);
    Py_DECREF((PyObject*) x);
    return NULL;
  }
  pdf((int)ndata,
      (double*) PyArray_DATA(data),
      (int) (PyArray_STRIDE(data, 0) / sizeof(double)),
      weight ? (double*) PyArray_DATA(weight) : NULL,
      weight ? (int) (PyArray_STRIDE(weight, 0) / sizeof(double)) : 0,
      h,
      kernel,
      (int) n,
      (double*) PyArray_DATA(x),
      (int) (PyArray_STRIDE(x, 0) / sizeof(double)),
      (double*) (PyArray_DATA(y)));
  Py_DECREF((PyObject*) data);
  if (weight)
  { Py_DECREF((PyObject*) weight);
    /* These braces are just to make the compiler shut up */
  }
  if (X)
  { Py_DECREF((PyObject*)x);
    return PyArray_Return(y);
  }
  return Py_BuildValue("NN", PyArray_Return(y), PyArray_Return(x));
} 
/* end of wrapper for pdf */

/* cpdf */
static char cpdf__doc__[] =
"y, x = cpdf(data, h = None, kernel = 'Epanechnikov', n = 100)\n"
"or\n"
"y = cpdf(data, x, h = None, kernel = 'Epanechnikov')\n"
"\n"
"This function estimates the cumulative probability density function from\n"
"the random numbers in the array data, using the bandwidth h and the\n"
"specified kernel function.\n"
"\n"
"You can either explicitly specify the points x at which the cumulative\n"
"probability density function is to be calculated, or let the function do\n"
"it for you. In this case, the function will estimate the cumulative\n"
"probability density function at n equally spaced points (defaulting\n"
"to 100).\n"
"\n"
"If you wish to estimate the cumulative probability density function at\n"
"specific values of x, you can pass an array of numbers via the argument x.\n"
"\n"
"Arguments:\n"
"o) The one-dimensional array data contains the sample data from which the\n"
"   cumulative probability density function is to be estimated.\n"
"o) Use the keyword argument 'x' to specify for which value of x you wish\n"
"   to estimate the cumulative probability density function. You can either\n"
"   pass a single value for x, or a 1D array of values. If you don't specify\n"
"   x, the function will create x as a 1D array of n values for you and\n"
"   return it together with the estimated cumulative probability density\n"
"   function.\n"
"o) The keyword argument 'h' specifies the bandwidth to be used for the\n"
"   estimation. If h is None (and also if the user specifies a zero or\n"
"   negative h), the optimal bandwidth is used (which can be calculated\n"
"   explicitly by the function 'bandwidth').\n"
"o) The keyword argument 'kernel' specifies the kernel function:\n"
"   -'E' or 'Epanechnikov' : Epanechnikov kernel (default)\n"
"   -'U' or 'Uniform'      : Uniform kernel\n"
"   -'T' or 'Triangle'     : Triangle kernel\n"
"   -'G' or 'Gaussian'     : Gaussian kernel\n"
"   -'B' or 'Biweight'     : Quartic/biweight kernel\n"
"   -'3' or 'Triweight'    : Triweight kernel\n"
"   -'C' or 'Cosine'       : Cosine kernel\n"
"   (case is ignored)\n"
"o) The keyword argument 'n' specifies the number of points for which the\n"
"   cumulative probability density function is to be estimated. This argument\n"
"   is meaningful only if you don't specify x explicitly; passing both x and\n"
"   n raises an error. Default value of n is 100.\n"
"\n"
"Return values:\n"
"o) If you specified x explicitly: the cumulative probability density,\n"
"   estimated at the values in x.\n"
"o) If you did not specify x explicitly: the estimated cumulative probability\n"
"   density, as well as the corresponding values of x.\n";

static PyObject*
py_cpdf(PyObject* self, PyObject* args, PyObject* keywords)
{ PyObject *DATA = NULL;
  PyArrayObject* data = NULL;
  npy_intp ndata = 0;
  PyObject *X = NULL;
  PyArrayObject* x;
  PyArrayObject* y;
  double h = -1.0;
  char* KERNEL = NULL;
  char kernel = 'e';
  npy_intp n = 100;
  PyObject* N = NULL;

  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = {"data", "x", "h", "kernel", "n", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|OdsO", kwlist,
                                  &DATA, &X, &h, &KERNEL, &N)) return NULL;

  /* Set the function name for error messages */
  strcpy (buffer, "cpdf: ");

  /* -- Check the kernel variable ---------------------------------------- */
  if (KERNEL) 
  { kernel = parse_kernel(KERNEL);
    if(!kernel)
    { strcat(buffer, "unknown kernel specified");
      PyErr_SetString(PyExc_ValueError, buffer);
      return NULL;
    }
#ifndef HAVE_ERF
    if (kernel=='g')
    { PyErr_SetString(PyExc_RuntimeError, "Gaussian kernel unavailable in cpdf due to missing erf function in this computer's math library. Try reinstalling Statistics for Python; locate the erf function by running 'python setup.py config' before build, install");
      return NULL;
    }
#endif
  }
  /* -- Check the data input array --------------------------------------- */
  data = parse_data(DATA);
  if (!data) return NULL;
  ndata = PyArray_DIM(data, 0);
  if (ndata < 2)
  { char* message = strchr(buffer, '\0');
    sprintf(message,
            "insufficient data (need at least two data points, received %" NPY_INTP_FMT ")", PyArray_DIM(data, 0));
    PyErr_SetString (PyExc_ValueError, buffer);
    Py_DECREF((PyObject*) data);
    return NULL;
  }

  /* -- Check the bandwidth argument ------------------------------------- */
  if (h <= 0.0)
  { /* Need to choose the bandwidth ourselves */
    h = bandwidth(ndata,
                  (double*) PyArray_DATA(data),
                  (int) (PyArray_STRIDE(data, 0) / sizeof(double)),
                  NULL,
                  0,
                  kernel);
    if (h < 0.0)
    { PyErr_SetString(PyExc_ValueError, "cpdf: not enough data to calculate the bandwidth");
      Py_DECREF((PyObject*) data);
      return NULL;
    }
  }
 
  /* -- Check if the user passed x ------------------------------------------ */
  if (X)
  { if (N)
    { PyErr_SetString(PyExc_TypeError, "cpdf: don't specify n if you specify x");
      Py_DECREF((PyObject*) data);
      return NULL;
    }
    x = parse_data(X);
    if (!x)
    { Py_DECREF((PyObject*) data);
      return NULL;
    }
    n = PyArray_DIM(x, 0);
  }
  else /* Create the x array */
  { npy_intp i;
    double xmin = DBL_MAX;
    double xmax = DBL_MIN;
    double dx;
    /* -- Create the x array ------------------------------------------------ */
    double* p = (double*) PyArray_DATA(data);
    int stride = (int) (PyArray_STRIDE(data, 0)/sizeof(double));
    for (i = 0; i < ndata; i++, p+=stride)
    { const double value = *p;
      if (value < xmin) xmin = value; 
      if (value > xmax) xmax = value; 
    }
    if (kernel=='g')
    { xmax += 3*h;
      xmin -= 3*h;
    }
    else
    { xmax += h;
      xmin -= h;
    }

    if (N)
    { n = parse_integer(N);
      if (n < 0)
      { Py_DECREF((PyObject*) data);
        return NULL;
      }
    }

    x = (PyArrayObject*) PyArray_SimpleNew(1, &n, NPY_DOUBLE);
    if (!x)
    { Py_DECREF((PyObject*) data);
      return NULL;
    }

    p = (double*) PyArray_DATA(x);
    dx = (xmax-xmin)/(n-1);
    for (i = 0; i < n; i++) p[i] = xmin + i*dx;
  }
  
  y = (PyArrayObject*) PyArray_SimpleNew(1, &n, NPY_DOUBLE);
  if (!y)
  { Py_DECREF((PyObject*) data);
    Py_DECREF((PyObject*) x);
    return NULL;
  }
  cpdf((int)ndata,
       (double*) PyArray_DATA(data),
       (int) (PyArray_STRIDE(data, 0) / sizeof(double)),
       h,
       kernel,
       (int) n,
       (double*) PyArray_DATA(x),
       (int) (PyArray_STRIDE(x, 0)/ sizeof(double)),
       (double*) PyArray_DATA(y));
  Py_DECREF((PyObject*) data);
  if (X)
  { Py_DECREF((PyObject*)x);
    return PyArray_Return(y);
  }
  return Py_BuildValue("NN", PyArray_Return(y), PyArray_Return(x));
} 
/* end of wrapper for cpdf */

/* cpdfc */
static char cpdfc__doc__[] =
"y, x = cpdfc(data, h = None, kernel = 'Epanechnikov', n = 100)\n"
"or\n"
"y = cpdfc(data, x, h = None, kernel = 'Epanechnikov')\n"
"\n"
"This function estimates the complement of the cumulative probability\n"
"density function from the random numbers in the array data, using the\n"
"bandwidth h and the specified kernel function.\n"
"\n"
"You can either explicitly specify the points x at which the complement\n"
"of the cumulative probability density function is to be calculated, or let\n"
"the function do it for you. In this case, the function will estimate the\n"
"complement of the cumulative probability density function at n equally\n"
"spaced points (defaulting to 100).\n"
"\n"
"If you wish to estimate the complement of the cumulative probability\n"
"density function at specific values of x, you can pass an array of numbers\n"
"via the argument x.\n"
"\n"
"Arguments:\n"
"o) The one-dimensional array data contains the sample data from which the\n"
"   complement of the cumulative probability density function is to be\n"
"   estimated.\n"
"o) Use the keyword argument 'x' to specify for which value of x you wish\n"
"   to estimate the complement of the cumulative probability density\n"
"   function. You can either pass a single value for x, or a 1D array of\n"
"   values. If you don't specify x, the function will create x as a 1D\n"
"   array of n values for you and return it together with the estimated\n"
"   complement of the cumulative probability density function.\n"
"o) The keyword argument 'h' specifies the bandwidth to be used for the\n"
"   estimation. If h is None (and also if the user specifies a zero or\n"
"   negative h), the optimal bandwidth is used (which can be calculated\n"
"   explicitly by the function 'bandwidth').\n"
"o) The keyword argument 'kernel' specifies the kernel function:\n"
"   -'E' or 'Epanechnikov' : Epanechnikov kernel (default)\n"
"   -'U' or 'Uniform'      : Uniform kernel\n"
"   -'T' or 'Triangle'     : Triangle kernel\n"
"   -'G' or 'Gaussian'     : Gaussian kernel\n"
"   -'B' or 'Biweight'     : Quartic/biweight kernel\n"
"   -'3' or 'Triweight'    : Triweight kernel\n"
"   -'C' or 'Cosine'       : Cosine kernel\n"
"   (case is ignored)\n"
"o) The keyword argument 'n' specifies the number of points for which the\n"
"   complement of the cumulative probability density function is to be\n"
"   estimated. This argument is meaningful only if you don't specify x\n"
"   explicitly; passing both x and n raises an error. Default value of n\n"
"   is 100.\n"
"\n"
"Return values:\n"
"o) If you specified x explicitly: the complement of the cumulative\n"
"   probability density, estimated at the values in x.\n"
"o) If you did not specify x explicitly: the estimated complement of the\n"
"   cumulative probability density, as well as the corresponding values\n"
"   of x.\n";

static PyObject*
py_cpdfc(PyObject* self, PyObject* args, PyObject* keywords)
{ PyObject *DATA = NULL;
  PyArrayObject* data = NULL;
  npy_intp ndata = 0;
  PyObject *X = NULL;
  PyArrayObject* x;
  PyArrayObject* y;
  double h = -1.0;
  char* KERNEL = NULL;
  char kernel = 'e';
  npy_intp n = 100;
  PyObject* N = NULL;

  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = {"data", "x", "h", "kernel", "n", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|OdsO", kwlist,
                                  &DATA, &X, &h, &KERNEL, &N)) return NULL;

  /* Set the function name for error messages */
  strcpy (buffer, "cpdfc: ");

  /* -- Check the kernel variable ---------------------------------------- */
  if (KERNEL) 
  { kernel = parse_kernel(KERNEL);
    if(!kernel)
    { strcat(buffer, "unknown kernel specified");
      PyErr_SetString(PyExc_ValueError, buffer);
      return NULL;
    }
#ifndef HAVE_ERFC
    if (kernel=='g')
    { PyErr_SetString(PyExc_RuntimeError, "cpdfc: Gaussian kernel unavailable in cpdf due to missing erfc function in this computer's math library. Try reinstalling Statistics for Python; locate the erfc function by running 'python setup.py config' before build, install");
      return NULL;
    }
#endif
  }
  /* -- Check the data input array --------------------------------------- */
  data = parse_data(DATA);
  if (!data) return NULL;
  ndata = PyArray_DIM(data, 0);
  if (ndata < 2)
  { sprintf(buffer, "insufficient data (need at least two data points, received %" NPY_INTP_FMT ")", PyArray_DIM(data, 0));
    PyErr_SetString(PyExc_ValueError, buffer);
    Py_DECREF((PyObject*) data);
    return NULL;
  }

  /* -- Check the bandwidth argument ------------------------------------- */
  if (h <= 0.0)
  { /* Need to choose the bandwidth ourselves */
    h = bandwidth((int)ndata,
                  (double*) PyArray_DATA(data),
                  (int) (PyArray_STRIDE(data, 0) / sizeof(double)),
                  NULL,
                  0,
                  kernel);
    if (h < 0.0)
    { PyErr_SetString(PyExc_ValueError, "cpdfc: not enough data to calculate the bandwidth");
      Py_DECREF((PyObject*) data);
      return NULL;
    }
  }
 
  /* -- Check if the user passed x ------------------------------------------ */
  if (X)
  { if (N)
    { PyErr_SetString(PyExc_TypeError, "cpdfc: don't specify n if you specify x");
      Py_DECREF((PyObject*) data);
      return NULL;
    }
    x = parse_data(X);
    if (!x)
    { Py_DECREF((PyObject*) data);
      return NULL;
    }
    n = PyArray_DIM(x, 0);
  }
  else /* Create the x array */
  { npy_intp i;
    double xmin = DBL_MAX;
    double xmax = DBL_MIN;
    double dx;
    /* -- Create the x array ------------------------------------------------ */
    double* p = (double*) PyArray_DATA(data);
    int stride = (int) (PyArray_STRIDE(data, 0)/sizeof(double));
    for (i = 0; i < ndata; i++, p+=stride)
    { const double value = *p;
      if (value < xmin) xmin = value; 
      if (value > xmax) xmax = value; 
    }
    if (kernel=='g')
    { xmax += 3*h;
      xmin -= 3*h;
    }
    else
    { xmax += h;
      xmin -= h;
    }

    if (N)
    { n = parse_integer(N);
      if (n < 0)
      { Py_DECREF((PyObject*) data);
        return NULL;
      }
    }

    x = (PyArrayObject*) PyArray_SimpleNew(1, &n, NPY_DOUBLE);
    if (!x)
    { Py_DECREF((PyObject*) data);
      return NULL;
    }

    p = (double*) PyArray_DATA(x);
    dx = (xmax-xmin)/(n-1);
    for (i = 0; i < n; i++) p[i] = xmin + i*dx;
  }
  
  y = (PyArrayObject*) PyArray_SimpleNew(1, &n, NPY_DOUBLE);
  if (!y)
  { Py_DECREF((PyObject*) data);
    Py_DECREF((PyObject*) x);
    return NULL;
  }
  cpdfc(ndata,
        (double*) PyArray_DATA(data),
        (int) (PyArray_STRIDE(data, 0) / sizeof(double)),
        h,
        kernel,
        (int) n,
        (double*) PyArray_DATA(x),
        (int) (PyArray_STRIDE(x, 0) / sizeof(double)),
        (double*) PyArray_DATA(y));
  Py_DECREF((PyObject*) data);
  if (X)
  { Py_DECREF((PyObject*)x);
    return PyArray_Return(y);
  }
  return Py_BuildValue("NN", PyArray_Return(y), PyArray_Return(x));
} 
/* end of wrapper for cpdfc */

/* localfit */
static char localfit__doc__[] =
"y, x = localfit(data, h = None, kernel = 'Epanechnikov', n = 100)\n"
"or\n"
"y = localfit(data, x, h = None, kernel = 'Epanechnikov')\n"
"\n"
"This function estimates the probability density function from the random\n"
"numbers in the array data, using the bandwidth h and the specified kernel\n"
"function.\n"
"\n"
"You can either explicitly specify the points x at which the probability\n"
"density function is to be calculated, or let the function do it for you.\n"
"In this case, the function will estimate the probability density function\n"
"at n equally spaced points (defaulting to 100).\n"
"\n"
"If you wish to estimate the probability density function at specific\n"
"values of x, you can either pass a single number or an array of numbers\n"
"via the argument x.\n"
"\n"
"Arguments:\n"
"o) The one-dimensional array data contains the sample data from which the\n"
"   probability density function is estimated.\n"
"o) Use the keyword argument 'x' to specify for which value of x you wish\n"
"   to estimate the probability density function. You can either pass a\n"
"   single value for x, or a 1D array of values. If you don't specify x,\n"
"   the function will create x as a 1D array of n values for you and return\n"
"   it together with the estimated probability density function.\n"
"o) The keyword argument 'h' specifies the bandwidth to be used for the\n"
"   estimation. If h is None (and also if the user specifies a zero or\n"
"   negative h), the optimal bandwidth is used (which can be calculated\n"
"   explicitly by the function 'bandwidth').\n"
"o) The keyword argument 'kernel' specifies the kernel function:\n"
"   -'E' or 'Epanechnikov' : Epanechnikov kernel (default)\n"
"   -'U' or 'Uniform'      : Uniform kernel\n"
"   -'T' or 'Triangle'     : Triangle kernel\n"
"   -'G' or 'Gaussian'     : Gaussian kernel\n"
"   -'B' or 'Biweight'     : Quartic/biweight kernel\n"
"   -'3' or 'Triweight'    : Triweight kernel\n"
"   -'C' or 'Cosine'       : Cosine kernel\n"
"   (case is ignored)\n"
"o) The keyword argument 'n' specifies the number of points between\n"
"   for which the probability density function is to be estimated. This\n"
"   argument is meaningful only if you don't specify x explicitly; passing\n"
"   both x and n raises an error. Default value of n is 100.\n"
"\n"
"Return values:\n"
"o) If you specified x explicitly: the probability density, estimated at\n"
"   the value(s) of x.\n"
"o) If you did not specify x explicitly: the estimated probability density,\n"
"   as well as the corresponding values of x.\n";

static PyObject*
py_localfit(PyObject* self, PyObject* args, PyObject* keywords)
{ PyObject *DATA = NULL;
  PyArrayObject* data = NULL;
  npy_intp ndata = 0;
  int strides[2];
  PyObject *X = NULL;
  PyArrayObject* x;
  PyArrayObject* y;
  double h = +1.0;
  char* KERNEL = NULL;
  char kernel = 'e';
  npy_intp n = 100;
  PyObject* N = NULL;

  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = {"data", "x", "h", "kernel", "n", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|OdsO", kwlist,
                                  &DATA, &X, &h, &KERNEL, &N)) return NULL;

  /* Set the function name for error messages */
  strcpy (buffer, "pdf: ");

  /* -- Check the kernel variable ---------------------------------------- */
  if (KERNEL) 
  { kernel = parse_kernel(KERNEL);
    if(!kernel)
    { PyErr_SetString(PyExc_ValueError,
                      "localfit: unknown kernel specified");
      return NULL;
    }
  }
  /* -- Check the data input array --------------------------------------- */
  data = parse_dataxy(DATA);
  if (!data) return NULL;
  ndata = PyArray_DIM(data, 0);

  /* -- Check if the user passed x ------------------------------------------ */
  if (X)
  { if (N)
    { PyErr_SetString(PyExc_TypeError, "Don't specify n if you specify x");
      Py_DECREF((PyObject*) data);
      return NULL;
    }
    x = parse_data(X);
    if (!x)
    { Py_DECREF((PyObject*) data);
      return NULL;
    }
    n = PyArray_DIM(x, 0);
  }
  else /* Create the x array */
  { npy_intp i;
    double xmin = DBL_MAX;
    double xmax = DBL_MIN;
    double dx;
    /* -- Create the x array ------------------------------------------------ */
    double* p = (double*) PyArray_DATA(data);
    int stride = (int) (PyArray_STRIDE(data, 0)/sizeof(double));
    for (i = 0; i < ndata; i++, p+=stride)
    { const double value = *p;
      if (value < xmin) xmin = value; 
      if (value > xmax) xmax = value; 
    }
    if (kernel=='g')
    { xmax += 3*h;
      xmin -= 3*h;
    }
    else
    { xmax += h;
      xmin -= h;
    }

    if (N)
    { n = parse_integer(N);
      if (n < 0)
      { Py_DECREF((PyObject*) data);
        return NULL;
      }
    }

    x = (PyArrayObject*) PyArray_SimpleNew(1, &n, NPY_DOUBLE);
    if (!x)
    { Py_DECREF((PyObject*) data);
      return NULL;
    }

    p = (double*) PyArray_DATA(x);
    dx = (xmax-xmin)/(n-1);
    for (i = 0; i < n; i++) p[i] = xmin + i*dx;
  }
  
  y = (PyArrayObject*) PyArray_SimpleNew(1, &n, NPY_DOUBLE);
  if (!y)
  { Py_DECREF((PyObject*) data);
    Py_DECREF((PyObject*) x);
    return NULL;
  }
  strides[0] = (int) (PyArray_STRIDE(data, 0) / sizeof(double));
  strides[1] = (int) (PyArray_STRIDE(data, 1) / sizeof(double));
  localfit((int)ndata,
           (double*) PyArray_DATA(data),
           strides,
           h,
           kernel,
           (int) n,
           (double*) PyArray_DATA(x),
           (int) (PyArray_STRIDE(x, 0) / sizeof(double)),
           (double*) (PyArray_DATA(y)));
  Py_DECREF((PyObject*)data);
  if (X)
  { Py_DECREF((PyObject*)x);
    return PyArray_Return(y);
  }
  return Py_BuildValue("NN", PyArray_Return(y), PyArray_Return(x));
} 
/* end of wrapper for localfit */

/* bandwidth */
static char bandwidth__doc__[] =
"bandwidth(data, weight = None, kernel='Epanechnikov') -> Float"
"\n"
"This function calculates the optimal bandwidth to be used for the kernel-\n"
"based estimation of the probability density function, using the kernel\n"
"function specified by the keyword argument 'kernel'.\n"
"\n"
"Arguments:\n"
"o) The one-dimensional array data contains the sample data from which the\n"
"   probability density function is calculated\n"
"o) The one-dimensional array weight contains weights for the sample data.\n"
"   If weight==None, then each data point receives an equal weight of 1.\n"
"o) The keyword argument 'kernel' specifies the kernel function:\n"
"   -'E' or 'Epanechnikov' : Epanechnikov kernel (default)\n"
"   -'U' or 'Uniform'      : Uniform kernel\n"
"   -'T' or 'Triangle'     : Triangle kernel\n"
"   -'G' or 'Gaussian'     : Gaussian kernel\n"
"   -'B' or 'Biweight'     : Quartic/biweight kernel\n"
"   -'3' or 'Triweight'    : Triweight kernel\n"
"   -'C' or 'Cosine'       : Cosine kernel\n"
"   (case is ignored)\n"
"\n"
"Return values:\n"
"The optimal bandwidth for the data, using the given kernel.\n"
"This bandwidth can then be used when estimating the probability density\n"
"function with the function pdf, or the cumulative probability density\n"
"function with the function cpdf or cpdfc.\n";

static PyObject*
py_bandwidth(PyObject* self, PyObject* args, PyObject* keywords)
{ PyObject *DATA = NULL;
  PyObject *WEIGHT = NULL;
  PyArrayObject* weight = NULL;
  char* KERNEL = NULL;
  char kernel = 'e';
  npy_intp ndata;
  PyArrayObject* data = NULL;
  double result;

  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = {"data", "weight", "kernel", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "OO|s", kwlist,
                                  &DATA,
                                  &WEIGHT,
                                  &KERNEL)) return NULL;
  /* Set the function name for error messages */
  strcpy (buffer, "bandwidth: ");

  /* -- Check the kernel variable ---------------------------------------- */
  if (KERNEL) 
  { kernel = parse_kernel(KERNEL);
    if(!kernel)
    { strcat(buffer, "unknown kernel specified");
      PyErr_SetString(PyExc_ValueError, buffer);
      return NULL;
    }
  }
  /* -- Check the data input array --------------------------------------- */
  data = parse_data(DATA);
  if (!data) return NULL;
  ndata = PyArray_DIM(data, 0);
  if (ndata < 2)
  { char* message = strchr(buffer, '\0');
    sprintf(message, "insufficient data (need at least two data points, received %" NPY_INTP_FMT ")", ndata);
    PyErr_SetString (PyExc_MemoryError, buffer);
    Py_DECREF((PyObject*) data);
    return NULL;
  }

  /* -- Check the weight input array ------------------------------------- */
  if (WEIGHT && WEIGHT!=Py_None)
  { weight = parse_data(WEIGHT);
    if (!weight)
    { Py_DECREF((PyObject*) data);
      return NULL;
    }
    if (PyArray_DIM(weight, 0) != ndata)
    { sprintf(buffer,
             "pdf: the array weight has an incorrect size (expected %" NPY_INTP_FMT ", received %" NPY_INTP_FMT ")", ndata, PyArray_DIM(weight, 0));
      PyErr_SetString (PyExc_ValueError, buffer);
      Py_DECREF((PyObject*) data);
      return NULL;
    }
  }

  /* --------------------------------------------------------------------- */
  result = bandwidth(ndata,
                     (double*) PyArray_DATA(data),
                     (int) (PyArray_STRIDE(data, 0) / sizeof(double)),
                      weight ? (double*) PyArray_DATA(weight) : NULL,
                      weight ? PyArray_STRIDE(weight, 0) / sizeof(double) : 0,
                     kernel);
  Py_DECREF((PyObject*) data);

  if (result < 0.0)
  { PyErr_SetString(PyExc_ValueError, "Not enough data to calculate the bandwidth");
    Py_DECREF((PyObject*) data);
    return NULL;
  }
  return PyFloat_FromDouble(result);
} 
/* end of wrapper for bandwidth */

/* mean */
static char mean__doc__[] =
"mean (data)\n"
"This function returns the mean of the 1D array data.\n";

static PyObject*
py_mean (PyObject* unused, PyObject* args)
{ double result;
  PyObject* DATA = NULL;
  PyArrayObject* aDATA = NULL;

  /* -- Read the input variables ----------------------------------------- */
  if(!PyArg_ParseTuple(args, "O", &DATA)) return NULL;

  /* -- Check the input variable ----------------------------------------- */
#if PY_MAJOR_VERSION < 3
  if(PyFloat_Check(DATA) || PyInt_Check(DATA) || PyLong_Check(DATA))
#else
  if(PyFloat_Check(DATA) || PyLong_Check(DATA))
#endif
  { Py_INCREF(DATA);
    return DATA;
  }
  if(!PyArray_Check (DATA))
  { aDATA = (PyArrayObject *) PyArray_ContiguousFromObject(DATA, PyArray_NOTYPE, 0, 0);
    if (!aDATA) return NULL;
  }
  else
  { aDATA = (PyArrayObject*) DATA;
    Py_INCREF(DATA);
  }
  if (PyArray_TYPE(aDATA) != NPY_DOUBLE)
  { PyObject* av = PyArray_Cast (aDATA, NPY_DOUBLE);
    Py_DECREF((PyObject*) aDATA);
    aDATA = (PyArrayObject*) av;
    if (!aDATA) return NULL;
  } 
  if ((PyArray_NDIM(aDATA) != 1) &&
      (PyArray_NDIM(aDATA) > 0 || PyArray_DIM(aDATA, 0) != 1))
  { sprintf(buffer,
            "mean: Argument has incorrect rank (%d expected 1).",
            PyArray_NDIM(aDATA));
    PyErr_SetString(PyExc_ValueError, buffer);
    Py_DECREF((PyObject*) aDATA);
    return NULL;
  }
  if (PyArray_DIM(aDATA, 0) > INT_MAX)
  { sprintf(buffer,
            "mean: Received more data than this function can handle (%" NPY_INTP_FMT " data points received)",
            PyArray_DIM(aDATA, 0));
    PyErr_SetString (PyExc_ValueError, buffer);
    Py_DECREF((PyObject*) aDATA);
    return NULL;
  }

  if (!PyArray_ISCONTIGUOUS(aDATA))
  { PyObject* av =
      PyArray_ContiguousFromObject((PyObject*) aDATA, PyArray_TYPE(aDATA), 0, 0);
    Py_DECREF((PyObject*)aDATA);
    if(!av) return NULL;
    aDATA = (PyArrayObject*) av;
  }
  /* --------------------------------------------------------------------- */
  result = mean((int) PyArray_DIM(aDATA, 0),
                (double*) PyArray_DATA(aDATA));
  /* --------------------------------------------------------------------- */
  Py_DECREF((PyObject*) aDATA);
  /* --------------------------------------------------------------------- */
  return PyFloat_FromDouble(result);
} 
/* end of wrapper for mean */

/* median */
static char median__doc__[] =
"median (data)\n"
"This function returns the median of the 1D array data.\n"
"Note: data will be partially ordered upon return.\n";

static PyObject*
py_median (PyObject* unused, PyObject* args)
{ double result;
  PyObject* DATA = NULL;
  PyArrayObject* aDATA = NULL;

  /* -- Read the input variables ----------------------------------------- */
  if(!PyArg_ParseTuple(args, "O", &DATA)) return NULL;

  /* -- Check the input variable ----------------------------------------- */
#if PY_MAJOR_VERSION < 3
  if(PyFloat_Check(DATA) || PyInt_Check(DATA) || PyLong_Check(DATA))
#else
  if(PyFloat_Check(DATA) || PyLong_Check(DATA))
#endif
  { Py_INCREF(DATA);
    return DATA;
  }
  if(!PyArray_Check (DATA))
  { aDATA = (PyArrayObject *) PyArray_ContiguousFromObject(DATA, PyArray_NOTYPE, 0, 0);
    if (!aDATA) return NULL;
  }
  else
  { aDATA = (PyArrayObject*) DATA;
    Py_INCREF(DATA);
  }
  if (PyArray_TYPE(aDATA) != NPY_DOUBLE)
  { PyObject* av = PyArray_Cast (aDATA, NPY_DOUBLE);
    Py_DECREF((PyObject*) aDATA);
    aDATA = (PyArrayObject*) av;
    if (!aDATA) return NULL;
  } 
  if ((PyArray_NDIM(aDATA) != 1) &&
      (PyArray_NDIM(aDATA) > 0 || PyArray_DIM(aDATA, 0) != 1))
  { sprintf(buffer,
            "median: Argument has incorrect rank (%d expected 1).",
            PyArray_NDIM(aDATA));
    PyErr_SetString(PyExc_ValueError, buffer);
    Py_DECREF((PyObject*) aDATA);
    return NULL;
  }
  if (PyArray_DIM(aDATA, 0) > INT_MAX)
  { sprintf(buffer,
            "median: Received more data than this function can handle (%" NPY_INTP_FMT " data points received)",
            PyArray_DIM(aDATA, 0));
    PyErr_SetString (PyExc_ValueError, buffer);
    Py_DECREF((PyObject*) aDATA);
    return NULL;
  }

  if (!PyArray_ISCONTIGUOUS(aDATA))
  { PyObject* av =
      PyArray_ContiguousFromObject((PyObject*) aDATA, PyArray_TYPE(aDATA), 0, 0);
    Py_DECREF((PyObject*)aDATA);
    if(!av) return NULL;
    aDATA = (PyArrayObject*) av;
  }
  /* --------------------------------------------------------------------- */
  result = median((int)PyArray_DIM(aDATA, 0),
                  (double*)PyArray_DATA(aDATA));
  /* --------------------------------------------------------------------- */
  Py_DECREF((PyObject*) aDATA);
  /* --------------------------------------------------------------------- */
  return PyFloat_FromDouble(result);
} 
/* end of wrapper for median */

/* variance */
static char variance__doc__[] =
"variance(x, mode = 'Unbiased')\n"
"\n"
"This function estimates the variance in the one-dimensional array x.\n"
"\n"
"Arguments:\n"
"o) The array x is a one-dimensional array containing the data for which\n"
"   to calculate the variance.\n"
"o) For mode=='Unbiased', which is the default, this function calculates the\n"
"   unbiased estimate of the variance. If mode=='ML', the maximum-\n"
"   likelihood estimate is used, which is a biased estimate.\n"
"\n"
"Return values:\n"
"o) The variance in x\n";


static PyObject*
py_variance(PyObject* self, PyObject* args, PyObject* keywords)
{ int n;
  PyObject* X = NULL;
  PyArrayObject* aX = NULL;
  char* MODE = NULL;
  char mode = 'u';
  int stride;
  double result;

  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = {"x", "mode", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|s", kwlist,
                                  &X, &MODE)) return NULL;

  /* -- Check the input variables ---------------------------------------- */
  if (MODE) 
  { const char* names[] = {"unbiased", "ml", NULL};
    const char mnemonics[] = {'u', 'm', '\0'};
    mode = find_mnemonic(MODE, names, mnemonics);
    if(!mode)
    { PyErr_SetString(PyExc_ValueError,
                      "unknown mode specified (should be 'unbiased' or 'ml')");
      return NULL;
    }
  }
  /* --------------------------------------------------------------------- */
  if(!PyArray_Check(X))
  { X = PyArray_ContiguousFromObject(X, NPY_DOUBLE, 1, 1);
    if (!X) return NULL;
    aX = (PyArrayObject*)X;
  }
  else
  { aX = (PyArrayObject*)X;
    if (PyArray_TYPE(aX) != NPY_DOUBLE)
    { X = PyArray_Cast(aX, NPY_DOUBLE);
      if (!X) return NULL;
      aX = (PyArrayObject*)X;
    }
    else Py_INCREF(X);
  }
  if (PyArray_NDIM(aX) != 1)
  { sprintf(buffer,
            "variance: Argument x has incorrect rank (%d expected 1).",
            PyArray_NDIM(aX));
    PyErr_SetString(PyExc_ValueError, buffer);
    Py_DECREF(X);
    return NULL;
  }
  if (PyArray_DIM(aX, 0) > INT_MAX)
  { sprintf(buffer,
            "variance: Received more data than this function can handle (%" NPY_INTP_FMT " data points received)",
            PyArray_DIM(aX, 0));
    PyErr_SetString (PyExc_ValueError, buffer);
    Py_DECREF(X);
    return NULL;
  }

  n = (int) PyArray_DIM(aX, 0);
  stride = (int) (PyArray_STRIDE(aX, 0) / sizeof(double));
  result = variance(n, (double*)PyArray_DATA(aX), stride, mode);
  Py_DECREF(X);
  return PyFloat_FromDouble(result);
} 
/* end of wrapper for variance */

/* covariance */
static char covariance__doc__[] =
"covariance(x, y = None, mode = 'Unbiased')\n"
"\n"
"This function estimates the covariance between x and y (if both are one-\n"
"dimensional have have the same length), or between the columns of x (if x\n"
"is a two-dimensional array of data).\n"
"\n"
"Arguments:\n"
"o) The array x should be either one-dimensional or two-dimensional. If x\n"
"   is one-dimensional and y==None, then this function returns the variance\n"
"   in x. If both x and y are one-dimensional and have the same length,\n"
"   this function returns the covariance between x and y. If x is two-\n"
"   dimensional, then this function calculates the covariance matrix between\n"
"   the columns in x. In this case, y is ignored.\n"
"o) If passed, the array y should be one-dimensional and contain the same\n"
"   number of elements as x; the function will return the covariance between\n"
"   x and y in this case.\n"
"o) For mode=='Unbiased', which is the default, this function calculates the\n"
"   unbiased estimate of the covariance. If mode=='ML', the maximum-\n"
"   likelihood estimate is used, which is a biased estimate.\n"
"\n"
"Return values:\n"
"o) If x is a one-dimensional array and y==None: the variance in x\n"
"o) If both x and y are one-dimensional arrays of the same length: the\n"
"   covariance between the variables x and y.\n"
"o) If x is a two-dimensional array: the covariance matrix between the\n"
"   columns in x. Element [i,j] in the covariance matrix contains the\n"
"   covariance between the columns x[:,i] and x[:,j].\n";


static PyObject*
py_covariance(PyObject* self, PyObject* args, PyObject* keywords)
{ npy_intp n, m;
  PyObject* X = NULL;
  PyArrayObject* aX = NULL;
  PyObject* Y = NULL;
  PyArrayObject* aY = NULL;
  char* MODE = NULL;
  char mode = 'u';
  int ok;

  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = {"x", "y", "mode", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|Os", kwlist,
                                  &X, &Y, &MODE)) return NULL;

  /* -- Check the input variables ---------------------------------------- */
  if (MODE) 
  { const char* names[] = {"unbiased", "ml", NULL};
    const char mnemonics[] = {'u', 'm', '\0'};
    mode = find_mnemonic(MODE, names, mnemonics);
    if(!mode)
    { strcat(buffer, "unknown mode specified (should be 'unbiased' or 'ml')");
      PyErr_SetString(PyExc_ValueError, buffer);
      return NULL;
    }
  }
  /* --------------------------------------------------------------------- */
  if(!PyArray_Check(X))
  { X = PyArray_ContiguousFromObject(X, NPY_DOUBLE, 1, 2);
    if (!X) return NULL;
    aX = (PyArrayObject*)X;
  }
  else
  { aX = (PyArrayObject*)X;
    if (PyArray_TYPE(aX) != NPY_DOUBLE)
    { X = PyArray_Cast(aX, NPY_DOUBLE);
      if (!X) return NULL;
      aX = (PyArrayObject*)X;
    }
    else Py_INCREF(X);
  }
  switch (PyArray_NDIM(aX))
  { case 2:
      if(!Y)
      { PyArrayObject* matrix = NULL;
        npy_intp dims[2];
        int strides[2];
        strides[0] = (int) (PyArray_STRIDE(aX, 0) / sizeof(double));
        strides[1] = (int) (PyArray_STRIDE(aX, 1) / sizeof(double));
        n = PyArray_DIM(aX, 0);
        m = PyArray_DIM(aX, 1);
        if (n > INT_MAX || m > INT_MAX)
        { sprintf(buffer,
                  "covariance: Received more data than this function can handle (received a %" NPY_INTP_FMT " x %" NPY_INTP_FMT "matrix)",
                  n, m);
          PyErr_SetString (PyExc_ValueError, buffer);
          Py_DECREF(X);
          return NULL;
        }
        dims[0] = m;
        dims[1] = m;
        matrix = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        if(!matrix)
        { Py_DECREF(X);
          return NULL;
        }
        ok = covariance((int) n,
                        (int) m,
                        (double*)PyArray_DATA(aX),
                        strides,
                        mode,
                        (double*)PyArray_DATA(matrix));
        Py_DECREF(X);
        if(!ok)
        { PyErr_SetString(PyExc_MemoryError, "covariance: failed to calculate the covariance matrix.");
          return NULL;
        }
        return PyArray_Return(matrix);
      }
      else
      { PyErr_SetString(PyExc_TypeError,
                        "covariance: Argument y passed, but argument x is two-dimensional.");
        Py_DECREF(X);
        return NULL;
      }
    case 1:
      if(!Y)
      { double value;
        const int stride = (int) (PyArray_STRIDE(aX, 0) / sizeof(double));
        m = PyArray_DIM(aX, 0);
        if (m > INT_MAX)
        { sprintf(buffer,
                  "covariance: Received more data than this function can handle (received %" NPY_INTP_FMT " data points)",
                  m);
          PyErr_SetString (PyExc_ValueError, buffer);
          Py_DECREF(X);
          return NULL;
        }
        value = variance((int)m, (double*)PyArray_DATA(aX), stride, mode);
        Py_DECREF(X);
        return PyFloat_FromDouble(value);
      }
      else
      { double* data;
        double matrix[4];
        int strides[2];
        n = PyArray_DIM(aX, 0);
        m = 2;
        if(!PyArray_Check(Y))
        { Y = PyArray_ContiguousFromObject(Y, NPY_DOUBLE, 1, 1);
          if (!Y)
          { Py_DECREF(X);
            return NULL;
          }
          aY = (PyArrayObject*)Y;
        }
        else
        { aY = (PyArrayObject*)Y;
          if (PyArray_NDIM(aY) != 1)
          { sprintf(buffer,
                    "covariance: Argument y has incorrect rank (%d expected 1).",
                    PyArray_NDIM(aY));
            PyErr_SetString(PyExc_ValueError, buffer);
            Py_DECREF(X);
            return NULL;
          }
          if (PyArray_TYPE(aY) != NPY_DOUBLE)
          { Y = PyArray_Cast(aY, NPY_DOUBLE);
            if (!Y)
            { Py_DECREF(X);
              return NULL;
            }
            aY = (PyArrayObject*)Y;
          } 
          else Py_INCREF(Y);
          /* X and Y should have equal strides. Make both strides equal to one. */
          if (!PyArray_ISCONTIGUOUS(aX))
          { X = PyArray_ContiguousFromObject(X, NPY_DOUBLE, 1, 1);
            Py_DECREF((PyObject*)aX);
            if(!X)
            { Py_DECREF(Y);
              return NULL;
            }
            aX = (PyArrayObject*)X;
          }
          if (!PyArray_ISCONTIGUOUS(aY))
          { Y = PyArray_ContiguousFromObject(Y, NPY_DOUBLE, 1, 1);
            Py_DECREF((PyObject*)aY);
            if(!Y)
            { Py_DECREF(X);
              return NULL;
            }
            aY = (PyArrayObject*)Y;
          }
        }
        if (PyArray_DIM(aY, 0) != n)
        { PyErr_SetString(PyExc_TypeError,
                          "covariance: arguments x and y have different lengths.");
          Py_DECREF(X);
          Py_DECREF(Y);
          return NULL;
        }
        if (PyArray_DATA(aX) < PyArray_DATA(aY))
        { data = (double*) PyArray_DATA(aX);
          strides[1] = (int) ((PyArray_BYTES(aY) - PyArray_BYTES(aX)) / sizeof(double));
        }
        else
        { data = (double*) PyArray_DATA(aY);
          strides[1] = (int) ((PyArray_BYTES(aX) - PyArray_BYTES(aY)) / sizeof(double));
        }
        strides[0] = 1;
        ok = covariance((int)n,
                        (int) m,
                        data,
                        strides,
                        mode,
                        matrix);
        Py_DECREF(X);
        Py_DECREF(Y);
        if(!ok)
        { strcpy(buffer, "covariance: Memory error.");
          PyErr_SetString (PyExc_MemoryError, buffer);
          return NULL;
        }
        return PyFloat_FromDouble(matrix[1]);
      }
    default:
      sprintf(buffer,
              "covariance: Argument x has incorrect rank (%d expected 1 or 2).",
              PyArray_NDIM(aX));
      PyErr_SetString(PyExc_ValueError, buffer);
      return NULL;
  }
  return NULL; /* Never get here */
} 
/* end of wrapper for covariance */

static PyObject* _py_pearson(PyObject* X, PyObject* Y)
{ npy_intp n, m;
  PyArrayObject* aX = NULL;
  PyArrayObject* aY = NULL;
  int ok;

  if(!PyArray_Check(X))
  { X = PyArray_ContiguousFromObject(X, NPY_DOUBLE, 1, 2);
    if (!X) return NULL;
    aX = (PyArrayObject*)X;
  }
  else
  { aX = (PyArrayObject*)X;
    if (PyArray_TYPE(aX) != NPY_DOUBLE)
    { X = PyArray_Cast(aX, NPY_DOUBLE);
      if(!X) return NULL;
      aX = (PyArrayObject*)X;
    }
    else Py_INCREF(X);
  }
  switch (PyArray_NDIM(aX))
  { case 2:
      if(!Y)
      { PyArrayObject* matrix = NULL;
        npy_intp dims[2];
        int strides[2];
        strides[0] = (int) (PyArray_STRIDE(aX, 0) / sizeof(double));
        strides[1] = (int) (PyArray_STRIDE(aX, 1) / sizeof(double));
        n = PyArray_DIM(aX, 0);
        m = PyArray_DIM(aX, 1);
        dims[0] = m;
        dims[1] = m;
        if (n > INT_MAX || m > INT_MAX)
        { sprintf(buffer,
                  "correlation: Received more data than this function can handle (received a %" NPY_INTP_FMT " x %" NPY_INTP_FMT "matrix)",
                  n, m);
          PyErr_SetString (PyExc_ValueError, buffer);
          Py_DECREF(X);
          return NULL;
        }
        matrix = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        if(!matrix)
        { Py_DECREF(X);
          return NULL;
        }
        ok = pearson((int)n,
                     (int)m,
                     (double*)PyArray_DATA(aX),
                     strides,
                     (double*)PyArray_DATA(matrix));
        Py_DECREF(X);
        if(!ok)
        { strcpy(buffer, "correlation: failed to calculate the correlation matrix.");
          PyErr_SetString(PyExc_MemoryError, buffer);
          return NULL;
        }
        return PyArray_Return(matrix);
      }
      else
      { PyErr_SetString(PyExc_ValueError,
                        "correlation: Argument y passed, but argument x is two-dimensional.");
        Py_DECREF(X);
        return NULL;
      }
    case 1:
      if(!Y)
        return PyFloat_FromDouble(1.0);
      else
      { double* data;
        double matrix[4];
        int strides[2];
        n = PyArray_DIM(aX, 0);
        m = 2;
        if(!PyArray_Check(Y))
        { Y = PyArray_ContiguousFromObject(Y, NPY_DOUBLE, 1, 1);
          if (!Y)
          { Py_DECREF(X);
            return NULL;
          }
          aY = (PyArrayObject*)Y;
        }
        else
        { aY = (PyArrayObject*)Y;
          if (PyArray_NDIM(aY) != 1)
          { sprintf(buffer,
                    "correlation: Argument y has incorrect rank (%d expected 1).",
                    PyArray_NDIM(aY));
            PyErr_SetString(PyExc_ValueError, buffer);
            Py_DECREF(X);
            return NULL;
          }
          if (PyArray_TYPE(aY) != NPY_DOUBLE)
          { Y = PyArray_Cast(aY, NPY_DOUBLE);
            if (!Y)
            { Py_DECREF(X);
              return NULL;
            }
            aY = (PyArrayObject*)Y;
          } 
          else Py_INCREF(Y);
          /* X and Y should have equal strides. Make both strides equal to one. */
          if (!PyArray_ISCONTIGUOUS(aX))
          { X = PyArray_ContiguousFromObject(X, NPY_DOUBLE, 1, 1);
            Py_DECREF((PyObject*)aX);
            if(!X)
            { Py_DECREF(Y);
              return NULL;
            }
            aX = (PyArrayObject*)X;
          }
          if (!PyArray_ISCONTIGUOUS(aY))
          { Y = PyArray_ContiguousFromObject(Y, NPY_DOUBLE, 1, 1);
            Py_DECREF((PyObject*)aY);
            if(!Y)
            { Py_DECREF(X);
              return NULL;
            }
            aY = (PyArrayObject*)Y;
          }
        }
        if (PyArray_DIM(aY, 0) != n)
        { PyErr_SetString(PyExc_ValueError,
                          "correlation: arguments x and y have different lengths.");
          Py_DECREF(X);
          Py_DECREF(Y);
          return NULL;
        }
        if (PyArray_DATA(aX) < PyArray_DATA(aY))
        { data = (double*) PyArray_DATA(aX);
          strides[1] = (int) ((PyArray_BYTES(aY) - PyArray_BYTES(aX)) / sizeof(double));
        }
        else
        { data = (double*) PyArray_DATA(aY);
          strides[1] = (int) ((PyArray_BYTES(aX) - PyArray_BYTES(aY)) / sizeof(double));
        }
        strides[0] = 1;
        ok = pearson((int)n,
                     (int)m,
                     data,
                     strides,
                     matrix);
        Py_DECREF(X);
        Py_DECREF(Y);
        if(!ok)
        { PyErr_SetString(PyExc_MemoryError,
                          "correlation: Memory error.");
          return NULL;
        }
        return PyFloat_FromDouble(matrix[1]);
      }
    default:
      sprintf(buffer,
              "correlation: Argument x has incorrect rank (%d expected 1 or 2).",
              PyArray_NDIM(aX));
      PyErr_SetString(PyExc_ValueError, buffer);
      return NULL;
  }
  return NULL;
}

static PyObject* _py_spearman(PyObject* X, PyObject* Y)
{ npy_intp n, m;
  PyArrayObject* aX = NULL;
  PyArrayObject* aY = NULL;
  int ok;

  if(!PyArray_Check(X))
  { X = PyArray_ContiguousFromObject(X, NPY_DOUBLE, 1, 2);
    if (!X) return NULL;
    aX = (PyArrayObject*)X;
  }
  else
  { aX = (PyArrayObject*)X;
    if (PyArray_TYPE(aX) != NPY_DOUBLE)
    { X = PyArray_Cast(aX, NPY_DOUBLE);
      if(!X) return NULL;
      aX = (PyArrayObject*)X;
    }
    else Py_INCREF(X);
  }

  switch (PyArray_NDIM(aX))
  { case 2:
    {
      int i, j;
      PyArrayObject* matrix = NULL;
      npy_intp dims[2];
      double* data;
      const int rowstride = (const int) (PyArray_STRIDE(aX, 0)/sizeof(double));
      const int colstride = (const int) (PyArray_STRIDE(aX, 1)/sizeof(double));
      const double* p;
      double* q;
      if(Y)
      { PyErr_SetString(PyExc_ValueError,
                        "correlation: Argument y passed, but argument x is two-dimensional.");
        Py_DECREF(X);
        return NULL;
      }
      n = PyArray_DIM(aX, 0);
      m = PyArray_DIM(aX, 1);
      if (n > INT_MAX || m > INT_MAX)
      { sprintf(buffer,
                "correlation: Received more data than this function can handle (received a %" NPY_INTP_FMT " x %" NPY_INTP_FMT "matrix)",
                n, m);
        PyErr_SetString (PyExc_ValueError, buffer);
        Py_DECREF(X);
        return NULL;
      }
      dims[0] = m;
      dims[1] = m;
      matrix = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
      if(!matrix)
      { Py_DECREF(X);
        strcpy(buffer, "correlation: failed to calculate the correlation matrix.");
        PyErr_SetString(PyExc_MemoryError, buffer);
        return NULL;
      }
      data = malloc(n*m*sizeof(double));
      if(!data)
      { Py_DECREF(X);
        Py_DECREF(matrix);
        strcpy(buffer, "correlation: failed to calculate the correlation matrix.");
        PyErr_SetString(PyExc_MemoryError, buffer);
        return NULL;
      }
      q = data;
      for (j = 0; j < m; j++)
      { p = ((double*)PyArray_DATA(aX)) + j*colstride;
        for (i = 0; i < n; i++, q++, p+=rowstride) *q = *p;
      }
      /* -------------------------------------------------------------------- */
      ok = spearman((int)m,
                    (int)n,
                    data,
                    (double*)PyArray_DATA(matrix));
      /* -------------------------------------------------------------------- */
      Py_DECREF(X);
      free(data);
      if(!ok)
      { strcpy(buffer, "correlation: failed to calculate the correlation matrix.");
        PyErr_SetString(PyExc_MemoryError, buffer);
        return NULL;
      }
      return PyArray_Return(matrix);
    }
    case 1:
    {
      double* data;
      double matrix[4];
      int stridex;
      int stridey;
      const double* p;
      double* q;
      int j;

      if(!Y) return PyFloat_FromDouble(1.0);
      n = PyArray_DIM(aX, 0);
      m = 2;
      if(!PyArray_Check(Y))
      { Y = PyArray_ContiguousFromObject(Y, NPY_DOUBLE, 1, 1);
        if (!Y)
        { Py_DECREF(X);
          return NULL;
        }
        aY = (PyArrayObject*)Y;
      }
      else
      { aY = (PyArrayObject*)Y;
        if (PyArray_NDIM(aY) != 1)
        { sprintf(buffer,
                  "correlation: Argument y has incorrect rank (%d expected 1).",
                  PyArray_NDIM(aY));
          PyErr_SetString(PyExc_ValueError, buffer);
          Py_DECREF(X);
          return NULL;
        }
        if (PyArray_TYPE(aY) != NPY_DOUBLE)
        { Y = PyArray_Cast(aY, NPY_DOUBLE);
          if (!Y)
          { Py_DECREF(X);
            return NULL;
          }
          aY = (PyArrayObject*)Y;
        } 
        else Py_INCREF(Y);
      }
      if (PyArray_DIM(aY, 0) != n)
      { PyErr_SetString(PyExc_ValueError,
                        "correlation: arguments x and y have different lengths.");
        Py_DECREF(X);
        Py_DECREF(Y);
        return NULL;
      }
      data = malloc(n*m*sizeof(double));
      if(data==NULL)
      { Py_DECREF(X);
        Py_DECREF(Y);
        PyErr_SetString(PyExc_MemoryError,
                        "correlation: failed to calculate the correlation.");
        return NULL;
      }
      stridex = (int) (PyArray_STRIDE(aX, 0) / sizeof(double));
      stridey = (int) (PyArray_STRIDE(aY, 0) / sizeof(double));
      q = data;
      p = (double*)PyArray_DATA(aX);
      for (j = 0; j < n; j++, p+=stridex, q++) *q = *p;
      p = (double*)PyArray_DATA(aY);
      for (j = 0; j < n; j++, p+=stridey, q++) *q = *p;

      /* -------------------------------------------------------------------- */
      ok = spearman(m, n, data, matrix);
      /* -------------------------------------------------------------------- */
      Py_DECREF(X);
      Py_DECREF(Y);
      free(data);
      if(!ok)
      { PyErr_SetString(PyExc_MemoryError, "correlation: Memory error.");
        return NULL;
      }
      return PyFloat_FromDouble(matrix[1]);
    }
    default:
      sprintf(buffer,
              "correlation: Argument x has incorrect rank (%d expected 1 or 2).",
              PyArray_NDIM(aX));
      PyErr_SetString(PyExc_ValueError, buffer);
      return NULL;
  }
  return NULL;
}

static PyObject* _py_intraclass(PyObject* X, PyObject* Y)
{ npy_intp n, m;
  PyArrayObject* aX = NULL;
  PyArrayObject* aY = NULL;
  int ok;

  if(!PyArray_Check(X))
  { X = PyArray_ContiguousFromObject(X, NPY_DOUBLE, 1, 2);
    if (!X) return NULL;
    aX = (PyArrayObject*)X;
  }
  else
  { aX = (PyArrayObject*)X;
    if (PyArray_TYPE(aX) != NPY_DOUBLE)
    { X = PyArray_Cast(aX, NPY_DOUBLE);
      if(!X) return NULL;
      aX = (PyArrayObject*)X;
    }
    else Py_INCREF(X);
  }
  switch (PyArray_NDIM(aX))
  { case 2:
      if(!Y)
      { PyArrayObject* matrix = NULL;
        npy_intp dims[2];
        int strides[2];
        strides[0] = (int) (PyArray_STRIDE(aX, 0) / sizeof(double));
        strides[1] = (int) (PyArray_STRIDE(aX, 1) / sizeof(double));
        n = PyArray_DIM(aX, 0);
        m = PyArray_DIM(aX, 1);
        if (n > INT_MAX || m > INT_MAX)
        { sprintf(buffer,
                  "correlation: Received more data than this function can handle (received a %" NPY_INTP_FMT " x %" NPY_INTP_FMT "matrix)",
                  n, m);
          PyErr_SetString (PyExc_ValueError, buffer);
          Py_DECREF(X);
          return NULL;
        }
        dims[0] = m;
        dims[1] = m;
        matrix = (PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        if(!matrix)
        { Py_DECREF(X);
          return NULL;
        }
        ok = intraclass((int)n,
                        (int)m,
                        (double*)PyArray_DATA(aX),
                        strides,
                        (double*)PyArray_DATA(matrix));
        Py_DECREF(X);
        if(!ok)
        { strcpy(buffer, "correlation: failed to calculate the correlation matrix.");
          PyErr_SetString(PyExc_MemoryError, buffer);
          return NULL;
        }
        return PyArray_Return(matrix);
      }
      else
      { PyErr_SetString(PyExc_ValueError,
                        "correlation: Argument y passed, but argument x is two-dimensional.");
        Py_DECREF(X);
        return NULL;
      }
    case 1:
      if(!Y)
        return PyFloat_FromDouble(1.0);
      else
      { double* data;
        double matrix[4];
        int strides[2];
        n = PyArray_DIM(aX, 0);
        m = 2;
        if(!PyArray_Check(Y))
        { Y = PyArray_ContiguousFromObject(Y, NPY_DOUBLE, 1, 1);
          if (!Y)
          { Py_DECREF(X);
            return NULL;
          }
          aY = (PyArrayObject*)Y;
        }
        else
        { aY = (PyArrayObject*)Y;
          if (PyArray_NDIM(aY) != 1)
          { sprintf(buffer,
                    "correlation: Argument y has incorrect rank (%d expected 1).",
                    PyArray_NDIM(aY));
            PyErr_SetString(PyExc_ValueError, buffer);
            Py_DECREF(X);
            return NULL;
          }
          if (PyArray_TYPE(aY) != NPY_DOUBLE)
          { Y = PyArray_Cast(aY, NPY_DOUBLE);
            if (!Y)
            { Py_DECREF(X);
              return NULL;
            }
            aY = (PyArrayObject*)Y;
          } 
          else Py_INCREF(Y);
          /* X and Y should have equal strides. Make both strides equal to one. */
          if (!PyArray_ISCONTIGUOUS(aX))
          { X = PyArray_ContiguousFromObject(X, NPY_DOUBLE, 1, 1);
            Py_DECREF((PyObject*)aX);
            if(!X)
            { Py_DECREF(Y);
              return NULL;
            }
            aX = (PyArrayObject*)X;
          }
          if (!PyArray_ISCONTIGUOUS(aY))
          { Y = PyArray_ContiguousFromObject(Y, NPY_DOUBLE, 1, 1);
            Py_DECREF((PyObject*)aY);
            if(!Y)
            { Py_DECREF(X);
              return NULL;
            }
            aY = (PyArrayObject*)Y;
          }
        }
        if (PyArray_DIM(aY, 0) != n)
        { PyErr_SetString(PyExc_ValueError,
                          "correlation: arguments x and y have different lengths.");
          Py_DECREF(X);
          Py_DECREF(Y);
          return NULL;
        }
        if (PyArray_DATA(aX) < PyArray_DATA(aY))
        { data = (double*) PyArray_DATA(aX);
          strides[1] = (int) ((PyArray_BYTES(aY) - PyArray_BYTES(aX)) / sizeof(double));
        }
        else
        { data = (double*) PyArray_DATA(aY);
          strides[1] = (int) ((PyArray_BYTES(aX) - PyArray_BYTES(aY)) / sizeof(double));
        }
        strides[0] = 1;
        ok = intraclass((int)n,
                        (int)m,
                        data,
                        strides,
                        matrix);
        Py_DECREF(X);
        Py_DECREF(Y);
        if(!ok)
        { PyErr_SetString(PyExc_MemoryError,
                          "correlation: Memory error.");
          return NULL;
        }
        return PyFloat_FromDouble(matrix[1]);
      }
    default:
      sprintf(buffer,
              "correlation: Argument x has incorrect rank (%d expected 1 or 2).",
              PyArray_NDIM(aX));
      PyErr_SetString(PyExc_ValueError, buffer);
      return NULL;
  }
  return NULL;
}

/* correlation */
static char correlation__doc__[] =
"correlation(x, y = None, method = 'Pearson')\n"
"\n"
"This function estimates the correlation between x and y (if both are one-\n"
"dimensional have have the same length), or between the columns of x (if x\n"
"is a two-dimensional array of data).\n"
"\n"
"Arguments:\n"
"o) The array x should be either one-dimensional or two-dimensional. If x\n"
"   is one-dimensional and y==None, then this function returns the variance\n"
"   in x. If both x and y are one-dimensional and have the same length,\n"
"   this function returns the correlation between x and y. If x is two-\n"
"   dimensional, then this function calculates the correlation matrix\n"
"   between the columns in x. In this case, y is ignored.\n"
"o) If passed, the array y should be one-dimensional and contain the same\n"
"   number of elements as x; the function will return the correlation\n"
"   between x and y in this case.\n"
"o) For method=='Pearson', which is the default, this function calculates\n"
"   the Pearson correlation. If mode=='Spearman', the Spearman rank\n"
"   correlation is calculated instead. For mode=='Intraclass', the function\n"
"   calculates the intraclass correlation.\n"
"\n"
"Return values:\n"
"o) If x is a one-dimensional array and y==None: the function returns 1.0\n"
"o) If both x and y are one-dimensional arrays of the same length: the\n"
"   correlation between the variables x and y.\n"
"o) If x is a two-dimensional array: the correlation matrix between the\n"
"   columns in x. Element [i,j] in the correlation matrix contains the\n"
"   correlation between the columns x[:,i] and x[:,j].\n";


static PyObject*
py_correlation(PyObject* self, PyObject* args, PyObject* keywords)
{ PyObject* X = NULL;
  PyObject* Y = NULL;
  char* METHOD = NULL;
  char method = 'p';

  /* -- Read the input variables ----------------------------------------- */
  static char* kwlist[] = {"x", "y", "method", NULL};
  if(!PyArg_ParseTupleAndKeywords(args, keywords, "O|Os", kwlist,
                                  &X, &Y, &METHOD)) return NULL;

  /* -- Check the input variables ---------------------------------------- */
  if (METHOD)
  { const char* names[] = {"pearson", "spearman", "intraclass", NULL};
    const char mnemonics[] = {'p', 's', 'i', '\0'};
    method = find_mnemonic(METHOD, names, mnemonics);
    if(!method)
    { PyErr_SetString(PyExc_ValueError,
                      "correlation: Unknown method specified (should be 'pearson' or 'spearman')");
      return NULL;
    }
  }

  switch (method)
  { case 'p':
      return _py_pearson(X, Y);
    case 's':
      return _py_spearman(X, Y);
    case 'i':
      return _py_intraclass(X, Y);
  }

  /* --------------------------------------------------------------------- */
  return NULL; /* Never get here */
} 
/* end of wrapper for correlation */

/* regression */
static char regression__doc__[] =
"regression(x, y)\n"
"This function returns the intercept and the slope of the linear regression\n"
"line through the data in x and y.\n";

static PyObject*
py_regression (PyObject* unused, PyObject* args)
{ double a, b;
  int result;
  npy_intp n;
  PyObject* X = NULL;
  PyArrayObject* aX = NULL;
  PyObject* Y = NULL;
  PyArrayObject* aY = NULL;

  /* -- Read the input variables ----------------------------------------- */
  if(!PyArg_ParseTuple(args, "OO", &X, &Y)) return NULL;

  /* -- Check the input variables ---------------------------------------- */
  if(!PyArray_Check(X))
  { X = PyArray_ContiguousFromObject(X, NPY_DOUBLE, 1, 1);
    if (!X) return NULL;
    aX = (PyArrayObject*)X;
  }
  else
  { aX = (PyArrayObject*)X;
    if ((PyArray_NDIM(aX) != 1) &&
        (PyArray_NDIM(aX) > 0 || PyArray_DIM(aX, 0) != 1))
    { sprintf(buffer,
              "regression: Argument x has incorrect rank (%d expected 1).",
              PyArray_NDIM(aX));
      PyErr_SetString(PyExc_ValueError, buffer);
      return NULL;
    }
    if (PyArray_TYPE(aX) != NPY_DOUBLE)
    { X = PyArray_Cast(aX, NPY_DOUBLE);
      if (!X) return NULL;
      aX = (PyArrayObject*)X;
    } 
    else Py_INCREF(X);
    if (!PyArray_ISCONTIGUOUS(aX))
    { X = PyArray_ContiguousFromObject(X, PyArray_TYPE(aX), 1, 1);
      Py_DECREF((PyObject*)aX);
      if(!X) return NULL;
      aX = (PyArrayObject*)X;
    }
  }
  /* --------------------------------------------------------------------- */
  n = PyArray_DIM(aX, 0);
  /* --------------------------------------------------------------------- */
  if(!PyArray_Check(Y))
  { Y = PyArray_ContiguousFromObject(Y, NPY_DOUBLE, 1, 1);
    if (!Y)
    { Py_DECREF(X);
      return NULL;
    }
    aY = (PyArrayObject*)Y;
  }
  else
  { aY = (PyArrayObject*)Y;
    if ((PyArray_NDIM(aY) != 1) &&
        (PyArray_NDIM(aY) > 0 || PyArray_DIM(aY, 0) != 1))
    { sprintf(buffer,
              "regression: Argument y has incorrect rank (%d expected 1).",
              PyArray_NDIM(aY));
      PyErr_SetString(PyExc_ValueError, buffer);
      Py_DECREF(X);
      return NULL;
    }
    if (PyArray_TYPE(aY) != NPY_DOUBLE)
    { Y = PyArray_Cast(aY, NPY_DOUBLE);
      if (!Y)
      { Py_DECREF(X);
        return NULL;
      }
      aY = (PyArrayObject*)Y;
    } 
    else Py_INCREF(Y);
    if (!PyArray_ISCONTIGUOUS(aY))
    { Y = PyArray_ContiguousFromObject(Y, PyArray_TYPE(aY), 1, 1);
      Py_DECREF((PyObject*)aY);
      if(!Y)
      { Py_DECREF(X);
        return NULL;
      }
      aY = (PyArrayObject*)Y;
    }
  }
  /* --------------------------------------------------------------------- */
  if (PyArray_DIM(aY, 0) != n)
  { PyErr_SetString(PyExc_TypeError,
                    "regression: arguments x and y should have an equal size.");
    Py_DECREF(X);
    Py_DECREF(Y);
    return NULL;
  }
  if (n > INT_MAX)
  { sprintf(buffer,
            "received more data than this function can handle (%" NPY_INTP_FMT " data points received)",
            n);
    PyErr_SetString (PyExc_ValueError, buffer);
    Py_DECREF(X);
    Py_DECREF(Y);
    return NULL;
  }
  /* --------------------------------------------------------------------- */
  result = regression((int)n,
                      (double*) PyArray_DATA(aX),
                      (double*) PyArray_DATA(aY),
                      &a, &b);
  /* --------------------------------------------------------------------- */
  Py_DECREF(X);
  Py_DECREF(Y);
  /* --------------------------------------------------------------------- */
  if (!result) /* Error occured due to zero denominator */
  { PyErr_SetString(PyExc_RuntimeError,
                    "regression: failed to calculate the slope of the regression line due to a zero denominator.");
    return NULL;
  }
  /* --------------------------------------------------------------------- */
  return Py_BuildValue("(ff)", a, b);
} 
/* end of wrapper for regression */

/* ========================================================================== */
/* -- The methods table ----------------------------------------------------- */
/* ========================================================================== */


static struct PyMethodDef methods[] = {
   {"mean", (PyCFunction) py_mean, METH_VARARGS, mean__doc__},
   {"median", (PyCFunction) py_median, METH_VARARGS, median__doc__},
   {"variance", (PyCFunction) py_variance, METH_VARARGS | METH_KEYWORDS, variance__doc__},
   {"covariance", (PyCFunction) py_covariance, METH_VARARGS | METH_KEYWORDS, covariance__doc__},
   {"correlation", (PyCFunction) py_correlation, METH_VARARGS | METH_KEYWORDS, correlation__doc__},
   {"regression", (PyCFunction) py_regression, METH_VARARGS, regression__doc__},
   {"bandwidth", (PyCFunction) py_bandwidth, METH_VARARGS | METH_KEYWORDS, bandwidth__doc__},
   {"pdf", (PyCFunction) py_pdf, METH_VARARGS | METH_KEYWORDS, pdf__doc__},
   {"cpdf", (PyCFunction) py_cpdf, METH_VARARGS | METH_KEYWORDS, cpdf__doc__},
   {"cpdfc", (PyCFunction) py_cpdfc, METH_VARARGS | METH_KEYWORDS, cpdfc__doc__},
   /* {"localfit", (PyCFunction) py_localfit, METH_VARARGS | METH_KEYWORDS, localfit__doc__}, */
   {NULL,          NULL, 0, NULL} /* sentinel */
};

/* ========================================================================== */
/* -- Initialization -------------------------------------------------------- */
/* ========================================================================== */

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef module = {
   PyModuleDef_HEAD_INIT,
   "statistics",   /* name of module */
   "Statistics for Python", /* module documentation, may be NULL */
   -1,       /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
   methods
};



PyMODINIT_FUNC
PyInit_statistics(void)

#else

void initstatistics(void)
#endif
{
  import_array ();
#if PY_MAJOR_VERSION < 3
  Py_InitModule3("statistics", methods, "Statistics for Python");
  if (PyErr_Occurred()) Py_FatalError("can't initialize module statistics");
#else
  return PyModule_Create(&module);
#endif
}
