cdef extern from "<armadillo>" namespace "arma":
    cdef cppclass Mat[T]:
        unsigned int n_rows
        unsigned int n_cols
        unsigned int n_elem
        Mat()
        Mat(T*, unsigned int, unsigned int, bool, bool)
        Mat(T*, unsigned int, unsigned int)
        T *memptr()

import numpy as np

cimport cython
cimport numpy as np
from libcpp cimport bool
from cython.operator cimport dereference as deref

# datatypes:
# double, cython.double, np.double_t, np.float64_t
# float, cython.float, np.float32_t
# dtype:
# np.float, np.double, np.float64 # python float is 64 bit

cdef Mat[double] dToArmaMat(np.ndarray[double, ndim=2] X):
    if not X.flags.f_contiguous:
        X = X.copy(order="F")
    return Mat[double](<double*> X.data, X.shape[0], X.shape[1], 0, 1)

cdef Mat[float] fToArmaMat(np.ndarray[float, ndim=2] X):
    if not X.flags.f_contiguous:
        X = X.copy(order="F")
    return Mat[float](<float*> X.data, X.shape[0], X.shape[1], 0, 1)

@cython.boundscheck(False)
cdef np.ndarray[double, ndim=2] dToNdarray(Mat[double] &m):
    cdef np.ndarray[double, ndim=2] arr
    cdef double *pArr
    cdef double *pM
    arr = np.ndarray((m.n_rows, m.n_cols), dtype=np.float64, order='F')
    pArr = <double *>arr.data
    pM = m.memptr()
    for i in range(m.n_rows*m.n_cols):
        pArr[i] = pM[i]
    return arr

@cython.boundscheck(False)
cdef np.ndarray[float, ndim=2] fToNdarray(Mat[float] &m):
    cdef np.ndarray[float, ndim=2] arr
    cdef float *pArr
    cdef float *pM
    arr = np.ndarray((m.n_rows, m.n_cols), dtype=np.float32, order='F')
    pArr = <float *>arr.data
    pM = m.memptr()
    for i in range(m.n_rows*m.n_cols):
        pArr[i] = pM[i]
    return arr

      
cdef extern from "vlfeat.hpp":
    cdef void c_vl_sift "vl_sift" (Mat[float] &, Mat[double] &, Mat[float] &, \
                      int, int, int, double, double, double, double, double, \
                      bool, bool, int)

def vl_sift(np.ndarray[np.float32_t, ndim=2] data, 
            frames=None,
            int octaves = -1,
            int levels = -1,
            int first_octave = -1,
            double peak_thresh = -1,
            double edge_thresh = -1,
            double norm_thresh = -1,
            double magnif = -1,
            double window_size = -1,
            bool orientations = False,
            bool float_descriptors = False,
            int verbose = 0):

    cdef np.ndarray[double, ndim=2] f
    cdef np.ndarray[float, ndim=2] d
    cdef Mat[float] _I
    cdef Mat[double] _f
    cdef Mat[float] _d
    
    _I = fToArmaMat(data)
    if frames == None:
        _f = Mat[double]()
    else:
        _f = dToArmaMat(frames)
    _d = Mat[float]()

    c_vl_sift(<const Mat[float] &>_I, _f, _d, 
              octaves, levels, first_octave, peak_thresh,
              edge_thresh, norm_thresh, magnif, window_size, orientations,
              float_descriptors, verbose)
    if frames == None:
        f = dToNdarray(_f)
    else:
        f = frames
    d = fToNdarray(_d)
    
    return f, d
