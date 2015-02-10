import numpy as np
cimport numpy as np
cimport cython
from libcpp cimport bool

# datatypes:
# 64 bit: double, cython.double, np.double_t, np.float64_t
# 32 bit: float, cython.float, np.float32_t
# dtype:
# 64 bit: np.float, np.double, np.float64 # python float is 64 bit
# 32 bit: np.float32

cdef extern from "<armadillo>" namespace "arma":
    cdef cppclass Mat[T]:
        unsigned int n_rows
        unsigned int n_cols
        unsigned int n_elem
        Mat()
        Mat(T*, unsigned int, unsigned int, bool, bool)
        Mat(T*, unsigned int, unsigned int)
        T *memptr()

cdef inline Mat[double] pyarma_to_double(np.ndarray[double, ndim=2] X):
    if not X.flags.f_contiguous:
        X = X.copy(order="F")
    return Mat[double](<double*> X.data, X.shape[0], X.shape[1], 0, 1)

cdef inline Mat[float] pyarma_to_float(np.ndarray[float, ndim=2] X):
    if not X.flags.f_contiguous:
        X = X.copy(order="F")
    return Mat[float](<float*> X.data, X.shape[0], X.shape[1], 0, 1)

@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline np.ndarray[double, ndim=2] pyarma_from_double(Mat[double] &m):
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
@cython.wraparound(False)
cdef inline np.ndarray[float, ndim=2] pyarma_from_float(Mat[float] &m):
    cdef np.ndarray[float, ndim=2] arr
    cdef float *pArr
    cdef float *pM
    arr = np.ndarray((m.n_rows, m.n_cols), dtype=np.float32, order='F')
    pArr = <float *>arr.data
    pM = m.memptr()
    for i in range(m.n_rows*m.n_cols):
        pArr[i] = pM[i]
    return arr
