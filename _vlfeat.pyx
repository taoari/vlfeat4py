import numpy as np
cimport numpy as np
cimport cython
from libcpp cimport bool

from libcpp.string cimport string
from cython.operator cimport dereference as deref

from arma cimport Mat, pyarma_from_double, pyarma_to_double
from arma cimport pyarma_from_float, pyarma_to_float

cdef extern from "vlfeat.hpp":
    cdef void c_vl_sift "vl_sift" (Mat[float] &, Mat[double] &, Mat[float] &, \
                      int, int, int, double, double, double, double, double, \
                      bool, bool, int)
    cdef void c_vl_dsift "vl_dsift" (Mat[float] &, Mat[double] &, Mat[float] &, \
                                     Mat[double] &, Mat[double] &, Mat[double] &, \
                                     Mat[double] &, bool, bool, double, bool, int)
    cdef void c_vl_imsmooth_f "vl_imsmooth" (Mat[float] &, Mat[float] &, \
            double, string, string, int, int)
    cdef void c_vl_imsmooth_d "vl_imsmooth" (Mat[double] &, Mat[double] &, \
            double, string, string, int, int)
    cdef void c_vl_homkermap "vl_homkermap" (Mat[float] &, Mat[float] &, int, \
            string, string, double, double) 
    cdef void c_vl_kmeans "vl_kmeans" (Mat[float] &, Mat[float] &, int, \
            string, string, string, double, int, int, int, int, int)
    cdef void c_vl_gmm "vl_gmm" (Mat[float] &, int, \
            Mat[float] &, Mat[float] &, Mat[float] &, Mat[float] &, \
            string, Mat[float] &, Mat[float] &, Mat[float] &, \
            int, int, Mat[double] &, int)
    cdef void c_vl_fisher "vl_fisher" (Mat[float] &, \
            Mat[float] &, Mat[float] &, Mat[float] &, \
            Mat[float] &, bool, bool, bool, bool, int)

def vl_sift(np.ndarray[np.float32_t, ndim=2] data, 
            frames=None,
            int octaves = -1,
            int levels = -1,
            int firstOctave = -1,
            double peakThresh = -1,
            double edgeThresh = -1,
            double normThresh = -1,
            double magnif = -1,
            double windowSize = -1,
            bool orientations = False,
            bool floatDescriptors = False,
            int verbose = 0):

    cdef np.ndarray[double, ndim=2] f
    cdef np.ndarray[float, ndim=2] d
    cdef Mat[float] _I
    cdef Mat[double] _f
    cdef Mat[float] _d
    
    _I = pyarma_to_float(data)
    if frames == None:
        _f = Mat[double]()
    else:
        _f = pyarma_to_double(frames)
    _d = Mat[float]()

    c_vl_sift(<const Mat[float] &>_I, _f, _d, 
              octaves, levels, firstOctave, peakThresh,
              edgeThresh, normThresh, magnif, windowSize, orientations,
              floatDescriptors, verbose)
    if frames == None:
        f = pyarma_from_double(_f)
    else:
        f = frames
    d = pyarma_from_float(_d)
    
    return f, d
    
def vl_dsift(np.ndarray[float, ndim=2] data,
             bounds = None,
             step = None,
             size = None,
             geometry = None,
             bool fast = True,
             bool norm = False,
             double windowSize = -1.0,
             bool floatDescriptors = False,
             int verbose = 0):

    cdef np.ndarray[double, ndim=2] f
    cdef np.ndarray[float, ndim=2] d
    
    cdef Mat[float] _I
    cdef Mat[double] _f
    cdef Mat[float] _d
    cdef Mat[double] _bounds
    cdef Mat[double] _step
    cdef Mat[double] _size
    cdef Mat[double] _geometry
    
    _I = pyarma_to_float(data)
    _f = Mat[double]()
    _d = Mat[float]()
    
    if bounds == None:
        _bounds = Mat[double]()
    else:
        _bounds = pyarma_to_double(np.float64(bounds, order='F').reshape((-1,1)))
        
    if step == None:
        _step = Mat[double]()
    else:
        _step = pyarma_to_double(np.float64(step, order='F').reshape((-1,1)))

    if size == None:
        _size = Mat[double]()
    else:
        _size = pyarma_to_double(np.float64(size, order='F').reshape((-1,1)))
        
    if geometry == None:
        _geometry = Mat[double]()
    else:
        _geometry = pyarma_to_double(np.float64(geometry, order='F').reshape((-1,1)))

    c_vl_dsift(<const Mat[float] &>_I, _f, _d,
               _bounds, _step, _size, _geometry, 
               fast, norm, windowSize, floatDescriptors, verbose)
               
    f = pyarma_from_double(_f)
    d = pyarma_from_float(_d)
    
    return f, d
    
def vl_imsmooth_f(np.ndarray[float, ndim=2] I,
        double sigma,
        string padding = "continuity",
        string kernel = "gaussian",
        int subsample = 1,
        int verbose = 0):
    
    cdef np.ndarray[float, ndim=2] Is
    
    cdef Mat[float] _I
    cdef Mat[float] _Is
    
    _I = pyarma_to_float(I)
    _Is = Mat[float]()
    
    c_vl_imsmooth_f(<const Mat[float] &>_I, _Is, 
            sigma, padding, kernel, subsample, verbose)
    
    Is = pyarma_from_float(_Is)
    
    return Is

def vl_imsmooth_d(np.ndarray[double, ndim=2] I,
        double sigma,
        string padding = "continuity",
        string kernel = "gaussian",
        int subsample = 1,
        int verbose = 0):
    
    cdef np.ndarray[double, ndim=2] Is
    
    cdef Mat[double] _I
    cdef Mat[double] _Is
    
    _I = pyarma_to_double(I)
    _Is = Mat[double]()
    
    c_vl_imsmooth_d(<const Mat[double] &>_I, _Is, 
            sigma, padding, kernel, subsample, verbose)
    
    Is = pyarma_from_double(_Is)
    
    return Is

def vl_imsmooth(I,
        double sigma,
        string padding = "continuity",
        string kernel = "gaussian",
        int subsample = 1,
        int verbose = 0):
    if I.dtype == np.float32:
        return vl_imsmooth_f(I, sigma, padding, kernel, subsample, verbose)
    elif I.dtype == np.float64:
        return vl_imsmooth_d(I, sigma, padding, kernel, subsample, verbose)
    else:
        return None

def vl_homkermap(np.ndarray[float, ndim=2] X, int n,
        string kernel = "kchi2",
        string window = "rectangular",
        double gamma = 1.0,
        double period = -1):

    cdef Mat[float] _X
    cdef Mat[float] _V
    
    _X = pyarma_to_float(X)
    _V = Mat[float]()
    
    c_vl_homkermap(<const Mat[float] &>_X, _V, n,
            kernel, window, gamma, period)
            
    cdef np.ndarray[float, ndim=2] V
    V = pyarma_from_float(_V)
    
    return V
    
def vl_kmeans(np.ndarray[float, ndim=2] X, int numCenters,
        string algorithm = "lloyd",
        string distance = "l2",
        string initialization = "plusplus",
        double minEnergyVariation = -1,
        int numRepetitions = 1,
        int numTrees = 3,
        int maxNumComparisons = 100,
        int maxNumIterations = 100,
        int verbose = 0):

    cdef Mat[float] _X
    cdef Mat[float] _Y
    
    _X = pyarma_to_float(X)
    _Y = Mat[float]()
    
    c_vl_kmeans(<const Mat[float] &>_X, _Y, numCenters,
            algorithm, distance, initialization, minEnergyVariation,
            numRepetitions, numTrees, maxNumComparisons, maxNumIterations,
            verbose)
            
    cdef np.ndarray[float, ndim=2] Y
    Y = pyarma_from_float(_Y)
    
    return Y
    
def vl_gmm(np.ndarray[float, ndim=2] X, int numClusters,
        string initialization = 'rand',
        initMeans = None,
        initCovariances = None,
        initPriors = None,
        int maxNumIterations = 100,
        int numRepetitions = 1,
        covarianceBound = None,
        int verbose = 0):
    cdef Mat[float] _X
    cdef Mat[float] _initMeans
    cdef Mat[float] _initCovariances
    cdef Mat[float] _initPriors
    cdef Mat[double] _covarianceBound
    
    _X = pyarma_to_float(X)
    if initMeans is not None:
        _initMeans = pyarma_to_float(initMeans)
    if initCovariances is not None:
        _initCovariances = pyarma_to_float(initCovariances)
    if initPriors is not None:
        _initPriors = pyarma_to_float(initPriors)
    if covarianceBound is not None:
        _covarianceBound = pyarma_to_double(np.float64(covarianceBound, order='F').reshape((-1,1)))
    
    cdef Mat[float] _means, _covariances, _priors, _posteriors
    c_vl_gmm(_X, numClusters, \
        _means, _covariances, _priors, _posteriors, \
        initialization, _initMeans, _initCovariances, _initPriors, \
        maxNumIterations, numRepetitions, _covarianceBound, verbose)
    cdef np.ndarray[float, ndim=2] means, covariances, priors, posteriors
    means = pyarma_from_float(_means)
    covariances = pyarma_from_float(_covariances)
    priors = pyarma_from_float(_priors)
    posteriors = pyarma_from_float(_posteriors)
    return means, covariances, priors, posteriors

def vl_fisher(np.ndarray[float, ndim=2] X,
        np.ndarray[float, ndim=2] means, 
        np.ndarray[float, ndim=2] covariances,
        np.ndarray[float, ndim=2] priors,
        bool normalized = False,
        bool squareRoot = False,
        bool improved = False,
        bool fast = False,
        int verbosity = 0):
    cdef Mat[float] _X, _means, _covariances, _priors
    _X = pyarma_to_float(X)
    _means = pyarma_to_float(means)
    _covariances = pyarma_to_float(covariances)
    _priors = pyarma_to_float(priors)
    
    cdef Mat[float] _enc
    c_vl_fisher(_X, _means, _covariances, _priors, _enc, 
        normalized, squareRoot, improved, fast, verbosity)
    
    cdef np.ndarray[float, ndim=2] enc
    enc = pyarma_from_float(_enc)
    return enc
