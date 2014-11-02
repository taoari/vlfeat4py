#!/usr/bin/python
from _vlfeat import *
from numpy import *

I = asfortranarray(genfromtxt('lena.txt'), dtype=float32)
Is = vl_imsmooth(I, 5.0, verbose=1)
savetxt('lena_smooth_python.txt', Is, fmt='%.6f')

I = float64(I)
Is = vl_imsmooth(I, 5.0, kernel="triangular", verbose=1)
