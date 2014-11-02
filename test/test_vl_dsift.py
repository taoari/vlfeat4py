#!/usr/bin/python
from _vlfeat import *
from numpy import *

I = asfortranarray(genfromtxt('lena.txt'), dtype=float32)
f, d = vl_dsift(I, fast=True, norm=True, step=100, floatDescriptors=True, verbose=True)
savetxt('dsift_f_python.txt', f, fmt='%.6f')
savetxt('dsift_d_python.txt', d, fmt='%.6f')

try:
    magic = get_ipython().magic
    magic(u'%timeit f, d = vl_sift(I)') # 182 ms
    import vlfeat # pyvlfeat
    magic(u'%timeit f, d = vlfeat.vl_sift(I)') # 216 ms
    # MATLAB 180 ms
    
except:
    print 'use ipython <script>.py to see speed'
