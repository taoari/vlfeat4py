#!/usr/bin/python
from _vlfeat import *
from numpy import *
from scipy.io import savemat, loadmat

I = loadmat('lena.mat')['lena_gray']
I = asfortranarray(I, dtype=float32)/255.0 
# SIFT frames and descriptors are independent of image scales

f, d = vl_sift(I, floatDescriptors=True, verbose=True)
savetxt('sift_f_python.txt', f, fmt='%.6f')
savetxt('sift_d_python.txt', d, fmt='%.6f')

try:
    magic = get_ipython().magic
    magic(u'%timeit f, d = vl_sift(I)') # 182 ms
    # MATLAB 180 ms
    
except:
    print 'use ipython <script>.py to see speed'
