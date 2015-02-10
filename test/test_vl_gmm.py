#!/usr/bin/python
from _vlfeat import *
from numpy import *

random.seed(0)

I = asfortranarray(genfromtxt('lena.txt'), dtype=float32)
f, d = vl_dsift(I, fast=True, norm=True, step=5, floatDescriptors=True, verbose=True)
means, covs, priors, _ = vl_gmm(d, 30)
enc = vl_fisher(d, means, covs, priors)

print means
print covs
print enc
print enc.shape
print linalg.norm(enc,2)

try:
    magic = get_ipython().magic
    magic(u'%timeit f, d = vl_sift(I)') # 182 ms
    import vlfeat # pyvlfeat
    magic(u'%timeit f, d = vlfeat.vl_sift(I)') # 216 ms
    # MATLAB 180 ms
    
except:
    print 'use ipython <script>.py to see speed'
