#!/usr/bin/python
from _vlfeat import *
from numpy import *
from scipy.io import loadmat

r = loadmat('rand.mat')['r']
r = asfortranarray(r, dtype=float32)

y = vl_homkermap(r, 1)

print r.shape
print y.shape
