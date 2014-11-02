#!/usr/bin/python
from _vlfeat import *
from numpy import *
from scipy.io import loadmat

r = loadmat('rand.mat')['r']
r = asfortranarray(r, dtype=float32)

# X = float32(random.rand(5,100))
X = r
Y = vl_kmeans(X, 10, verbose=1)

print X.shape
print Y.shape
