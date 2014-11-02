#!/usr/bin/python
from vl_phow import vl_phow
from numpy import *
from scipy.io import savemat, loadmat

I = loadmat('lena.mat')['lena_gray']
I = asfortranarray(I, dtype=float32)

f, d = vl_phow(I, step=50, floatDescriptors=True, verbose=True)
savetxt('phow_f_python.txt', f.T, fmt='%.2f')
savetxt('phow_d_python.txt', d.T, fmt='%.2f')
