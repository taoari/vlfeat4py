import numpy as np
import cv2

def vl_imsmooth(I, sigma):
    w = int(np.ceil(4.0*sigma))
    ksize = 2*w + 1
    I = cv2.GaussianBlur(I, (ksize, ksize), sigma)
    return np.asfortranarray(I, dtype=I.dtype)

#def vl_imsmooth2(I, sigma):
#    w = int(np.ceil(4.0*sigma))
#    ksize = 2*w + 1
#    J = np.arange(ksize, dtype=I.dtype)
#    z = (J-w)/(sigma)
#    k = np.exp(-0.5*z**2)
#    k = np.c_[k/sum(k)]
#    kern = k * k.T
#    I = cv2.filter2D(I, -1, kern)
#    I = cv2.filter2D(I, -1, kern)
#    return I
    
if __name__ == '__main__':
    I = np.asfortranarray(np.genfromtxt('lena.txt'), dtype=np.float32)
    Is = vl_imsmooth(I, 5.0)
    Is2 = np.genfromtxt('lena_smooth_matlab.txt', dtype=np.float32)
    np.savetxt('lena_smooth_python_opencv.txt', Is, fmt='%.6f')
