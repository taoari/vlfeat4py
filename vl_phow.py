# credit: https://github.com/shackenberg/phow_caltech101.py/blob/master/vl_phow.py
# author: Tao Wei <taowei@buffalo.edu>
# about 4 times faster than the original vl_phow.py

import numpy as np
import _vlfeat
from _vlfeat import vl_dsift
import cv2
from sys import maxint

"""
Python rewrite of https://github.com/vlfeat/vlfeat/blob/master/toolbox/sift/vl_phow.m
### notice no hsv support atm
### comments are largely copied from the code

"""

def vl_imsmooth(I, sigma):
    '''vl_imsmooth.
    
    Currently only comptabile with vl_imsmooth for grayscale images'''
    if I.ndim == 2:
        return _vlfeat.vl_imsmooth(I, sigma)
        
    w = int(np.ceil(4.0*sigma))
    ksize = 2*w + 1
    I = cv2.GaussianBlur(I, (ksize, ksize), sigma)
    return np.asfortranarray(I, dtype=I.dtype)
    
def vl_phow(I,
            verbose=False,
            fast=True,
            sizes=[4, 6, 8, 10],
            step=2,
            color='gray',
            floatDescriptors=False,
            magnif=6,
            windowSize=1.5,
            contrastThreshold=0.005):
    '''extracts PHOW features [1] from the image IM. 
    PHOW is simply dense SIFT applied at several resolutions.
    
    Parameters
    ----------
    
    data : float32 ndarray
        float32 image
        
    Returns
    -------
    
    frames : float64 ndarray with shape (4, n_descrs)
        Dense SIFT DoG keypoint frames
        
    descrs : float32 ndarray with shape (128, n_descrs)
        Dense SIFT descriptors at frames
    
    Options
    -------
    Verbose : false
      Set to true to turn on verbose output.
 
    Sizes : [4 6 8 10]
      Scales at which the dense SIFT features are extracted. Each
      value is used as bin size for the VL_DSIFT() function.
 
    Fast : true
      Set to false to turn off the fast SIFT features computation by
      VL_DSIFT().
 
    Step : 2
      Step (in pixels) of the grid at which the dense SIFT features
      are extracted.
 
    Color : 'gray'
      Choose between 'gray' (PHOW-gray), 'rgb', 'hsv', and 'opponent'
      (PHOW-color).
 
    ContrastThreshold : 0.005
      Contrast threshold below which SIFT features are mapped to
      zero. The input image is scaled to have intensity range in [0,1]
      (rather than [0,255]) and this value is compared to the
      descriptor norm as returned by VL_DSIFT().
 
    WindowSize : 1.5
      Size of the Gaussian window in units of spatial bins.
 
    Magnif : 6
      The image is smoothed by a Gaussian kernel of standard deviation
      SIZE / MAGNIF. Note that, in the standard SIFT descriptor, the
      magnification value is 3; here the default one is 6 as it seems
      to perform better in applications.
 
    FloatDescriptors : false
      If set to TRUE, the descriptors are returned in floating point
      format.
    '''

    opts = Options(verbose, fast, sizes, step, color, floatDescriptors,
                   magnif, windowSize, contrastThreshold)
    dsiftOpts = DSiftOptions(opts)

    # Extract the features
    imageSize = I.shape
    if I.ndim == 3:
        if imageSize[2] != 3:
            # "IndexError: tuple index out of range" if both if's are checked at the same time
            raise ValueError("Image data in unknown format/shape")
    if opts.color == 'gray':
        numChannels = 1
        if (I.ndim != 2):
            I = cv2.cvtColor(I, cv2.COLOR_RGB2GRAY)
    else:
        numChannels = 3
        if (I.ndim == 2):
            I = np.dstack([I, I, I])
        if opts.color == 'rgb':
            pass
        elif opts.color == 'hsv':
            I = cv2.cvtColor(I, cv2.COLOR_RGB2HSV)
        elif opts.color == 'opponent':
             # from https://github.com/vlfeat/vlfeat/blob/master/toolbox/sift/vl_phow.m
             # Note that the mean differs from the standard definition of opponent
             # space and is the regular intesity (for compatibility with
             # the contrast thresholding).
             # Note also that the mean is added pack to the other two
             # components with a small multipliers for monochromatic
             # regions.

            mu = 0.3 * I[:, :, 0] + 0.59 * I[:, :, 1] + 0.11 * I[:, :, 2]
            alpha = 0.01
            I = np.dstack([mu,
                         (I[:, :, 0] - I[:, :, 1]) / np.sqrt(2) + alpha * mu,
                         (I[:, :, 0] + I[:, :, 1] - 2 * I[:, :, 2]) / np.sqrt(6) + alpha * mu])
        else:
            raise ValueError('Color option ' + str(opts.color) + ' not recognized')
    if opts.verbose:
        print('{0}: color space: {1}'.format('vl_phow', opts.color))
        print('{0}: image size: {1} x {2}'.format('vl_phow', imageSize[0], imageSize[1]))
        print('{0}: sizes: [{1}]'.format('vl_phow', opts.sizes))
    
    # make sure image I is float32, f_order
    I = np.asfortranarray(I, dtype=np.float32)

    frames_all = []
    descrs_all = []
    for size_of_spatial_bins in opts.sizes:
        # from https://github.com/vlfeat/vlfeat/blob/master/toolbox/sift/vl_phow.m
        # Recall from VL_DSIFT() that the first descriptor for scale SIZE has
        # center located at XC = XMIN + 3/2 SIZE (the Y coordinate is
        # similar). It is convenient to align the descriptors at different
        # scales so that they have the same geometric centers. For the
        # maximum size we pick XMIN = 1 and we get centers starting from
        # XC = 1 + 3/2 MAX(OPTS.SIZES). For any other scale we pick XMIN so
        # that XMIN + 3/2 SIZE = 1 + 3/2 MAX(OPTS.SIZES).
        # In pracrice, the offset must be integer ('bounds'), so the
        # alignment works properly only if all OPTS.SZES are even or odd.

        off = np.floor(3.0 / 2 * (max(opts.sizes) - size_of_spatial_bins)) + 1

        # smooth the image to the appropriate scale based on the size
        # of the SIFT bins
        sigma = size_of_spatial_bins / float(opts.magnif)
        ims = vl_imsmooth(I, sigma)

        if opts.verbose:
            print('smooth sigma: %.2f' % sigma)
        if ims.ndim == 2:
            ims = ims[:,:,np.newaxis]

        # extract dense SIFT features from all channels
        frames = []
        descrs = []
        for k in range(numChannels):
            size_of_spatial_bins = int(size_of_spatial_bins)
            # vl_dsift does not accept numpy.int64 or similar
            f_temp, d_temp = vl_dsift(ims[:, :, k],
                                      step=dsiftOpts.step,
                                      size=size_of_spatial_bins,
                                      fast=dsiftOpts.fast,
                                      floatDescriptors=dsiftOpts.floatDescriptors,
                                      verbose=dsiftOpts.verbose,
                                      norm=dsiftOpts.norm,
                                      bounds=[off, off, maxint, maxint],
                                      windowSize=dsiftOpts.windowSize)
            if (not opts.floatDescriptors):
                d_temp = np.floor(d_temp);
            frames.append(f_temp)
            descrs.append(d_temp)
            
        if (opts.color == 'gray') or (opts.color == 'opponent'):
            contrast = frames[0][2,:]
        elif opts.color == 'rgb':
            contrast = np.mean([frames[0][2,:], frames[1][2,:], frames[2][2,:]], axis=0)
        elif opts.color == 'hsv':
            contrast = frames[2][2,:]
        else:
            raise ValueError('Color option ' + str(opts.color) + ' not recognized')
        
        for k in xrange(numChannels):
            descrs[k][:,contrast < opts.contrastThreshold] = 0
            

        _pos = frames[0][0:2,:]   # x, y
        _contrast = contrast
        _binSize = size_of_spatial_bins * np.ones(frames[0].shape[1])
        
        frames_all.append(np.vstack((_pos, _contrast, _binSize)))
        descrs_all.append(np.vstack((descrs)))

    frames_all = np.hstack(frames_all)
    descrs_all = np.hstack(descrs_all)
    return frames_all, descrs_all


class Options(object):
    def __init__(self, verbose, fast, sizes, step, color,
            floatDescriptors, magnif, windowSize,
            contrastThreshold):
        self.verbose = verbose
        self.fast = fast
        if (type(sizes) is not np.ndarray) & (type(sizes) is not list):
            sizes = np.asarray([sizes])
        self.sizes = sizes
        self.step = step
        self.color = color
        self.floatDescriptors = floatDescriptors
        self.magnif = magnif
        self.windowSize = windowSize
        self.contrastThreshold = contrastThreshold


class DSiftOptions(object):
    def __init__(self, opts):
        self.norm = True
        self.windowSize = opts.windowSize
        self.verbose = opts.verbose
        self.fast = opts.fast
        self.floatDescriptors = opts.floatDescriptors
        self.step = opts.step

if __name__ == "__main__":
    from scipy.misc import lena
    import time
    I = np.float32(lena())
    # I = cv2.imread('test/lena_color.jpg')
    __TIC = time.time()
    frames, descrs = vl_phow(I, color='hsv', verbose=1) 
    print time.time() - __TIC
