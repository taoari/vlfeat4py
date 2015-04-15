VLFeat4Py: A wrapper library for VLFeat
=======================================

This is a Python wrapper library for the open source computer vision library [VLFeat](http://www.vlfeat.org/)
using the [Armadillo](http://arma.sourceforge.net/) C++ Matrix Library.

VLFeat4Py is not targeted to be a complete wrapper for VLFeat, it is best used with the combination of 
[sklearn](http://scikit-learn.org/stable/) and [skimage](http://scikit-image.org/docs/dev/api/skimage.html).

Current wrapped functions are:

* vl_sift
* vl_dsift
* vl_imsmooth (obsolted, use `skimage.filters.gaussian_filter` instead)
* vl_phow
* vl_kmeans
* vl_gmm
* vl_fisher

INSTALL
-------

1. Install Armadillo
2. Install VLFeat (Properly setup include files and libvl.so)
3. Run `python setup.py build_ext --inplace`
4. Copy `_vlfeat.so` and `vl_phow.py` under folder `vlfeat/`, and move to PYTHONPATH directory

Usage
-----

```
from vlfeat import vl_dsift, vl_phow
```
