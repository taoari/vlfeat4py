from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext = Extension('_vlfeat',
                sources=['_vlfeat.pyx', 
                'vl_sift.cpp', 
                'vl_dsift.cpp', 
                'vl_imsmooth.cpp', 
                'vl_homkermap.cpp', 
                'vl_kmeans.cpp'],
                libraries = ['armadillo', 'vl'],
                language='c++',
                extra_compile_args=['-D NPY_NO_DEPRECATED_API'])

setup(ext_modules=[ext],
      cmdclass={'build_ext': build_ext})
