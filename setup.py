from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext = Extension('_vlfeat',
                include_dirs=['include/'],
                sources=['_vlfeat.pyx', 
                'src/vl_sift.cpp', 
                'src/vl_dsift.cpp', 
                'src/vl_imsmooth.cpp', 
                'src/vl_homkermap.cpp', 
                'src/vl_kmeans.cpp',
                'src/vl_gmm.cpp',
                'src/vl_fisher.cpp'],
                libraries = ['armadillo', 'vl'],
                language='c++',
                extra_compile_args=['-D NPY_NO_DEPRECATED_API'])

setup(ext_modules=[ext],
      cmdclass={'build_ext': build_ext})
