from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('pygrackle',parent_package,top_path)
    config.add_extension("grackle_wrapper",
                ["pygrackle/grackle_wrapper.pyx"],
                include_dirs=["../clib/"],
                library_dirs=["../clib/.libs/"],
                libraries=["grackle"],
                define_macros=[("CONFIG_BFLOAT_8", True),
                               ("LARGE_INTS", True)]
    )
    return config
