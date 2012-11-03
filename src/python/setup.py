from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("grackle",
        ["grackle.pyx"],
        language="c++",
        include_dirs=["../clib/"],
        library_dirs=["../clib/"],
        libraries=["grackle"],
        define_macros=[("BFLOAT_8", True),
                       ("PFLOAT_8", True)]
        )]
)
