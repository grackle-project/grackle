from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.extension import Extension

cython_extensions = [
    Extension(
        "pygrackle.grackle_wrapper",
        ["pygrackle/grackle_wrapper.pyx"],
        include_dirs=["../clib"],
        library_dirs=["../clib/.libs/"],
        libraries=["grackle"],
        define_macros=[
            ("CONFIG_BFLOAT_8", True),
            ("LARGE_INTS", True),
        ],
    ),
]


class build_ext(_build_ext):
    # subclass setuptools extension builder to avoid importing numpy
    # at top level in setup.py. See http://stackoverflow.com/a/21621689/1382869
    def finalize_options(self):
        _build_ext.finalize_options(self)
        # Prevent numpy from thinking it is still in its setup process
        # see http://stackoverflow.com/a/21621493/1382869
        __builtins__.__NUMPY_SETUP__ = False
        import numpy
        self.include_dirs.append(numpy.get_include())


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)
    config.add_subpackage('pygrackle', 'pygrackle')
    return config


setup(
    name="pygrackle",
    version="0.1",
    description="A wrapper for the Grackle chemistry library",
    packages=find_packages(),
    setup_requires=[
        'numpy',
        'cython',
    ],
    install_requires=[
        'setuptools',
        'numpy',
        'matplotlib',
    ],
    cmdclass={'build_ext': build_ext},
    license="BSD",
    ext_modules=cython_extensions
)
