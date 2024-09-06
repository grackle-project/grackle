import os
import sys
import warnings
from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.extension import Extension

# Process the PYGRACKLE_CMAKE_BUILD_DIR env-variable configuration option

grackle_build_dir = os.getenv("PYGRACKLE_CMAKE_BUILD_DIR", "")
if grackle_build_dir == "": # traditional in-source build
    TRADITIONAL_IN_SOURCE_BUILD_MACRO = '1'
    autogen_public_hdr_dir = "../clib/autogen"
    library_dir = "../clib/.libs/"
else: # CMAKE-based out-of-source build
    TRADITIONAL_IN_SOURCE_BUILD_MACRO = '0'
    autogen_public_hdr_dir = f"{grackle_build_dir}/generated_public_headers"
    if os.path.isdir(f"{grackle_build_dir}/grackle/lib"):
        library_dir = f"{grackle_build_dir}/grackle/lib"
    elif os.path.isdir(f"{grackle_build_dir}/grackle/lib64"):
        # preferred on Red Hat related Linux distributions (like Fedora)
        library_dir = f"{grackle_build_dir}/grackle/lib64"
    else:
        raise RuntimeError("could not find the cmake-based grackle build dir")

if not os.path.isfile(f'{library_dir}/libgrackle.so'):
    _message = f"""\
libgrackle.so not found in {library_dir!r}.
Unless the library is installed in a system search-path, problems will arise.

If you are using using cmake, did you remember to:
  1. configure the build with -DBUILD_SHARED_LIBS=ON AND
  2. actually execute the build (& installation) after configuring it?
  3. use the PYGRACKLE_CMAKE_BUILD_DIR variable to tell setup.py where
     libgrackle.so was built?"""

    warnings.warn(_message)

cython_extensions = [
    Extension(
        "pygrackle.grackle_wrapper",
        ["pygrackle/grackle_wrapper.pyx"],
        include_dirs=[autogen_public_hdr_dir, "../include"],
        library_dirs=[library_dir],
        libraries=["grackle"],
        define_macros=[
            ("CONFIG_BFLOAT_8", True),
            ("TRADITIONAL_IN_SOURCE_BUILD", TRADITIONAL_IN_SOURCE_BUILD_MACRO),
        ],
    ),
]

# on some platforms the cython bindings don't work unless the
# language_level matches the python version. To specify the level
# see https://stackoverflow.com/a/58116368
for e in cython_extensions:
    e.cython_directives = {'language_level': sys.version_info[0]}

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


dev_requirements = [
    'flake8',
    'packaging',
    'pytest',
    'sphinx',
    'sphinx-tabs',
    'sphinx_rtd_theme',
]

setup(
    name="pygrackle",
    version="1.1.dev",
    description="A wrapper for the Grackle chemistry library",
    keywords=["simulation", "chemistry", "cooling", "astronomy", "astrophysics"],
    url="https://github.com/grackle-project/grackle",
    project_urls={
        'Homepage': 'https://github.com/grackle-project/grackle',
        'Documentation': 'https://grackle.readthedocs.io/',
        'Source': 'https://github.com/grackle-project/grackle',
        'Tracker': 'https://github.com/grackle-project/grackle/issues'
    },
    include_package_data=True,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: BSD License",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: POSIX :: Linux",
        "Operating System :: Unix",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    packages=find_packages(),
    setup_requires=[
        'numpy',
        'cython',
    ],
    install_requires=[
        'cython',
        'h5py',
        'numpy',
        'matplotlib',
        'yt>=4.0.2',
    ],
    cmdclass={'build_ext': build_ext},
    license="BSD 3-clause",
    ext_modules=cython_extensions,
    extras_require={
        'dev': dev_requirements,
    },
    python_requires='>=3.7'
)
