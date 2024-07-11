import sys
from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.extension import Extension

cython_extensions = [
    Extension(
        "pygrackle.grackle_wrapper",
        ["pygrackle/grackle_wrapper.pyx"],
        include_dirs=["../clib/autogen", "../include"],
        library_dirs=["../clib/.libs/"],
        libraries=["grackle"],
        define_macros=[
            ("CONFIG_BFLOAT_8", True),
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
    'pytest',
    'sphinx',
    'packaging'
]

setup(
    name="pygrackle",
    version="1.0.1",
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
