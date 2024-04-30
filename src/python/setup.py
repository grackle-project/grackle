import glob
import os
import shutil
import subprocess
import sys
from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.extension import Extension

# helper functions:

def _getenv_bool(key, default = None):
    val = os.getenv(key, default = default)
    if val is None:
        return val
    elif val not in ['0','1']:
        raise RuntimeError(f"when set, {key!r} env variable must be 0 or 1")
    return bool(int(val))

def _find_solib_basename(dirpath, nominal_root = 'libgrackle',
                         return_single = True):
    # try to find the exact name of the library in the specified directory
    _paths = set()
    for cur_path in glob.iglob(f'{dirpath}/{nominal_root}*.so'):
        if os.path.islink(cur_path) and return_single:
            cur_path = os.path.abspath(os.path.realpath(cur_path))
            if not cur_path.startswith(os.path.abspath(dirpath)):
                raise ValueError("symlink links outside of {dirpath!r}")
        _paths.add(os.path.abspath(cur_path))
    if not return_single:
        return [os.path.basename(p) for p in _paths]
    elif len(_paths) == 0:
        return None
    elif len(_paths) == 1:
        return os.path.basename(_paths.pop())
    raise ValueError("Can't identify shared lib.")

def _config_kw(package_locallib_copy, library_dir):
    if package_locallib_copy:
        basename = _find_solib_basename(library_dir)
        assert basename is not None
        # copy the shared library into the pygrackle directory
        shutil.copy2(f'{library_dir}/{basename}', './pygrackle')

        extension_kw = {'library_dirs' : [os.path.abspath('./pygrackle')],
                        'libraries' : [ basename[3:].split('.so')[0] ] }
        # library_dirs tells the linker where to search for lib during linking,
        # next, we specify where to search at runtime (we use ORIGIN to
        # indicate that the library should be right next to the module)
        if sys.platform.startswith('darwin'): # untested
            extension_kw['extra_link_args'] = '-Wl,rpath,$ORIGIN'
        else:
            extension_kw['runtime_library_dirs'] = ['$ORIGIN']

        # this ensures that the shared library will be properly copied out of
        # the source-directory into the installed python module
        setup_pkgdata_kw = { 'package_data' : {"": [basename]} }
    else:
        # library_dirs used to find libgrackle while building the extension.
        # The system installation is linked at runtime
        extension_kw = {'library_dirs' : [library_dir],
                        'libraries' : ["grackle"]}
        setup_pkgdata_kw = {}

    return extension_kw, setup_pkgdata_kw


# Process the env-variable configuration options:
# - PYGRACKLE_PACKAGE_LOCALLIB:
#   -> should be set to a value of 0 or 1
#   -> in both cases, the copy of the libgrackle shared library from the
#      build-tree is used when linking the entension module. However, the
#      libgrakle shared library used at runtime differs.
#   -> when this is set to `1`, a copy of that shared library is packaged with
#      the python module. The python module is configured to always loads that
#      packaged library at runtime
#   -> when set to a value of 0, we use the historical behavior. The libgrackle
#      shared library is NOT packaged with the python module. At runtime when
#      extension module is loaded, we depend on the operating system's dynamic
#      linker to search for the libgrackle shared library (e.g. first in paths
#      given by LD_LIBRARY_PATH and then at other standard system locations)
# - PYGRACKLE_CMAKE_BUILD_DIR
#   -> this variable should be used when you want to indicate that the python
#      extension module should be built against a copy of libgrackle that was
#      constructed using cmake
#   -> in more detail, this should specify the cmake build directory. A common
#      choice is "../../build"
package_locallib_copy = _getenv_bool("PYGRACKLE_PACKAGE_LOCALLIB", None)
grackle_build_dir = os.getenv("PYGRACKLE_CMAKE_BUILD_DIR", "")

# handle different choices of grackle_build_dir
if grackle_build_dir == "": # traditional in-source build
    TRADITIONAL_IN_SOURCE_BUILD_MACRO = '1'
    library_dir="../clib/.libs/"
    include_dirs = ["../clib"]
    if package_locallib_copy is None: # fallback to historic behavior
        package_locallib_copy = False

else: # CMAKE-based out-of-source build
    abs_grackle_build_dir = os.path.abspath(
        path = os.path.expandvars(os.path.expanduser(grackle_build_dir))
    )
    TRADITIONAL_IN_SOURCE_BUILD_MACRO = '0'
    library_dir=os.path.join(abs_grackle_build_dir, "src/clib")

    if not os.path.isdir(abs_grackle_build_dir):
        raise RuntimeError(f"{grackle_build_dir} is not a directory")

    # unlike the traditional in-source build, the cmake approach puts the
    # the generated headers in a different location
    # -> therefore, we specify location where both sets of headers are found
    # -> to be safe, ensure that there are no autogened files from alt approach
    include_dirs = [
        "../clib",
        os.path.join(abs_grackle_build_dir, "generated_public_headers")
    ]
    # don't bother checking result code. This will fail if we've never built
    # grackle (if the depend file is missing
    subprocess.run(["make", "clean_autogen"], cwd = "../clib")

    if package_locallib_copy is None:
        package_locallib_copy = True
    elif package_locallib_copy == False:
        raise ValueError("we haven't tested this configuration.")

if package_locallib_copy and (_find_solib_basename(library_dir) is None):
    raise RuntimeError(
        f"no shared library found in {library_dir!r}. If using cmake, did you "
        "remember to: "
        "1. configure the build with -DBUILD_SHARED_LIBS=ON AND "
        "2. actually execute the build after configuring it?"
        "3. use the PYGRACKLE_CMAKE_BUILD_DIR variable to tell setup.py where "
        "   libgrackle was built"
    )

# for added safety, remove copies of the libgrackle shared library that were
# previously copied into the python directory
for so_basename in _find_solib_basename('./pygrackle', return_single = False):
    os.remove(f'./pygrackle/{so_basename}')

# handle different choices of package_locallib_copy
extension_kw, setup_pkgdata_kw = _config_kw(package_locallib_copy, library_dir)

cython_extensions = [
    Extension(
        "pygrackle.grackle_wrapper",
        ["pygrackle/grackle_wrapper.pyx"],
        include_dirs=include_dirs,
        define_macros=[
            ("CONFIG_BFLOAT_8", True),
            ("TRADITIONAL_IN_SOURCE_BUILD", TRADITIONAL_IN_SOURCE_BUILD_MACRO),
        ],
        **extension_kw
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
    version="1.0.0",
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
    python_requires='>=3.7',
    **setup_pkgdata_kw
)
