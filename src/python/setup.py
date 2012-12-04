import os
import os.path
import glob
import sys
import time

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(ignore_setup_xxx_py=True,
                       assume_default_configuration=True,
                       delegate_options_to_subpackages=True,
                       quiet=True)
    config.add_subpackage('pygrackle', 'pygrackle')
    return config


def setup_package():

    from numpy.distutils.core import setup

    setup(
        name="pygrackle",
        version="0.1",
        description="A wrapper for the Grackle chemistry library",
        configuration=configuration,
        )
    return

######
# This next bit comes from Matthew Brett, to get Cython working with NumPy
# distutils.  I added a bit to get C++ Cython working.
from os.path import join as pjoin, dirname
from distutils.dep_util import newer_group
from distutils.errors import DistutilsError
from numpy.distutils import log
from numpy.distutils.misc_util import appendpath

def generate_a_pyrex_source(self, base, ext_name, source, extension):
    ''' Monkey patch for numpy build_src.build_src method

    Uses Cython instead of Pyrex.

    Assumes Cython is present
    '''
    if self.inplace:
        target_dir = dirname(base)
    else:
        target_dir = appendpath(self.build_src, dirname(base))
    if extension.language == "c++":
        cplus = True
        file_ext = ".cpp"
    else:
        cplus = False
        file_ext = ".c"
    target_file = pjoin(target_dir, ext_name + file_ext)
    depends = [source] + extension.depends
    if self.force or newer_group(depends, target_file, 'newer'):
        import Cython.Compiler.Main
        log.info("cythonc:> %s" % (target_file))
        self.mkpath(target_dir)
        options = Cython.Compiler.Main.CompilationOptions(
            defaults=Cython.Compiler.Main.default_options,
            include_path=extension.include_dirs,
            language=extension.language, cplus=cplus,
            output_file=target_file)
        cython_result = Cython.Compiler.Main.compile(source,
                                                   options=options)
        if cython_result.num_errors != 0:
            raise DistutilsError("%d errors while compiling %r with Cython" \
                  % (cython_result.num_errors, source))
    return target_file

from numpy.distutils.command import build_src
build_src.build_src.generate_a_pyrex_source = generate_a_pyrex_source
# End snippet
######

if __name__ == '__main__':
    setup_package()
