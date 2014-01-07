#!/usr/bin/env python

import os
compiler = os.environ['FORTRAN_COMPILER']

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('fabm',parent_package,top_path)
    config.add_extension('fabm', ['fabm_python.F90'],
                         library_dirs=['../../../lib/python/%s' % compiler],
                         libraries=['fabm_prod'],
                         include_dirs=['../../../include','../../../modules/python/%s' % compiler])
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration,
          author='Jorn Bruggeman',
          author_email='jorn@bolding-burchard.com',
          version='0.2.0')


