#!/usr/bin/env python

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('fabm_fortran',parent_package,top_path)
    config.add_extension('fabm_fortran', ['fabm_python.F90'],
                         library_dirs=['../../../lib/python/GFORTRAN'],
                         libraries=['fabm_prod'],
                         include_dirs=['../../../include','../../../modules/python/GFORTRAN'])
    return config

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(configuration=configuration,
          author='Jorn Bruggeman',
          author_email='jbr@pml.ac.uk',
          version='0.2.0')


