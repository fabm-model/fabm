import sys
import glob
from setuptools import setup, find_packages

#def readme():
#    with open('README.rst') as f:
#        return f.read()

scripts = glob.glob('scripts/*.py')

# This part needs to be cleaned up :-)
if sys.platform.startswith('win32'):
    lib = 'python_fabm.dll'
    pkgpath = 'pyfabm'
    libpath = 'build/python_fabm.dll'
    platform = 'win-x86_64'
elif sys.platform.startswith('linux'):
    lib = 'libpython_fabm.so'
    pkgpath = 'lib/python2.7/site-packages/pyfabm'
#   some way we need to get the build path - here I just made a link in the source dir and named the link build
    libpath = 'build/libpython_fabm.so'
    platform = 'linux-x86_64'

setup(name='pyfabm',
      version='0.1',
      description='Python driver for FABM',
#      long_description=readme(),
      url='https://github.com/fabm-model/fabm/tree/master/src/drivers/python',
      author='Jorn Bruggeman',
      author_email='jorn@bolding-bruggeman.com',
      license='GPL',
      scripts=scripts,
#      packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
      packages=['pyfabm'],
#      package_data={'pyfabm': [lib]},
      include_package_data=True,
      data_files=[(pkgpath, [libpath])],
      platforms=platform,
      zip_safe=False)


