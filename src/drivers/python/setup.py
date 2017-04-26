import sys
import glob
from setuptools import setup, find_packages

def readme():
    with open('README.rst') as f:
        return f.read()

scripts = glob.glob('scripts/*.py')

# This part needs to be cleaned up :-)
# Some way we need to get the build path and the name of the .dll/.so
# Here I assume a link in the source/pyfabm to the library 

if sys.platform.startswith('win32'):
    lib = 'python_fabm.dll'
#    platform = 'win-x86_64'
elif sys.platform.startswith('linux'):
    lib = 'libpython_fabm.so'
#    platform = 'linux-x86_64'

setup(name='pyfabm',
      version='0.1',
      description='Python driver for FABM',
      long_description=readme(),
      url='https://github.com/fabm-model/fabm/tree/master/src/drivers/python',
      author='Jorn Bruggeman',
      author_email='jorn@bolding-bruggeman.com',
      license='GPL',
      scripts=scripts,
#      packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
      packages=['pyfabm'],
      package_data={'pyfabm': [lib]},
      include_package_data=True,
#      platforms = platform,
      zip_safe=False)


