import sys
import glob
import os
from setuptools import setup, find_packages

def readme():
    with open(os.path.join(os.path.dirname(__file__), 'README.rst'), 'rU') as f:
        return f.read()

setup(name='pyfabm',
      version='0.1',
      description='Python driver for FABM',
      long_description=readme(),
      url='https://github.com/fabm-model/fabm/tree/master/src/drivers/python',
      author='Jorn Bruggeman',
      author_email='jorn@bolding-bruggeman.com',
      license='GPL',
      entry_points={
          'console_scripts': [
                'fabm_complete_yaml=pyfabm.utils.fabm_complete_yaml:main',
                'fabm_configuration_gui=pyfabm.utils.fabm_configuration_gui:main',
                'fabm_describe_model=pyfabm.utils.fabm_describe_model:main',
#KB - script needs further checks                'fabm_evaluate=pyfabm.utils.fabm_evaluate:main',
                'fabm_stress_test=pyfabm.utils.fabm_stress_test:main',
          ]
      },
#      packages=find_packages(exclude=['contrib', 'docs', 'tests*']),
      packages=['pyfabm', 'pyfabm/utils'],
      #package_data={'pyfabm': [lib]},
      include_package_data=True,
#      platforms = platform,
      zip_safe=False)


