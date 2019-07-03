import os.path
from setuptools import setup

try:
    import wheel.bdist_wheel
    class bdist_wheel(wheel.bdist_wheel.bdist_wheel):
        def finalize_options(self):
            wheel.bdist_wheel.bdist_wheel.finalize_options(self)
            self.root_is_pure = False
        def get_tag(self):
            python, abi, plat = wheel.bdist_wheel.bdist_wheel.get_tag(self)
            python, abi = 'py2.py3', 'none'
            return python, abi, plat
except ImportError:
    bdist_wheel = None

def readme():
    with open(os.path.join(os.path.dirname(__file__), 'README.rst'), 'rU') as f:
        return f.read()

setup(
    name='pyfabm',
    version='0.1',
    description='Python driver for FABM',
    long_description=readme(),
    long_description_content_type='text/x-rst',
    url='https://github.com/fabm-model/fabm/tree/master/src/drivers/python',
    author='Jorn Bruggeman',
    author_email='jorn@bolding-bruggeman.com',
    license='GPL',
    entry_points={
        'console_scripts': [
            'fabm_complete_yaml=pyfabm.utils.fabm_complete_yaml:main',
            'fabm_configuration_gui=pyfabm.utils.fabm_configuration_gui:main',
            'fabm_describe_model=pyfabm.utils.fabm_describe_model:main',
            'fabm_evaluate=pyfabm.utils.fabm_evaluate:main',
            'fabm_stress_test=pyfabm.utils.fabm_stress_test:main',
        ]
    },
    packages=['pyfabm', 'pyfabm/utils'],
    package_data={'pyfabm': ['*.so', '*.dll', '*.dylib']},
    cmdclass={'bdist_wheel': bdist_wheel},
    zip_safe=False
)


