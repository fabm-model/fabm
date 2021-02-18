import os.path
import subprocess
import io
import shutil
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

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
    raise Exception('wheel must be installed to build pyfabm. Try "python -m pip install wheel".')

def readme():
    with io.open(os.path.join(os.path.dirname(__file__), 'README.rst'), 'r') as f:
        return f.read()

class CMakeExtension(Extension):
    def __init__(self, name, sourcedir='', *cmake_args):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)
        self.cmake_args = cmake_args

class CMakeBuild(build_ext):
    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        if not os.path.isdir(self.build_temp):
            os.makedirs(self.build_temp)

        ext_path = self.get_ext_fullpath(ext.name)
        libname = ext.name.rsplit('.', 1)[-1]

        # Directory where your build output should go
        install_prefix = os.path.abspath(os.path.dirname(ext_path))

        # Temporary directory where all intermediate build files should go.
        build_dir = os.path.join(self.build_temp, ext.name)
        if self.force and os.path.isdir(build_dir):
            shutil.rmtree(build_dir)
        if not os.path.isdir(build_dir):
            os.makedirs(build_dir)

        build_type = 'Debug' if self.debug else 'Release'
        cmake_args = list(ext.cmake_args) + ['-DCMAKE_BUILD_TYPE=%s' % build_type]
        if self.compiler is not None:
            cmake_args.append('-DCMAKE_Fortran_COMPILER=%s' % self.compiler)
        subprocess.check_call(['cmake', ext.sourcedir, '-DPYFABM_NAME=%s' % libname, '-DPYFABM_DIR=%s' % install_prefix] + cmake_args, cwd=build_dir)
        subprocess.check_call(['cmake', '--build', '.', '--config', build_type], cwd=build_dir)

setup(
    name='pyfabm',
    version='1.0.2',
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
    ext_modules=[CMakeExtension('pyfabm.fabm_0d', '.'), CMakeExtension('pyfabm.fabm_1d', '.', '-DPYFABM_DIM_COUNT=1')],
    cmdclass={'bdist_wheel': bdist_wheel, 'build_ext': CMakeBuild},
    zip_safe=False
)


