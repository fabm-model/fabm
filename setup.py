import os
import sys
from setuptools import setup

FABM_BASE = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(FABM_BASE, "src/drivers/python"))

from build_cmake import CMakeExtension, bdist_wheel, CMakeBuild

setup(
    packages=["pyfabm", "pyfabm.utils"],
    package_dir={"": "src"},
    ext_modules=[
        CMakeExtension(
            "pyfabm.fabm_0d", "-DPYFABM_DEFINITIONS=_FABM_DIMENSION_COUNT_=0"
        ),
        CMakeExtension(
            "pyfabm.fabm_1d", "-DPYFABM_DEFINITIONS=_FABM_DIMENSION_COUNT_=1"
        ),
    ],
    cmdclass={"bdist_wheel": bdist_wheel, "build_ext": CMakeBuild},
    zip_safe=False,
)
