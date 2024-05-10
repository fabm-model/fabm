import os
import subprocess
import shutil
from setuptools import Extension
from setuptools.command.build_ext import build_ext
import shlex


try:
    import wheel.bdist_wheel
except ImportError:
    raise Exception(
        "wheel must be installed to build pyfabm. Try 'python -m pip install wheel'."
    )


class bdist_wheel(wheel.bdist_wheel.bdist_wheel):
    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False

    def get_tag(self):
        python, abi, plat = super().get_tag()
        python, abi = "py3", "none"
        return python, abi, plat


class CMakeExtension(Extension):
    def __init__(self, name, *cmake_args):
        super().__init__(name, sources=[])
        self.cmake_args = cmake_args


class CMakeBuild(build_ext):
    user_options = build_ext.user_options + [
        ("cmake-opts=", None, "additional options to pass to cmake"),
    ]

    def run(self):
        for ext in self.extensions:
            self.build_extension(ext)

    def initialize_options(self):
        super().initialize_options()
        self.cmake_opts = None

    def build_extension(self, ext):
        if not os.path.isdir(self.build_temp):
            os.makedirs(self.build_temp)

        ext_path = self.get_ext_fullpath(ext.name)
        libname = ext.name.rsplit(".", 1)[-1]

        # Directory where your build output should go
        install_prefix = os.path.abspath(os.path.dirname(ext_path))

        # Temporary directory where all intermediate build files should go.
        build_dir = os.path.join(self.build_temp, ext.name)
        if self.force and os.path.isdir(build_dir):
            print(f"Emptying existing build directory {build_dir}")
            shutil.rmtree(build_dir)
        if not os.path.isdir(build_dir):
            os.makedirs(build_dir)

        build_type = "Debug" if self.debug else "Release"
        cmake_args = list(ext.cmake_args) + [f"-DCMAKE_BUILD_TYPE={build_type}"]
        if self.cmake_opts is not None:
            cmake_args += shlex.split(self.cmake_opts)
        if self.compiler is not None:
            cmake_args.append(f"-DCMAKE_Fortran_COMPILER={self.compiler}")
        subprocess.check_call(
            [
                "cmake",
                "-B",
                build_dir,
                "-S",
                os.path.dirname(__file__),
                f"-DPYFABM_NAME={libname}",
                f"-DPYFABM_DIR={install_prefix}",
            ]
            + cmake_args,
        )
        subprocess.check_call(["cmake", "--build", build_dir, "--config", build_type])
