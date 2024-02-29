from glob import glob
from setuptools import setup, Command
import os

try:
    from pybind11.setup_helpers import Pybind11Extension, build_ext
except ImportError:
    from setuptools import Extension as Pybind11Extension, build_ext

__version__ = "0.3.0"
name = "jetsimpy"

ext_modules = [
    Pybind11Extension(
        f"{name}.jetsimpy_extension",
        sorted(glob("jetsimpy/src/*.cpp")),
        define_macros = [('VERSION_INFO', __version__)],
        extra_compile_args=["-O3", "-std=c++11"],
        include_dirs = ["jetsimpy/src/"]
    ),
]

class CleanCommand(Command):
    """Custom clean command to tidy up the project root."""
    user_options = []
    def initialize_options(self):
        pass
    def finalize_options(self):
        pass
    def run(self):
        os.system('rm -vrf ./build ./dist ./*.pyc ./*.egg-info')


setup(
    name=name,
    version=__version__,
    author="Hao Wang",
    author_email="haowang.astro@gmail.com",
    description="Gamma-ray burst jet simulation & afterglow modeling",
    packages=[name],
    ext_modules=ext_modules,
    python_requires=">=3.6",
    install_requires=[
        'numpy',
    ],
    include_package_data=True,
    cmdclass={"clean": CleanCommand},
)

