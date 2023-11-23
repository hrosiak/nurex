from pathlib import Path

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

DIR = Path(__file__).parents[0]

SRC = [str((DIR/'nurex.cpp').resolve())]
example_module = Pybind11Extension(
    'pynurex',
    SRC,
    include_dirs=['../build/include'],
    library_dirs=['../build/lib','../build','../build/Release'],
    libraries=['nurex']
)

setup(
    name='pynurex',
    version=1.3,
    author='Andrej Prochazka',
    author_email='hrocho@vodacionline.sk',
    description='python interface to nurex library',
    long_description = (DIR/"README.md").read_text(),
    long_description_content_type='text/markdown',
    url='https://github.com/hrosiak/nurex',
    ext_modules=[example_module],
    cmdclass={"build_ext": build_ext},
)
