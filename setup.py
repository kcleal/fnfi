from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy


extensions = [
    Extension(
        "align.align_path_c",
        ["align/align_path_c.pyx"],
        library_dirs=[numpy.get_include()],
        language="c++",
        extra_compile_args=["-std=c++11", "-mmacosx-version-min=10.9"],
        extra_link_args=["-std=c++11"]
        ),
]

setup(
    name="fufi",
    version='0.2.0',
    packages=find_packages(),
    ext_modules=cythonize(extensions),
    include_dirs=[numpy.get_include()],
    include_package_data=True,
    install_requires=[
        'click',
        'numpy',
        'pandas',
        'pysam',
        'quicksect',
        'pybedtools',
        'natsort',
        'networkx>=2.0',
        'scikit-learn',
        'pytest'
    ],
    entry_points='''
        [console_scripts]
        fufi=src.fufi:cli
    ''',
)
