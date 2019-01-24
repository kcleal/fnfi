from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy


extensions = [
    Extension(
        "fnfi.align_path_c",
        ["fnfi/align_path_c.pyx"],
        library_dirs=[numpy.get_include()],
        language="c++",
        #extra_compile_args=["-std=c++11", "-mmacosx-version-min=10.9"],
        #extra_link_args=["-std=c++11"]
        ),

    Extension(
            "fnfi.c_io_funcs",
            ["fnfi/c_io_funcs.pyx"],
            library_dirs=[numpy.get_include()],
            language="c++",
            #extra_compile_args=["-std=c++11", "-mmacosx-version-min=10.9"],
            #extra_link_args=["-std=c++11"]
            ),

    Extension(
            "fnfi.c_samflags",
            ["fnfi/c_samflags.pyx"],
            library_dirs=[numpy.get_include()],
            language="c",
            ),
]
print("Found packages", find_packages(where="."))
setup(
    name="fnfi",
    version='0.8.0',
    packages=find_packages(where="."),
    ext_modules=cythonize(extensions),
    include_dirs=[numpy.get_include()],
    include_package_data=True,
    install_requires=[
        'click',
        'numpy',
        'pandas',
        'pysam',
        'pybedtools',
        'networkx',
        'scikit-learn',
        'ncls',
    ],
    entry_points='''
        [console_scripts]
        fnfi=fnfi.main:cli
    ''',
)
