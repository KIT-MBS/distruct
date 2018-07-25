#!/usr/bin/env python3
#####################################
#
# Filename : setup.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Thu 28 Jun 2018 12:50:34 PM CEST
#
# Last Modified : Wed 25 Jul 2018 05:08:13 PM CEST
#
#####################################

# TODO look at what kind of structure/files pypi expects
# TODO install requirements

from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize


USE_CYTHON = True
fext = '.pyx' if USE_CYTHON else '.cpp'

def networkit_get_include():
    # TODO just wget the needed nwk headers
    return

sources = ['MOBi/src/BioMaxentStress.cpp', "MOBi/_MOBi" + fext]

# TODO correct include path
from os import environ
includeDirs = [environ["HOME"] + "/Projects/networkit/networkit/", "MOBi/src/"]

# TODO determine lib path from networkit install path
libraryDirs = ["/usr/lib/python3.6/site-packages/"]

# TODO determine lib name at install time
libraries = [":_NetworKit.cpython-36m-x86_64-linux-gnu.so"]

compile_args = ["-fopenmp"]
link_args = ["-fopenmp"]

# TODO basically use the pythonpath as LD_LIBRARY_PATH

extensions = [
        Extension(
            "_MOBi",
            sources,
            language = "c++",
            extra_compile_args = compile_args,
            extra_link_args = link_args,
            include_dirs = includeDirs,
            libraries = libraries,
            library_dirs = libraryDirs)
]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions, language = 'c++')
    pass

setup(
        name = "MOBi",
        ext_modules = extensions,
        author = "Oskar Taubert",
        author_email = "oskar.taubert@kit.edu",
        url = "TODO",
        download_url = "TODO",
        packages = find_packages(),
        package_data = {},
        keywords = "",
        platforms = "",
        classifiers = "",
        zip_safe = False)
