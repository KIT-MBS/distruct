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
# Last Modified : Tue 03 Jul 2018 07:48:29 PM CEST
#
#####################################

# NOTE for now i use distutils. maybe have to switch to setuptools?
# TODO better way to find the networkit c++ headers?
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
# TODO download nwk c++ headers
from os import environ
includeDirs = [environ["HOME"] + "/Projects/networkit/networkit/", "MOBi/src/"]
libraryDirs = ["/usr/lib/python3.6/site-packages/"]
# TODO fix this
libraries = [libraryDirs[0] + "_NetworKit.cpython-36m-x86_64-linux-gnu.so"]
# TODO try global and local, use local first

# TODO optimization

extensions = [
        Extension(
            "_MOBi",
            sources,
            language = "c++",
            extra_compile_args = ["-fopenmp"],
            extra_link_args = ["-fopenmp"],
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
