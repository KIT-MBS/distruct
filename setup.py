#!/usr/bin/env python3
#####################################
#
# Filename : setup.py
#
# Projectname : diSTruct
#
# Author : Oskar Taubert
#
# Creation Date : Thu 28 Jun 2018 12:50:34 PM CEST
#
# Last Modified : Mon 30 Jul 2018 04:17:49 PM CEST
#
#####################################

# TODO look at what kind of structure/files pypi expects
# TODO install requirements

from setuptools import setup, Extension, find_packages
from Cython.Build import cythonize

import os


USE_CYTHON = True
fext = '.pyx' if USE_CYTHON else '.cpp'

# NOTE we need the cpp headers from networkit, which they do not redistribute
# NOTE so we do a shallow clone from their repo
# TODO this will hopefully be fixed at some point
def networkit_get_include(dir):
    import subprocess as sp
    # TODO get the downloadurl from the networkit package
    nwkurl = "https://github.com/kit-parco/networkit.git"
    sp.call(['git', 'clone', '--depth=1', nwkurl, dir])
    return

sources = ['distruct/src/BioMaxentStress.cpp', "distruct/_diSTruct" + fext]

nwkDir = "shallownwk"
print("downloading networkit header files...")
if not os.path.isdir(nwkDir):
    networkit_get_include(nwkDir)
    pass

includeDir = nwkDir + "/networkit"
includeDirs = [includeDir, "distruct/src/"]

import _NetworKit
libraryDir, library = tuple(os.path.split(_NetworKit.__file__))

libraryDirs = [libraryDir]
libraries = [':' + library]  # NOTE the : is to tell the linker to use the actual filename

import _NetworKit
nwkpath = os.path.split(_NetworKit.__file__)[0]

compile_args = ["-fopenmp"]
link_args = ["-fopenmp"]

extensions = [
        Extension(
            "_diSTruct",
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
        name = "distruct",
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
