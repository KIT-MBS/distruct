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
# Last Modified : Fri 22 Mar 2019 06:11:01 PM CET
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

includeDir = nwkDir + "/include/networkit"
includeDirs = [includeDir, "distruct/src/"]

import _NetworKit
libraryDir = os.path.split(_NetworKit.__file__)[0]

libraryDirs = [libraryDir]
libraries =['networkit']

nwkpath = os.path.split(_NetworKit.__file__)[0]

compile_args = ["-fopenmp", "-std=c++11"]
link_args = ["-fopenmp"]

extensions = [
        Extension(
            "_diSTruct",
            sources,
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

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
        name = "distruct",
        version = "0.0.7",
        ext_modules = extensions,
        author = "Oskar Taubert",
        author_email = "oskar.taubert@kit.edu",
        description = "a package to generate 3d molecular structures from distance constraints",
        long_description = long_description,
        long_description_content_type = "text/markdown",
        url = "https://github.com/KIT-MBS/distruct",
        download_url = "https://pypi.python.org/pypi/distruct",
        packages = find_packages(),
        package_data = {},
        keywords = ["biomolecules", "graph drawing"],
        classifiers = [
            "Programming Language :: Python :: 3",
            "Programming Language :: C++",
            "Environment :: Console",
            "Natural Language :: English",
            "License :: OSI Approved :: MIT License",
            "Operating System :: POSIX :: Linux",
            "Topic :: Scientific/Engineering :: Bio-Informatics"
            ],
        install_requires = ["numpy", "cython", "networkit", "biopython", "lxml"],
        zip_safe = False)
