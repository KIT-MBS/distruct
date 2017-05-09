#!/usr/bin/env python3
#####################################
#
# Filename : setup.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Tue 09 May 2017 13:40:25 CEST
#
# Last Modified : Tue 09 May 2017 13:50:41 CEST
#
#####################################

try:
    from setuptools import setup
except ImportError:
    from distutile.core import setup
    pass

config = {
        'description': 'Maxent-stress optimization of Biomolecular systems: a tool for molecular structure generation from sparse atomic distance information',
        'author': 'Oskar Taubert',
        'url': 'http://www.scc.kit.edu/en/research/8900.php',
        'download_url': 'https://github.com/KIT-MBS/MOBi',
        'author_email': 'oskar.taubert@partner.kit.edu',
        'version': '0.1',
        'install_requires': ['nose', 'networkit'],
        'packages': ['MOBi'],
        'scripts': [],
        'name': 'MOBi'}

setup(**config)
