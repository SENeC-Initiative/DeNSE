#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os, errno
from setuptools import setup, find_packages


# create directory
directory = 'PyNCulture/'
try:
    os.makedirs(directory)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise


# move important 
move = (
    '__init__.py',
    'LICENSE',
    'dxf_import',
    'examples',
    'backup_shape.py',
    'dxftools.py',
    'geom_utils.py',
    'plot.py',
    'pync_log.py',
    'shape.py',
    'svgtools.py',
)


for fname in move:
    os.rename(fname, directory + fname)


from PyNCulture import __version__


try:
    # install
    setup(
        name = 'PyNCulture',
        version = __version__,
        description = 'Python module to describe neuronal cultures as '+\
                      ' complex shapes.',
        package_dir = {'': '.'},
        packages = find_packages('.'),

        # Requirements
        install_requires = ['numpy', 'scipy>=0.11', 'matplotlib'],
        extras_require = {
            'dxfgrabber': 'dxfgrabber',
            'shapely': 'shapely',
            'svgtools': 'svgtools'
        },

        # Metadata
        url = 'https://github.com/Silmathoron/PyNCulture',
        author = 'Tanguy Fardet',
        author_email = 'tanguy.fardet@univ-paris-diderot.fr',
        license = 'GPL3',
        keywords = 'neuronal cultures geometry'
    )
finally:
    for fname in move:
        os.rename(directory + fname, fname)
