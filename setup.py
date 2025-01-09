# -*- coding: utf-8 -*-

from distutils.core import setup

setup(
    name = "infraga2",
    license='LANL-MIT',
    version = '1.1.0',
    description = "A tool for modeling the propagation of infrasound in the limit of geometric acoustics (Python interface).",
    keywords=['infrasound', 'geophysics', 'seismic', 'array'],
    author = "LANL Seismoacoustics Infrasound (LANL-SA) Team",
    author_email = 'pblom@lanl.gov',

    packages = ['infraga',
                'infraga.cli',
                'infraga.doc',
                'infraga.resources',
                'infraga.src',
                'infraga.src.atmo',
                'infraga.src.geoac',
                'infraga.src.util'],

    entry_points = {'console_scripts':['infraga=infraga.cli.__main__:main']},

    install_requires = ['cartopy',
                        'click',
                        'netCDF4',
                        'numpy',
                        'matplotlib',
                        'pyproj',
                        'scipy',
                        'ipython',
                        'pip',
                        'wget']
)
