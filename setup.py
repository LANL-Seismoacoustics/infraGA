# -*- coding: utf-8 -*-
try:
    import setuptools
except ImportError:
    pass

from distutils.core import setup

setup(
    name = "infraga",
    license='LANL-MIT',
    version = '1.0.3',
    description = "A tool for modeling the propagation of infrasound in the limit of geometric acosutics (Python interface).",
    keywords=['infrasound', 'geophysics', 'seismic', 'array'],
    author = "LANL Seismoacoustics Infrasound (LANL-SA) Team",
    author_email = 'pblom@lanl.gov',

    packages = ['infraga'],

    entry_points = {'console_scripts':['infraga=infraga.__main__:main']},

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
