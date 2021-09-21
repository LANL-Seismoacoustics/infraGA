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

    packages = ['scripts'],

    entry_points = {'console_scripts':['infraga=scripts.__main__.__main__']},

    install_requires = ['click',
                        'numpy',
                        'scipy',
                        'ipython',
                        'matplotlib']
)
