#!/usr/bin/env python
from email.policy import default
import sys

import os
import click

import numpy as np 

from . import run_infraga
from . import multipath_wvfrm as mltwvfrm

from . import plotting
from . import utils as infraga_utils


@click.group(context_settings={'help_option_names': ['-h', '--help']})
def main(args=None):
    '''

    infraga - Python interface for using the infraGA/GeoAc software tools
    '''
    pass


@click.group('2d', short_help="Run 2d (effective sound speed) ray tracing", context_settings={'help_option_names': ['-h', '--help']})
def run_2d():
    '''
    infraga 2d - run ray tracing in an azimuthal plane (r vs. z) using the effective sound speed 
    
    '''
    pass 


@click.group('3d', short_help="Run 3d (moving medium) ray tracing", context_settings={'help_option_names': ['-h', '--help']})
def run_3d():
    '''
    infraga 3d - run Cartesian (x, y, z) geometry ray tracing assuming a moving medium
    
    '''
    pass 


@click.group('sph', short_help="Run spherical layer (moving medium) ray tracing", context_settings={'help_option_names': ['-h', '--help']})
def run_sph():
    '''
    infraga sph - run spherical atmosphere layer geometry (latitude, longitude, altitude) ray tracing assuming a moving medium
    
    '''
    pass 


@click.group('plot', short_help="Various visualization functions", context_settings={'help_option_names': ['-h', '--help']})
def plot():
    '''
    infraga plot - visualization functions for atmospheric data and infraga results
    
    '''
    pass 


@click.group('utils', short_help="Various utility functions", context_settings={'help_option_names': ['-h', '--help']})
def utils():
    '''
    infraga utils - utility functions for infraga usage
    
    '''
    pass 


main.add_command(run_2d)
main.add_command(run_3d)
main.add_command(run_sph)
main.add_command(plot)
main.add_command(utils)

run_2d.add_command(run_infraga.run_2d_prop)
run_2d.add_command(run_infraga.run_2d_wvfrm)

run_3d.add_command(run_infraga.run_3d_prop)
run_3d.add_command(run_infraga.run_3d_eig)
run_3d.add_command(run_infraga.run_3d_wvfrm)
run_3d.add_command(run_infraga.run_3d_eig_wvfrm)

run_sph.add_command(run_infraga.run_sph_prop)
run_sph.add_command(run_infraga.run_sph_eig)
run_sph.add_command(run_infraga.run_sph_wvfrm)
run_sph.add_command(run_infraga.run_sph_eig_wvfrm)

plot.add_command(plotting.plot_atmo)
plot.add_command(plotting.plot_azimuthal)
plot.add_command(plotting.plot_eigenrays)
plot.add_command(plotting.plot_eig_wvfrm)
plot.add_command(plotting.plot_map)
    
utils.add_command(infraga_utils.build_g2s_grid)
utils.add_command(infraga_utils.extract_ecmwf)
utils.add_command(infraga_utils.extract_terrain)
utils.add_command(infraga_utils.nearby_arrivals)

if __name__ == '__main__':
    main()
