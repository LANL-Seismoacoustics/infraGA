#!/usr/bin/env python
from email.policy import default
import sys

import os
import click

import numpy as np 

from . import plot_on_map as map
from . import multipath_wvfrm as mltwvfrm
from . import topo_extractor as terrain 
from . import ecmwf_extractor as ecmwf
from . import g2s_grid
from . import analyze_atmo


@click.group(context_settings={'help_option_names': ['-h', '--help']})
def main(args=None):
    '''

    infraga - Python interface for various infraGA/GeoAc features
    
    Methods are included here to enable visualization of infraga-sph output on a map, identification of all eigenrays and computation of the combined waveform for a given source-receiver geometry, extraction of the terrain along a propagation path or within a grid, and extraction of usable atmospheric specifications from ECMWF netCDF format files

    -------------------------------------------------------------------------

    Note: the methods here are Python supplements to the C++ implementation of propagation simulation methods that are the core of infraGA/GeoAc.  To view usage of the propagation simulation methods, run one of the following on the command line:

    \b
    \t infraga-2d  
    \t infraga-3d  
    \t infraga-3d-rngdep
    \t infraga-sph
    \t infraga-sph-rngdep

    -------------------------------------------------------------------------
    '''
    pass


@main.command('map-results', short_help="Visualize results on a cartopy map")
@click.option("--arrivals", help="Arrivals file from an infraga-sph simulation", default=None)
@click.option("--ray-paths", help="Ray path file from an infraga-sph simulation", default=None)
@click.option("--plot-option", help="Parameter to visualize for arrivals ('amplitude', 'turning-height', 'celerity', or 'none')", default='amplitude')
@click.option("--figure-name", help="Name of output figure", default="arrivals.png")
@click.option("--rcvrs-file", help="File containing receiver locations (optional)", default=None)
@click.option("--title", help="Title for the figure", default="infraga-sph predictions")
@click.option("--start-time", help="Propagation time [hours] for plotting sub-set of data", default=None, type=float)
@click.option("--end-time", help="Propagation time [hours] for plotting sub-set of data", default=None, type=float)
@click.option("--include-absorption", help="Include Sutherland & Bass losses", default=True)
def run_map(arrivals, ray_paths, plot_option, figure_name, rcvrs_file, title, start_time, end_time, include_absorption):
    '''
    Visualize arrivals or ray paths computed using infraga-sph methods on a Cartopy map

    \b
    Examples:
    \t infraga map-results --arrivals ToyAtmo.arrivals.dat --plot-option amplitude --title 'Toy Atmo arrival amplitudes' --figure-name ToyAtmo.arrivals.png
    \t infraga map-results --arrivals ToyAtmo.arrivals.dat --plot-option celerity --title 'Toy Atmo arrival celerity' --figure-name ToyAtmo.celerities.png
    \t infraga map-results --ray-paths ToyAtmo.raypaths.dat --title 'Toy Atmo ray paths' --figure-name ToyAtmo.raypaths.png

    '''
    if arrivals is None and ray_paths is None:
        click.echo("Mapping requires either arrivals or ray paths file")
    elif arrivals is not None and ray_paths is not None:
        click.echo("Mapping cannot use both arrivals and ray paths (only include one)")
    else:
        if arrivals:
            click.echo("Plotting arrival information in '" + arrivals + "' with option: '" + plot_option + "'...")
        else:
            click.echo("Plotting ray path information in '" + ray_paths + "'...")

        map.run(arrivals, ray_paths, plot_option, figure_name, rcvrs_file=rcvrs_file, title_text=title, time1=start_time, time2=end_time, include_absorp=include_absorption)


@main.command('multi-wvfrm', short_help="Identify eigenrays and compute a combined waveform ")
@click.option("--specification", help="Atmospheric specification", prompt = "Specify atmospheric specification file: ")
@click.option("--geom", help="Geometry('3d' or 'sph')", default='sph', prompt = "Specify geometry ('3d' or 'sph'): ")
@click.option("--src-lat", help="Source latitude (only for 'sph' methods)", default=30.0)
@click.option("--src-lon", help="Source longitude (only for 'sph' methods)", default=-110.0)
@click.option("--src-alt", help="Source altitude", default=0.0)
@click.option("--rcvr-x", help="Receiver x location [km] (relative to source, only for '3d' methods)", default=-250.0)
@click.option("--rcvr-y", help="Receiver y location [km] (relative to source, only for '3d' methods)", default=0.0)
@click.option("--rcvr-lat", help="Source latitude [deg] (only for 'sph' methods)", default=30.0)
@click.option("--rcvr-lon", help="Source latitude [deg] (only for 'sph' methods)", default=-113.0)
@click.option("--z-grnd", help="Ground elevation [km]", default=0.0)
@click.option("--bnc-max", help="Number of reflections to consider", default=0)
@click.option("--incl-min", help="Minimum inclination angle for eigenray search", default=0.5)
@click.option("--incl-max", help="Maximum inclination angle for eigenray search", default=45.0)
@click.option("--incl-step-max", help="Maximum inclination angle step during for eigenray search", default=0.1)
@click.option("--rng-max", help="Maximum range for eigenray search calculation", default=1000.0)
@click.option("--verbose", help="Option to output verbose information about eigenray search", default=True)
@click.option("--wvfrm-ref", help="Reference distance from source for initiating waveform calculation (used with wvfrm-file)", default=1.0)
@click.option("--wvfrm-yld", help="Source yield [kg eq. TNT] (default: 10 tons)", default=10.0e3)
@click.option("--wvfrm-file", help="File containing a starting waveform for reference", default=None)
@click.option("--cpu-cnt", help="CPUs to use in calculation (requires OpenMPI methods to be installed)", default=None)
def run_wvfrm(specification, geom, src_lat, src_lon, src_alt, rcvr_x, rcvr_y, rcvr_lat, rcvr_lon, z_grnd, bnc_max, incl_min, incl_max, incl_step_max, rng_max, verbose, wvfrm_ref, wvfrm_yld, wvfrm_file, cpu_cnt):
    '''
    Identify eigenrays and compute combined waveform by running the -eig_search and -wnl_wvfrm methods in combination

    \b
    Examples:
    \t infraga multi-wvfrm --specification ToyAtmo.met --geom 3d --rcvr-x 175.0 --rcvr-y 75.0 --bnc-max 0 --wvfrm-yld 10.0e3
    \t infraga multi-wvfrm --specification ToyAtmo.met --geom sph --src-lat 30.0 --src-lon -110.0 --rcvr-lat 30.5 --rcvr-lon 115.0 --bnc-max 1 --wvfrm-yld 1.0e3

    '''
    click.echo("Running multi-waveform methods with specification '" + specification + "' and option: '" + geom + "'...")
    mltwvfrm.run(specification, geom, src_lat, src_lon, src_alt, rcvr_x, rcvr_y, rcvr_lat, rcvr_lon, z_grnd, bnc_max, incl_min=incl_min, incl_max=incl_max, incl_step_max=incl_step_max, rng_max=rng_max, verbose_output=verbose, wvfrm_ref=wvfrm_ref, wvfrm_yld=wvfrm_yld, wvfrm_file=wvfrm_file, cpu_cnt=cpu_cnt)


@main.command('extract-terrain', short_help="Extract a line or grid of terrain information")
@click.option("--geom", help="Geometry option ('line', 'pnt2pnt', 'xy-grid' or 'latlon-grid')", prompt="Enter terrain option  ('line', 'pnt2pnt', 'xy-grid' or 'latlon-grid')")
@click.option("--lat1", help="Latitude of first point (starting point for 'pnt2pnt', lower-left corner for grids)", default=30.0)
@click.option("--lon1", help="Longitude of first point (starting point for 'pnt2pnt', lower-left corner for grids)", default=-110.0)
@click.option("--lat2", help="Latitude of second point (end point for 'pnt2pnt', upper-right corner for grids)", default=30.0)
@click.option("--lon2", help="Longitude of second point (end point for 'pnt2pnt', upper-right corner for grids)", default=-114.0)
@click.option("--lat-ref", help="Reference latitude of second point (0.0 for xy-grid option)", default=30.0)
@click.option("--lon-ref", help="Reference longitude of second point (0.0 for xy-grid option)", default=-110.0)
@click.option("--azimuth", help="Azimuth of great circle path for line option", default=-90.0)
@click.option("--range", help="Great circle distance for line option", default=1000.0)
@click.option("--output-file", help="Output file", prompt="Specify output file: ")
@click.option("--show-terrain", help="Visualize terrain results", default=True)
def run_terrain(geom, lat1, lat2, lon1, lon2, lat_ref, lon_ref, azimuth, range, output_file, show_terrain):
    '''
    Extract lines or grids of terrain information from an ETOPO1 file

    \b
    Examples:
    \t infraga extract-terrain --geom line --lat1 40.0 --lon1 -102.5 --azimuth -90.0 --range 750.0 --output-file line_topo.dat
    \t infraga extract-terrain --geom pnt2pnt --lat1 40.0 --lon1 -102.5 --lat2 40.0 --lon2 -110.0 --output-file line_topo.dat
    \t infraga extract-terrain --geom xy-grid --lat1 35.0 --lon1 -110.0 --lat2 45.0 --lon2 -100.0 --lat-ref 40.0 --lon-ref -105.0 --output-file xy_topo.dat
    \t infraga extract-terrain --geom latlon-grid --lat1 35.0 --lon1 -110.0 --lat2 45.0 --lon2 -100.0 --output-file sph_topo.dat

    '''
    terrain.run(geom, lat1, lat2, lon1, lon2, lat_ref, lon_ref, azimuth, range, output_file, show_terrain)



@main.command('extract-ecmwf', short_help="Extract atmospheric information from an ECMWF netCDF file")
@click.option("--ecmwf-file", help="ECMWF netCDF file")
@click.option("--option", help="Extraction option ('single' or 'grid')", prompt="Enter terrain option  ('single' or 'grid')")
@click.option("--lat1", help="Latitude of first point (starting point for 'pnt2pnt', lower-left corner for grids)", default=30.0)
@click.option("--lon1", help="Longitude of first point (starting point for 'pnt2pnt', lower-left corner for grids)", default=-110.0)
@click.option("--lat2", help="Latitude of second point (end point for 'pnt2pnt', upper-right corner for grids)", default=30.0)
@click.option("--lon2", help="Longitude of second point (end point for 'pnt2pnt', upper-right corner for grids)", default=-114.0)
@click.option("--sample_skips", help="Frequency of samples in the grid option (defaults to 1 to keep every node)", default=1)
@click.option("--output-path", help="Output file", prompt="Specify output path: ")

def run_ecmwf(ecmwf_file, option, lat1, lon1, lat2, lon2, output_path, sample_skips):
    '''
    Extract vertical profiles from an ECMWF NetCDF format file.  Able to extract a single file or a range dependent grid.

    \b
    Examples:
    \t infraga extract-ecmwf --option single --ecmwf-file ECMWF.nc --lat1 40.0 --lon1 -100.0 --output-path ECMWF-40.0_-100.0.dat
    \t infraga extract-ecmwf --option grid --ecmwf-file ECMWF.nc --lat1 35.0 --lon1 -110.0 --lat2 45.0 --lon2 -100.0 --output-path grid/ECMWF-

    '''
    ecmwf.run(ecmwf_file, option, lat1, lon1, lat2, lon2, output_path, sample_skips)



@main.command('build-g2s-grid', short_help="Build grid for -rngdep analysis from G2S grid")
@click.option("--g2s-path", help="Path to G2S specifications", prompt="Specify directory containing G2S specifications: ")
@click.option("--output-path", help="Output dir + label", prompt="Specify output directory and label: ")
@click.option("--src-info", help="Source info (lat, lon, time) (optional)", default=None)
@click.option("--celerity-est", help="Celerity estimate [km/s] (optional)", default=0.29)
def run_ecmwf(g2s_path, output_path, src_info, celerity_est):
    '''
    Construct the numbered specifications and grid files needed to run -rngdep methods.
    Assumes file format from https://g2s.ncpa.olemiss.edu/ (e.g., g2stxt_2022011506_-3.0000_-54.0000.dat)

    Inclusion of source info (location and time), an estimated celerity, and profiles at multiple reference times enables construction of a temporally varying grid where a node at a distance, r, from the source has an estimated time delay of, dt = r / cel, and uses the appropriate atmospheric information

    \b
    Examples:
    \t infraga build-g2s-grid --g2s-path g2s_dir/ --output-path grid/g2s_grid
    \t infraga build-g2s-grid --g2s-path g2s_dir/ --output-path grid/g2s_grid --src-info '[-20.56989, -175.379975, 2022-01-15T04:14:45]' --celerity-est 0.29

    '''

    g2s_grid.run(g2s_path, output_path, src_info, celerity_est)



@main.command('visualize-atmo', short_help="Visualize information about an atmospheric specification")
@click.option("--specification", help="Atmospheric specification file")
@click.option("--max-alt", help="Maximum altitude for analysis (default: 120 km)", default=120)
@click.option("--format", help="Atmospheric specification format (default: 'zTuvdp')", default='zTuvdp')
def run_visualize_atmo(specification, format, max_alt):
    '''
    Visualize the sound speed, wind fields, and effective sound speed ratio ducting information for an atmospheric specification


    \b
    Examples:
    \t infraga visualize-atmo --specification examples/ToyAtmo.met
    \t infraga visualize-atmo --specification examples/G2S_example.met

    '''

    analyze_atmo.run(specification, format, max_alt)
    

if __name__ == '__main__':
    main()
