#!which python
"""
run_infraga.py

Run the various infraga methods using a simplified Python Click interface

Author: pblom@lanl.gov    
"""

import os 
import sys

import click
import subprocess
import configparser as cnfg

from importlib.util import find_spec

# Set up default configuation
defaults = cnfg.ConfigParser()
defaults.read(find_spec('infraga').submodule_search_locations[0] + '/default.config')
bin_path =  find_spec('infraga').submodule_search_locations[0][:-8] + "/bin/"

def set_param(user_config, section, param, cli_val, format='float'):
    if cli_val is not None:
        # return CLI value if entered
        return cli_val
    elif user_config is not None:
        # check if parameter is in user config and return default if it's not
        try:
            if user_config[section][param] == "None":
                return None
            else:
                if format == 'float':
                    return float(user_config[section][param])
                elif format == 'int':
                    return int(user_config[section][param])
                elif format == 'bool':
                    return user_config[section].getboolean(param)
                else:
                    return user_config[section][param]
        except:
            try:
                if defaults[section][param] == "None":
                    return None
                else:
                    if format == 'float':
                        return float(defaults[section][param])
                    elif format == 'int':
                        return int(defaults[section][param])
                    elif format == 'bool':
                        return defaults[section].getboolean(param)
                    else:
                        return defaults[section][param]
            except:
                return None
    else:
        # use default values if no CLI and no user config
        try:
            if defaults[section][param] == "None":
                return None
            else:
                if format == 'float':
                    return float(defaults[section][param])
                elif format == 'int':
                    return int(defaults[section][param])
                elif format == 'bool':
                    return defaults[section].getboolean(param)
                else:
                    return defaults[section][param]
        except:
            return None


@click.command('prop', short_help="Run infraga-sph methods")
@click.option("--option", help="Simulation option (e.g., prop, eig-search)", default="prop")
@click.option("--atmo-file", help="Atmosphere file", default=None)
def run_2d_prop(option, atmo_file):
    command = bin_path + "infraga-2d"

    subprocess.call(command, shell=True)

    



@click.command('prop', short_help="Run infraga-sph methods")
@click.option("--option", help="Simulation option (e.g., prop, eig-search)", default="prop")
@click.option("--atmo-file", help="Atmosphere file", default=None, prompt="Enter atmosphere file:")
def run_3d_prop(option, atmo_file):
    command = bin_path + "infraga-3d"

    subprocess.call(command, shell=True)


@click.command('prop', short_help="Run infraga-sph methods")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)
@click.option("--atmo-prefix", help="Atmosphere file prefix (range dependent)", default=None)
@click.option("--grid-lats", help="Atmosphere grid latitudes file", default=None)
@click.option("--grid-lons", help="Atmosphere grid longitudes file", default=None)
@click.option("--incl-min", help="Minimum inclination angle (rel. horizontal)", default=None, type=float)
@click.option("--incl-max", help="Maximum inclination angle (rel. horizontal)", default=None, type=float)
@click.option("--incl-step", help="Inclination angle step resolution", default=None, type=float)
@click.option("--inclination", help="Single inclination angle", default=None, type=float)
@click.option("--az-min", help="Minimum azimuth angle (clockwise rel. N)", default=None, type=float)
@click.option("--az-max", help="Maximum azimuth angle (clockwise rel. N)", default=None, type=float)
@click.option("--az-step", help="Azimuth angle step resolution", default=None, type=float)
@click.option("--azimuth", help="Single azimuth angle", default=None, type=float)
@click.option("--bounces", help="Number of ground reflections (bounces) to consider", default=None, type=int)
@click.option("--src-lat", help="Source latitude", default=None, type=float)
@click.option("--src-lon", help="Source lontidue", default=None, type=float)
@click.option("--src-alt", help="Source altitude", default=None, type=float)
@click.option("--write-rays", help="Option to write [...].raypaths.dat output", default=None, type=bool)
@click.option("--write-topo", help="Option to write terrain info under first ray", default=None, type=bool)

@click.option("--freq", help="Frequency for Sutherland-Bass losses", default=None, type=float)
@click.option("--abs-coeff", help="Scaling coefficient for Sutherland-Bass losses", default=None, type=float)
@click.option("--z-grnd", help="Ground elevation (flat ground, km rel. sea level)", default=None, type=float)
@click.option("--write-atmo", help="Option to write atmosphere data (for QC)", default=None, type=float)
@click.option("--output-id", help="User specified output file path", default=None, type=str)
@click.option("--calc-amp", help="Option to turn off transport coefficient calculation", default=None, type=bool)
@click.option("--max-alt", help="Maximum altitude for ray calculation", default=None, type=float)
@click.option("--max-rng", help="Maximum range for ray calculation", default=None, type=float)

@click.option("--min-lat", help="Minimum latitude for user defined bounds", default=None, type=float)
@click.option("--max-lat", help="Maximum latitude for user defined bounds", default=None, type=float)
@click.option("--min-lon", help="Minimum latitude for user defined bounds", default=None, type=float)
@click.option("--max-lon", help="Maximum latitude for user defined bounds", default=None, type=float)

@click.option("--min-ds", help="Minimum step size (near-ground) in RK4 solver", default=None, type=float)
@click.option("--max-ds", help="Maximum step size in RK4 solver", default=None, type=float)
@click.option("--max-s", help="Maximum ray length between bounces", default=None, type=float)

@click.option("--topo-file", help="Terrain file", default=None, type=str)
@click.option("--topo-BL-wind", help="Use terrain corrected boundary layer winds", default=None, type=bool)


@click.option("--cpu-cnt", help="Number of CPUs to use in analysis", default=None, type=int)
def run_sph_prop(config_file, atmo_file, atmo_prefix, grid_lats, grid_lons, incl_min, incl_max, incl_step, inclination, 
                az_min, az_max, az_step, azimuth, bounces, src_lat, src_lon, src_alt, write_rays, write_topo, freq, 
                abs_coeff, z_grnd, write_atmo, output_id, calc_amp, max_alt, max_rng, max_lat, min_lat, max_lon, min_lon,
                min_ds, max_ds, max_s, topo_file, topo_bl_wind, cpu_cnt):
    '''
    Run spherical atmosphere geometry ray tracing analysis (e.g., infraga-sph)

    \b
    Examples:
    \t infraga sph prop --atmo-file ToyAtmo.met
    \t infraga sph prop --atmo-file ToyAtmo.met --azimuth -90.0

    '''

    if config_file:
        click.echo('\n' + "Loading configuration info from: " + config_file)
        if os.path.isfile(config_file):
            user_config = cnfg.ConfigParser()
            user_config.read(config_file)
        else:
            click.echo('\n' + "Invalid configuration file (file not found)")
            return 0
    else:
        user_config = None

    incl_min = set_param(user_config, 'PROP', 'incl_min', incl_min, 'str')
    incl_max = set_param(user_config, 'PROP', 'incl_max', incl_max, 'str')
    incl_step = set_param(user_config, 'PROP', 'incl_step', incl_step, 'str')

    az_min = set_param(user_config, 'PROP', 'az_min', az_min, 'str')
    az_max = set_param(user_config, 'PROP', 'az_max', az_max, 'str')
    az_step = set_param(user_config, 'PROP', 'az_step', az_step, 'str')

    bounces = set_param(user_config, 'PROP', 'bounces', bounces, 'str')


    cpu_cnt = int(cpu_cnt)


    # Build command
    if cpu_cnt < 2:
        command = bin_path + "infraga-sph"
    else:
        command = "mpirun -np " + str(cpu_cnt) + " " + bin_path + "infraga-accel-sph"

    if atmo_prefix is not None:
        command = command + "-rng-dep"

    if atmo_file is not None:
        command = command + " -prop " + atmo_file
    elif atmo_prefix is not None and grid_lats is not None and grid_lons is not None:
        command = command + " -prop " + atmo_prefix + " " + grid_lats + " " + grid_lons
    else:
        click.echo("Simulation requires either an '--atmo-file' or --atmo-prefix' with grid info (--grid-lats and --grid-lons)")
        return 0




    print(command)

    subprocess.call(command, shell=True)
