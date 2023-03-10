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

from scipy.interpolate import interp1d 

from importlib.util import find_spec

import numpy as np

bin_path =  find_spec('infraga').submodule_search_locations[0][:-8] + "/bin/"

###################################
#                                 #
#   Parameter ingestion and use   #
#                                 #
###################################
def define_param(user_config, section, param, cli_val, format='str'):
    if cli_val is not None:
        # return CLI value if entered
        return cli_val
    elif user_config is not None:
        # check if parameter is in user config
        try:
            if user_config[section][param] == "None":
                return None
            else:
                if format == 'bool':
                    return user_config[section].getboolean(param)
                else:
                    return user_config[section][param]
        except:
            return None
    else:
        # If not in CLI or user config, set to None and use default value in C++ implementation
        return None

def set_param(command, param, param_label):
    if param is not None:
        command = command + " " + param_label + "=" + param
    return command 

########################
#                      #
#    Kinney & Graham   #
#     Scaling Laws     #
#                      #
########################

def kg_op(W, r, p_amb=101.325, T_amb=288.15, type="chemical"):
    """
        Kinney & Graham scaling law peak overpressure model
                
        Parameters
        ----------
        W : float
            Explosive yield of the source [kg eq. TNT]
        r : float
            Propagation distance [km]
        p_amb : float
            Ambient atmospheric pressure [kPa]
        T_amb : float
            Ambient atmospheric temperature [deg K]
        type : string
            Type of explosion modeled, options are "chemical" or "nuclear"
        
        Returns
        -------
        p0 : float
            Peak overpressure [Pa]    
    """
    
    fd = (p_amb / 101.325)**(1.0 / 3.0) * (T_amb / 288.15)**(1.0 / 3.0)
    sc_rng = fd / W**(1.0 / 3.0) * r * 1000.0
    
    if type=="chemical":
        term1 = 1.0 + (sc_rng / 4.5)**2
        term2 = np.sqrt(1.0 + (sc_rng / 0.048)**2)
        term3 = np.sqrt(1.0 + (sc_rng / 0.32)**2)
        term4 = np.sqrt(1.0 + (sc_rng / 1.35)**2)
        
        result = 808.0 * term1 / (term2 * term3 * term4)
    else:
        term1 = (1.0 + sc_rng / 800.0)
        term2 = np.sqrt(1.0 + (sc_rng / 87.0)**2)
        
        result = 3.2e6 / sc_rng**3 * term1 * term2

    return p_amb * 1.0e3 * result


def kg_ppd(W, r, p_amb=101.325, T_amb=288.15, type="chemical"):
    """
        Kinney & Graham scaling law positive phase duration model
                
        Parameters
        ----------
        W : float
            Explosive yield of the source [kg eq. TNT]
        r : float
            Propagation distance [km]
        p_amb : float
            Ambient atmospheric pressure [kPa]
        T_amb : float
            Ambient atmospheric temperature [deg K]
        type : string
            Type of explosion modeled, options are "chemical" or "nuclear"
            
        Returns
        -------
        t0 : float
            Positive phase duration [s]
    """

    fd = (p_amb / 101.325)**(1.0 / 3.0) * (T_amb / 288.15)**(1.0 / 3.0)
    sc_rng = fd / W**(1.0 / 3.0) * r * 1000.0
    
    if type=="chemical":
        term1 = 1.0 + (sc_rng / 0.54)**10
        term2 = 1.0 + (sc_rng / 0.02)**3
        term3 = 1.0 + (sc_rng / 0.74)**6
        term4 = np.sqrt(1.0 + (sc_rng / 6.9)**2)
        
        result = 980.0 * term1 / (term2 * term3 * term4)
    else:
        term1 = np.sqrt(1.0 + (sc_rng / 100.0)**3)
        term2 = np.sqrt(1.0 + (sc_rng / 40.0))
        term3 = (1.0 + (sc_rng / 285.0)**5)**(1.0 / 6.0)
        term4 = (1.0 + (sc_rng / 50000.0))**(1.0 / 6.0)
        
        result = 180.0 * term1 / (term2 * term3 * term4)

    return result * W**(1.0 / 3.0) / 1e3



###################################
##                               ##
##   Python CLI for infraga-2d   ##
##                               ##
###################################
@click.command('prop', short_help="Run a point source propagation simulation")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)

@click.option("--incl-min", help="Minimum inclination angle (rel. horizontal)", default=None)
@click.option("--incl-max", help="Maximum inclination angle (rel. horizontal)", default=None)
@click.option("--incl-step", help="Inclination angle step resolution", default=None)
@click.option("--inclination", help="Single inclination angle", default=None)
@click.option("--azimuth", help="Single azimuth angle", default=None)
@click.option("--bounces", help="Number of ground reflections (bounces) to consider", default=None)
@click.option("--src-alt", help="Source altitude", default=None)

@click.option("--freq", help="Frequency for Sutherland-Bass losses", default=None)
@click.option("--abs-coeff", help="Scaling coefficient for Sutherland-Bass losses", default=None)
@click.option("--z-grnd", help="Ground elevation (flat ground, km rel. sea level)", default=None)
@click.option("--write-atmo", help="Option to write atmosphere data (for QC)", default=None)
@click.option("--output-id", help="User specified output file path", default=None)
@click.option("--calc-amp", help="Option to turn off transport coefficient calculation", default=None, type=bool)
@click.option("--max-alt", help="Maximum altitude for ray calculation", default=None)
@click.option("--max-rng", help="Maximum range for ray calculation", default=None)

@click.option("--min-ds", help="Minimum step size (near-ground) in RK4 solver", default=None)
@click.option("--max-ds", help="Maximum step size in RK4 solver", default=None)
@click.option("--max-s", help="Maximum ray length between bounces", default=None)

@click.option("--topo-file", help="Terrain file", default=None)
@click.option("--topo-BL-wind", help="Use terrain corrected boundary layer winds", default=None, type=bool)
def run_2d_prop(config_file, atmo_file, incl_min, incl_max, incl_step, inclination, azimuth, bounces, src_alt,
                freq, abs_coeff, z_grnd, write_atmo, output_id, calc_amp, max_alt, max_rng, min_ds, max_ds, max_s,
                topo_file, topo_bl_wind):
    '''
    Run 2D geometry ray tracing analysis for a point source using the effective sound speed approximation

    \b
    Examples:
    \t infraga 2d prop --atmo-file ToyAtmo.met
    \t infraga 2d prop --atmo-file ToyAtmo.met --config-file example.cnfg

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

    # Set propagation specific parameters
    incl_min = define_param(user_config, 'PROP', 'incl_min', incl_min)
    incl_max = define_param(user_config, 'PROP', 'incl_max', incl_max)
    incl_step = define_param(user_config, 'PROP', 'incl_step', incl_step)

    inclination = define_param(user_config, 'PROP', 'inclination', inclination)
    azimuth = define_param(user_config, 'PROP', 'azimuth', azimuth)

    bounces = define_param(user_config, 'PROP', 'bounces', bounces)

    src_alt = define_param(user_config, 'PROP', 'src_alt', src_alt)

    # Set general parameters
    freq = define_param(user_config, 'GENERAL', 'freq', freq)
    abs_coeff = define_param(user_config, 'GENERAL', 'abs_coeff', abs_coeff)

    write_atmo = define_param(user_config, 'GENERAL', 'write_atmo', write_atmo)
    calc_amp = define_param(user_config, 'GENERAL', 'calc_amp', calc_amp)

    output_id = define_param(user_config, 'GENERAL', 'output_id', output_id)

    max_alt = define_param(user_config, 'GENERAL', 'max_alt', max_alt)
    max_rng = define_param(user_config, 'GENERAL', 'max_rng', max_rng)

    min_ds = define_param(user_config, 'GENERAL', 'min_ds', min_ds)
    max_ds = define_param(user_config, 'GENERAL', 'max_ds', max_ds)
    max_s = define_param(user_config, 'GENERAL', 'max_s', max_s)

    z_grnd = define_param(user_config, 'GENERAL', 'z_grnd', z_grnd)
    topo_file = define_param(user_config, 'GENERAL', 'topo_file', topo_file)
    topo_bl_wind = define_param(user_config, 'GENERAL', 'topo_bl_wind', topo_bl_wind)

    # Build run command
    if atmo_file is not None:
        command = bin_path + "infraga-2d -prop " + atmo_file
    else:
        click.echo("Simulation requires an '--atmo-file'")
        return 0

    # Set parameters
    command = set_param(command, incl_min, "incl_min")
    command = set_param(command, incl_max, "incl_max")
    command = set_param(command, incl_step, "incl_step")

    command = set_param(command, inclination, "inclination")
    command = set_param(command, azimuth, "azimuth")

    command = set_param(command, bounces, "bounces")

    command = set_param(command, src_alt, "src_alt")


    command = set_param(command, freq, "freq")
    command = set_param(command, abs_coeff, "abs_coeff")

    if write_atmo is not None:
        command = set_param(command, str(write_atmo), "write_atmo")

    if calc_amp is not None:
        command = set_param(command, str(calc_amp), "calc_amp")

    command = set_param(command, output_id, "output_id")

    command = set_param(command, max_alt, "max_alt")
    command = set_param(command, max_rng, "max_rng")

    command = set_param(command, min_ds, "min_ds")
    command = set_param(command, max_ds, "max_ds")
    command = set_param(command, max_s, "max_s")

    command = set_param(command, z_grnd, "z_grnd")
    command = set_param(command, topo_file, "topo_file")
    if topo_file is not None:
        if topo_bl_wind is not None:
            command = set_param(command, str(topo_bl_wind), "topo_bl_wind")

    # print(command)
    subprocess.call(command, shell=True)



@click.command('wnl-wvfrm', short_help="Run weakly non-linear waveform calculation")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)

@click.option("--inclination", help="Inclination angle (rel. horizontal)", default=None)
@click.option("--azimuth", help="Azimuth angle (clockwise rel. N)", default=None)
@click.option("--bounces", help="Number of ground reflections (bounces) to consider", default=None)

@click.option("--src-alt", help="Source altitude", default=None)

@click.option("--write-ray", help="Option to write ray path info", default=None, type=bool)

@click.option("--wvfrm-file", help="File containing reference waveform", default=None)
@click.option("--wvfrm-opt", help="Waveform option ('impulse', 'Uwave', 'Nwave')", default=None)

@click.option("--wvfrm-p0", help="Waveform peak overpressure [Pa]", default=None)
@click.option("--wvfrm-t0", help="Waveform positive phase duration [s]", default=None)
@click.option("--wvfrm-alpha", help="Waveform shaping parameter", default=None)
@click.option("--wvfrm-ref", help="Waveform reference distance", default=None)
@click.option("--wvfrm-out-step", help="Waveform output step along ray path", default=None)
@click.option("--wvfrm-yield", help="Eq. TNT yield (kg) for an explosive source", default=None)

@click.option("--wvfrm-ds", help="Burgers equation solver resolution [km]", default=None)
@click.option("--wvfrm-len", help="Waveform sample count (default = 2e13)", default=None)

@click.option("--freq", help="Frequency for Sutherland-Bass losses", default=None)
@click.option("--abs-coeff", help="Scaling coefficient for Sutherland-Bass losses", default=None)
@click.option("--z-grnd", help="Ground elevation (flat ground, km rel. sea level)", default=None)
@click.option("--write-atmo", help="Option to write atmosphere data (for QC)", default=None)
@click.option("--output-id", help="User specified output file path", default=None)
@click.option("--max-alt", help="Maximum altitude for ray calculation", default=None)
@click.option("--max-rng", help="Maximum range for ray calculation", default=None)

@click.option("--min-ds", help="Minimum step size (near-ground) in RK4 solver", default=None)
@click.option("--max-ds", help="Maximum step size in RK4 solver", default=None)
@click.option("--max-s", help="Maximum ray length between bounces", default=None)

@click.option("--topo-file", help="Terrain file", default=None)
@click.option("--topo-BL-wind", help="Use terrain corrected boundary layer winds", default=None, type=bool)
def run_2d_wvfrm(config_file, atmo_file, inclination, azimuth, bounces, src_alt, write_ray, wvfrm_file, wvfrm_opt,
                wvfrm_p0, wvfrm_t0, wvfrm_alpha, wvfrm_ref, wvfrm_out_step, wvfrm_yield, wvfrm_ds, wvfrm_len, 
                freq, abs_coeff, z_grnd, write_atmo, output_id, max_alt, max_rng, min_ds, max_ds, max_s, topo_file, 
                topo_bl_wind):
    '''
    Run weakly non-linear waveform analysis along a 2D ray path using the effective sound speed approximation.

    \b
    Examples:
    \t infraga 2d wnl-wvfrm --atmo-file ToyAtmo.met --inclination 12.0 --azimuth -90.0 --wvfrm-p0 500.0
    \t infraga 2d wnl-wvfrm --atmo-file ToyAtmo.met --config-file example.cnfg
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

    # Set propagation specific parameters
    inclination = define_param(user_config, 'WAVEFORM', 'inclination', inclination)
    azimuth = define_param(user_config, 'WAVEFORM', 'azimuth', azimuth)
    bounces = define_param(user_config, 'WAVEFORM', 'bounces', bounces)

    src_alt = define_param(user_config, 'WAVEFORM', 'src_alt', src_alt)

    write_ray = define_param(user_config, 'WAVEFORM', 'write_ray', write_ray)

    wvfrm_file = define_param(user_config, 'WAVEFORM', 'wvfrm_file', wvfrm_file)
    wvfrm_opt = define_param(user_config, 'WAVEFORM', 'wvfrm_opt', wvfrm_opt)
    wvfrm_p0 = define_param(user_config, 'WAVEFORM', 'wvfrm_p0', wvfrm_p0)
    wvfrm_t0 = define_param(user_config, 'WAVEFORM', 'wvfrm_t0', wvfrm_t0)
    wvfrm_alpha = define_param(user_config, 'WAVEFORM', 'wvfrm_alpha', wvfrm_alpha)
    wvfrm_ref = define_param(user_config, 'WAVEFORM', 'wvfrm_ref', wvfrm_ref)
    wvfrm_out_step = define_param(user_config, 'WAVEFORM', 'wvfrm_out_step', wvfrm_out_step)
    wvfrm_yield = define_param(user_config, 'WAVEFORM', 'wvfrm_yield', wvfrm_yield)

    wvfrm_ds = define_param(user_config, 'WAVEFORM', 'wvfrm_ds', wvfrm_ds)
    wvfrm_len = define_param(user_config, 'WAVEFORM', 'wvfrm_len', wvfrm_len)

    # Set general parameters
    freq = define_param(user_config, 'GENERAL', 'freq', freq)
    abs_coeff = define_param(user_config, 'GENERAL', 'abs_coeff', abs_coeff)

    write_atmo = define_param(user_config, 'GENERAL', 'write_atmo', write_atmo)
    output_id = define_param(user_config, 'GENERAL', 'output_id', output_id)

    max_alt = define_param(user_config, 'GENERAL', 'max_alt', max_alt)
    max_rng = define_param(user_config, 'GENERAL', 'max_rng', max_rng)

    min_ds = define_param(user_config, 'GENERAL', 'min_ds', min_ds)
    max_ds = define_param(user_config, 'GENERAL', 'max_ds', max_ds)
    max_s = define_param(user_config, 'GENERAL', 'max_s', max_s)

    z_grnd = define_param(user_config, 'GENERAL', 'z_grnd', z_grnd)
    topo_file = define_param(user_config, 'GENERAL', 'topo_file', topo_file)
    topo_bl_wind = define_param(user_config, 'GENERAL', 'topo_bl_wind', topo_bl_wind)

    # Build run command
    if atmo_file is not None:
        command = bin_path + "infraga-2d -wnl_wvfrm " + atmo_file
    else:
        click.echo("Simulation requires an '--atmo-file'")
        return 0

    # Set parameters
    command = set_param(command, inclination, "inclination")
    command = set_param(command, azimuth, "azimuth")
    command = set_param(command, bounces, "bounces")

    command = set_param(command, src_alt, "src_alt")

    if write_ray is not None:
        command = set_param(command, str(write_ray), "write_ray")

    command = set_param(command, wvfrm_file, "wvfrm_file")
    if wvfrm_yield is not None:
        # use atmosphere to define ambient pressure (1 mbar = 0.1 kPa)
        atmo = np.loadtxt(atmo_file)

        if src_alt is not None:
            src_ht = float(src_alt)
        elif z_grnd is not None:
            src_ht = float(z_grnd)
        else:
            src_ht = 0.0

        p_ambient = atmo[np.argmin(abs(atmo[:, 0] - src_ht)), 5] * 0.1
        T_ambient = atmo[np.argmin(abs(atmo[:, 0] - src_ht)), 1]

        wvfrm_opt = 'impulse'
        wvfrm_ref = str(0.035 * float(wvfrm_yield)**(1.0 / 3.0))
        wvfrm_p0 = str(kg_op(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient))
        wvfrm_t0 = str(kg_ppd(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient))
        wvfrm_alpha = '0.01'

    command = set_param(command, wvfrm_opt, "wvfrm_opt")
    command = set_param(command, wvfrm_p0, "wvfrm_p0")
    command = set_param(command, wvfrm_t0, "wvfrm_t0")
    command = set_param(command, wvfrm_alpha, "wvfrm_alpha")
    command = set_param(command, wvfrm_ref, "wvfrm_ref")
    command = set_param(command, wvfrm_out_step, "wvfrm_out_step")

    command = set_param(command, wvfrm_ds, "wvfrm_ds")
    command = set_param(command, wvfrm_len, "wvfrm_len")

    command = set_param(command, freq, "freq")
    command = set_param(command, abs_coeff, "abs_coeff")

    if write_atmo is not None:
        command = set_param(command, str(write_atmo), "write_atmo")

    command = set_param(command, output_id, "output_id")

    command = set_param(command, max_alt, "max_alt")
    command = set_param(command, max_rng, "max_rng")

    command = set_param(command, min_ds, "min_ds")
    command = set_param(command, max_ds, "max_ds")
    command = set_param(command, max_s, "max_s")

    command = set_param(command, z_grnd, "z_grnd")
    command = set_param(command, topo_file, "topo_file")
    if topo_file is not None:
        if topo_bl_wind is not None:
            command = set_param(command, str(topo_bl_wind), "topo_bl_wind")

    # print(command)
    subprocess.call(command, shell=True)


###################################
##                               ##
##   Python CLI for infraga-3d   ##
##                               ##
###################################
@click.command('prop', short_help="Run a point source propagation simulation")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)
@click.option("--atmo-prefix", help="Atmosphere file prefix (range dependent)", default=None)
@click.option("--grid-lats", help="Atmosphere grid latitudes file", default=None)
@click.option("--grid-lons", help="Atmosphere grid longitudes file", default=None)

@click.option("--incl-min", help="Minimum inclination angle (rel. horizontal)", default=None)
@click.option("--incl-max", help="Maximum inclination angle (rel. horizontal)", default=None)
@click.option("--incl-step", help="Inclination angle step resolution", default=None)
@click.option("--inclination", help="Single inclination angle", default=None)
@click.option("--az-min", help="Minimum azimuth angle (clockwise rel. N)", default=None)
@click.option("--az-max", help="Maximum azimuth angle (clockwise rel. N)", default=None)
@click.option("--az-step", help="Azimuth angle step resolution", default=None)
@click.option("--azimuth", help="Single azimuth angle", default=None)
@click.option("--bounces", help="Number of ground reflections (bounces) to consider", default=None)
@click.option("--src-x", help="Source location (E/W)", default=None)
@click.option("--src-y", help="Source location (N/S)", default=None)
@click.option("--src-alt", help="Source altitude", default=None)
@click.option("--write-rays", help="Option to write [...].raypaths.dat output", default=None, type=bool)
@click.option("--write-topo", help="Option to write terrain info under first ray", default=None, type=bool)

@click.option("--freq", help="Frequency for Sutherland-Bass losses", default=None)
@click.option("--abs-coeff", help="Scaling coefficient for Sutherland-Bass losses", default=None)
@click.option("--z-grnd", help="Ground elevation (flat ground, km rel. sea level)", default=None)
@click.option("--write-atmo", help="Option to write atmosphere data (for QC)", default=None)
@click.option("--output-id", help="User specified output file path", default=None)
@click.option("--calc-amp", help="Option to turn off transport coefficient calculation", default=None, type=bool)
@click.option("--max-alt", help="Maximum altitude for ray calculation", default=None)
@click.option("--max-rng", help="Maximum range for ray calculation", default=None)

@click.option("--min-x", help="Minimum E/W offset for user defined bounds", default=None)
@click.option("--max-x", help="Maximum E/W offset for user defined bounds", default=None)
@click.option("--min-y", help="Minimum N/S offset for user defined bounds", default=None)
@click.option("--max-y", help="Maximum N/S offset for user defined bounds", default=None)

@click.option("--min-ds", help="Minimum step size (near-ground) in RK4 solver", default=None)
@click.option("--max-ds", help="Maximum step size in RK4 solver", default=None)
@click.option("--max-s", help="Maximum ray length between bounces", default=None)

@click.option("--topo-file", help="Terrain file", default=None)
@click.option("--topo-BL-wind", help="Use terrain corrected boundary layer winds", default=None, type=bool)
@click.option("--cpu-cnt", help="Number of CPUs to use in analysis", default=1)
def run_3d_prop(config_file, atmo_file, atmo_prefix, grid_lats, grid_lons, incl_min, incl_max, incl_step, inclination, 
                az_min, az_max, az_step, azimuth, bounces, src_x, src_y, src_alt, write_rays, write_topo, freq, 
                abs_coeff, z_grnd, write_atmo, output_id, calc_amp, max_alt, max_rng, min_x, max_x, min_y, max_y, 
                min_ds, max_ds, max_s, topo_file, topo_bl_wind, cpu_cnt):
    '''
    Run 3D geometry ray tracing analysis for a point source

    \b
    Examples:
    \t infraga 3d prop --atmo-file ToyAtmo.met
    \t infraga 3d prop --atmo-file ToyAtmo.met --config-file example.cnfg

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

    # Set propagation specific parameters
    incl_min = define_param(user_config, 'PROP', 'incl_min', incl_min)
    incl_max = define_param(user_config, 'PROP', 'incl_max', incl_max)
    incl_step = define_param(user_config, 'PROP', 'incl_step', incl_step)

    az_min = define_param(user_config, 'PROP', 'az_min', az_min)
    az_max = define_param(user_config, 'PROP', 'az_max', az_max)
    az_step = define_param(user_config, 'PROP', 'az_step', az_step)

    inclination = define_param(user_config, 'PROP', 'inclination', inclination)
    azimuth = define_param(user_config, 'PROP', 'azimuth', azimuth)

    bounces = define_param(user_config, 'PROP', 'bounces', bounces)

    src_x = define_param(user_config, 'PROP', 'src_x', src_x)
    src_y = define_param(user_config, 'PROP', 'src_y', src_y)
    src_alt = define_param(user_config, 'PROP', 'src_alt', src_alt)

    write_rays = define_param(user_config, 'PROP', 'write_rays', write_rays)
    write_topo = define_param(user_config, 'PROP', 'write_topo', write_topo)

    # Set general parameters
    freq = define_param(user_config, 'GENERAL', 'freq', freq)
    abs_coeff = define_param(user_config, 'GENERAL', 'abs_coeff', abs_coeff)

    write_atmo = define_param(user_config, 'GENERAL', 'write_atmo', write_atmo)
    calc_amp = define_param(user_config, 'GENERAL', 'calc_amp', calc_amp)

    output_id = define_param(user_config, 'GENERAL', 'output_id', output_id)

    max_alt = define_param(user_config, 'GENERAL', 'max_alt', max_alt)
    max_rng = define_param(user_config, 'GENERAL', 'max_rng', max_rng)

    min_x = define_param(user_config, 'GENERAL', 'min_x', min_x)
    max_x = define_param(user_config, 'GENERAL', 'max_x', max_x)
    min_y = define_param(user_config, 'GENERAL', 'min_y', min_y)
    max_y = define_param(user_config, 'GENERAL', 'max_y', max_y)

    min_ds = define_param(user_config, 'GENERAL', 'min_ds', min_ds)
    max_ds = define_param(user_config, 'GENERAL', 'max_ds', max_ds)
    max_s = define_param(user_config, 'GENERAL', 'max_s', max_s)

    z_grnd = define_param(user_config, 'GENERAL', 'z_grnd', z_grnd)
    topo_file = define_param(user_config, 'GENERAL', 'topo_file', topo_file)
    topo_bl_wind = define_param(user_config, 'GENERAL', 'topo_bl_wind', topo_bl_wind)

    cpu_cnt = define_param(user_config, 'GENERAL', 'cpu_cnt', cpu_cnt)
    cpu_cnt = int(cpu_cnt)

    # Build run command
    if cpu_cnt < 2:
        command = bin_path + "infraga-3d"
    else:
        command = "mpirun -np " + str(cpu_cnt) + " " + bin_path + "infraga-accel-3d"

    if atmo_prefix is not None:
        command = command + "-rngdep"

    if atmo_file is not None:
        command = command + " -prop " + atmo_file
    elif atmo_prefix is not None and grid_lats is not None and grid_lons is not None:
        command = command + " -prop " + atmo_prefix + " " + grid_lats + " " + grid_lons
    else:
        click.echo("Simulation requires either an '--atmo-file' or --atmo-prefix' with grid info (--grid-lats and --grid-lons)")
        return 0

    # Set parameters
    command = set_param(command, incl_min, "incl_min")
    command = set_param(command, incl_max, "incl_max")
    command = set_param(command, incl_step, "incl_step")

    command = set_param(command, az_min, "az_min")
    command = set_param(command, az_max, "az_max")
    command = set_param(command, az_step, "az_step")

    command = set_param(command, inclination, "inclination")
    command = set_param(command, azimuth, "azimuth")

    command = set_param(command, bounces, "bounces")

    command = set_param(command, src_x, "src_x")
    command = set_param(command, src_y, "src_y")
    command = set_param(command, src_alt, "src_alt")

    command = set_param(command, az_step, "az_step")

    if write_rays is not None:
        command = set_param(command, str(write_rays), "write_rays")

    command = set_param(command, freq, "freq")
    command = set_param(command, abs_coeff, "abs_coeff")

    if write_atmo is not None:
        command = set_param(command, str(write_atmo), "write_atmo")

    if calc_amp is not None:
        command = set_param(command, str(calc_amp), "calc_amp")

    command = set_param(command, output_id, "output_id")

    command = set_param(command, max_alt, "max_alt")
    command = set_param(command, max_rng, "max_rng")

    command = set_param(command, min_x, "min_x")
    command = set_param(command, max_x, "max_x")
    command = set_param(command, min_y, "min_y")
    command = set_param(command, max_y, "max_y")

    command = set_param(command, min_ds, "min_ds")
    command = set_param(command, max_ds, "max_ds")
    command = set_param(command, max_s, "max_s")

    command = set_param(command, z_grnd, "z_grnd")
    command = set_param(command, topo_file, "topo_file")
    if topo_file is not None:
        if write_topo is not None:
            command = set_param(command, str(write_topo), "write_topo")
        if topo_bl_wind is not None:
            command = set_param(command, str(topo_bl_wind), "topo_bl_wind")

    # print(command)
    subprocess.call(command, shell=True)


@click.command('eigenray', short_help="Run eigenray search for specific source-receiver")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)
@click.option("--atmo-prefix", help="Atmosphere file prefix (range dependent)", default=None)
@click.option("--grid-lats", help="Atmosphere grid latitudes file", default=None)
@click.option("--grid-lons", help="Atmosphere grid longitudes file", default=None)

@click.option("--incl-min", help="Minimum inclination angle (rel. horizontal)", default=None)
@click.option("--incl-max", help="Maximum inclination angle (rel. horizontal)", default=None)

@click.option("--bnc-min", help="Minimum number of ground reflections (bounces)", default=None)
@click.option("--bnc-max", help="Maximum number of ground reflections (bounces)", default=None)
@click.option("--bounces", help="Number of ground reflections (bounces) to consider", default=None)

@click.option("--src-x", help="Source location (E/W)", default=None)
@click.option("--src-y", help="Source location (N/S)", default=None)
@click.option("--src-alt", help="Source altitude", default=None)
@click.option("--rcvr-x", help="Receiver latitude", default=None)
@click.option("--rcvr-y", help="Receiver longitude", default=None)

@click.option("--verbose", help="Option to print output to screen", default=None, type=bool)
@click.option("--iterations", help="Maximum Levenberg-Marquardt iterations", default=None)
@click.option("--damping", help="Damping factor in Levenberg-Marquardt algorithm", default=None)
@click.option("--tolerance", help="Receiver accuracy threshold", default=None)
@click.option("--az-dev-lim", help="Azimuth deviation limit", default=None)
@click.option("--incl_step-min", help="Minimum inclination step", default=None)
@click.option("--incl_step-max", help="Maximum inclination step", default=None)

@click.option("--freq", help="Frequency for Sutherland-Bass losses", default=None)
@click.option("--abs-coeff", help="Scaling coefficient for Sutherland-Bass losses", default=None)
@click.option("--z-grnd", help="Ground elevation (flat ground, km rel. sea level)", default=None)
@click.option("--write-atmo", help="Option to write atmosphere data (for QC)", default=None)
@click.option("--output-id", help="User specified output file path", default=None)
@click.option("--max-alt", help="Maximum altitude for ray calculation", default=None)
@click.option("--max-rng", help="Maximum range for ray calculation", default=None)

@click.option("--min-x", help="Minimum E/W offset for user defined bounds", default=None)
@click.option("--max-x", help="Maximum E/W offset for user defined bounds", default=None)
@click.option("--min-y", help="Minimum N/S offset for user defined bounds", default=None)
@click.option("--max-y", help="Maximum N/S offset for user defined bounds", default=None)

@click.option("--min-ds", help="Minimum step size (near-ground) in RK4 solver", default=None)
@click.option("--max-ds", help="Maximum step size in RK4 solver", default=None)
@click.option("--max-s", help="Maximum ray length between bounces", default=None)

@click.option("--topo-file", help="Terrain file", default=None)
@click.option("--topo-BL-wind", help="Use terrain corrected boundary layer winds", default=None, type=bool)
@click.option("--cpu-cnt", help="Number of CPUs to use in analysis", default=1)
def run_3d_eig(config_file, atmo_file, atmo_prefix, grid_lats, grid_lons, incl_min, incl_max, bnc_min, bnc_max, bounces, 
                src_x, src_y, src_alt, rcvr_x, rcvr_y, verbose, iterations, damping, tolerance, az_dev_lim, 
                incl_step_min, incl_step_max, freq, abs_coeff, z_grnd, write_atmo, output_id, max_alt, max_rng, 
                min_x, max_x, min_y, max_y, min_ds, max_ds, max_s, topo_file, topo_bl_wind, cpu_cnt):
    '''
    Run 3D Cartesian eigenray analysis to identify propagation paths connecting a specific source-receiver geometry

    \b
    Examples:
    \t infraga 3d eigenray --atmo-file ToyAtmo.met --rcvr-x -175.0 --rcvr-y 75.0 --bnc-max 1 --verbose True
    \t infraga 3d eigenray --atmo-file ToyAtmo.met --config-file example.cnfg

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

    # Set propagation specific parameters
    incl_min = define_param(user_config, 'EIGENRAY', 'incl_min', incl_min)
    incl_max = define_param(user_config, 'EIGENRAY', 'incl_max', incl_max)

    bnc_min = define_param(user_config, 'EIGENRAY', 'bnc_min', bnc_min)
    bnc_max = define_param(user_config, 'EIGENRAY', 'bnc_max', bnc_max)
    bounces = define_param(user_config, 'EIGENRAY', 'bounces', bounces)

    src_x = define_param(user_config, 'EIGENRAY', 'src_x', src_x)
    src_y = define_param(user_config, 'EIGENRAY', 'src_y', src_y)
    src_alt = define_param(user_config, 'EIGENRAY', 'src_alt', src_alt)

    rcvr_x = define_param(user_config, 'EIGENRAY', 'rcvr_x', rcvr_x)
    rcvr_y = define_param(user_config, 'EIGENRAY', 'rcvr_y', rcvr_y)

    verbose = define_param(user_config, 'EIGENRAY', 'verbose', verbose)
    iterations = define_param(user_config, 'EIGENRAY', 'iterations', iterations)
    damping = define_param(user_config, 'EIGENRAY', 'damping', damping)
    tolerance = define_param(user_config, 'EIGENRAY', 'tolerance', tolerance)
    az_dev_lim = define_param(user_config, 'EIGENRAY', 'az_dev_lim', az_dev_lim)
    incl_step_min = define_param(user_config, 'EIGENRAY', 'incl_step_min', incl_step_min)
    incl_step_max = define_param(user_config, 'EIGENRAY', 'incl_step_max', incl_step_max)

    # Set general parameters
    freq = define_param(user_config, 'GENERAL', 'freq', freq)
    abs_coeff = define_param(user_config, 'GENERAL', 'abs_coeff', abs_coeff)

    write_atmo = define_param(user_config, 'GENERAL', 'write_atmo', write_atmo)
    output_id = define_param(user_config, 'GENERAL', 'output_id', output_id)

    max_alt = define_param(user_config, 'GENERAL', 'max_alt', max_alt)
    max_rng = define_param(user_config, 'GENERAL', 'max_rng', max_rng)

    min_x = define_param(user_config, 'GENERAL', 'min_x', min_x)
    max_x = define_param(user_config, 'GENERAL', 'max_x', max_x)
    min_y = define_param(user_config, 'GENERAL', 'min_y', min_y)
    max_y = define_param(user_config, 'GENERAL', 'max_y', max_y)

    min_ds = define_param(user_config, 'GENERAL', 'min_ds', min_ds)
    max_ds = define_param(user_config, 'GENERAL', 'max_ds', max_ds)
    max_s = define_param(user_config, 'GENERAL', 'max_s', max_s)

    z_grnd = define_param(user_config, 'GENERAL', 'z_grnd', z_grnd)
    topo_file = define_param(user_config, 'GENERAL', 'topo_file', topo_file)
    topo_bl_wind = define_param(user_config, 'GENERAL', 'topo_bl_wind', topo_bl_wind)

    cpu_cnt = define_param(user_config, 'GENERAL', 'cpu_cnt', cpu_cnt)
    cpu_cnt = int(cpu_cnt)

    # Build run command
    if cpu_cnt < 2:
        command = bin_path + "infraga-3d"
    else:
        command = "mpirun -np " + str(cpu_cnt) + " " + bin_path + "infraga-accel-3d"

    if atmo_prefix is not None:
        command = command + "-rngdep"

    if atmo_file is not None:
        command = command + " -eig_search " + atmo_file
    elif atmo_prefix is not None and grid_lats is not None and grid_lons is not None:
        command = command + " -eig_search " + atmo_prefix + " " + grid_lats + " " + grid_lons
    else:
        click.echo("Simulation requires either an '--atmo-file' or --atmo-prefix' with grid info (--grid-lats and --grid-lons)")
        return 0

    # Set parameters
    command = set_param(command, incl_min, "incl_min")
    command = set_param(command, incl_max, "incl_max")

    command = set_param(command, bnc_min, "bnc_min")
    command = set_param(command, bnc_max, "bnc_max")
    command = set_param(command, bounces, "bounces")

    command = set_param(command, src_x, "src_x")
    command = set_param(command, src_y, "src_y")
    command = set_param(command, src_alt, "src_alt")

    command = set_param(command, rcvr_x, "rcvr_x")
    command = set_param(command, rcvr_y, "rcvr_y")

    if verbose is not None:
        command = set_param(command, str(verbose), "verbose")

    command = set_param(command, iterations, "iterations")
    command = set_param(command, damping, "damping")
    command = set_param(command, tolerance, "tolerance")
    command = set_param(command, az_dev_lim, "az_dev_lim")
    command = set_param(command, incl_step_min, "incl_step_min")
    command = set_param(command, incl_step_max, "incl_step_max")

    command = set_param(command, freq, "freq")
    command = set_param(command, abs_coeff, "abs_coeff")

    if write_atmo is not None:
        command = set_param(command, str(write_atmo), "write_atmo")

    command = set_param(command, output_id, "output_id")

    command = set_param(command, max_alt, "max_alt")
    command = set_param(command, max_rng, "max_rng")

    command = set_param(command, min_x, "min_x")
    command = set_param(command, max_x, "max_x")
    command = set_param(command, min_y, "min_y")
    command = set_param(command, max_y, "max_y")

    command = set_param(command, min_ds, "min_ds")
    command = set_param(command, max_ds, "max_ds")
    command = set_param(command, max_s, "max_s")

    command = set_param(command, z_grnd, "z_grnd")
    command = set_param(command, topo_file, "topo_file")
    if topo_file is not None:
        if topo_bl_wind is not None:
            command = set_param(command, str(topo_bl_wind), "topo_bl_wind")

    # print(command)
    subprocess.call(command, shell=True)


@click.command('wnl-wvfrm', short_help="Run weakly non-linear waveform calculation")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)
@click.option("--atmo-prefix", help="Atmosphere file prefix (range dependent)", default=None)
@click.option("--grid-lats", help="Atmosphere grid latitudes file", default=None)
@click.option("--grid-lons", help="Atmosphere grid longitudes file", default=None)

@click.option("--inclination", help="Inclination angle (rel. horizontal)", default=None)
@click.option("--azimuth", help="Azimuth angle (clockwise rel. N)", default=None)
@click.option("--bounces", help="Number of ground reflections (bounces) to consider", default=None)

@click.option("--src-x", help="Source location (E/W)", default=None)
@click.option("--src-y", help="Source location (N/S)", default=None)
@click.option("--src-alt", help="Source altitude", default=None)

@click.option("--write-ray", help="Option to write ray path info", default=None, type=bool)

@click.option("--wvfrm-file", help="File containing reference waveform", default=None)
@click.option("--wvfrm-opt", help="Waveform option ('impulse', 'Uwave', 'Nwave')", default=None)

@click.option("--wvfrm-p0", help="Waveform peak overpressure [Pa]", default=None)
@click.option("--wvfrm-t0", help="Waveform positive phase duration [s]", default=None)
@click.option("--wvfrm-alpha", help="Waveform shaping parameter", default=None)
@click.option("--wvfrm-ref", help="Waveform reference distance", default=None)
@click.option("--wvfrm-out-step", help="Waveform output step along ray path", default=None)
@click.option("--wvfrm-yield", help="Eq. TNT yield (kg) for an explosive source", default=None)

@click.option("--wvfrm-ds", help="Burgers equation solver resolution [km]", default=None)
@click.option("--wvfrm-len", help="Waveform sample count (default = 2e13)", default=None)

@click.option("--freq", help="Frequency for Sutherland-Bass losses", default=None)
@click.option("--abs-coeff", help="Scaling coefficient for Sutherland-Bass losses", default=None)
@click.option("--z-grnd", help="Ground elevation (flat ground, km rel. sea level)", default=None)
@click.option("--write-atmo", help="Option to write atmosphere data (for QC)", default=None)
@click.option("--output-id", help="User specified output file path", default=None)
@click.option("--max-alt", help="Maximum altitude for ray calculation", default=None)
@click.option("--max-rng", help="Maximum range for ray calculation", default=None)

@click.option("--min-x", help="Minimum E/W offset for user defined bounds", default=None)
@click.option("--max-x", help="Maximum E/W offset for user defined bounds", default=None)
@click.option("--min-y", help="Minimum N/S offset for user defined bounds", default=None)
@click.option("--max-y", help="Maximum N/S offset for user defined bounds", default=None)

@click.option("--min-ds", help="Minimum step size (near-ground) in RK4 solver", default=None)
@click.option("--max-ds", help="Maximum step size in RK4 solver", default=None)
@click.option("--max-s", help="Maximum ray length between bounces", default=None)

@click.option("--topo-file", help="Terrain file", default=None)
@click.option("--topo-BL-wind", help="Use terrain corrected boundary layer winds", default=None, type=bool)
def run_3d_wvfrm(config_file, atmo_file, atmo_prefix, grid_lats, grid_lons, inclination, azimuth, bounces, src_x, src_y, 
                src_alt, write_ray, wvfrm_file, wvfrm_opt, wvfrm_p0, wvfrm_t0, wvfrm_alpha, wvfrm_ref, wvfrm_out_step,
                wvfrm_yield, wvfrm_ds, wvfrm_len, freq, abs_coeff, z_grnd, write_atmo, output_id, max_alt, max_rng, min_x, max_x, min_y, max_y, 
                min_ds, max_ds, max_s, topo_file, topo_bl_wind):
    '''
    Run weakly non-linear waveform analysis along a 3D Cartesian ray path.

    \b
    Examples:
    \t infraga 3d wnl-wvfrm --atmo-file ToyAtmo.met --inclination 4.1314876 --azimuth -84.969455 --bounces 1 --wvfrm-p0 500.0
    \t infraga 3d wnl-wvfrm --atmo-file ToyAtmo.met --config-file example.cnfg
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

    # Set propagation specific parameters
    inclination = define_param(user_config, 'WAVEFORM', 'inclination', inclination)
    azimuth = define_param(user_config, 'WAVEFORM', 'azimuth', azimuth)
    bounces = define_param(user_config, 'WAVEFORM', 'bounces', bounces)

    src_x = define_param(user_config, 'EIGENRAY', 'src_x', src_x)
    src_y = define_param(user_config, 'EIGENRAY', 'src_y', src_y)
    src_alt = define_param(user_config, 'WAVEFORM', 'src_alt', src_alt)

    write_ray = define_param(user_config, 'WAVEFORM', 'write_ray', write_ray)

    wvfrm_file = define_param(user_config, 'WAVEFORM', 'wvfrm_file', wvfrm_file)
    wvfrm_opt = define_param(user_config, 'WAVEFORM', 'wvfrm_opt', wvfrm_opt)
    wvfrm_p0 = define_param(user_config, 'WAVEFORM', 'wvfrm_p0', wvfrm_p0)
    wvfrm_t0 = define_param(user_config, 'WAVEFORM', 'wvfrm_t0', wvfrm_t0)
    wvfrm_alpha = define_param(user_config, 'WAVEFORM', 'wvfrm_alpha', wvfrm_alpha)
    wvfrm_ref = define_param(user_config, 'WAVEFORM', 'wvfrm_ref', wvfrm_ref)
    wvfrm_out_step = define_param(user_config, 'WAVEFORM', 'wvfrm_out_step', wvfrm_out_step)
    wvfrm_yield = define_param(user_config, 'WAVEFORM', 'wvfrm_yield', wvfrm_yield)

    wvfrm_ds = define_param(user_config, 'WAVEFORM', 'wvfrm_ds', wvfrm_ds)
    wvfrm_len = define_param(user_config, 'WAVEFORM', 'wvfrm_len', wvfrm_len)

    # Set general parameters
    freq = define_param(user_config, 'GENERAL', 'freq', freq)
    abs_coeff = define_param(user_config, 'GENERAL', 'abs_coeff', abs_coeff)

    write_atmo = define_param(user_config, 'GENERAL', 'write_atmo', write_atmo)
    output_id = define_param(user_config, 'GENERAL', 'output_id', output_id)

    max_alt = define_param(user_config, 'GENERAL', 'max_alt', max_alt)
    max_rng = define_param(user_config, 'GENERAL', 'max_rng', max_rng)

    min_x = define_param(user_config, 'GENERAL', 'min_x', min_x)
    max_x = define_param(user_config, 'GENERAL', 'max_x', max_x)
    min_y = define_param(user_config, 'GENERAL', 'min_y', min_y)
    max_y = define_param(user_config, 'GENERAL', 'max_y', max_y)

    min_ds = define_param(user_config, 'GENERAL', 'min_ds', min_ds)
    max_ds = define_param(user_config, 'GENERAL', 'max_ds', max_ds)
    max_s = define_param(user_config, 'GENERAL', 'max_s', max_s)

    z_grnd = define_param(user_config, 'GENERAL', 'z_grnd', z_grnd)
    topo_file = define_param(user_config, 'GENERAL', 'topo_file', topo_file)
    topo_bl_wind = define_param(user_config, 'GENERAL', 'topo_bl_wind', topo_bl_wind)

    # Build run command
    command = bin_path + "infraga-3d"
    
    if atmo_prefix is not None:
        command = command + "-rngdep"

    if atmo_file is not None:
        command = command + " -wnl_wvfrm " + atmo_file
    elif atmo_prefix is not None and grid_lats is not None and grid_lons is not None:
        command = command + " -wnl_wvfrm " + atmo_prefix + " " + grid_lats + " " + grid_lons
    else:
        click.echo("Simulation requires either an '--atmo-file' or --atmo-prefix' with grid info (--grid-lats and --grid-lons)")
        return 0

    # Set parameters
    command = set_param(command, inclination, "inclination")
    command = set_param(command, azimuth, "azimuth")
    command = set_param(command, bounces, "bounces")

    command = set_param(command, src_x, "src_x")
    command = set_param(command, src_y, "src_y")
    command = set_param(command, src_alt, "src_alt")

    if write_ray is not None:
        command = set_param(command, str(write_ray), "write_ray")

    command = set_param(command, wvfrm_file, "wvfrm_file")
    if wvfrm_yield is not None:
        # use atmosphere to define ambient pressure (1 mbar = 0.1 kPa)
        atmo = np.loadtxt(atmo_file)

        if src_alt is not None:
            src_ht = float(src_alt)
        elif z_grnd is not None:
            src_ht = float(z_grnd)
        else:
            src_ht = 0.0

        p_ambient = atmo[np.argmin(abs(atmo[:, 0] - src_ht)), 5] * 0.1
        T_ambient = atmo[np.argmin(abs(atmo[:, 0] - src_ht)), 1]

        wvfrm_opt = 'impulse'
        wvfrm_ref = str(0.035 * float(wvfrm_yield)**(1.0 / 3.0))
        wvfrm_p0 = str(kg_op(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient))
        wvfrm_t0 = str(kg_ppd(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient))
        wvfrm_alpha = '0.01'

    command = set_param(command, wvfrm_opt, "wvfrm_opt")
    command = set_param(command, wvfrm_p0, "wvfrm_p0")
    command = set_param(command, wvfrm_t0, "wvfrm_t0")
    command = set_param(command, wvfrm_alpha, "wvfrm_alpha")
    command = set_param(command, wvfrm_ref, "wvfrm_ref")
    command = set_param(command, wvfrm_out_step, "wvfrm_out_step")

    command = set_param(command, wvfrm_ds, "wvfrm_ds")
    command = set_param(command, wvfrm_len, "wvfrm_len")
    
    command = set_param(command, freq, "freq")
    command = set_param(command, abs_coeff, "abs_coeff")

    if write_atmo is not None:
        command = set_param(command, str(write_atmo), "write_atmo")

    command = set_param(command, output_id, "output_id")

    command = set_param(command, max_alt, "max_alt")
    command = set_param(command, max_rng, "max_rng")

    command = set_param(command, min_x, "min_x")
    command = set_param(command, max_x, "max_x")
    command = set_param(command, min_y, "min_y")
    command = set_param(command, max_y, "max_y")

    command = set_param(command, min_ds, "min_ds")
    command = set_param(command, max_ds, "max_ds")
    command = set_param(command, max_s, "max_s")

    command = set_param(command, z_grnd, "z_grnd")
    command = set_param(command, topo_file, "topo_file")
    if topo_file is not None:
        if topo_bl_wind is not None:
            command = set_param(command, str(topo_bl_wind), "topo_bl_wind")

    # print(command)
    subprocess.call(command, shell=True)


@click.command('eig_wvfrm', short_help="Run eigenray search and compute waveform contributions")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)
@click.option("--atmo-prefix", help="Atmosphere file prefix (range dependent)", default=None)
@click.option("--grid-lats", help="Atmosphere grid latitudes file", default=None)
@click.option("--grid-lons", help="Atmosphere grid longitudes file", default=None)

@click.option("--incl-min", help="Minimum inclination angle (rel. horizontal)", default=None)
@click.option("--incl-max", help="Maximum inclination angle (rel. horizontal)", default=None)

@click.option("--bnc-min", help="Minimum number of ground reflections (bounces)", default=None)
@click.option("--bnc-max", help="Maximum number of ground reflections (bounces)", default=None)
@click.option("--bounces", help="Specific number of ground reflections (bounces)", default=None)

@click.option("--src-x", help="Source latitude", default=None)
@click.option("--src-y", help="Source longitude", default=None)
@click.option("--src-alt", help="Source altitude", default=None)
@click.option("--rcvr-x", help="Receiver latitude", default=None)
@click.option("--rcvr-y", help="Receiver longitude", default=None)

@click.option("--wvfrm-file", help="File containing reference waveform", default=None)
@click.option("--wvfrm-opt", help="Waveform option ('impulse', 'Uwave', 'Nwave')", default=None)

@click.option("--verbose", help="Option to print output to screen", default=None, type=bool)
@click.option("--iterations", help="Maximum Levenberg-Marquardt iterations", default=None)
@click.option("--damping", help="Damping factor in Levenberg-Marquardt algorithm", default=None)
@click.option("--tolerance", help="Receiver accuracy threshold", default=None)
@click.option("--az-dev-lim", help="Azimuth deviation limit", default=None)
@click.option("--incl_step-min", help="Minimum inclination step", default=None)
@click.option("--incl_step-max", help="Maximum inclination step", default=None)

@click.option("--wvfrm-p0", help="Waveform peak overpressure [Pa]", default=None)
@click.option("--wvfrm-t0", help="Waveform positive phase duration [s]", default=None)
@click.option("--wvfrm-alpha", help="Waveform shaping parameter", default=None)
@click.option("--wvfrm-ref", help="Waveform reference distance", default=None)
@click.option("--wvfrm-out-step", help="Waveform output step along ray path", default=None)
@click.option("--wvfrm-yield", help="Eq. TNT yield (kg) for an explosive source", default=None)

@click.option("--wvfrm-ds", help="Burgers equation solver resolution [km]", default=None)
@click.option("--wvfrm-len", help="Waveform sample count (default = 2e13)", default=None)

@click.option("--freq", help="Frequency for Sutherland-Bass losses", default=None)
@click.option("--abs-coeff", help="Scaling coefficient for Sutherland-Bass losses", default=None)
@click.option("--z-grnd", help="Ground elevation (km rel. sea level)", default=None)
@click.option("--write-atmo", help="Option to write atmosphere data (for QC)", default=None)
@click.option("--output-id", help="User specified output file path", default=None)
@click.option("--max-alt", help="Maximum altitude for ray calculation", default=None)
@click.option("--max-rng", help="Maximum range for ray calculation", default=None)

@click.option("--min-x", help="Minimum latitude for user defined bounds", default=None)
@click.option("--max-x", help="Maximum latitude for user defined bounds", default=None)
@click.option("--min-y", help="Minimum longitude for user defined bounds", default=None)
@click.option("--max-y", help="Maximum longitude for user defined bounds", default=None)

@click.option("--min-ds", help="Minimum step size (near-ground) in RK4 solver", default=None)
@click.option("--max-ds", help="Maximum step size in RK4 solver", default=None)
@click.option("--max-s", help="Maximum ray length between bounces", default=None)

@click.option("--topo-file", help="Terrain file", default=None)
@click.option("--topo-BL-wind", help="Use terrain corrected boundary layer winds", default=None, type=bool)
@click.option("--cpu-cnt", help="Number of CPUs to use in analysis", default=1)

@click.option("--keep-eig-arrivals", help="Keep eigenray arrivals", default=False)
def run_3d_eig_wvfrm(config_file, atmo_file, atmo_prefix, grid_lats, grid_lons, incl_min, incl_max, bnc_min, bnc_max, bounces, 
                src_x, src_y, src_alt, rcvr_x, rcvr_y, verbose, iterations, damping, tolerance, az_dev_lim, incl_step_min, 
                incl_step_max, wvfrm_file, wvfrm_opt, wvfrm_p0, wvfrm_t0, wvfrm_alpha, wvfrm_ref, wvfrm_out_step, wvfrm_yield, 
                wvfrm_ds, wvfrm_len, freq, abs_coeff, z_grnd, write_atmo, output_id, max_alt, max_rng, min_x, max_x, min_y, 
                max_y, min_ds, max_ds, max_s, topo_file, topo_bl_wind, cpu_cnt, keep_eig_arrivals):
    '''
    Run 3d geometry eigenray analysis to identify propagation paths connecting a specific source-receiver geometry and then compute weakly-nonlinear waveform predictions for each eigenray

    \b
    Examples:
    \t infraga 3d eig_wvfrm --atmo-file ToyAtmo.met --rcvr-x -175.0 --rcvr-y 75.0 --bnc-max 1 --verbose True --wvfrm-yield 10e3
    \t infraga 3d eig_wvfrm --atmo-file ToyAtmo.met --config-file example.cnfg

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

    # Set eigenray specific parameters
    incl_min = define_param(user_config, 'EIGENRAY', 'incl_min', incl_min)
    incl_max = define_param(user_config, 'EIGENRAY', 'incl_max', incl_max)

    bnc_min = define_param(user_config, 'EIGENRAY', 'bnc_min', bnc_min)
    bnc_max = define_param(user_config, 'EIGENRAY', 'bnc_max', bnc_max)
    bounces = define_param(user_config, 'EIGENRAY', 'bounces', bounces)

    src_x = define_param(user_config, 'EIGENRAY', 'src_x', src_x)
    src_y = define_param(user_config, 'EIGENRAY', 'src_y', src_y)
    src_alt = define_param(user_config, 'EIGENRAY', 'src_alt', src_alt)

    rcvr_x = define_param(user_config, 'EIGENRAY', 'rcvr_x', rcvr_x)
    rcvr_y = define_param(user_config, 'EIGENRAY', 'rcvr_y', rcvr_y)

    verbose = define_param(user_config, 'EIGENRAY', 'verbose', verbose)
    iterations = define_param(user_config, 'EIGENRAY', 'iterations', iterations)
    damping = define_param(user_config, 'EIGENRAY', 'damping', damping)
    tolerance = define_param(user_config, 'EIGENRAY', 'tolerance', tolerance)
    az_dev_lim = define_param(user_config, 'EIGENRAY', 'az_dev_lim', az_dev_lim)
    incl_step_min = define_param(user_config, 'EIGENRAY', 'incl_step_min', incl_step_min)
    incl_step_max = define_param(user_config, 'EIGENRAY', 'incl_step_max', incl_step_max)

    # Set waveform specific parameters 
    wvfrm_file = define_param(user_config, 'WAVEFORM', 'wvfrm_file', wvfrm_file)
    wvfrm_opt = define_param(user_config, 'WAVEFORM', 'wvfrm_opt', wvfrm_opt)
    wvfrm_p0 = define_param(user_config, 'WAVEFORM', 'wvfrm_p0', wvfrm_p0)
    wvfrm_t0 = define_param(user_config, 'WAVEFORM', 'wvfrm_t0', wvfrm_t0)
    wvfrm_alpha = define_param(user_config, 'WAVEFORM', 'wvfrm_alpha', wvfrm_alpha)
    wvfrm_ref = define_param(user_config, 'WAVEFORM', 'wvfrm_ref', wvfrm_ref)
    wvfrm_out_step = define_param(user_config, 'WAVEFORM', 'wvfrm_out_step', wvfrm_out_step)
    wvfrm_yield = define_param(user_config, 'WAVEFORM', 'wvfrm_yield', wvfrm_yield)

    wvfrm_ds = define_param(user_config, 'WAVEFORM', 'wvfrm_ds', wvfrm_ds)
    wvfrm_len = define_param(user_config, 'WAVEFORM', 'wvfrm_len', wvfrm_len)
    
    # Set general parameters
    freq = define_param(user_config, 'GENERAL', 'freq', freq)
    abs_coeff = define_param(user_config, 'GENERAL', 'abs_coeff', abs_coeff)

    write_atmo = define_param(user_config, 'GENERAL', 'write_atmo', write_atmo)
    output_id = define_param(user_config, 'GENERAL', 'output_id', output_id)

    max_alt = define_param(user_config, 'GENERAL', 'max_alt', max_alt)
    max_rng = define_param(user_config, 'GENERAL', 'max_rng', max_rng)

    min_x = define_param(user_config, 'GENERAL', 'min_x', min_x)
    max_x = define_param(user_config, 'GENERAL', 'max_x', max_x)
    min_y = define_param(user_config, 'GENERAL', 'min_y', min_y)
    max_y = define_param(user_config, 'GENERAL', 'max_y', max_y)

    min_ds = define_param(user_config, 'GENERAL', 'min_ds', min_ds)
    max_ds = define_param(user_config, 'GENERAL', 'max_ds', max_ds)
    max_s = define_param(user_config, 'GENERAL', 'max_s', max_s)

    z_grnd = define_param(user_config, 'GENERAL', 'z_grnd', z_grnd)
    topo_file = define_param(user_config, 'GENERAL', 'topo_file', topo_file)
    topo_bl_wind = define_param(user_config, 'GENERAL', 'topo_bl_wind', topo_bl_wind)

    cpu_cnt = define_param(user_config, 'GENERAL', 'cpu_cnt', cpu_cnt)
    cpu_cnt = int(cpu_cnt)

    # Check if eigenray analysis is already done
    if output_id is not None:
        eig_arrivals_file = output_id + ".arrivals.dat"
        profile_id = output_id
    else:
        if atmo_prefix is not None:
            eig_arrivals_file = atmo_prefix + ".arrivals.dat"
            profile_id = atmo_prefix
        else:
            eig_arrivals_file = os.path.splitext(atmo_file)[0] + ".arrivals.dat"
            profile_id = os.path.splitext(atmo_file)[0]
    
    if os.path.isfile(eig_arrivals_file):
        print("Eigenray results found.  Skipping to waveform calculation...")
    else:
        print("Running eigenray analysis...")
        # Build run command
        if cpu_cnt < 2:
            command = bin_path + "infraga-3d"
        else:
            command = "mpirun -np " + str(cpu_cnt) + " " + bin_path + "infraga-accel-3d"

        if atmo_prefix is not None:
            command = command + "-rngdep"

        if atmo_file is not None:
            command = command + " -eig_search " + atmo_file
        elif atmo_prefix is not None and grid_lats is not None and grid_lons is not None:
            command = command + " -eig_search " + atmo_prefix + " " + grid_lats + " " + grid_lons
        else:
            click.echo("Simulation requires either an '--atmo-file' or --atmo-prefix' with grid info (--grid-lats and --grid-lons)")
            return 0

        # Set parameters
        command = set_param(command, incl_min, "incl_min")
        command = set_param(command, incl_max, "incl_max")

        command = set_param(command, bnc_min, "bnc_min")
        command = set_param(command, bnc_max, "bnc_max")
        command = set_param(command, bounces, "bounces")

        command = set_param(command, src_x, "src_x")
        command = set_param(command, src_y, "src_y")
        command = set_param(command, src_alt, "src_alt")

        command = set_param(command, rcvr_x, "rcvr_x")
        command = set_param(command, rcvr_y, "rcvr_y")

        if verbose is not None:
            command = set_param(command, str(verbose), "verbose")

        command = set_param(command, iterations, "iterations")
        command = set_param(command, damping, "damping")
        command = set_param(command, tolerance, "tolerance")
        command = set_param(command, az_dev_lim, "az_dev_lim")
        command = set_param(command, incl_step_min, "incl_step_min")
        command = set_param(command, incl_step_max, "incl_step_max")

        command = set_param(command, freq, "freq")
        command = set_param(command, abs_coeff, "abs_coeff")

        if write_atmo is not None:
            command = set_param(command, str(write_atmo), "write_atmo")

        command = set_param(command, output_id, "output_id")

        command = set_param(command, max_alt, "max_alt")
        command = set_param(command, max_rng, "max_rng")

        command = set_param(command, min_x, "min_x")
        command = set_param(command, max_x, "max_x")
        command = set_param(command, min_y, "min_y")
        command = set_param(command, max_y, "max_y")

        command = set_param(command, min_ds, "min_ds")
        command = set_param(command, max_ds, "max_ds")
        command = set_param(command, max_s, "max_s")

        command = set_param(command, z_grnd, "z_grnd")
        command = set_param(command, topo_file, "topo_file")
        if topo_file is not None:
            if topo_bl_wind is not None:
                command = set_param(command, str(topo_bl_wind), "topo_bl_wind")

        # print(command)
        subprocess.call(command, shell=True)


    # Analyze eigenray results and compute waveform contributions...
    print("Computing waveform...")
    results_id = eig_arrivals_file[:-13]

    eig_results = np.loadtxt(eig_arrivals_file)
    if len(eig_results) > 0:
        eig_results = np.atleast_2d(eig_results)
        os.system("rm " + results_id + ".eigenray*")

        mask = np.ones(len(eig_results), dtype=bool)
        for n, line in enumerate(eig_results):
            if mask[n]:   
                for m in range(n + 1, len(eig_results)):   
                    arrivals_diff  = (line[0] - eig_results[m][0])**2
                    arrivals_diff += (line[1] - eig_results[m][1])**2
                    arrivals_diff += (line[5] - eig_results[m][5])**2
                    arrivals_diff += (line[7] - eig_results[m][7])**2
                    if arrivals_diff < 0.1:
                        mask[m] = False
        eig_results = eig_results[mask]

        # Create waveform interpolations
        file_out = open(results_id + ".eigenrays.dat", 'w')
        print("# 'infraga 3d eig_wvfrm' eigenray results", '\n#', file=file_out)
        if atmo_prefix is not None:
            print("# 	profile: " + atmo_prefix, file=file_out)
        else:
            print("# 	profile: " + atmo_file, file=file_out)
        print("# 	source location: " + str(src_x) + ", " + str(src_y) + ", " + str(src_alt) , file=file_out)
        print("# 	receiver location: " + str(rcvr_x) + ", " + str(rcvr_y) + ", 0.0", file=file_out)
        print("# 	inclination range: " + str(incl_min) + ", " + str(incl_max), file=file_out)
        print("#    inclination step max:", incl_step_max, file=file_out)
        print("# 	bounces: " + str(bnc_min) + ", " + str(bnc_max), file=file_out)
        if topo_file is not None:
            print("# 	terrain file: " + str(topo_file), file=file_out)        
        else:
            print("# 	ground elevation: " + str(z_grnd), file=file_out)
        print("# 	damping:", damping, file=file_out)
        print("# 	range max:", max_rng, file=file_out)

        wvfrms = []
        t_lims = [np.inf, 0.0]
        for n, line in enumerate(eig_results):
            command = bin_path + "infraga-3d"
            
            if atmo_prefix is not None:
                command = command + "-rngdep"

            if atmo_file is not None:
                command = command + " -wnl_wvfrm " + atmo_file
            elif atmo_prefix is not None and grid_lats is not None and grid_lons is not None:
                command = command + " -wnl_wvfrm " + atmo_prefix + " " + grid_lats + " " + grid_lons
            else:
                click.echo("Simulation requires either an '--atmo-file' or --atmo-prefix' with grid info (--grid-lats and --grid-lons)")
                return 0

            # Set parameters
            command = set_param(command, str(line[0]), "inclination")
            command = set_param(command, str(line[1]), "azimuth")
            command = set_param(command, str(line[2]), "bounces")

            command = set_param(command, src_x, "src_x")
            command = set_param(command, src_y, "src_y")
            command = set_param(command, src_alt, "src_alt")
            command = command + " write_ray=true"

            command = set_param(command, wvfrm_file, "wvfrm_file")
            if wvfrm_yield is not None:
                # use atmosphere to define ambient pressure (1 mbar = 0.1 kPa)
                atmo = np.loadtxt(atmo_file)

                if src_alt is not None:
                    src_ht = float(src_alt)
                elif z_grnd is not None:
                    src_ht = float(z_grnd)
                else:
                    src_ht = 0.0

                p_ambient = atmo[np.argmin(abs(atmo[:, 0] - src_ht)), 5] * 0.1
                T_ambient = atmo[np.argmin(abs(atmo[:, 0] - src_ht)), 1]

                wvfrm_opt = 'impulse'
                wvfrm_ref = str(0.035 * float(wvfrm_yield)**(1.0 / 3.0))
                wvfrm_p0 = str(kg_op(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient))
                wvfrm_t0 = str(kg_ppd(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient))
                wvfrm_alpha = '0.01'

            command = set_param(command, wvfrm_opt, "wvfrm_opt")
            command = set_param(command, wvfrm_p0, "wvfrm_p0")
            command = set_param(command, wvfrm_t0, "wvfrm_t0")
            command = set_param(command, wvfrm_alpha, "wvfrm_alpha")
            command = set_param(command, wvfrm_ref, "wvfrm_ref")
            command = set_param(command, wvfrm_out_step, "wvfrm_out_step")

            command = set_param(command, wvfrm_ds, "wvfrm_ds")
            command = set_param(command, wvfrm_len, "wvfrm_len")
            
            command = set_param(command, freq, "freq")
            command = set_param(command, abs_coeff, "abs_coeff")

            if write_atmo is not None:
                command = set_param(command, str(write_atmo), "write_atmo")

            command = set_param(command, output_id, "output_id")

            command = set_param(command, max_alt, "max_alt")
            command = set_param(command, max_rng, "max_rng")

            command = set_param(command, min_x, "min_x")
            command = set_param(command, max_x, "max_x")
            command = set_param(command, min_y, "min_y")
            command = set_param(command, max_y, "max_y")

            command = set_param(command, min_ds, "min_ds")
            command = set_param(command, max_ds, "max_ds")
            command = set_param(command, max_s, "max_s")

            command = set_param(command, z_grnd, "z_grnd")
            command = set_param(command, topo_file, "topo_file")
            if topo_file is not None:
                if topo_bl_wind is not None:
                    command = set_param(command, str(topo_bl_wind), "topo_bl_wind")

            # print(command)
            subprocess.call(command, shell=True)

            temp = np.loadtxt(profile_id + ".wvfrm_out.dat")
            wvfrms += [interp1d(temp[:, 0], temp[:, 1], bounds_error=False, fill_value=0.0,  kind='cubic')]
            t_lims[0] = min(t_lims[0], temp[0][0])
            t_lims[1] = max(t_lims[1], temp[-1][0])
            dt = abs(temp[1][0] - temp[0][0])

            temp = np.loadtxt(profile_id + ".raypaths.dat")
            for line in temp:
                print(*line, file=file_out)
            print('\n', file=file_out)

            command = "rm " + profile_id + ".wvfrm_out.dat"
            command = command + " " + profile_id + ".raypaths.dat"
            command = command + " " + profile_id + ".wvfrm_init.dat"
            os.system(command)

        file_out.close()

        # combine waveforms and write to file
        print("Interpolating and merging waveforms...")
        print('\t' + "Eigenrays written into " + profile_id + ".eigenrays.dat")
        print('\t' + "Arrival waveform written into " + profile_id + ".wvfrms.dat")

        t_vals = np.arange(t_lims[0], t_lims[1], dt)
        file_out = open(profile_id + ".wvfrms.dat", 'w')
        print("# 'infraga 3f eig_wvfrm' waveform results", '\n#', file=file_out)
        if atmo_prefix is not None:
            print("# 	profile: " + atmo_prefix, file=file_out)        
        else:
            print("# 	profile: " + atmo_file, file=file_out)
        print("# 	source location: " + str(src_x) + ", " + str(src_y) + ", " + str(src_alt) , file=file_out)
        print("# 	receiver location: " + str(rcvr_x) + ", " + str(rcvr_y) + ", 0.0", file=file_out)
        print("# 	inclination range: " + str(incl_min) + ", " + str(incl_max), file=file_out)
        print("#    inclination step max:", incl_step_max, file=file_out)
        print("# 	bounces: " + str(bnc_min) + ", " + str(bnc_max), file=file_out)
        if topo_file is not None:
            print("# 	terrain file: " + str(topo_file), file=file_out)        
        else:
            print("# 	ground elevation: " + str(z_grnd), file=file_out)
        print("#    damping:", damping, file=file_out)
        print("#    range max:" + '\n#', max_rng, file=file_out)

        # Waveform calculation parameters
        if wvfrm_ref is not None:
            print("#    waveform reference distance:", wvfrm_ref, file=file_out)
        else:
            print("#    waveform reference distance: 1.0", file=file_out)

        if wvfrm_len is not None:
            print("#    waveform length:", wvfrm_len, file=file_out)
        else:
            print("#    waveform length: 0.05", wvfrm_len, file=file_out)

        if wvfrm_yield:
            print("#    waveform source yield:", wvfrm_yield, file=file_out)
        else:
            print("#    waveform option:", wvfrm_opt, file=file_out)
            print("#    waveform peak op:", wvfrm_p0, file=file_out)
            print("#    waveform time sc.:", wvfrm_t0, file=file_out)
            print("#    waveform shaping param.:", wvfrm_alpha, file=file_out)

        print("#", file=file_out)
        print("# Eigenray arrivals:", file=file_out)
        print("# incl [deg]	az [deg]	n_b	x_0 [deg]	y_0 [deg]	time [s]	cel [km/s]	turning ht [km]	inclination [deg]	back azimuth [deg]	trans. coeff. [dB]	absorption [dB]", file=file_out)
        for line in eig_results:
            print("#", *line, file=file_out) 

        print('\n#', "t [s]" + '\t' + "p1 [Pa]", '\t' + "p2 [Pa] ...", file=file_out)
        for n in range(len(t_vals)):
            print(t_vals[n], end='\t', file=file_out)
            for wvfrm in wvfrms:
                print(wvfrm(t_vals[n]), end='\t', file=file_out)
            print('', file=file_out)
        file_out.close()
    else:
        print('\n' + "No waveforms to compute.")
    
    if not keep_eig_arrivals:
        os.system("rm " + profile_id + ".arrivals.dat")




####################################
##                                ##
##   Python CLI for infraga-sph   ##
##                                ##
####################################
@click.command('prop', short_help="Run a point source propagation simulation")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)
@click.option("--atmo-prefix", help="Atmosphere file prefix (range dependent)", default=None)
@click.option("--grid-lats", help="Atmosphere grid latitudes file", default=None)
@click.option("--grid-lons", help="Atmosphere grid longitudes file", default=None)
@click.option("--incl-min", help="Minimum inclination angle (rel. horizontal)", default=None)
@click.option("--incl-max", help="Maximum inclination angle (rel. horizontal)", default=None)
@click.option("--incl-step", help="Inclination angle step resolution", default=None)
@click.option("--inclination", help="Single inclination angle", default=None)
@click.option("--az-min", help="Minimum azimuth angle (clockwise rel. N)", default=None)
@click.option("--az-max", help="Maximum azimuth angle (clockwise rel. N)", default=None)
@click.option("--az-step", help="Azimuth angle step resolution", default=None)
@click.option("--azimuth", help="Single azimuth angle", default=None)
@click.option("--bounces", help="Number of ground reflections (bounces) to consider", default=None)
@click.option("--src-lat", help="Source latitude", default=None)
@click.option("--src-lon", help="Source longitude", default=None)
@click.option("--src-alt", help="Source altitude", default=None)
@click.option("--write-rays", help="Option to write [...].raypaths.dat output", default=None, type=bool)
@click.option("--write-topo", help="Option to write terrain info under first ray", default=None, type=bool)

@click.option("--freq", help="Frequency for Sutherland-Bass losses", default=None)
@click.option("--abs-coeff", help="Scaling coefficient for Sutherland-Bass losses", default=None)
@click.option("--z-grnd", help="Ground elevation (flat ground, km rel. sea level)", default=None)
@click.option("--write-atmo", help="Option to write atmosphere data (for QC)", default=None)
@click.option("--output-id", help="User specified output file path", default=None)
@click.option("--calc-amp", help="Option to turn off transport coefficient calculation", default=None, type=bool)
@click.option("--max-alt", help="Maximum altitude for ray calculation", default=None)
@click.option("--max-rng", help="Maximum range for ray calculation", default=None)

@click.option("--min-lat", help="Minimum latitude for user defined bounds", default=None)
@click.option("--max-lat", help="Maximum latitude for user defined bounds", default=None)
@click.option("--min-lon", help="Minimum longitude for user defined bounds", default=None)
@click.option("--max-lon", help="Maximum longitude for user defined bounds", default=None)

@click.option("--min-ds", help="Minimum step size (near-ground) in RK4 solver", default=None)
@click.option("--max-ds", help="Maximum step size in RK4 solver", default=None)
@click.option("--max-s", help="Maximum ray length between bounces", default=None)

@click.option("--topo-file", help="Terrain file", default=None)
@click.option("--topo-BL-wind", help="Use terrain corrected boundary layer winds", default=None, type=bool)
@click.option("--cpu-cnt", help="Number of CPUs to use in analysis", default=1)
def run_sph_prop(config_file, atmo_file, atmo_prefix, grid_lats, grid_lons, incl_min, incl_max, incl_step, inclination, 
                az_min, az_max, az_step, azimuth, bounces, src_lat, src_lon, src_alt, write_rays, write_topo, freq, 
                abs_coeff, z_grnd, write_atmo, output_id, calc_amp, max_alt, max_rng, max_lat, min_lat, max_lon, min_lon,
                min_ds, max_ds, max_s, topo_file, topo_bl_wind, cpu_cnt):
    '''
    Run spherical atmosphere geometry ray tracing analysis for a point source

    \b
    Examples:
    \t infraga sph prop --atmo-file ToyAtmo.met
    \t infraga sph prop --atmo-file ToyAtmo.met --config-file example.cnfg

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

    # Set propagation specific parameters
    incl_min = define_param(user_config, 'PROP', 'incl_min', incl_min)
    incl_max = define_param(user_config, 'PROP', 'incl_max', incl_max)
    incl_step = define_param(user_config, 'PROP', 'incl_step', incl_step)

    az_min = define_param(user_config, 'PROP', 'az_min', az_min)
    az_max = define_param(user_config, 'PROP', 'az_max', az_max)
    az_step = define_param(user_config, 'PROP', 'az_step', az_step)

    inclination = define_param(user_config, 'PROP', 'inclination', inclination)
    azimuth = define_param(user_config, 'PROP', 'azimuth', azimuth)

    bounces = define_param(user_config, 'PROP', 'bounces', bounces)

    src_lat = define_param(user_config, 'PROP', 'src_lat', src_lat)
    src_lon = define_param(user_config, 'PROP', 'src_lon', src_lon)
    src_alt = define_param(user_config, 'PROP', 'src_alt', src_alt)

    write_rays = define_param(user_config, 'PROP', 'write_rays', write_rays)
    write_topo = define_param(user_config, 'PROP', 'write_topo', write_topo)

    # Set general parameters
    freq = define_param(user_config, 'GENERAL', 'freq', freq)
    abs_coeff = define_param(user_config, 'GENERAL', 'abs_coeff', abs_coeff)

    write_atmo = define_param(user_config, 'GENERAL', 'write_atmo', write_atmo)
    calc_amp = define_param(user_config, 'GENERAL', 'calc_amp', calc_amp)

    output_id = define_param(user_config, 'GENERAL', 'output_id', output_id)

    max_alt = define_param(user_config, 'GENERAL', 'max_alt', max_alt)
    max_rng = define_param(user_config, 'GENERAL', 'max_rng', max_rng)

    min_lat = define_param(user_config, 'GENERAL', 'min_lat', min_lat)
    max_lat = define_param(user_config, 'GENERAL', 'max_lat', max_lat)
    min_lon = define_param(user_config, 'GENERAL', 'min_lon', min_lon)
    max_lon = define_param(user_config, 'GENERAL', 'max_lon', max_lon)

    min_ds = define_param(user_config, 'GENERAL', 'min_ds', min_ds)
    max_ds = define_param(user_config, 'GENERAL', 'max_ds', max_ds)
    max_s = define_param(user_config, 'GENERAL', 'max_s', max_s)

    z_grnd = define_param(user_config, 'GENERAL', 'z_grnd', z_grnd)
    topo_file = define_param(user_config, 'GENERAL', 'topo_file', topo_file)
    topo_bl_wind = define_param(user_config, 'GENERAL', 'topo_bl_wind', topo_bl_wind)

    cpu_cnt = define_param(user_config, 'GENERAL', 'cpu_cnt', cpu_cnt)
    cpu_cnt = int(cpu_cnt)

    # Build run command
    if cpu_cnt < 2:
        command = bin_path + "infraga-sph"
    else:
        command = "mpirun -np " + str(cpu_cnt) + " " + bin_path + "infraga-accel-sph"

    if atmo_prefix is not None:
        command = command + "-rngdep"

    if atmo_file is not None:
        command = command + " -prop " + atmo_file
    elif atmo_prefix is not None and grid_lats is not None and grid_lons is not None:
        command = command + " -prop " + atmo_prefix + " " + grid_lats + " " + grid_lons
    else:
        click.echo("Simulation requires either an '--atmo-file' or --atmo-prefix' with grid info (--grid-lats and --grid-lons)")
        return 0

    # Set parameters
    command = set_param(command, incl_min, "incl_min")
    command = set_param(command, incl_max, "incl_max")
    command = set_param(command, incl_step, "incl_step")

    command = set_param(command, az_min, "az_min")
    command = set_param(command, az_max, "az_max")
    command = set_param(command, az_step, "az_step")

    command = set_param(command, inclination, "inclination")
    command = set_param(command, azimuth, "azimuth")

    command = set_param(command, bounces, "bounces")

    command = set_param(command, src_lat, "src_lat")
    command = set_param(command, src_lon, "src_lon")
    command = set_param(command, src_alt, "src_alt")

    command = set_param(command, az_step, "az_step")

    if write_rays is not None:
        command = set_param(command, str(write_rays), "write_rays")

    command = set_param(command, freq, "freq")
    command = set_param(command, abs_coeff, "abs_coeff")

    if write_atmo is not None:
        command = set_param(command, str(write_atmo), "write_atmo")

    if calc_amp is not None:
        command = set_param(command, str(calc_amp), "calc_amp")

    command = set_param(command, output_id, "output_id")

    command = set_param(command, max_alt, "max_alt")
    command = set_param(command, max_rng, "max_rng")

    command = set_param(command, min_lat, "min_lat")
    command = set_param(command, max_lat, "max_lat")
    command = set_param(command, min_lon, "min_lon")
    command = set_param(command, max_lon, "max_lon")

    command = set_param(command, min_ds, "min_ds")
    command = set_param(command, max_ds, "max_ds")
    command = set_param(command, max_s, "max_s")

    command = set_param(command, z_grnd, "z_grnd")
    command = set_param(command, topo_file, "topo_file")
    if topo_file is not None:
        if write_topo is not None:
            command = set_param(command, str(write_topo), "write_topo")
        if topo_bl_wind is not None:
            command = set_param(command, str(topo_bl_wind), "topo_bl_wind")

    # print(command)
    subprocess.call(command, shell=True)




@click.command('eigenray', short_help="Run eigenray search for specific source-receiver")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)
@click.option("--atmo-prefix", help="Atmosphere file prefix (range dependent)", default=None)
@click.option("--grid-lats", help="Atmosphere grid latitudes file", default=None)
@click.option("--grid-lons", help="Atmosphere grid longitudes file", default=None)

@click.option("--incl-min", help="Minimum inclination angle (rel. horizontal)", default=None)
@click.option("--incl-max", help="Maximum inclination angle (rel. horizontal)", default=None)

@click.option("--bnc-min", help="Minimum number of ground reflections (bounces)", default=None)
@click.option("--bnc-max", help="Maximum number of ground reflections (bounces)", default=None)
@click.option("--bounces", help="Number of ground reflections (bounces) to consider", default=None)

@click.option("--src-lat", help="Source latitude", default=None)
@click.option("--src-lon", help="Source longitude", default=None)
@click.option("--src-alt", help="Source altitude", default=None)
@click.option("--rcvr-lat", help="Receiver latitude", default=None)
@click.option("--rcvr-lon", help="Receiver longitude", default=None)

@click.option("--verbose", help="Option to print output to screen", default=None, type=bool)
@click.option("--iterations", help="Maximum Levenberg-Marquardt iterations", default=None)
@click.option("--damping", help="Damping factor in Levenberg-Marquardt algorithm", default=None)
@click.option("--tolerance", help="Receiver accuracy threshold", default=None)
@click.option("--az-dev-lim", help="Azimuth deviation limit", default=None)
@click.option("--incl_step-min", help="Minimum inclination step", default=None)
@click.option("--incl_step-max", help="Maximum inclination step", default=None)

@click.option("--freq", help="Frequency for Sutherland-Bass losses", default=None)
@click.option("--abs-coeff", help="Scaling coefficient for Sutherland-Bass losses", default=None)
@click.option("--z-grnd", help="Ground elevation (flat ground, km rel. sea level)", default=None)
@click.option("--write-atmo", help="Option to write atmosphere data (for QC)", default=None)
@click.option("--output-id", help="User specified output file path", default=None)
@click.option("--max-alt", help="Maximum altitude for ray calculation", default=None)
@click.option("--max-rng", help="Maximum range for ray calculation", default=None)

@click.option("--min-lat", help="Minimum latitude for user defined bounds", default=None)
@click.option("--max-lat", help="Maximum latitude for user defined bounds", default=None)
@click.option("--min-lon", help="Minimum longitude for user defined bounds", default=None)
@click.option("--max-lon", help="Maximum longitude for user defined bounds", default=None)

@click.option("--min-ds", help="Minimum step size (near-ground) in RK4 solver", default=None)
@click.option("--max-ds", help="Maximum step size in RK4 solver", default=None)
@click.option("--max-s", help="Maximum ray length between bounces", default=None)

@click.option("--topo-file", help="Terrain file", default=None)
@click.option("--topo-BL-wind", help="Use terrain corrected boundary layer winds", default=None, type=bool)
@click.option("--cpu-cnt", help="Number of CPUs to use in analysis", default=1)
def run_sph_eig(config_file, atmo_file, atmo_prefix, grid_lats, grid_lons, incl_min, incl_max, bnc_min, bnc_max, bounces, 
                src_lat, src_lon, src_alt, rcvr_lat, rcvr_lon, verbose, iterations, damping, tolerance, az_dev_lim, 
                incl_step_min, incl_step_max, freq, abs_coeff, z_grnd, write_atmo, output_id, max_alt, max_rng, 
                min_lat, max_lat, min_lon, max_lon, min_ds, max_ds, max_s, topo_file, topo_bl_wind, cpu_cnt):
    '''
    Run spherical atmospheric layer eigenray analysis to identify propagation paths connecting a specific source-receiver geometry

    \b
    Examples:
    \t infraga sph eigenray --atmo-file ToyAtmo.met --src-lat 30.0 --src-lon -100.0 --rcvr-lat 30.25 --rcvr-lon -104.25 --bnc-max 1 --verbose True
    \t infraga sph eigenray --atmo-file ToyAtmo.met --config-file example.cnfg

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

    # Set propagation specific parameters
    incl_min = define_param(user_config, 'EIGENRAY', 'incl_min', incl_min)
    incl_max = define_param(user_config, 'EIGENRAY', 'incl_max', incl_max)

    bnc_min = define_param(user_config, 'EIGENRAY', 'bnc_min', bnc_min)
    bnc_max = define_param(user_config, 'EIGENRAY', 'bnc_max', bnc_max)
    bounces = define_param(user_config, 'EIGENRAY', 'bounces', bounces)

    src_lat = define_param(user_config, 'EIGENRAY', 'src_lat', src_lat)
    src_lon = define_param(user_config, 'EIGENRAY', 'src_lon', src_lon)
    src_alt = define_param(user_config, 'EIGENRAY', 'src_alt', src_alt)

    rcvr_lat = define_param(user_config, 'EIGENRAY', 'rcvr_lat', rcvr_lat)
    rcvr_lon = define_param(user_config, 'EIGENRAY', 'rcvr_lon', rcvr_lon)

    verbose = define_param(user_config, 'EIGENRAY', 'verbose', verbose)
    iterations = define_param(user_config, 'EIGENRAY', 'iterations', iterations)
    damping = define_param(user_config, 'EIGENRAY', 'damping', damping)
    tolerance = define_param(user_config, 'EIGENRAY', 'tolerance', tolerance)
    az_dev_lim = define_param(user_config, 'EIGENRAY', 'az_dev_lim', az_dev_lim)
    incl_step_min = define_param(user_config, 'EIGENRAY', 'incl_step_min', incl_step_min)
    incl_step_max = define_param(user_config, 'EIGENRAY', 'incl_step_max', incl_step_max)


    # Set general parameters
    freq = define_param(user_config, 'GENERAL', 'freq', freq)
    abs_coeff = define_param(user_config, 'GENERAL', 'abs_coeff', abs_coeff)

    write_atmo = define_param(user_config, 'GENERAL', 'write_atmo', write_atmo)
    output_id = define_param(user_config, 'GENERAL', 'output_id', output_id)

    max_alt = define_param(user_config, 'GENERAL', 'max_alt', max_alt)
    max_rng = define_param(user_config, 'GENERAL', 'max_rng', max_rng)

    min_lat = define_param(user_config, 'GENERAL', 'min_lat', min_lat)
    max_lat = define_param(user_config, 'GENERAL', 'max_lat', max_lat)
    min_lon = define_param(user_config, 'GENERAL', 'min_lon', min_lon)
    max_lon = define_param(user_config, 'GENERAL', 'max_lon', max_lon)

    min_ds = define_param(user_config, 'GENERAL', 'min_ds', min_ds)
    max_ds = define_param(user_config, 'GENERAL', 'max_ds', max_ds)
    max_s = define_param(user_config, 'GENERAL', 'max_s', max_s)

    z_grnd = define_param(user_config, 'GENERAL', 'z_grnd', z_grnd)
    topo_file = define_param(user_config, 'GENERAL', 'topo_file', topo_file)
    topo_bl_wind = define_param(user_config, 'GENERAL', 'topo_bl_wind', topo_bl_wind)

    cpu_cnt = define_param(user_config, 'GENERAL', 'cpu_cnt', cpu_cnt)
    cpu_cnt = int(cpu_cnt)

    # Build run command
    if cpu_cnt < 2:
        command = bin_path + "infraga-sph"
    else:
        command = "mpirun -np " + str(cpu_cnt) + " " + bin_path + "infraga-accel-sph"

    if atmo_prefix is not None:
        command = command + "-rngdep"

    if atmo_file is not None:
        command = command + " -eig_search " + atmo_file
    elif atmo_prefix is not None and grid_lats is not None and grid_lons is not None:
        command = command + " -eig_search " + atmo_prefix + " " + grid_lats + " " + grid_lons
    else:
        click.echo("Simulation requires either an '--atmo-file' or --atmo-prefix' with grid info (--grid-lats and --grid-lons)")
        return 0

    # Set parameters
    command = set_param(command, incl_min, "incl_min")
    command = set_param(command, incl_max, "incl_max")

    command = set_param(command, bnc_min, "bnc_min")
    command = set_param(command, bnc_max, "bnc_max")
    command = set_param(command, bounces, "bounces")

    command = set_param(command, src_lat, "src_lat")
    command = set_param(command, src_lon, "src_lon")
    command = set_param(command, src_alt, "src_alt")

    command = set_param(command, rcvr_lat, "rcvr_lat")
    command = set_param(command, rcvr_lon, "rcvr_lon")

    if verbose is not None:
        command = set_param(command, str(verbose), "verbose")

    command = set_param(command, iterations, "iterations")
    command = set_param(command, damping, "damping")
    command = set_param(command, tolerance, "tolerance")
    command = set_param(command, az_dev_lim, "az_dev_lim")
    command = set_param(command, incl_step_min, "incl_step_min")
    command = set_param(command, incl_step_max, "incl_step_max")

    command = set_param(command, freq, "freq")
    command = set_param(command, abs_coeff, "abs_coeff")

    if write_atmo is not None:
        command = set_param(command, str(write_atmo), "write_atmo")

    command = set_param(command, output_id, "output_id")

    command = set_param(command, max_alt, "max_alt")
    command = set_param(command, max_rng, "max_rng")

    command = set_param(command, min_lat, "min_lat")
    command = set_param(command, max_lat, "max_lat")
    command = set_param(command, min_lon, "min_lon")
    command = set_param(command, max_lon, "max_lon")

    command = set_param(command, min_ds, "min_ds")
    command = set_param(command, max_ds, "max_ds")
    command = set_param(command, max_s, "max_s")

    command = set_param(command, z_grnd, "z_grnd")
    command = set_param(command, topo_file, "topo_file")
    if topo_file is not None:
        if topo_bl_wind is not None:
            command = set_param(command, str(topo_bl_wind), "topo_bl_wind")

    # print(command)
    subprocess.call(command, shell=True)




@click.command('wnl-wvfrm', short_help="Run weakly non-linear waveform calculation")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)
@click.option("--atmo-prefix", help="Atmosphere file prefix (range dependent)", default=None)
@click.option("--grid-lats", help="Atmosphere grid latitudes file", default=None)
@click.option("--grid-lons", help="Atmosphere grid longitudes file", default=None)

@click.option("--inclination", help="Inclination angle (rel. horizontal)", default=None)
@click.option("--azimuth", help="Azimuth angle (clockwise rel. N)", default=None)
@click.option("--bounces", help="Number of ground reflections (bounces) to consider", default=None)

@click.option("--src-lat", help="Source latitude", default=None)
@click.option("--src-lon", help="Source longitude", default=None)
@click.option("--src-alt", help="Source altitude", default=None)

@click.option("--write-ray", help="Option to write ray path info", default=None, type=bool)

@click.option("--wvfrm-file", help="File containing reference waveform", default=None)
@click.option("--wvfrm-opt", help="Waveform option ('impulse', 'Uwave', 'Nwave')", default=None)

@click.option("--wvfrm-p0", help="Waveform peak overpressure [Pa]", default=None)
@click.option("--wvfrm-t0", help="Waveform positive phase duration [s]", default=None)
@click.option("--wvfrm-alpha", help="Waveform shaping parameter", default=None)
@click.option("--wvfrm-ref", help="Waveform reference distance", default=None)
@click.option("--wvfrm-out-step", help="Waveform output step along ray path", default=None)
@click.option("--wvfrm-yield", help="Eq. TNT yield (kg) for an explosive source", default=None)

@click.option("--wvfrm-ds", help="Burgers equation solver resolution [km]", default=None)
@click.option("--wvfrm-len", help="Waveform sample count (default = 2e13)", default=None)

@click.option("--freq", help="Frequency for Sutherland-Bass losses", default=None)
@click.option("--abs-coeff", help="Scaling coefficient for Sutherland-Bass losses", default=None)
@click.option("--z-grnd", help="Ground elevation (flat ground, km rel. sea level)", default=None)
@click.option("--write-atmo", help="Option to write atmosphere data (for QC)", default=None)
@click.option("--output-id", help="User specified output file path", default=None)
@click.option("--max-alt", help="Maximum altitude for ray calculation", default=None)
@click.option("--max-rng", help="Maximum range for ray calculation", default=None)

@click.option("--min-lat", help="Minimum latitude for user defined bounds", default=None)
@click.option("--max-lat", help="Maximum latitude for user defined bounds", default=None)
@click.option("--min-lon", help="Minimum longitude for user defined bounds", default=None)
@click.option("--max-lon", help="Maximum longitude for user defined bounds", default=None)

@click.option("--min-ds", help="Minimum step size (near-ground) in RK4 solver", default=None)
@click.option("--max-ds", help="Maximum step size in RK4 solver", default=None)
@click.option("--max-s", help="Maximum ray length between bounces", default=None)

@click.option("--topo-file", help="Terrain file", default=None)
@click.option("--topo-BL-wind", help="Use terrain corrected boundary layer winds", default=None, type=bool)
def run_sph_wvfrm(config_file, atmo_file, atmo_prefix, grid_lats, grid_lons, inclination, azimuth, bounces, src_lat, src_lon, 
                src_alt, write_ray, wvfrm_file, wvfrm_opt, wvfrm_p0, wvfrm_t0, wvfrm_alpha, wvfrm_ref, wvfrm_out_step,
                wvfrm_yield, wvfrm_ds, wvfrm_len, freq, abs_coeff, z_grnd, write_atmo, output_id, max_alt, max_rng, max_lat, min_lat, max_lon, 
                min_lon, min_ds, max_ds, max_s, topo_file, topo_bl_wind):
    '''
        Run weakly non-linear waveform analysis along a spherical atmospheric layer ray path.

    \b
    Examples:
    \t infraga sph wnl-wvfrm --atmo-file ToyAtmo.met --src-lat 30.0 --src-lon -100.0 --inclination 4.1314876 --azimuth -84.969455 --bounces 1 --wvfrm-p0 500.0
    \t infraga sph wnl-wvfrm --atmo-file ToyAtmo.met --config-file example.cnfg
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

    # Set propagation specific parameters
    inclination = define_param(user_config, 'WAVEFORM', 'inclination', inclination)
    azimuth = define_param(user_config, 'WAVEFORM', 'azimuth', azimuth)
    bounces = define_param(user_config, 'WAVEFORM', 'bounces', bounces)

    src_lat = define_param(user_config, 'WAVEFORM', 'src_lat', src_lat)
    src_lon = define_param(user_config, 'WAVEFORM', 'src_lon', src_lon)
    src_alt = define_param(user_config, 'WAVEFORM', 'src_alt', src_alt)

    write_ray = define_param(user_config, 'WAVEFORM', 'write_ray', write_ray)

    wvfrm_file = define_param(user_config, 'WAVEFORM', 'wvfrm_file', wvfrm_file)
    wvfrm_opt = define_param(user_config, 'WAVEFORM', 'wvfrm_opt', wvfrm_opt)
    wvfrm_p0 = define_param(user_config, 'WAVEFORM', 'wvfrm_p0', wvfrm_p0)
    wvfrm_t0 = define_param(user_config, 'WAVEFORM', 'wvfrm_t0', wvfrm_t0)
    wvfrm_alpha = define_param(user_config, 'WAVEFORM', 'wvfrm_alpha', wvfrm_alpha)
    wvfrm_ref = define_param(user_config, 'WAVEFORM', 'wvfrm_ref', wvfrm_ref)
    wvfrm_out_step = define_param(user_config, 'WAVEFORM', 'wvfrm_out_step', wvfrm_out_step)
    wvfrm_yield = define_param(user_config, 'WAVEFORM', 'wvfrm_yield', wvfrm_yield)

    wvfrm_ds = define_param(user_config, 'WAVEFORM', 'wvfrm_ds', wvfrm_ds)
    wvfrm_len = define_param(user_config, 'WAVEFORM', 'wvfrm_len', wvfrm_len)

    # Set general parameters
    freq = define_param(user_config, 'GENERAL', 'freq', freq)
    abs_coeff = define_param(user_config, 'GENERAL', 'abs_coeff', abs_coeff)

    write_atmo = define_param(user_config, 'GENERAL', 'write_atmo', write_atmo)
    output_id = define_param(user_config, 'GENERAL', 'output_id', output_id)

    max_alt = define_param(user_config, 'GENERAL', 'max_alt', max_alt)
    max_rng = define_param(user_config, 'GENERAL', 'max_rng', max_rng)

    min_lat = define_param(user_config, 'GENERAL', 'min_lat', min_lat)
    max_lat = define_param(user_config, 'GENERAL', 'max_lat', max_lat)
    min_lon = define_param(user_config, 'GENERAL', 'min_lon', min_lon)
    max_lon = define_param(user_config, 'GENERAL', 'max_lon', max_lon)

    min_ds = define_param(user_config, 'GENERAL', 'min_ds', min_ds)
    max_ds = define_param(user_config, 'GENERAL', 'max_ds', max_ds)
    max_s = define_param(user_config, 'GENERAL', 'max_s', max_s)

    z_grnd = define_param(user_config, 'GENERAL', 'z_grnd', z_grnd)
    topo_file = define_param(user_config, 'GENERAL', 'topo_file', topo_file)
    topo_bl_wind = define_param(user_config, 'GENERAL', 'topo_bl_wind', topo_bl_wind)
    
    # Build run command
    command = bin_path + "infraga-sph"
    
    if atmo_prefix is not None:
        command = command + "-rngdep"

    if atmo_file is not None:
        command = command + " -wnl_wvfrm " + atmo_file
    elif atmo_prefix is not None and grid_lats is not None and grid_lons is not None:
        command = command + " -wnl_wvfrm " + atmo_prefix + " " + grid_lats + " " + grid_lons
    else:
        click.echo("Simulation requires either an '--atmo-file' or --atmo-prefix' with grid info (--grid-lats and --grid-lons)")
        return 0

    # Set parameters
    command = set_param(command, inclination, "inclination")
    command = set_param(command, azimuth, "azimuth")
    command = set_param(command, bounces, "bounces")

    command = set_param(command, src_lat, "src_lat")
    command = set_param(command, src_lon, "src_lon")
    command = set_param(command, src_alt, "src_alt")

    if write_ray is not None:
        command = set_param(command, str(write_ray), "write_ray")

    command = set_param(command, wvfrm_file, "wvfrm_file")
    if wvfrm_yield is not None:
        # use atmosphere to define ambient pressure (1 mbar = 0.1 kPa)
        atmo = np.loadtxt(atmo_file)

        if src_alt is not None:
            src_ht = float(src_alt)
        elif z_grnd is not None:
            src_ht = float(z_grnd)
        else:
            src_ht = 0.0

        p_ambient = atmo[np.argmin(abs(atmo[:, 0] - src_ht)), 5] * 0.1
        T_ambient = atmo[np.argmin(abs(atmo[:, 0] - src_ht)), 1]

        wvfrm_opt = 'impulse'
        wvfrm_ref = str(0.035 * float(wvfrm_yield)**(1.0 / 3.0))
        wvfrm_p0 = str(kg_op(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient))
        wvfrm_t0 = str(kg_ppd(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient))
        wvfrm_alpha = '0.01'

    command = set_param(command, wvfrm_opt, "wvfrm_opt")
    command = set_param(command, wvfrm_p0, "wvfrm_p0")
    command = set_param(command, wvfrm_t0, "wvfrm_t0")
    command = set_param(command, wvfrm_alpha, "wvfrm_alpha")
    command = set_param(command, wvfrm_ref, "wvfrm_ref")
    command = set_param(command, wvfrm_out_step, "wvfrm_out_step")
        
    command = set_param(command, wvfrm_ds, "wvfrm_ds")
    command = set_param(command, wvfrm_len, "wvfrm_len")

    command = set_param(command, freq, "freq")
    command = set_param(command, abs_coeff, "abs_coeff")

    if write_atmo is not None:
        command = set_param(command, str(write_atmo), "write_atmo")

    command = set_param(command, output_id, "output_id")

    command = set_param(command, max_alt, "max_alt")
    command = set_param(command, max_rng, "max_rng")

    command = set_param(command, min_lat, "min_lat")
    command = set_param(command, max_lat, "max_lat")
    command = set_param(command, min_lon, "min_lon")
    command = set_param(command, max_lon, "max_lon")

    command = set_param(command, min_ds, "min_ds")
    command = set_param(command, max_ds, "max_ds")
    command = set_param(command, max_s, "max_s")

    command = set_param(command, z_grnd, "z_grnd")
    command = set_param(command, topo_file, "topo_file")
    if topo_file is not None:
        if topo_bl_wind is not None:
            command = set_param(command, str(topo_bl_wind), "topo_bl_wind")

    # print(command)
    subprocess.call(command, shell=True)



@click.command('eig_wvfrm', short_help="Run eigenray search and compute waveform contributions")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)
@click.option("--atmo-prefix", help="Atmosphere file prefix (range dependent)", default=None)
@click.option("--grid-lats", help="Atmosphere grid latitudes file", default=None)
@click.option("--grid-lons", help="Atmosphere grid longitudes file", default=None)

@click.option("--incl-min", help="Minimum inclination angle (rel. horizontal)", default=None)
@click.option("--incl-max", help="Maximum inclination angle (rel. horizontal)", default=None)

@click.option("--bnc-min", help="Minimum number of ground reflections (bounces)", default=None)
@click.option("--bnc-max", help="Maximum number of ground reflections (bounces)", default=None)
@click.option("--bounces", help="Specific number of ground reflections (bounces)", default=None)

@click.option("--src-lat", help="Source latitude", default=None)
@click.option("--src-lon", help="Source longitude", default=None)
@click.option("--src-alt", help="Source altitude", default=None)
@click.option("--rcvr-lat", help="Receiver latitude", default=None)
@click.option("--rcvr-lon", help="Receiver longitude", default=None)

@click.option("--wvfrm-file", help="File containing reference waveform", default=None)
@click.option("--wvfrm-opt", help="Waveform option ('impulse', 'Uwave', 'Nwave')", default=None)

@click.option("--verbose", help="Option to print output to screen", default=None, type=bool)
@click.option("--iterations", help="Maximum Levenberg-Marquardt iterations", default=None)
@click.option("--damping", help="Damping factor in Levenberg-Marquardt algorithm", default=None)
@click.option("--tolerance", help="Receiver accuracy threshold", default=None)
@click.option("--az-dev-lim", help="Azimuth deviation limit", default=None)
@click.option("--incl_step-min", help="Minimum inclination step", default=None)
@click.option("--incl_step-max", help="Maximum inclination step", default=None)

@click.option("--wvfrm-p0", help="Waveform peak overpressure [Pa]", default=None)
@click.option("--wvfrm-t0", help="Waveform positive phase duration [s]", default=None)
@click.option("--wvfrm-alpha", help="Waveform shaping parameter", default=None)
@click.option("--wvfrm-ref", help="Waveform reference distance", default=None)
@click.option("--wvfrm-out-step", help="Waveform output step along ray path", default=None)
@click.option("--wvfrm-yield", help="Eq. TNT yield (kg) for an explosive source", default=None)

@click.option("--wvfrm-ds", help="Burgers equation solver resolution [km]", default=None)
@click.option("--wvfrm-len", help="Waveform sample count (default = 2e13)", default=None)

@click.option("--freq", help="Frequency for Sutherland-Bass losses", default=None)
@click.option("--abs-coeff", help="Scaling coefficient for Sutherland-Bass losses", default=None)
@click.option("--z-grnd", help="Ground elevation (km rel. sea level)", default=None)
@click.option("--write-atmo", help="Option to write atmosphere data (for QC)", default=None)
@click.option("--output-id", help="User specified output file path", default=None)
@click.option("--max-alt", help="Maximum altitude for ray calculation", default=None)
@click.option("--max-rng", help="Maximum range for ray calculation", default=None)

@click.option("--min-lat", help="Minimum latitude for user defined bounds", default=None)
@click.option("--max-lat", help="Maximum latitude for user defined bounds", default=None)
@click.option("--min-lon", help="Minimum longitude for user defined bounds", default=None)
@click.option("--max-lon", help="Maximum longitude for user defined bounds", default=None)

@click.option("--min-ds", help="Minimum step size (near-ground) in RK4 solver", default=None)
@click.option("--max-ds", help="Maximum step size in RK4 solver", default=None)
@click.option("--max-s", help="Maximum ray length between bounces", default=None)

@click.option("--topo-file", help="Terrain file", default=None)
@click.option("--topo-BL-wind", help="Use terrain corrected boundary layer winds", default=None, type=bool)
@click.option("--cpu-cnt", help="Number of CPUs to use in analysis", default=1)

@click.option("--keep-eig-arrivals", help="Keep eigenray arrivals", default=False)
def run_sph_eig_wvfrm(config_file, atmo_file, atmo_prefix, grid_lats, grid_lons, incl_min, incl_max, bnc_min, bnc_max, bounces, 
                src_lat, src_lon, src_alt, rcvr_lat, rcvr_lon, verbose, iterations, damping, tolerance, az_dev_lim, 
                incl_step_min, incl_step_max, wvfrm_file, wvfrm_opt, wvfrm_p0, wvfrm_t0, wvfrm_alpha, wvfrm_ref, wvfrm_out_step,
                wvfrm_yield, wvfrm_ds, wvfrm_len, freq, abs_coeff, z_grnd, write_atmo, output_id, max_alt, max_rng, min_lat, 
                max_lat, min_lon, max_lon, min_ds, max_ds, max_s, topo_file, topo_bl_wind, cpu_cnt, keep_eig_arrivals):
    '''
    Run spherical atmospheric layer eigenray analysis to identify propagation paths connecting a specific source-receiver geometr yand then compute weakly-nonlinear waveform predictions for each eigenray

    \b
    Examples:
    \t infraga sph eig_wvfrm --atmo-file ToyAtmo.met --src-lat 30.0 --src-lon -100.0 --rcvr-lat 30.25 --rcvr-lon -104.25 --bnc-max 1 --verbose True --keep-eig-arrivals True --wvfrm-yield 10e3
    \t infraga sph eig_wvfrm --atmo-file ToyAtmo.met --config-file example.cnfg

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

    # Set eigenray specific parameters
    incl_min = define_param(user_config, 'EIGENRAY', 'incl_min', incl_min)
    incl_max = define_param(user_config, 'EIGENRAY', 'incl_max', incl_max)

    bnc_min = define_param(user_config, 'EIGENRAY', 'bnc_min', bnc_min)
    bnc_max = define_param(user_config, 'EIGENRAY', 'bnc_max', bnc_max)
    bounces = define_param(user_config, 'EIGENRAY', 'bounces', bounces)

    src_lat = define_param(user_config, 'EIGENRAY', 'src_lat', src_lat)
    src_lon = define_param(user_config, 'EIGENRAY', 'src_lon', src_lon)
    src_alt = define_param(user_config, 'EIGENRAY', 'src_alt', src_alt)

    rcvr_lat = define_param(user_config, 'EIGENRAY', 'rcvr_lat', rcvr_lat)
    rcvr_lon = define_param(user_config, 'EIGENRAY', 'rcvr_lon', rcvr_lon)

    verbose = define_param(user_config, 'EIGENRAY', 'verbose', verbose)
    iterations = define_param(user_config, 'EIGENRAY', 'iterations', iterations)
    damping = define_param(user_config, 'EIGENRAY', 'damping', damping)
    tolerance = define_param(user_config, 'EIGENRAY', 'tolerance', tolerance)
    az_dev_lim = define_param(user_config, 'EIGENRAY', 'az_dev_lim', az_dev_lim)
    incl_step_min = define_param(user_config, 'EIGENRAY', 'incl_step_min', incl_step_min)
    incl_step_max = define_param(user_config, 'EIGENRAY', 'incl_step_max', incl_step_max)

    # Set waveform specific parameters 
    wvfrm_file = define_param(user_config, 'WAVEFORM', 'wvfrm_file', wvfrm_file)
    wvfrm_opt = define_param(user_config, 'WAVEFORM', 'wvfrm_opt', wvfrm_opt)
    wvfrm_p0 = define_param(user_config, 'WAVEFORM', 'wvfrm_p0', wvfrm_p0)
    wvfrm_t0 = define_param(user_config, 'WAVEFORM', 'wvfrm_t0', wvfrm_t0)
    wvfrm_alpha = define_param(user_config, 'WAVEFORM', 'wvfrm_alpha', wvfrm_alpha)
    wvfrm_ref = define_param(user_config, 'WAVEFORM', 'wvfrm_ref', wvfrm_ref)
    wvfrm_out_step = define_param(user_config, 'WAVEFORM', 'wvfrm_out_step', wvfrm_out_step)
    wvfrm_yield = define_param(user_config, 'WAVEFORM', 'wvfrm_yield', wvfrm_yield)

    wvfrm_ds = define_param(user_config, 'WAVEFORM', 'wvfrm_ds', wvfrm_ds)
    wvfrm_len = define_param(user_config, 'WAVEFORM', 'wvfrm_len', wvfrm_len)
    
    # Set general parameters
    freq = define_param(user_config, 'GENERAL', 'freq', freq)
    abs_coeff = define_param(user_config, 'GENERAL', 'abs_coeff', abs_coeff)

    write_atmo = define_param(user_config, 'GENERAL', 'write_atmo', write_atmo)
    output_id = define_param(user_config, 'GENERAL', 'output_id', output_id)

    max_alt = define_param(user_config, 'GENERAL', 'max_alt', max_alt)
    max_rng = define_param(user_config, 'GENERAL', 'max_rng', max_rng)

    min_lat = define_param(user_config, 'GENERAL', 'min_lat', min_lat)
    max_lat = define_param(user_config, 'GENERAL', 'max_lat', max_lat)
    min_lon = define_param(user_config, 'GENERAL', 'min_lon', min_lon)
    max_lon = define_param(user_config, 'GENERAL', 'max_lon', max_lon)

    min_ds = define_param(user_config, 'GENERAL', 'min_ds', min_ds)
    max_ds = define_param(user_config, 'GENERAL', 'max_ds', max_ds)
    max_s = define_param(user_config, 'GENERAL', 'max_s', max_s)

    z_grnd = define_param(user_config, 'GENERAL', 'z_grnd', z_grnd)
    topo_file = define_param(user_config, 'GENERAL', 'topo_file', topo_file)
    topo_bl_wind = define_param(user_config, 'GENERAL', 'topo_bl_wind', topo_bl_wind)

    cpu_cnt = define_param(user_config, 'GENERAL', 'cpu_cnt', cpu_cnt)
    cpu_cnt = int(cpu_cnt)

    # Check if eigenray analysis is already done
    if output_id is not None:
        eig_arrivals_file = output_id + ".arrivals.dat"
        profile_id = output_id
    else:
        if atmo_prefix is not None:
            eig_arrivals_file = atmo_prefix + ".arrivals.dat"
            profile_id = atmo_prefix
        else:
            eig_arrivals_file = os.path.splitext(atmo_file)[0] + ".arrivals.dat"
            profile_id = os.path.splitext(atmo_file)[0]
    
    if os.path.isfile(eig_arrivals_file):
        print("Eigenray results found.  Skipping to waveform calculation...")
    else:
        print("Running eigenray analysis...")

        if cpu_cnt < 2:
            command = bin_path + "infraga-sph"
        else:
            command = "mpirun -np " + str(cpu_cnt) + " " + bin_path + "infraga-accel-sph"

        if atmo_prefix is not None:
            command = command + "-rngdep"

        if atmo_file is not None:
            command = command + " -eig_search " + atmo_file
        elif atmo_prefix is not None and grid_lats is not None and grid_lons is not None:
            command = command + " -eig_search " + atmo_prefix + " " + grid_lats + " " + grid_lons
        else:
            click.echo("Simulation requires either an '--atmo-file' or --atmo-prefix' with grid info (--grid-lats and --grid-lons)")
            return 0

        # Set parameters
        command = set_param(command, incl_min, "incl_min")
        command = set_param(command, incl_max, "incl_max")

        command = set_param(command, bnc_min, "bnc_min")
        command = set_param(command, bnc_max, "bnc_max")
        command = set_param(command, bounces, "bounces")

        command = set_param(command, src_lat, "src_lat")
        command = set_param(command, src_lon, "src_lon")
        command = set_param(command, src_alt, "src_alt")

        command = set_param(command, rcvr_lat, "rcvr_lat")
        command = set_param(command, rcvr_lon, "rcvr_lon")

        if verbose is not None:
            command = set_param(command, str(verbose), "verbose")

        command = set_param(command, iterations, "iterations")
        command = set_param(command, damping, "damping")
        command = set_param(command, tolerance, "tolerance")
        command = set_param(command, az_dev_lim, "az_dev_lim")
        command = set_param(command, incl_step_min, "incl_step_min")
        command = set_param(command, incl_step_max, "incl_step_max")

        command = set_param(command, freq, "freq")
        command = set_param(command, abs_coeff, "abs_coeff")

        if write_atmo is not None:
            command = set_param(command, str(write_atmo), "write_atmo")

        command = set_param(command, output_id, "output_id")

        command = set_param(command, max_alt, "max_alt")
        command = set_param(command, max_rng, "max_rng")

        command = set_param(command, min_lat, "min_lat")
        command = set_param(command, max_lat, "max_lat")
        command = set_param(command, min_lon, "min_lon")
        command = set_param(command, max_lon, "max_lon")

        command = set_param(command, min_ds, "min_ds")
        command = set_param(command, max_ds, "max_ds")
        command = set_param(command, max_s, "max_s")

        command = set_param(command, z_grnd, "z_grnd")
        command = set_param(command, topo_file, "topo_file")
        if topo_file is not None:
            if topo_bl_wind is not None:
                command = set_param(command, str(topo_bl_wind), "topo_bl_wind")

        # print(command)
        subprocess.call(command, shell=True)


    # Analyze eigenray results and compute waveform contributions...
    print("Computing waveform...")
    results_id = eig_arrivals_file[:-13]

    eig_results = np.loadtxt(eig_arrivals_file)
    if len(eig_results) > 0:
        eig_results = np.atleast_2d(eig_results)
        os.system("rm " + results_id + ".eigenray*")

        mask = np.ones(len(eig_results), dtype=bool)
        for n, line in enumerate(eig_results):
            if mask[n]:   
                for m in range(n + 1, len(eig_results)):   
                    arrivals_diff  = (line[0] - eig_results[m][0])**2
                    arrivals_diff += (line[1] - eig_results[m][1])**2
                    arrivals_diff += (line[5] - eig_results[m][5])**2
                    arrivals_diff += (line[7] - eig_results[m][7])**2
                    if arrivals_diff < 0.1:
                        mask[m] = False
        eig_results = eig_results[mask]

        # Create waveform interpolations
        file_out = open(results_id + ".eigenrays.dat", 'w')
        print("# 'infraga sph eig_wvfrm' eigenray results", '\n#', file=file_out)
        if atmo_prefix is not None:
            print("# 	profile: " + atmo_prefix, file=file_out)
        else:
            print("# 	profile: " + atmo_file, file=file_out)
        print("# 	source location (lat, lon, alt): " + str(src_lat) + ", " + str(src_lon) + ", " + str(src_alt) , file=file_out)
        print("# 	receiver location (lat, lon, alt): " + str(rcvr_lat) + ", " + str(rcvr_lon) + ", 0.0", file=file_out)
        print("# 	inclination range: " + str(incl_min) + ", " + str(incl_max), file=file_out)
        print("#    inclination step max:", incl_step_max, file=file_out)
        print("# 	bounces: " + str(bnc_min) + ", " + str(bnc_max), file=file_out)
        if topo_file is not None:
            print("# 	terrain file: " + str(topo_file), file=file_out)        
        else:
            print("# 	ground elevation: " + str(z_grnd), file=file_out)
        print("# 	damping:", damping, file=file_out)
        print("# 	range max:", max_rng, file=file_out)

        wvfrms = []
        t_lims = [np.inf, 0.0]
        for n, line in enumerate(eig_results):
            # Build waveform simulation command
            command = bin_path + "infraga-sph"
            
            if atmo_prefix is not None:
                command = command + "-rngdep"

            if atmo_file is not None:
                command = command + " -wnl_wvfrm " + atmo_file
            elif atmo_prefix is not None and grid_lats is not None and grid_lons is not None:
                command = command + " -wnl_wvfrm " + atmo_prefix + " " + grid_lats + " " + grid_lons
            else:
                click.echo("Simulation requires either an '--atmo-file' or --atmo-prefix' with grid info (--grid-lats and --grid-lons)")
                return 0

            # Set parameters
            command = set_param(command, str(line[0]), "inclination")
            command = set_param(command, str(line[1]), "azimuth")
            command = set_param(command, str(line[2]), "bounces")

            command = set_param(command, src_lat, "src_lat")
            command = set_param(command, src_lon, "src_lon")
            command = set_param(command, src_alt, "src_alt")
            command = command + " write_ray=true"

            command = set_param(command, wvfrm_file, "wvfrm_file")
            if wvfrm_yield is not None:
                # use atmosphere to define ambient pressure (1 mbar = 0.1 kPa)
                atmo = np.loadtxt(atmo_file)

                if src_alt is not None:
                    src_ht = float(src_alt)
                elif z_grnd is not None:
                    src_ht = float(z_grnd)
                else:
                    src_ht = 0.0

                p_ambient = atmo[np.argmin(abs(atmo[:, 0] - src_ht)), 5] * 0.1
                T_ambient = atmo[np.argmin(abs(atmo[:, 0] - src_ht)), 1]

                wvfrm_opt = 'impulse'
                wvfrm_ref = str(0.035 * float(wvfrm_yield)**(1.0 / 3.0))
                wvfrm_p0 = str(kg_op(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient))
                wvfrm_t0 = str(kg_ppd(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient))
                wvfrm_alpha = '0.01'

            command = set_param(command, wvfrm_opt, "wvfrm_opt")
            command = set_param(command, wvfrm_p0, "wvfrm_p0")
            command = set_param(command, wvfrm_t0, "wvfrm_t0")
            command = set_param(command, wvfrm_alpha, "wvfrm_alpha")
            command = set_param(command, wvfrm_ref, "wvfrm_ref")
            command = set_param(command, wvfrm_out_step, "wvfrm_out_step")
               
            command = set_param(command, wvfrm_ds, "wvfrm_ds")
            command = set_param(command, wvfrm_len, "wvfrm_len")
    
            command = set_param(command, freq, "freq")
            command = set_param(command, abs_coeff, "abs_coeff")
            
            command = set_param(command, output_id, "output_id")

            command = set_param(command, max_alt, "max_alt")
            command = set_param(command, max_rng, "max_rng")

            command = set_param(command, min_lat, "min_lat")
            command = set_param(command, max_lat, "max_lat")
            command = set_param(command, min_lon, "min_lon")
            command = set_param(command, max_lon, "max_lon")

            command = set_param(command, min_ds, "min_ds")
            command = set_param(command, max_ds, "max_ds")
            command = set_param(command, max_s, "max_s")

            command = set_param(command, z_grnd, "z_grnd")
            command = set_param(command, topo_file, "topo_file")
            if topo_file is not None:
                if topo_bl_wind is not None:
                    command = set_param(command, str(topo_bl_wind), "topo_bl_wind")

            subprocess.call(command, shell=True)

            temp = np.loadtxt(profile_id + ".wvfrm_out.dat")
            wvfrms += [interp1d(temp[:, 0], temp[:, 1], bounds_error=False, fill_value=0.0,  kind='cubic')]
            t_lims[0] = min(t_lims[0], temp[0][0])
            t_lims[1] = max(t_lims[1], temp[-1][0])
            dt = abs(temp[1][0] - temp[0][0])

            temp = np.loadtxt(profile_id + ".raypaths.dat")
            for line in temp:
                print(*line, file=file_out)
            print('\n', file=file_out)

            command = "rm " + profile_id + ".wvfrm_out.dat"
            command = command + " " + profile_id + ".raypaths.dat"
            command = command + " " + profile_id + ".wvfrm_init.dat"
            os.system(command)
        
        file_out.close()

        # combine waveforms and write to file
        print("Interpolating and merging waveforms...")
        print('\t' + "Eigenrays written into " + profile_id + ".eigenrays.dat")
        print('\t' + "Arrival waveform written into " + profile_id + ".wvfrms.dat")

        t_vals = np.arange(t_lims[0], t_lims[1], dt)
        file_out = open(profile_id + ".wvfrms.dat", 'w')
        print("# 'infraga sph eig_wvfrm' waveform results", '\n#', file=file_out)
        if atmo_prefix is not None:
            print("# 	profile: " + atmo_prefix, file=file_out)        
        else:
            print("# 	profile: " + atmo_file, file=file_out)
        print("# 	source location (lat, lon, alt): " + str(src_lat) + ", " + str(src_lon) + ", " + str(src_alt) , file=file_out)
        print("# 	receiver location (lat, lon, alt): " + str(rcvr_lat) + ", " + str(rcvr_lon) + ", 0.0", file=file_out)
        print("# 	inclination range: " + str(incl_min) + ", " + str(incl_max), file=file_out)
        print("#    inclination step max:", incl_step_max, file=file_out)
        print("# 	bounces: " + str(bnc_min) + ", " + str(bnc_max), file=file_out)
        if topo_file is not None:
            print("# 	terrain file: " + str(topo_file), file=file_out)        
        else:
            print("# 	ground elevation: " + str(z_grnd), file=file_out)
        print("#    damping:", damping, file=file_out)
        print("#    range max:" + '\n#', max_rng, file=file_out)

        # Waveform calculation parameters
        if wvfrm_ref is not None:
            print("#    waveform reference distance:", wvfrm_ref, file=file_out)
        else:
            print("#    waveform reference distance: 1.0", file=file_out)

        if wvfrm_len is not None:
            print("#    waveform length:", wvfrm_len, file=file_out)
        else:
            print("#    waveform length: 0.05", wvfrm_len, file=file_out)

        if wvfrm_yield:
            print("#    waveform source yield:", wvfrm_yield, file=file_out)
        else:
            print("#    waveform option:", wvfrm_opt, file=file_out)
            print("#    waveform peak op:", wvfrm_p0, file=file_out)
            print("#    waveform time sc.:", wvfrm_t0, file=file_out)
            print("#    waveform shaping param.:", wvfrm_alpha, file=file_out)

        print("#", file=file_out)
        print("# Eigenray arrivals:", file=file_out)
        print("# incl [deg]	az [deg]	n_b	x_0 [deg]	y_0 [deg]	time [s]	cel [km/s]	turning ht [km]	inclination [deg]	back azimuth [deg]	trans. coeff. [dB]	absorption [dB]", file=file_out)
        for line in eig_results:
            print("#", *line, file=file_out) 

        print('\n#', "t [s]" + '\t' + "p1 [Pa]", '\t' + "p2 [Pa] ...", file=file_out)
        for n in range(len(t_vals)):
            print(t_vals[n], end='\t', file=file_out)
            for wvfrm in wvfrms:
                print(wvfrm(t_vals[n]), end='\t', file=file_out)
            print('', file=file_out)
        file_out.close()
    else:
        print('\n' + "No waveforms to compute.")
    
    if not keep_eig_arrivals:
        os.system("rm " + profile_id + ".arrivals.dat")