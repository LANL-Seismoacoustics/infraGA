#!which python
"""
run_infraga.py

Run the various infraga methods using a simplified Python Click interface

Author: pblom@lanl.gov    
"""

import os 
import warnings
import glob
import fnmatch
import shlex
import click
import subprocess
import tempfile 

import configparser as cnfg
import numpy as np
import matplotlib.pyplot as plt 

from importlib.util import find_spec
from scipy.interpolate import interp1d

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

def kg_op(W, r, p_amb=101.325, T_amb=288.15, exp_type="chemical"):
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
    
    fd = (p_amb / 101.325)**(1.0 / 3.0) * (288.15 / T_amb)**(1.0 / 3.0)
    
    if exp_type=="chemical":
        sc_rng = fd / W**(1.0 / 3.0) * r * 1000.0
    
        term1 = 1.0 + (sc_rng / 4.5)**2
        term2 = np.sqrt(1.0 + (sc_rng / 0.048)**2)
        term3 = np.sqrt(1.0 + (sc_rng / 0.32)**2)
        term4 = np.sqrt(1.0 + (sc_rng / 1.35)**2)
        
        result = 808.0 * term1 / (term2 * term3 * term4)
    else:
        sc_rng = fd / (W * 1.0e-6)**(1.0 / 3.0) * r * 1000.0

        term1 = (1.0 + sc_rng / 800.0)
        term2 = np.sqrt(1.0 + (sc_rng / 87.0)**2)
        
        result = 3.2e6 / sc_rng**3 * term1 * term2

    return p_amb * 1.0e3 * result


def kg_ppd(W, r, p_amb=101.325, T_amb=288.15, exp_type="chemical"):
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

    fd = (p_amb / 101.325)**(1.0 / 3.0) * (288.15 / T_amb)**(1.0 / 3.0)
    
    if exp_type=="chemical":
        sc_rng = fd / W**(1.0 / 3.0) * r * 1000.0

        term1 = 1.0 + (sc_rng / 0.54)**10
        term2 = 1.0 + (sc_rng / 0.02)**3
        term3 = 1.0 + (sc_rng / 0.74)**6
        term4 = np.sqrt(1.0 + (sc_rng / 6.9)**2)
        
        result = 980.0 * term1 / (term2 * term3 * term4) * W**(1.0 / 3.0)
    else:
        sc_rng = fd / (W * 1.0e-6)**(1.0 / 3.0) * r * 1000.0

        term1 = np.sqrt(1.0 + (sc_rng / 100.0)**3)
        term2 = np.sqrt(1.0 + (sc_rng / 40.0))
        term3 = (1.0 + (sc_rng / 285.0)**5)**(1.0 / 6.0)
        term4 = (1.0 + (sc_rng / 50000.0))**(1.0 / 6.0)
        
        result = 180.0 * term1 / (term2 * term3 * term4) * (W * 1.0e-6) **(1.0 / 3.0)

    return result / 1e3


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
@click.option("--prof-format", help="Option to specify the atmospheric file format", default=None)
@click.option("--reverse-winds", help="Option to reverse wind directions for back projection", default=None, type=bool)
@click.option("--output-id", help="User specified output file path", default=None)
@click.option("--calc-amp", help="Option to turn off transport coefficient calculation", default=None, type=bool)
@click.option("--max-alt", help="Maximum altitude for ray calculation", default=None)
@click.option("--max-rng", help="Maximum range for ray calculation", default=None)

@click.option("--min-ds", help="Minimum step size (near-ground) in RK4 solver", default=None)
@click.option("--max-ds", help="Maximum step size in RK4 solver", default=None)
@click.option("--max-s", help="Maximum ray length between bounces", default=None)
@click.option("--refl-alt", help="Partial reflection altitude", default=None)

@click.option("--topo-file", help="Terrain file", default=None)
@click.option("--topo-BL-wind", help="Use terrain corrected boundary layer winds", default=None, type=bool)
def run_2d_prop(config_file, atmo_file, incl_min, incl_max, incl_step, inclination, azimuth, bounces, src_alt,
                freq, abs_coeff, z_grnd, write_atmo, prof_format, reverse_winds, output_id, calc_amp, max_alt,
                max_rng, min_ds, max_ds, max_s, refl_alt, topo_file, topo_bl_wind):
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

    prof_format = define_param(user_config, 'GENERAL', 'prof_format', prof_format)
    reverse_winds = define_param(user_config, 'GENERAL', 'reverse_winds', reverse_winds)
    output_id = define_param(user_config, 'GENERAL', 'output_id', output_id)

    max_alt = define_param(user_config, 'GENERAL', 'max_alt', max_alt)
    max_rng = define_param(user_config, 'GENERAL', 'max_rng', max_rng)

    min_ds = define_param(user_config, 'GENERAL', 'min_ds', min_ds)
    max_ds = define_param(user_config, 'GENERAL', 'max_ds', max_ds)
    max_s = define_param(user_config, 'GENERAL', 'max_s', max_s)
    refl_alt = define_param(user_config, 'GENERAL', 'refl_alt', refl_alt)

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

    command = set_param(command, prof_format, "prof_format")

    if reverse_winds is not None:
        command = set_param(command, str(reverse_winds), "reverse_winds")

    command = set_param(command, output_id, "output_id")

    command = set_param(command, max_alt, "max_alt")
    command = set_param(command, max_rng, "max_rng")

    command = set_param(command, min_ds, "min_ds")
    command = set_param(command, max_ds, "max_ds")
    command = set_param(command, max_s, "max_s")

    command = set_param(command, refl_alt, "refl_alt")

    command = set_param(command, z_grnd, "z_grnd")
    command = set_param(command, topo_file, "topo_file")
    if topo_file is not None:
        if topo_bl_wind is not None:
            command = set_param(command, str(topo_bl_wind), "topo_bl_wind")

    click.echo(command)
    subprocess.run(shlex.split(command), shell=False)


@click.command('wnl_wvfrm', short_help="Run weakly non-linear waveform calculation")
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
@click.option("--prof-format", help="Option to specify the atmospheric file format", default=None)
@click.option("--reverse-winds", help="Option to reverse wind directions for back projection", default=None, type=bool)
@click.option("--output-id", help="User specified output file path", default=None)
@click.option("--max-alt", help="Maximum altitude for ray calculation", default=None)
@click.option("--max-rng", help="Maximum range for ray calculation", default=None)

@click.option("--min-ds", help="Minimum step size (near-ground) in RK4 solver", default=None)
@click.option("--max-ds", help="Maximum step size in RK4 solver", default=None)
@click.option("--max-s", help="Maximum ray length between bounces", default=None)
@click.option("--refl-alt", help="Partial reflection altitude", default=None)

@click.option("--topo-file", help="Terrain file", default=None)
@click.option("--topo-BL-wind", help="Use terrain corrected boundary layer winds", default=None, type=bool)
def run_2d_wvfrm(config_file, atmo_file, inclination, azimuth, bounces, src_alt, write_ray, wvfrm_file, wvfrm_opt,
                wvfrm_p0, wvfrm_t0, wvfrm_alpha, wvfrm_ref, wvfrm_out_step, wvfrm_yield, wvfrm_ds, wvfrm_len, 
                freq, abs_coeff, z_grnd, write_atmo, prof_format, reverse_winds, output_id, max_alt, max_rng, 
                min_ds, max_ds, max_s, refl_alt, topo_file, topo_bl_wind):
    '''
    Run weakly non-linear waveform analysis along a 2D ray path using the effective sound speed approximation.

    \b
    Examples:
    \t infraga 2d wnl_wvfrm --atmo-file ToyAtmo.met --inclination 12.0 --azimuth -90.0 --wvfrm-p0 500.0
    \t infraga 2d wnl_wvfrm --atmo-file ToyAtmo.met --config-file example.cnfg
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
    prof_format = define_param(user_config, 'GENERAL', 'prof_format', prof_format)
    reverse_winds = define_param(user_config, 'GENERAL', 'reverse_winds', reverse_winds)
    output_id = define_param(user_config, 'GENERAL', 'output_id', output_id)

    max_alt = define_param(user_config, 'GENERAL', 'max_alt', max_alt)
    max_rng = define_param(user_config, 'GENERAL', 'max_rng', max_rng)

    min_ds = define_param(user_config, 'GENERAL', 'min_ds', min_ds)
    max_ds = define_param(user_config, 'GENERAL', 'max_ds', max_ds)
    max_s = define_param(user_config, 'GENERAL', 'max_s', max_s)
    refl_alt = define_param(user_config, 'GENERAL', 'refl_alt', refl_alt)

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

    command = set_param(command, prof_format, "prof_format")
    if reverse_winds is not None:
        command = set_param(command, str(reverse_winds), "reverse_winds")

    command = set_param(command, output_id, "output_id")

    command = set_param(command, max_alt, "max_alt")
    command = set_param(command, max_rng, "max_rng")

    command = set_param(command, min_ds, "min_ds")
    command = set_param(command, max_ds, "max_ds")
    command = set_param(command, max_s, "max_s")

    command = set_param(command, refl_alt, "refl_alt")

    command = set_param(command, z_grnd, "z_grnd")
    command = set_param(command, topo_file, "topo_file")
    if topo_file is not None:
        if topo_bl_wind is not None:
            command = set_param(command, str(topo_bl_wind), "topo_bl_wind")

    click.echo(command)
    subprocess.run(shlex.split(command), shell=False)


@click.command('refl_eigs', short_help="Compute partial reflection eigenrays")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)

@click.option("--refl-alt-min", help="Lower partial reflection altitude [km]", default=25.0)
@click.option("--refl-alt-max", help="Upper partial reflection altitude [km]", default=40.0)
@click.option("--refl-alt-step", help="Partial reflection resolution [km]", default=0.5)

@click.option("--rcvr-rng", help="Receiver range [km]", default=100.0)
@click.option("--rcvr-az", help="REceiver azimuth (clockwise rel. N)", default=90.0)
@click.option("--local-temp-dir", help="Local temporary directory for results", default=None)
@click.option("--output-id", help="Output file path", default=None)
@click.option("--verbose", help="Verbose output of simulations", default=False)

@click.option("--incl-min", help="Minimum inclination angle (rel. horizontal)", default=None)
@click.option("--incl-max", help="Maximum inclination angle (rel. horizontal)", default=None)
@click.option("--incl-step", help="Inclination angle step resolution", default=None)
@click.option("--bounces", help="Maximum number of ground reflections (bounces)", default=None)
@click.option("--src-alt", help="Source altitude", default=None)

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
@click.option("--prof-format", help="Option to specify the atmospheric file format", default=None)
@click.option("--reverse-winds", help="Option to reverse wind directions for back projection", default=None, type=bool)
@click.option("--max-alt", help="Maximum altitude for ray calculation", default=None)
@click.option("--max-rng", help="Maximum range for ray calculation", default=None)

@click.option("--min-ds", help="Minimum step size (near-ground) in RK4 solver", default=None)
@click.option("--max-ds", help="Maximum step size in RK4 solver", default=None)
@click.option("--max-s", help="Maximum ray length between bounces", default=None)

@click.option("--topo-file", help="Terrain file", default=None)
@click.option("--topo-BL-wind", help="Use terrain corrected boundary layer winds", default=None, type=bool)
def run_2d_refl_eigs(config_file, atmo_file, refl_alt_min, refl_alt_max, refl_alt_step, rcvr_rng, rcvr_az, local_temp_dir, output_id, verbose,
                     incl_min, incl_max, incl_step, bounces, src_alt, wvfrm_file, wvfrm_opt, wvfrm_p0, wvfrm_t0, wvfrm_alpha, wvfrm_ref,
                     wvfrm_out_step, wvfrm_yield, wvfrm_ds, wvfrm_len, freq, abs_coeff, z_grnd, write_atmo, prof_format, reverse_winds, 
                     max_alt, max_rng, min_ds, max_ds, max_s,topo_file, topo_bl_wind):
    '''
    Run partial reflected paths to a specified range via 2D rays using the effective sound speed approximation.

    \b
    Examples:
    \t infraga 2d refl_eigs --atmo-file ToyAtmo.met --rcvr-rng 140.0 --rcvr-az -90.0 --wvfrm-yield 1.0e3 --local-temp-dir refl_temp
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

    bounces = define_param(user_config, 'PROP', 'bounces', bounces)
    src_alt = define_param(user_config, 'PROP', 'src_alt', src_alt)

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
    prof_format = define_param(user_config, 'GENERAL', 'prof_format', prof_format)
    reverse_winds = define_param(user_config, 'GENERAL', 'reverse_winds', reverse_winds)
    output_id = define_param(user_config, 'GENERAL', 'output_id', output_id)

    max_alt = define_param(user_config, 'GENERAL', 'max_alt', max_alt)
    max_rng = define_param(user_config, 'GENERAL', 'max_rng', max_rng)

    min_ds = define_param(user_config, 'GENERAL', 'min_ds', min_ds)
    max_ds = define_param(user_config, 'GENERAL', 'max_ds', max_ds)
    max_s = define_param(user_config, 'GENERAL', 'max_s', max_s)
 
    z_grnd = define_param(user_config, 'GENERAL', 'z_grnd', z_grnd)
    topo_file = define_param(user_config, 'GENERAL', 'topo_file', topo_file)
    topo_bl_wind = define_param(user_config, 'GENERAL', 'topo_bl_wind', topo_bl_wind)

    with tempfile.TemporaryDirectory(prefix='infraga_') as tmpdirname:
        if local_temp_dir is not None:
            click.echo("Writing individual reflection altitude results into " + local_temp_dir)
            if not os.path.isdir(local_temp_dir):
                os.mkdir(local_temp_dir)
            tmpdirname = local_temp_dir
        else:
            click.echo('Created temp directory: ' + tmpdirname)

        # Set parameter values for header output
        def set_header_val(param, default):
            if param is None:
                return default
            else:
                return param
            
        src_alt_out = set_header_val(src_alt, 0.0)
        incl_min_out = set_header_val(incl_min, 0.5)
        incl_max_out = set_header_val(incl_max, 45.0)
        incl_step_out = set_header_val(incl_step, 0.5)
        z_grnd_out = set_header_val(z_grnd, 0.0)

        wvfrm_opt_out = set_header_val(wvfrm_opt, 'impulse')
        wvfrm_p0_out = set_header_val(wvfrm_p0, 10.0)
        wvfrm_t0_out = set_header_val(wvfrm_t0, 1.0)
        wvfrm_alpha_out = set_header_val(wvfrm_alpha, 1.0)

        if output_id is None:
            output_id = os.path.splitext(atmo_file)[0]
                
        eigenrays_out = open(output_id + ".eigenrays.dat", 'w')
        wvfrms_out = open(output_id + ".wvfrms.dat", 'w')
        
        print("# 'infraga 2d refl_eigs' eigenray paths", '\n#', file=eigenrays_out)
        print("# 	profile: " + atmo_file, file=eigenrays_out)
        print("# 	source location: 0.0, " + str(src_alt_out), file=eigenrays_out)
        print("# 	receiver range [km], azimuth [deg]: " + str(rcvr_rng) + ", " + str(rcvr_az), file=eigenrays_out)
        print("# 	inclination range: " + str(incl_min_out) + ", " + str(incl_max_out) + ", " + str(incl_step_out), file=eigenrays_out)
        if topo_file is not None:
            print("# 	terrain file: " + str(topo_file), file=eigenrays_out)        
        else:
            print("# 	ground elevation: " + str(z_grnd_out), file=eigenrays_out)
        print('\n# r [km]	z [km]	trans. coeff. [dB]	absorption [dB]	time [s]', file=eigenrays_out)

        print("# 'infraga 2d refl_eigs' waveform results", '\n#', file=wvfrms_out)
        print("# 	profile: " + atmo_file, file=wvfrms_out)
        print("# 	source location: 0.0, " + str(src_alt_out), file=wvfrms_out)
        print("# 	receiver range [km], azimuth [deg]: " + str(rcvr_rng) + ", " + str(rcvr_az), file=wvfrms_out)
        print("# 	inclination range: " + str(incl_min_out) + ", " + str(incl_max_out) + ", " + str(incl_step_out), file=wvfrms_out)
        if topo_file is not None:
            print("# 	terrain file: " + str(topo_file), file=wvfrms_out)        
        else:
            print("# 	ground elevation: " + str(z_grnd_out), file=wvfrms_out)

        # Waveform calculation parameters
        if wvfrm_ref is not None:
            print("#    waveform reference distance:", wvfrm_ref, file=wvfrms_out)
        else:
            print("#    waveform reference distance: 1.0", file=wvfrms_out)

        if wvfrm_yield:
            print("#    waveform source yield:", wvfrm_yield, file=wvfrms_out)
        else:
            print("#    waveform option:", wvfrm_opt_out, file=wvfrms_out)
            print("#    waveform peak op:", wvfrm_p0_out, file=wvfrms_out)
            print("#    waveform time sc.:", wvfrm_t0_out, file=wvfrms_out)
            print("#    waveform shaping param.:", wvfrm_alpha_out, file=wvfrms_out)

        print('#\n# Eigenray arrivals:', file=wvfrms_out)
        print('\n# incl [deg]	az [deg]	n_b	r_0 [km]	time [s]	celerity [km/s]	turning ht [km]	arrival incl. [deg]	trans. coeff. [dB]	absorption [dB]	perp. dist. [km]', file=wvfrms_out)

        wvfrms = []
        t_lims = [np.inf, 0.0]
        for refl_alt in np.arange(refl_alt_min, refl_alt_max, refl_alt_step):
            click.echo("Identifying propagation paths with partial reflection altitude at " + str(refl_alt) + "km...")

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

            command = set_param(command, str(rcvr_az), "azimuth")
            command = set_param(command, bounces, "bounces")
            command = set_param(command, src_alt, "src_alt")

            command = set_param(command, freq, "freq")
            command = set_param(command, abs_coeff, "abs_coeff")

            if write_atmo is not None:
                command = set_param(command, str(write_atmo), "write_atmo")

            command = set_param(command, prof_format, "prof_format")

            if reverse_winds is not None:
                command = set_param(command, str(reverse_winds), "reverse_winds")

            command = set_param(command, max_alt, "max_alt")
            command = set_param(command, str(1.1 * rcvr_rng), "max_rng")

            command = set_param(command, min_ds, "min_ds")
            command = set_param(command, max_ds, "max_ds")
            command = set_param(command, max_s, "max_s")

            command = set_param(command, str(refl_alt), "refl_alt")

            command = set_param(command, z_grnd, "z_grnd")
            command = set_param(command, topo_file, "topo_file")
            if topo_file is not None:
                if topo_bl_wind is not None:
                    command = set_param(command, str(topo_bl_wind), "topo_bl_wind")

            command = command + " output_id=" + tmpdirname + "/temp_" + str(refl_alt) + "km"

            if verbose:
                click.echo(command)
                subprocess.run(shlex.split(command), shell=False)
            else:
                subprocess.run(shlex.split(command), shell=False, stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)

            # interpolate inclination vs. range for each bounce count to get arrival(s)
            arr_info = np.loadtxt( tmpdirname + "/temp_" + str(refl_alt) + "km.arrivals.dat")

            arr_incl = arr_info[:, 0]
            arr_bncs = arr_info[:, 2].astype(int)
            arr_rng = arr_info[:, 3]

            for bn in range(max(arr_bncs)):
                print('\tEstimating intercept with ' + str(bn) + " bounces...")
                bnc_mask = arr_bncs == bn
                
                k = np.argmin(abs(arr_rng[bnc_mask] - rcvr_rng))
                incl_est = arr_incl[bnc_mask][k] + (arr_incl[bnc_mask][k + 1] - arr_incl[bnc_mask][k - 1]) / (arr_rng[bnc_mask][k + 1] - arr_rng[bnc_mask][k - 1]) * (rcvr_rng - arr_rng[bnc_mask][k])

                # Write arrival into file
                print("# " + str(incl_est), file=wvfrms_out, end=' ')
                print(rcvr_az, file=wvfrms_out, end=' ')
                print(bn, file=wvfrms_out, end=' ')
                print(rcvr_rng, file=wvfrms_out, end=' ')
                for N in range(4, 11):
                    print(arr_info[:, N][bnc_mask][k] + ( arr_info[:, N][bnc_mask][k + 1] -  arr_info[:, N][bnc_mask][k - 1]) / (arr_rng[bnc_mask][k + 1] - arr_rng[bnc_mask][k - 1]) * (rcvr_rng - arr_rng[bnc_mask][k]), file=wvfrms_out, end=' ')
                print('',file=wvfrms_out)

                command = bin_path + "infraga-2d -wnl_wvfrm " + atmo_file
                command = set_param(command, str(incl_est), "inclination")
                command = set_param(command, str(rcvr_az), "azimuth")
                command = set_param(command, str(bn), "bounces")

                command = set_param(command, src_alt, "src_alt")

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

                command = set_param(command, "True", "write_ray")

                command = set_param(command, prof_format, "prof_format")
                if reverse_winds is not None:
                    command = set_param(command, str(reverse_winds), "reverse_winds")

                command = set_param(command, max_alt, "max_alt")
                command = set_param(command, max_rng, "max_rng")

                command = set_param(command, min_ds, "min_ds")
                command = set_param(command, max_ds, "max_ds")
                command = set_param(command, max_s, "max_s")

                command = set_param(command, str(refl_alt), "refl_alt")

                command = set_param(command, z_grnd, "z_grnd")
                command = set_param(command, topo_file, "topo_file")
                if topo_file is not None:
                    if topo_bl_wind is not None:
                        command = set_param(command, str(topo_bl_wind), "topo_bl_wind")

                command = command + " output_id=" + tmpdirname + "/temp_" + str(refl_alt) + "km-" + str(bn)

                if verbose:
                    click.echo(command)
                    subprocess.run(shlex.split(command), shell=False)
                else:
                    subprocess.run(shlex.split(command), shell=False, stdout=open(os.devnull, 'w'), stderr=subprocess.STDOUT)

                temp = np.loadtxt(tmpdirname + "/temp_" + str(refl_alt) + "km-" + str(bn) + ".wvfrm_out.dat")
                wvfrms += [interp1d(temp[:, 0], temp[:, 1], bounds_error=False, fill_value=0.0,  kind='cubic')]
                t_lims[0] = min(t_lims[0], temp[0][0])
                t_lims[1] = max(t_lims[1], temp[-1][0])
                dt = abs(temp[1][0] - temp[0][0])

                temp = np.loadtxt(tmpdirname + "/temp_" + str(refl_alt) + "km-" + str(bn) + ".raypaths.dat")
                for line in temp:
                    print(*line, file=eigenrays_out)
                print('\n', file=eigenrays_out)
                
        click.echo('\nInterpolating and merging waveforms...')
        print('\n#', "t [s]" + '\t' + "p1 [Pa]", '\t' + "p2 [Pa] ...", file=wvfrms_out)
        t_vals = np.arange(t_lims[0], t_lims[1], dt)
        for n in range(len(t_vals)):
            print(t_vals[n], end='\t', file=wvfrms_out)
            for wvfrm in wvfrms:
                print(wvfrm(t_vals[n]), end='\t', file=wvfrms_out)
            print('', file=wvfrms_out)

        eigenrays_out.close()
        wvfrms_out.close()



###################################
##                               ##
##   Python CLI for infraga-3d   ##
##                               ##
###################################
@click.command('prop', short_help="Run a point source propagation simulation")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)
@click.option("--atmo-prefix", help="Atmosphere file prefix (range dependent)", default=None)
@click.option("--grid-x", help="Atmosphere grid x (E/W) file", default=None)
@click.option("--grid-y", help="Atmosphere grid y (N/S) file", default=None)

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
@click.option("--prof-format", help="Option to specify the atmospheric file format", default=None)
@click.option("--reverse-winds", help="Option to reverse wind directions for back projection", default=None, type=bool)
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
@click.option("--refl-alt", help="Partial reflection altitude", default=None)


@click.option("--topo-file", help="Terrain file", default=None)
@click.option("--topo-BL-wind", help="Use terrain corrected boundary layer winds", default=None, type=bool)
@click.option("--cpu-cnt", help="Number of CPUs to use in analysis", default=1)
def run_3d_prop(config_file, atmo_file, atmo_prefix, grid_x, grid_y, incl_min, incl_max, incl_step, inclination, 
                az_min, az_max, az_step, azimuth, bounces, src_x, src_y, src_alt, write_rays, write_topo, freq, 
                abs_coeff, z_grnd, write_atmo, prof_format, reverse_winds, output_id, calc_amp, max_alt, max_rng,
                min_x, max_x, min_y, max_y, min_ds, max_ds, max_s, topo_file, topo_bl_wind, cpu_cnt):
    '''
    Run 3D geometry ray tracing analysis for a point source

    \b
    Examples:
    \t infraga 3d prop --atmo-file ToyAtmo.met
    \t infraga 3d prop --atmo-file ToyAtmo.met --config-file example.cnfg
    \t infraga 3d prop --atmo-prefix profs/example --grid-x profs/example_x.dat --grid-y profs/example_x.dat --azimuth -90.0

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

    prof_format = define_param(user_config, 'GENERAL', 'prof_format', prof_format)
    reverse_winds = define_param(user_config, 'GENERAL', 'reverse_winds', reverse_winds)
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
    if cpu_cnt is None:
        cpu_cnt = 1 
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
    elif atmo_prefix is not None and grid_x is not None and grid_y is not None:
        command = command + " -prop " + atmo_prefix + " " + grid_x + " " + grid_y
    else:
        click.echo("Simulation requires either an '--atmo-file' or --atmo-prefix' with grid info (--grid-x and --grid-y)")
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

    command = set_param(command, prof_format, "prof_format")
    if reverse_winds is not None:
        command = set_param(command, str(reverse_winds), "reverse_winds")
        
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

    click.echo(command)
    subprocess.run(shlex.split(command), shell=False)



@click.command('eigenray', short_help="Run eigenray search for specific source-receiver")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)
@click.option("--atmo-prefix", help="Atmosphere file prefix (range dependent)", default=None)
@click.option("--grid-x", help="Atmosphere grid x (E/W) file", default=None)
@click.option("--grid-y", help="Atmosphere grid y (N/S) file", default=None)

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
@click.option("--prof-format", help="Option to specify the atmospheric file format", default=None)
@click.option("--reverse-winds", help="Option to reverse wind directions for back projection", default=None, type=bool)
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
@click.option("--cpu-cnt", help="Number of CPUs to use in analysis", default=None)
def run_3d_eig(config_file, atmo_file, atmo_prefix, grid_x, grid_y, incl_min, incl_max, bnc_min, bnc_max, bounces, 
                src_x, src_y, src_alt, rcvr_x, rcvr_y, verbose, iterations, damping, tolerance, az_dev_lim, 
                incl_step_min, incl_step_max, freq, abs_coeff, z_grnd, write_atmo, prof_format, reverse_winds, 
                output_id, max_alt, max_rng, min_x, max_x, min_y, max_y, min_ds, max_ds, max_s, topo_file, topo_bl_wind, cpu_cnt):
    '''
    Run 3D Cartesian eigenray analysis to identify propagation paths connecting a specific source-receiver geometry

    \b
    Examples:
    \t infraga 3d eigenray --atmo-file ToyAtmo.met --rcvr-x -175.0 --rcvr-y 75.0 --bnc-max 1 --verbose True
    \t infraga 3d eigenray --atmo-file ToyAtmo.met --config-file example.cnfg
    \t infraga 3d eigenray --atmo-prefix profs/example --grid-x profs/example_x.dat --grid-y profs/example_x.dat --src-x 0.0 --src-y 0.0 --rcvr-x -500.0 --rcvr-y -100.0 --bnc-max 1 --incl-min 10.0 --incl-max 20.0 --verbose true

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
    prof_format = define_param(user_config, 'GENERAL', 'prof_format', prof_format)
    reverse_winds = define_param(user_config, 'GENERAL', 'reverse_winds', reverse_winds)
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
    if cpu_cnt is None:
        cpu_cnt = 1 
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
    elif atmo_prefix is not None and grid_x is not None and grid_y is not None:
        command = command + " -eig_search " + atmo_prefix + " " + grid_x + " " + grid_y
    else:
        click.echo("Simulation requires either an '--atmo-file' or --atmo-prefix' with grid info (--grid-x and --grid-y)")
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

    command = set_param(command, prof_format, "prof_format")
    if reverse_winds is not None:
        command = set_param(command, str(reverse_winds), "reverse_winds")

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

    click.echo(command)
    subprocess.run(shlex.split(command), shell=False)



@click.command('wnl_wvfrm', short_help="Run weakly non-linear waveform calculation")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)
@click.option("--atmo-prefix", help="Atmosphere file prefix (range dependent)", default=None)
@click.option("--grid-x", help="Atmosphere grid x (E/W) file", default=None)
@click.option("--grid-y", help="Atmosphere grid y (N/S) file", default=None)

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
@click.option("--expl-type", help="Explosion type ('chemical' vs. 'nuclear')", default="chemical")

@click.option("--wvfrm-ds", help="Burgers equation solver resolution [km]", default=None)
@click.option("--wvfrm-len", help="Waveform sample count (default = 2e13)", default=None)

@click.option("--freq", help="Frequency for Sutherland-Bass losses", default=None)
@click.option("--abs-coeff", help="Scaling coefficient for Sutherland-Bass losses", default=None)
@click.option("--z-grnd", help="Ground elevation (flat ground, km rel. sea level)", default=None)
@click.option("--write-atmo", help="Option to write atmosphere data (for QC)", default=None)
@click.option("--prof-format", help="Option to specify the atmospheric file format", default=None)
@click.option("--reverse-winds", help="Option to reverse wind directions for back projection", default=None, type=bool)
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
def run_3d_wvfrm(config_file, atmo_file, atmo_prefix, grid_x, grid_y, inclination, azimuth, bounces, src_x, src_y, 
                src_alt, write_ray, wvfrm_file, wvfrm_opt, wvfrm_p0, wvfrm_t0, wvfrm_alpha, wvfrm_ref, wvfrm_out_step,
                wvfrm_yield, expl_type, wvfrm_ds, wvfrm_len, freq, abs_coeff, z_grnd, write_atmo, prof_format, reverse_winds,
                output_id, max_alt, max_rng, min_x, max_x, min_y, max_y, min_ds, max_ds, max_s, topo_file, topo_bl_wind):
    '''
    Run weakly non-linear waveform analysis along a 3D Cartesian ray path.

    \b
    Examples:
    \t infraga 3d wnl_wvfrm --atmo-file ToyAtmo.met --inclination 4.1314876 --azimuth -84.969455 --bounces 1 --wvfrm-p0 500.0
    \t infraga 3d wnl_wvfrm --atmo-file ToyAtmo.met --config-file example.cnfg
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
    prof_format = define_param(user_config, 'GENERAL', 'prof_format', prof_format)
    reverse_winds = define_param(user_config, 'GENERAL', 'reverse_winds', reverse_winds)
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
    elif atmo_prefix is not None and grid_x is not None and grid_y is not None:
        command = command + " -wnl_wvfrm " + atmo_prefix + " " + grid_x + " " + grid_y
    else:
        click.echo("Simulation requires either an '--atmo-file' or --atmo-prefix' with grid info (--grid-x and --grid-y)")
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
        wvfrm_p0 = str(kg_op(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient, exp_type=expl_type))
        wvfrm_t0 = str(kg_ppd(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient, exp_type=expl_type))
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

    command = set_param(command, prof_format, "prof_format")
    if reverse_winds is not None:
        command = set_param(command, str(reverse_winds), "reverse_winds")

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

    click.echo(command)
    subprocess.run(shlex.split(command), shell=False)



@click.command('eig_wvfrm', short_help="Run eigenray search and compute waveform contributions")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)
@click.option("--atmo-prefix", help="Atmosphere file prefix (range dependent)", default=None)
@click.option("--grid-x", help="Atmosphere grid x (E/W) file", default=None)
@click.option("--grid-y", help="Atmosphere grid y (N/S) file", default=None)

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
@click.option("--expl-type", help="Explosion type ('chemical' vs. 'nuclear')", default="chemical")

@click.option("--wvfrm-ds", help="Burgers equation solver resolution [km]", default=None)
@click.option("--wvfrm-len", help="Waveform sample count (default = 2e13)", default=None)

@click.option("--freq", help="Frequency for Sutherland-Bass losses", default=None)
@click.option("--abs-coeff", help="Scaling coefficient for Sutherland-Bass losses", default=None)
@click.option("--z-grnd", help="Ground elevation (km rel. sea level)", default=None)
@click.option("--write-atmo", help="Option to write atmosphere data (for QC)", default=None)
@click.option("--prof-format", help="Option to specify the atmospheric file format", default=None)
@click.option("--reverse-winds", help="Option to reverse wind directions for back projection", default=None, type=bool)
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
@click.option("--cpu-cnt", help="Number of CPUs to use in analysis", default=None)

@click.option("--keep-eig-results", help="Keep eigenray arrivals", default=False)
def run_3d_eig_wvfrm(config_file, atmo_file, atmo_prefix, grid_x, grid_y, incl_min, incl_max, bnc_min, bnc_max, bounces, 
                src_x, src_y, src_alt, rcvr_x, rcvr_y, verbose, iterations, damping, tolerance, az_dev_lim, incl_step_min, 
                incl_step_max, wvfrm_file, wvfrm_opt, wvfrm_p0, wvfrm_t0, wvfrm_alpha, wvfrm_ref, wvfrm_out_step, wvfrm_yield, 
                expl_type, wvfrm_ds, wvfrm_len, freq, abs_coeff, z_grnd, write_atmo, prof_format, reverse_winds, output_id,
                max_alt, max_rng, min_x, max_x, min_y, max_y, min_ds, max_ds, max_s, topo_file, topo_bl_wind, cpu_cnt, 
                keep_eig_results):
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
    prof_format = define_param(user_config, 'GENERAL', 'prof_format', prof_format)
    reverse_winds = define_param(user_config, 'GENERAL', 'reverse_winds', reverse_winds)
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
    if cpu_cnt is None:
        cpu_cnt = 1 
    cpu_cnt = int(cpu_cnt)

    # Set parameter values for header output
    def set_header_val(param, default):
        if param is None:
            return default
        else:
            return param

    src_x_out = set_header_val(src_x, 0.0)
    src_y_out = set_header_val(src_y, 0.0)
    src_alt_out = set_header_val(src_alt, 0.0)

    rcvr_x_out = set_header_val(rcvr_x, 0.0)
    rcvr_y_out = set_header_val(rcvr_y, 0.0)

    incl_min_out = set_header_val(incl_min, 0.5)
    incl_max_out = set_header_val(incl_max, 45.0)
    incl_step_max_out = set_header_val(incl_step_max, 0.1)

    bnc_min_out = set_header_val(bnc_min, 0)
    bnc_max_out = set_header_val(bnc_max, 0)
    if bounces is not None:
        bnc_min_out = bounces
        bnc_max_out = bounces 

    z_grnd_out = set_header_val(z_grnd, 0.0)
    damping_out = set_header_val(damping, 1.0e-3)
    max_rng_out = set_header_val(max_rng, 2500.0)

    wvfrm_opt_out = set_header_val(wvfrm_opt, 'impulse')
    wvfrm_p0_out = set_header_val(wvfrm_p0, 10.0)
    wvfrm_t0_out = set_header_val(wvfrm_t0, 1.0)
    wvfrm_alpha_out = set_header_val(wvfrm_alpha, 1.0)

    # Check if eigenray analysis is already done
    if output_id is not None:
        eig_arrivals_file = output_id + ".arrivals.dat"
    else:
        if atmo_prefix is not None:
            eig_arrivals_file = atmo_prefix + ".arrivals.dat"
            output_id = atmo_prefix
        else:
            eig_arrivals_file = os.path.splitext(atmo_file)[0] + ".arrivals.dat"
            output_id = os.path.splitext(atmo_file)[0]
    
    eig_results_check = False
    if os.path.isfile(eig_arrivals_file):
        eig_results_option_check = False 
        eig_results_atmo_check = False
        eig_results_src_check = False
        eig_results_rcvr_check = False 

        arrivals_file = open(eig_arrivals_file, 'r')
        for line in arrivals_file:
            if "-eig_search" in line:
                eig_results_option_check = True
            if "profile" in line:
                eig_results_atmo_check = line.split(' ')[-1].strip('\n') == atmo_file
            if "source location" in line:
                ref_loc = np.array([float(item) for item in line.replace(",","").split(' ')[-3:]])
                eig_results_src_check = np.allclose(np.array([src_x_out, src_y_out, src_alt_out], dtype=float), ref_loc)
            if "receiver location" in line:
                ref_loc = [float(item) for item in line.replace(",","").split(' ')[-3:]]
                eig_results_rcvr_check = np.allclose(np.array([rcvr_x_out, rcvr_y_out, 0.0], dtype=float), ref_loc)

        if not eig_results_option_check:
            click.echo("Arrivals file detected, but not from eigenray analysis.")
        if not eig_results_atmo_check:
            click.echo("Arrivals file detected, but atmosphere file doesn't agree.")
        if not eig_results_src_check:
            click.echo("Arrivals file detected, but source location doesn't agree.")
        if not eig_results_rcvr_check:
            click.echo("Arrivals file detected, but receiver location doesn't agree.")

        eig_results_check = eig_results_option_check and eig_results_atmo_check and eig_results_src_check and eig_results_rcvr_check

    with tempfile.TemporaryDirectory(prefix='infraga_') as tmpdirname:
        print('Created temp directory:', tmpdirname)
        temp_path = tmpdirname + "/temp"
            
        if eig_results_check:
            click.echo('\nEigenray results found.\nSkipping to waveform calculation...')
        else:
            click.echo("Running eigenray analysis...")
            eig_arrivals_file = temp_path + ".arrivals.dat"

            # Build run command
            if cpu_cnt < 2:
                command = bin_path + "infraga-3d"
            else:
                command = "mpirun -np " + str(cpu_cnt) + " " + bin_path + "infraga-accel-3d"

            if atmo_prefix is not None:
                command = command + "-rngdep"

            if atmo_file is not None:
                command = command + " -eig_search " + atmo_file
            elif atmo_prefix is not None and grid_x is not None and grid_y is not None:
                command = command + " -eig_search " + atmo_prefix + " " + grid_x + " " + grid_y
            else:
                click.echo("Simulation requires either an '--atmo-file' or --atmo-prefix' with grid info (--grid-x and --grid-y)")
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

            command = set_param(command, prof_format, "prof_format")

            if reverse_winds is not None:
                command = set_param(command, str(reverse_winds), "reverse_winds")
                
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

            command = command + " output_id=" + temp_path

            click.echo(command)
            subprocess.run(shlex.split(command), shell=False)

        # Analyze eigenray results and compute waveform contributions...
        click.echo('\n' + "Loading eigenray results from " + eig_arrivals_file)
        eig_results = np.loadtxt(eig_arrivals_file)

        if len(eig_results) > 0:
            eig_results = np.atleast_2d(eig_results)

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
            command_list = []
            for n, line in enumerate(eig_results):
                command = bin_path + "infraga-3d"
                
                if atmo_prefix is not None:
                    command = command + "-rngdep"

                if atmo_file is not None:
                    command = command + " -wnl_wvfrm " + atmo_file
                elif atmo_prefix is not None and grid_x is not None and grid_y is not None:
                    command = command + " -wnl_wvfrm " + atmo_prefix + " " + grid_x + " " + grid_y
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
                    wvfrm_p0 = str(kg_op(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient, exp_type=expl_type))
                    wvfrm_t0 = str(kg_ppd(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient, exp_type=expl_type))
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

                command = set_param(command, prof_format, "prof_format")

                if reverse_winds is not None:
                    command = set_param(command, str(reverse_winds), "reverse_winds")
                    
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

                command = command + " output_id=" + temp_path + "-" + str(n)
                command = command + " > /dev/null"
                
                command_list = command_list + [command]
                
            click.echo('\n' + "Computing waveforms...")
            for j in range(0, len(command_list), cpu_cnt):              
                procs_list = [subprocess.Popen(cmd, shell=True) for cmd in command_list[j:j + cpu_cnt]]
                for proc in procs_list:
                    proc.communicate()
                    proc.wait()

            file_out = open(output_id + ".eigenrays.dat", 'w')
            print("# 'infraga 3d eig_wvfrm' eigenray results", '\n#', file=file_out)
            if atmo_prefix is not None:
                print("# 	profile: " + atmo_prefix, file=file_out)
            else:
                print("# 	profile: " + atmo_file, file=file_out)
            print("#    source location: " + str(src_x_out) + ", " + str(src_y_out) + ", " + str(src_alt_out), file=file_out)
            print("#    receiver location: " + str(rcvr_x_out) + ", " + str(rcvr_y_out) + ", 0.0", file=file_out)
            print("#    inclination range: " + str(incl_min_out) + ", " + str(incl_max_out), file=file_out)
            print("#    inclination step max:", incl_step_max_out, file=file_out)
            print("# 	bounces: " + str(bnc_min_out) + ", " + str(bnc_max_out), file=file_out)
            if topo_file is not None:
                print("# 	terrain file: " + str(topo_file), file=file_out)        
            else:
                print("# 	ground elevation: " + str(z_grnd_out), file=file_out)
            print("# 	damping:", damping_out, file=file_out)
            print("# 	range max:", max_rng_out, '\n#', file=file_out)
            print("# x [km]	y [km]	z [km]	trans. coeff. [dB]	absorption [dB]	time [s]", file=file_out)

            wvfrms = []
            t_lims = [np.inf, 0.0]
            for n, line in enumerate(eig_results):

                temp = np.loadtxt(temp_path + "-" + str(n) + ".wvfrm_out.dat")
                wvfrms += [interp1d(temp[:, 0], temp[:, 1], bounds_error=False, fill_value=0.0,  kind='cubic')]
                t_lims[0] = min(t_lims[0], temp[0][0])
                t_lims[1] = max(t_lims[1], temp[-1][0])
                dt = abs(temp[1][0] - temp[0][0])

                temp = np.loadtxt(temp_path + "-" + str(n) + ".raypaths.dat")
                for line in temp:
                    print(*line, file=file_out)
                print('\n', file=file_out)

            file_out.close()

            # combine waveforms and write to file
            click.echo("Interpolating and merging waveforms...")
            click.echo('\t' + "Eigenrays written into " + output_id + ".eigenrays.dat")
            click.echo('\t' + "Arrival waveform written into " + output_id + ".wvfrms.dat")

            t_vals = np.arange(t_lims[0], t_lims[1], dt)
            file_out = open(output_id + ".wvfrms.dat", 'w')
            print("# 'infraga 3d eig_wvfrm' waveform results", '\n#', file=file_out)
            if atmo_prefix is not None:
                print("# 	profile: " + atmo_prefix, file=file_out)        
            else:
                print("# 	profile: " + atmo_file, file=file_out)
            print("# 	source location: " + str(src_x_out) + ", " + str(src_y_out) + ", " + str(src_alt_out) , file=file_out)
            print("# 	receiver location: " + str(rcvr_x_out) + ", " + str(rcvr_y_out) + ", 0.0", file=file_out)
            print("# 	inclination range: " + str(incl_min_out) + ", " + str(incl_max_out), file=file_out)
            print("#    inclination step max:", incl_step_max_out, file=file_out)
            print("# 	bounces: " + str(bnc_min_out) + ", " + str(bnc_max_out), file=file_out)
            if topo_file is not None:
                print("# 	terrain file: " + str(topo_file), file=file_out)        
            else:
                print("# 	ground elevation: " + str(z_grnd_out), file=file_out)
            print("#    damping:", damping_out, file=file_out)
            print("#    range max:", max_rng_out, '\n#', file=file_out)

            # Waveform calculation parameters
            if wvfrm_ref is not None:
                print("#    waveform reference distance:", wvfrm_ref, file=file_out)
            else:
                print("#    waveform reference distance: 1.0", file=file_out)

            if wvfrm_len is not None:
                print("#    waveform length:", wvfrm_len, file=file_out)
            else:
                print("#    waveform length: 2e13", wvfrm_len, file=file_out)

            if wvfrm_yield is not None:
                print("#    waveform source yield:", wvfrm_yield, file=file_out)
            else:
                print("#    waveform option:", wvfrm_opt_out, file=file_out)
                print("#    waveform peak op:", wvfrm_p0_out, file=file_out)
                print("#    waveform time sc.:", wvfrm_t0_out, file=file_out)
                print("#    waveform shaping param.:", wvfrm_alpha_out, file=file_out)

            print("#", file=file_out)
            print("# Eigenray arrivals:", file=file_out)
            print("# incl [deg]	az [deg]	n_b	x_0 [km]	y_0 [km]	time [s]	cel [km/s]	turning ht [km]	inclination [deg]	back azimuth [deg]	trans. coeff. [dB]	absorption [dB]", file=file_out)
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
        
        if keep_eig_results and not os.path.isfile(output_id + ".arrivals.dat"):
            click.echo('\t' + "Copying eigenray results into " + output_id + ".arrivals.dat")
            os.system("cp " + eig_arrivals_file + " " + output_id + ".arrivals.dat")


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
@click.option("--prof-format", help="Option to specify the atmospheric file format", default=None)
@click.option("--reverse-winds", help="Option to reverse wind directions for back projection", default=None, type=bool)
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
@click.option("--cpu-cnt", help="Number of CPUs to use in analysis", default=None)
def run_sph_prop(config_file, atmo_file, atmo_prefix, grid_lats, grid_lons, incl_min, incl_max, incl_step, inclination, 
                az_min, az_max, az_step, azimuth, bounces, src_lat, src_lon, src_alt, write_rays, write_topo, freq, 
                abs_coeff, z_grnd, write_atmo, prof_format, reverse_winds, output_id, calc_amp, max_alt, max_rng,
                max_lat, min_lat, max_lon, min_lon, min_ds, max_ds, max_s, topo_file, topo_bl_wind, cpu_cnt):
    '''
    Run spherical atmosphere geometry ray tracing analysis for a point source

    \b
    Examples:
    \t infraga sph prop --atmo-file ToyAtmo.met
    \t infraga sph prop --atmo-file ToyAtmo.met --config-file example.cnfg
    \t infraga sph prop --atmo-prefix profs/example --grid-lats profs/example_lat.dat --grid-lons profs/example_lon.dat --src-lat 40.0 --src-lon -102.5 --azimuth -90.0 --z-grnd 1.0
    
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
    prof_format = define_param(user_config, 'GENERAL', 'prof_format', prof_format)
    reverse_winds = define_param(user_config, 'GENERAL', 'reverse_winds', reverse_winds)
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
    if cpu_cnt is None:
        cpu_cnt = 1 
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

    command = set_param(command, prof_format, "prof_format")
    if reverse_winds is not None:
        command = set_param(command, str(reverse_winds), "reverse_winds")

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

    click.echo(command)
    subprocess.run(shlex.split(command), shell=False)




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
@click.option("--prof-format", help="Option to specify the atmospheric file format", default=None)
@click.option("--reverse-winds", help="Option to reverse wind directions for back projection", default=None, type=bool)
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
@click.option("--cpu-cnt", help="Number of CPUs to use in analysis", default=None)
def run_sph_eig(config_file, atmo_file, atmo_prefix, grid_lats, grid_lons, incl_min, incl_max, bnc_min, bnc_max, bounces, 
                src_lat, src_lon, src_alt, rcvr_lat, rcvr_lon, verbose, iterations, damping, tolerance, az_dev_lim, 
                incl_step_min, incl_step_max, freq, abs_coeff, z_grnd, write_atmo, prof_format, reverse_winds, output_id,
                max_alt, max_rng, min_lat, max_lat, min_lon, max_lon, min_ds, max_ds, max_s, topo_file, topo_bl_wind, cpu_cnt):
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
    prof_format = define_param(user_config, 'GENERAL', 'prof_format', prof_format)
    reverse_winds = define_param(user_config, 'GENERAL', 'reverse_winds', reverse_winds)
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
    if cpu_cnt is None:
        cpu_cnt = 1 
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

    command = set_param(command, prof_format, "prof_format")
    if reverse_winds is not None:
        command = set_param(command, str(reverse_winds), "reverse_winds")

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

    click.echo(command)
    subprocess.run(shlex.split(command), shell=False)



@click.command('wnl_wvfrm', short_help="Run weakly non-linear waveform calculation")
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
@click.option("--expl-type", help="Explosion type ('chemical' vs. 'nuclear')", default="chemical")

@click.option("--wvfrm-ds", help="Burgers equation solver resolution [km]", default=None)
@click.option("--wvfrm-len", help="Waveform sample count (default = 2e13)", default=None)

@click.option("--freq", help="Frequency for Sutherland-Bass losses", default=None)
@click.option("--abs-coeff", help="Scaling coefficient for Sutherland-Bass losses", default=None)
@click.option("--z-grnd", help="Ground elevation (flat ground, km rel. sea level)", default=None)
@click.option("--write-atmo", help="Option to write atmosphere data (for QC)", default=None)
@click.option("--prof-format", help="Option to specify the atmospheric file format", default=None)
@click.option("--reverse-winds", help="Option to reverse wind directions for back projection", default=None, type=bool)
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
                wvfrm_yield, expl_type, wvfrm_ds, wvfrm_len, freq, abs_coeff, z_grnd, write_atmo, prof_format, reverse_winds, 
                output_id, max_alt, max_rng, max_lat, min_lat, max_lon, min_lon, min_ds, max_ds, max_s, topo_file, topo_bl_wind):
    '''
        Run weakly non-linear waveform analysis along a spherical atmospheric layer ray path.

    \b
    Examples:
    \t infraga sph wnl_wvfrm --atmo-file ToyAtmo.met --src-lat 30.0 --src-lon -100.0 --inclination 4.1314876 --azimuth -84.969455 --bounces 1 --wvfrm-yield 10e3
    \t infraga sph wnl_wvfrm --atmo-file ToyAtmo.met --config-file example.cnfg
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
    prof_format = define_param(user_config, 'GENERAL', 'prof_format', prof_format)
    reverse_winds = define_param(user_config, 'GENERAL', 'reverse_winds', reverse_winds)
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
        wvfrm_p0 = str(kg_op(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient, exp_type=expl_type))
        wvfrm_t0 = str(kg_ppd(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient, exp_type=expl_type))
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

    command = set_param(command, prof_format, "prof_format")
    if reverse_winds is not None:
        command = set_param(command, str(reverse_winds), "reverse_winds")
        
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

    click.echo(command)
    subprocess.run(shlex.split(command), shell=False)



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
@click.option("--expl-type", help="Explosion type ('chemical' vs. 'nuclear')", default="chemical")

@click.option("--wvfrm-ds", help="Burgers equation solver resolution [km]", default=None)
@click.option("--wvfrm-len", help="Waveform sample count (default = 2e13)", default=None)

@click.option("--freq", help="Frequency for Sutherland-Bass losses", default=None)
@click.option("--abs-coeff", help="Scaling coefficient for Sutherland-Bass losses", default=None)
@click.option("--z-grnd", help="Ground elevation (km rel. sea level)", default=None)
@click.option("--write-atmo", help="Option to write atmosphere data (for QC)", default=None)
@click.option("--prof-format", help="Option to specify the atmospheric file format", default=None)
@click.option("--reverse-winds", help="Option to reverse wind directions for back projection", default=None, type=bool)
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
@click.option("--cpu-cnt", help="Number of CPUs to use in analysis", default=None)

@click.option("--keep-eig-results", help="Keep eigenray arrivals", default=False)
def run_sph_eig_wvfrm(config_file, atmo_file, atmo_prefix, grid_lats, grid_lons, incl_min, incl_max, bnc_min, bnc_max, bounces, 
                src_lat, src_lon, src_alt, rcvr_lat, rcvr_lon, verbose, iterations, damping, tolerance, az_dev_lim, 
                incl_step_min, incl_step_max, wvfrm_file, wvfrm_opt, wvfrm_p0, wvfrm_t0, wvfrm_alpha, wvfrm_ref, wvfrm_out_step,
                wvfrm_yield, expl_type, wvfrm_ds, wvfrm_len, freq, abs_coeff, z_grnd, write_atmo, prof_format, reverse_winds, 
                output_id, max_alt, max_rng, min_lat, max_lat, min_lon, max_lon, min_ds, max_ds, max_s, topo_file, topo_bl_wind, 
                cpu_cnt, keep_eig_results):
    '''
    Run spherical atmospheric layer eigenray analysis to identify propagation paths connecting a specific source-receiver geometry and then compute weakly-nonlinear waveform predictions for each eigenray

    \b
    Examples:
    \t infraga sph eig_wvfrm --atmo-file ToyAtmo.met --src-lat 30.0 --src-lon -100.0 --rcvr-lat 30.25 --rcvr-lon -104.25 --bnc-max 1 --keep-eig-results True --wvfrm-yield 10e3
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
    atmo_file = define_param(user_config, 'GENERAL', 'atmo_file', atmo_file)

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
    prof_format = define_param(user_config, 'GENERAL', 'prof_format', prof_format)
    reverse_winds = define_param(user_config, 'GENERAL', 'reverse_winds', reverse_winds)
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
    if cpu_cnt is None:
        cpu_cnt = 1 
    cpu_cnt = int(cpu_cnt)

    # Set parameter values for header output
    def set_header_val(param, default):
        if param is None:
            return default
        else:
            return param

    src_lat_out = set_header_val(src_lat, 30.0)
    src_lon_out = set_header_val(src_lon, 0.0)
    src_alt_out = set_header_val(src_alt, 0.0)

    rcvr_lat_out = set_header_val(rcvr_lat, 30.0)
    rcvr_lon_out = set_header_val(rcvr_lon, 10.0)

    incl_min_out = set_header_val(incl_min, 0.5)
    incl_max_out = set_header_val(incl_max, 45.0)
    incl_step_max_out = set_header_val(incl_step_max, 0.1)

    bnc_min_out = set_header_val(bnc_min, 0)
    bnc_max_out = set_header_val(bnc_max, 0)
    if bounces is not None:
        bnc_min_out = bounces
        bnc_max_out = bounces 

    z_grnd_out = set_header_val(z_grnd, 0.0)
    damping_out = set_header_val(damping, 1.0e-3)
    max_rng_out = set_header_val(max_rng, 2500.0)

    wvfrm_opt_out = set_header_val(wvfrm_opt, 'impulse')
    wvfrm_p0_out = set_header_val(wvfrm_p0, 10.0)
    wvfrm_t0_out = set_header_val(wvfrm_t0, 1.0)
    wvfrm_alpha_out = set_header_val(wvfrm_alpha, 1.0)


    # Check if eigenray analysis is already done
    if output_id is not None:
        eig_arrivals_file = output_id + ".arrivals.dat"
    else:
        if atmo_prefix is not None:
            eig_arrivals_file = atmo_prefix + ".arrivals.dat"
            output_id = atmo_prefix
        else:
            eig_arrivals_file = os.path.splitext(atmo_file)[0] + ".arrivals.dat"
            output_id = os.path.splitext(atmo_file)[0]
    
    eig_results_check = False   
    if os.path.isfile(eig_arrivals_file):
        eig_results_option_check = False 
        eig_results_atmo_check = False
        eig_results_src_check = False
        eig_results_rcvr_check = False 

        arrivals_file = open(eig_arrivals_file, 'r')
        for line in arrivals_file:
            if "-eig_search" in line:
                eig_results_option_check = True
            if "profile" in line:
                eig_results_atmo_check = line.split(' ')[-1].strip('\n') == atmo_file
            if "source location" in line:
                ref_loc = np.array([float(item) for item in line.replace(",","").split(' ')[-3:]])
                eig_results_src_check = np.allclose(np.array([src_lat_out, src_lon_out, src_alt_out], dtype=float), ref_loc)
            if "receiver location" in line:
                ref_loc = [float(item) for item in line.replace(",","").split(' ')[-3:]]
                eig_results_rcvr_check = np.allclose(np.array([rcvr_lat_out, rcvr_lon_out, 0.0], dtype=float), ref_loc)

        if not eig_results_option_check:
            click.echo("Arrivals file detected, but not from eigenray analysis.")
        if not eig_results_atmo_check:
            click.echo("Arrivals file detected, but atmosphere file doesn't agree.")
        if not eig_results_src_check:
            click.echo("Arrivals file detected, but source location doesn't agree.")
        if not eig_results_rcvr_check:
            click.echo("Arrivals file detected, but receiver location doesn't agree.")

        eig_results_check = eig_results_option_check and eig_results_atmo_check and eig_results_src_check and eig_results_rcvr_check

    with tempfile.TemporaryDirectory(prefix='infraga_') as tmpdirname:
        print('Created temp directory:', tmpdirname)
        temp_path = tmpdirname + "/temp"

        if eig_results_check:
            click.echo('\nEigenray results found.\nSkipping to waveform calculation...')
        else:
            click.echo("Running eigenray analysis...")
            eig_arrivals_file = temp_path + ".arrivals.dat"

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

            command = set_param(command, prof_format, "prof_format")

            if reverse_winds is not None:
                command = set_param(command, str(reverse_winds), "reverse_winds")
                
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

            command = command + " output_id=" + temp_path

            click.echo(command)
            subprocess.run(shlex.split(command), shell=False)

        # Analyze eigenray results and compute waveform contributions...
        click.echo('\n' + "Loading eigenray results from " + eig_arrivals_file)
        eig_results = np.loadtxt(eig_arrivals_file)

        if len(eig_results) > 0:
            eig_results = np.atleast_2d(eig_results)

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
            command_list = []
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
                    wvfrm_p0 = str(kg_op(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient, exp_type=expl_type))
                    wvfrm_t0 = str(kg_ppd(float(wvfrm_yield), float(wvfrm_ref), p_amb=p_ambient, T_amb=T_ambient, exp_type=expl_type))
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

                command = set_param(command, prof_format, "prof_format")

                if reverse_winds is not None:
                    command = set_param(command, str(reverse_winds), "reverse_winds")
                                
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

                command = command + " output_id=" + temp_path + "-" + str(n)
                command = command + " > /dev/null"

                command_list = command_list + [command]

            click.echo('\n' + "Computing waveforms...")
            for j in range(0, len(command_list), cpu_cnt):              
                procs_list = [subprocess.Popen(cmd, shell=True) for cmd in command_list[j:j + cpu_cnt]]
                for proc in procs_list:
                    proc.communicate()
                    proc.wait()

            file_out = open(output_id + ".eigenrays.dat", 'w')
            print("# 'infraga sph eig_wvfrm' eigenray results", '\n#', file=file_out)
            if atmo_prefix is not None:
                print("# 	profile: " + atmo_prefix, file=file_out)
            else:
                print("# 	profile: " + atmo_file, file=file_out)
            print("# 	source location (lat, lon, alt): " + str(src_lat_out) + ", " + str(src_lon_out) + ", " + str(src_alt_out) , file=file_out)
            print("# 	receiver location (lat, lon, alt): " + str(rcvr_lat_out) + ", " + str(rcvr_lon_out) + ", 0.0", file=file_out)
            print("# 	inclination range: " + str(incl_min_out) + ", " + str(incl_max_out), file=file_out)
            print("#    inclination step max:", incl_step_max_out, file=file_out)
            print("# 	bounces: " + str(bnc_min_out) + ", " + str(bnc_max_out), file=file_out)
            if topo_file is not None:
                print("# 	terrain file: " + str(topo_file), file=file_out)        
            else:
                print("# 	ground elevation: " + str(z_grnd_out), file=file_out)
            print("# 	damping:", damping_out, file=file_out)
            print("# 	range max:", max_rng_out, '\n#', file=file_out)
            print("# lat [deg]	lon [deg]	z [km]	trans. coeff. [dB]	absorption [dB]	time [s]", file=file_out)
            
            wvfrms = []
            t_lims = [np.inf, 0.0]
            for n, line in enumerate(eig_results):

                temp = np.loadtxt(temp_path + "-" + str(n) + ".wvfrm_out.dat")
                wvfrms += [interp1d(temp[:, 0], temp[:, 1], bounds_error=False, fill_value=0.0,  kind='cubic')]
                t_lims[0] = min(t_lims[0], temp[0][0])
                t_lims[1] = max(t_lims[1], temp[-1][0])
                dt = abs(temp[1][0] - temp[0][0])

                temp = np.loadtxt(temp_path + "-" + str(n) + ".raypaths.dat")
                for line in temp:
                    print(*line, file=file_out)
                print('\n', file=file_out)
            
            file_out.close()

            # combine waveforms and write to file
            click.echo("Interpolating and merging waveforms...")
            click.echo('\t' + "Eigenrays written into " + output_id + ".eigenrays.dat")
            click.echo('\t' + "Arrival waveform written into " + output_id + ".wvfrms.dat")

            t_vals = np.arange(t_lims[0], t_lims[1], dt)
            file_out = open(output_id + ".wvfrms.dat", 'w')
            print("# 'infraga sph eig_wvfrm' waveform results", '\n#', file=file_out)
            if atmo_prefix is not None:
                print("# 	profile: " + atmo_prefix, file=file_out)        
            else:
                print("# 	profile: " + atmo_file, file=file_out)
            print("# 	source location (lat, lon, alt): " + str(src_lat_out) + ", " + str(src_lon_out) + ", " + str(src_alt_out) , file=file_out)
            print("# 	receiver location (lat, lon, alt): " + str(rcvr_lat_out) + ", " + str(rcvr_lon_out) + ", 0.0", file=file_out)
            print("# 	inclination range: " + str(incl_min_out) + ", " + str(incl_max_out), file=file_out)
            print("#    inclination step max:", incl_step_max_out, file=file_out)
            print("#    bounces: " + str(bnc_min_out) + ", " + str(bnc_max_out), file=file_out)
            if topo_file is not None:
                print("#    terrain file: " + str(topo_file), file=file_out)        
            else:
                print("#    ground elevation: " + str(z_grnd_out), file=file_out)
            print("#    damping:", damping_out, file=file_out)
            print("#    range max:", max_rng_out, '\n#', file=file_out)

            # Waveform calculation parameters
            if wvfrm_ref is not None:
                print("#    waveform reference distance:", wvfrm_ref, file=file_out)
            else:
                print("#    waveform reference distance: 1.0", file=file_out)

            if wvfrm_yield:
                print("#    waveform source yield:", wvfrm_yield, file=file_out)
            else:
                print("#    waveform option:", wvfrm_opt_out, file=file_out)
                print("#    waveform peak op:", wvfrm_p0_out, file=file_out)
                print("#    waveform time sc.:", wvfrm_t0_out, file=file_out)
                print("#    waveform shaping param.:", wvfrm_alpha_out, file=file_out)

            print("#", file=file_out)
            print("# Eigenray arrivals:", file=file_out)
            print("# incl [deg]	az [deg]	n_b	lat_0 [deg]	lon_0 [deg]	time [s]	cel [km/s]	turning ht [km]	inclination [deg]	back azimuth [deg]	trans. coeff. [dB]	absorption [dB]", file=file_out)
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
            click.echo('\n' + "No waveforms to compute.")
        
        if keep_eig_results and not os.path.isfile(output_id + ".arrivals.dat"):
            click.echo('\t' + "Copying eigenray results into " + output_id + ".arrivals.dat")
            os.system("cp " + eig_arrivals_file + " " + output_id + ".arrivals.dat")


def mach_cone_arrrival_header(atmo_file, traj_file, cone_resol, traj_resol, max_rng, lat_bnds=None, lon_bnds=None):
    header_text = "# infraga supersonic_source summary:" + '\n'
    header_text = header_text + "    profile: " + atmo_file + '\n'
    header_text = header_text + "    trajectory: " + traj_file + '\n'
    header_text = header_text + "    cone resolution: " + str(cone_resol) + '\n'
    header_text = header_text + "    trajectory resolution: " + str(traj_resol) + '\n'
    header_text = header_text + "    max_rng: " + str(max_rng) + '\n'
    if lat_bnds is not None:
        header_text = header_text + "    lat_bnds: (" + str(lat_bnds[0]) + ", " + str(lat_bnds[1]) + ')\n'
    if lon_bnds is not None:
        header_text = header_text + "    lon_bnds: (" + str(lon_bnds[0]) + ", " + str(lon_bnds[1]) + ')\n'

    header_text = header_text + '\n\n' + "incl [deg]	az [deg]	n_b	lat_0 [deg]	lon_0 [deg]	time [s]	cel [km/s]	turning ht [km]	inclination [deg]	back azimuth [deg]	trans. coeff. [dB]"
    header_text = header_text + " src_lat [deg]   src_lon [deg]   src_alt [deg]"

    return header_text


@click.command('supersonic', short_help="Run Mach cone source for a supersonic source trajectory")
@click.option("--config-file", help="Configuration file for simulation", default=None)
@click.option("--atmo-file", help="Atmosphere file", default=None)
@click.option("--atmo-prefix", help="Atmosphere file prefix (range dependent)", default=None)
@click.option("--grid-lats", help="Atmosphere grid latitudes file", default=None)
@click.option("--grid-lons", help="Atmosphere grid longitudes file", default=None)

@click.option("--trajectory", help="Trajectory file (time : lat : lon : alt)", default=None)
@click.option("--cone-resol", help="Mach cone angular resolution", default=None)
@click.option("--bounces", help="Number of ground reflections (bounces) to consider", default=None)
@click.option("--traj-step", help="Trajectory stepping factor (traj[::j]), default=1", default=1)

@click.option("--freq", help="Frequency for Sutherland-Bass losses", default=None)
@click.option("--abs-coeff", help="Scaling coefficient for Sutherland-Bass losses", default=None)
@click.option("--z-grnd", help="Ground elevation (km rel. sea level)", default=None)
@click.option("--write-atmo", help="Option to write atmosphere data (for QC)", default=None, type=bool)
@click.option("--write-rays", help="Option to write ray paths to file", default=None, type=bool)
@click.option("--prof-format", help="Option to specify the atmospheric file format", default=None)
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
@click.option("--cpu-cnt", help="Number of CPUs to use in analysis", default=None)
@click.option("--local-temp-dir", help="Local temporary directory for results", default=None)
def run_sph_supersonic(config_file, atmo_file, atmo_prefix, grid_lats, grid_lons, trajectory, cone_resol, traj_step,
                        bounces, freq, abs_coeff, z_grnd, write_atmo, write_rays, prof_format, output_id, max_alt, 
                        max_rng, min_lat, max_lat, min_lon, max_lon, min_ds, max_ds, max_s, topo_file, topo_bl_wind, 
                        cpu_cnt, local_temp_dir):
    
    '''
    Run spherical atmospheric layer eigenray analysis to identify propagation paths connecting a specific source-receiver geometr yand then compute weakly-nonlinear waveform predictions for each eigenray

    \b
    Examples:
    \t infraga sph supersonic --atmo-file G2S_example.met --trajectory trajectories/ballistic_traj.dat --traj-step 6 --cpu-cnt 12 --output-id ballistic --local-temp-dir ballistic_temp

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

    # Set supersonic source specific parameters
    trajectory = define_param(user_config, 'SUPERSONIC', 'trajectory', trajectory)
    cone_resol = define_param(user_config, 'SUPERSONIC', 'cone_resol', cone_resol)
    traj_step = define_param(user_config, 'SUPERSONIC', 'traj_step', traj_step)
    bounces = define_param(user_config, 'SUPERSONIC', 'bounces', bounces)
    local_temp_dir = define_param(user_config, 'SUPERSONIC', 'local_temp_dir', local_temp_dir)
   
    # Set general parameters
    freq = define_param(user_config, 'GENERAL', 'freq', freq)
    abs_coeff = define_param(user_config, 'GENERAL', 'abs_coeff', abs_coeff)

    write_atmo = define_param(user_config, 'GENERAL', 'write_atmo', write_atmo)
    write_rays = define_param(user_config, 'GENERAL', 'write_rays', write_rays)
    prof_format = define_param(user_config, 'GENERAL', 'prof_format', prof_format)
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
    if cpu_cnt is None:
        cpu_cnt = 1 
    cpu_cnt = int(cpu_cnt)

    if traj_step is not None:
        traj_step = int(max(1, traj_step))

    if output_id is None:
        if atmo_file is not None:
            output_id = atmo_file[:-4]
        else:
            output_id = atmo_prefix
        
    with tempfile.TemporaryDirectory(prefix='infraga_') as tmpdirname:
        if local_temp_dir is not None:
            click.echo("Writing individual Mach cone source results into " + local_temp_dir)
            if not os.path.isdir(local_temp_dir):
                os.mkdir(local_temp_dir)
            tmpdirname = local_temp_dir
        else:
            click.echo('Created temp directory:', tmpdirname)

        # read in atmospheric file and define sound speed
        if atmo_file is not None:
            g2s = np.loadtxt(atmo_file)
        else:
            g2s = np.loadtxt(atmo_prefix + "0.met")

        sndspd = interp1d(g2s[:, 0], np.sqrt(0.14 * g2s[:, 5] / g2s[:, 4]) * 1.0e-3) 

        # read in trajectory and interpolate
        traj = np.loadtxt(trajectory)
        time_vals = traj[:, 0]
        lat_vals = traj[:, 1]
        lon_vals = traj[:, 2]
        alt_vals = traj[:, 3]

        # compute the mach number, attack angle, and azimuth angle
        dz_dt = np.gradient(alt_vals) / np.gradient(time_vals)
        dlat_dt = np.gradient(np.radians(lat_vals)) / np.gradient(time_vals)
        dlon_dt = np.gradient(np.radians(lon_vals)) / np.gradient(time_vals)

        r0 = 6370.0
        dx_dt = (r0 + alt_vals) * dlat_dt
        dy_dt = (r0 + alt_vals) * np.cos(np.radians(lat_vals)) * dlon_dt
        ds_dt = np.sqrt(dx_dt**2 + dy_dt**2 + dz_dt**2)

        mach = ds_dt / sndspd(alt_vals)
        attack = np.degrees(np.arcsin(dz_dt / ds_dt))
        azimuth = np.degrees(np.arctan2(np.cos(np.radians(dlat_dt)) * np.radians(dlon_dt), np.radians(dlat_dt)))

        # Plot trajectory info and step through mach cone source instances
        _, ax = plt.subplots(3, 2, figsize=(15, 6), sharex=True)

        for n in range(2):
            ax[2][n].set_xlabel("Time [s]")

        ax[0][0].set_ylabel("Altitude [km]")
        ax[1][0].set_ylabel("Latitude [km]")
        ax[2][0].set_ylabel("Longitude [km]")

        ax[0][1].set_ylabel("Mach [-]")
        ax[1][1].set_ylabel("Attack [deg]")
        ax[2][1].set_ylabel("Azimuth [deg]")

        ax[1][1].set_ylim([-90, 90])
        ax[2][1].set_ylim([-180.0, 180.0])

        ax[0][0].plot(time_vals, alt_vals, '-k', linewidth=4)
        ax[1][0].plot(time_vals, lat_vals, '-k', linewidth=4)
        ax[2][0].plot(time_vals, lon_vals, '-k', linewidth=4)

        ax[0][1].plot(time_vals, mach, '--k', linewidth=2)
        ax[1][1].plot(time_vals, attack, '--k', linewidth=2)
        ax[2][1].plot(time_vals, azimuth, '--k', linewidth=2)

        ax[0][1].plot(time_vals[mach > 1], mach[mach > 1], '-k', linewidth=4)
        ax[1][1].plot(time_vals[mach > 1], attack[mach > 1], '-k', linewidth=4)
        ax[2][1].plot(time_vals[mach > 1], azimuth[mach > 1], '-k', linewidth=4)

        plt.show(block=False)

        # cycle through trajectory and 
        for jj in range(0, len(time_vals), traj_step):
            ax[0][0].plot([time_vals[jj]], [alt_vals[jj]], 'or', markersize=3)
            ax[1][0].plot([time_vals[jj]], [lat_vals[jj]], 'or', markersize=3)
            ax[2][0].plot([time_vals[jj]], [lon_vals[jj]], 'or', markersize=3)

            ax[0][1].plot([time_vals[jj]], [mach[jj]],  'or', markersize=3)
            ax[1][1].plot([time_vals[jj]], [attack[jj]],  'or', markersize=3)
            ax[2][1].plot([time_vals[jj]], [azimuth[jj]],  'or', markersize=3)

            plt.pause(0.001)

            if mach[jj] > 1.0:
                if not os.path.isfile(tmpdirname +"/t0_" + "%03f" % (time_vals[jj]) + ".arrivals.dat"):

                    # Build run command
                    if cpu_cnt < 2:
                        command = bin_path + "infraga-sph"
                    else:
                        command = "mpirun -np " + str(cpu_cnt) + " " + bin_path + "infraga-accel-sph"

                    if atmo_prefix is not None:
                        command = command + "-rngdep"

                    if atmo_file is not None:
                        command = command + " -mach_cone " + atmo_file
                    elif atmo_prefix is not None and grid_lats is not None and grid_lons is not None:
                        command = command + " -mach_cone " + atmo_prefix + " " + grid_lats + " " + grid_lons
                    else:
                        click.echo("Simulation requires either an '--atmo-file' or --atmo-prefix' with grid info (--grid-lats and --grid-lons)")
                        return 0              

                    command = command + " output_id=" + tmpdirname + "/t0_" + "%03f" % (time_vals[jj])

                    command = set_param(command, str(mach[jj]), "src_mach")
                    command = set_param(command, str(attack[jj]), "src_attack")
                    command = set_param(command, str(azimuth[jj]), "src_az")
                    command = set_param(command, cone_resol, "cone_resol")
                    command = set_param(command, bounces, "bounces")
                    
                    command = set_param(command, str(lat_vals[jj]), "src_lat")
                    command = set_param(command, str(lon_vals[jj]), "src_lon")
                    command = set_param(command, str(alt_vals[jj]), "src_alt")

                    if write_rays is not None:
                        command = set_param(command, str(write_rays), "write_rays")

                    command = set_param(command, freq, "freq")
                    command = set_param(command, abs_coeff, "abs_coeff")

                    if write_atmo is not None:
                        command = set_param(command, str(write_atmo), "write_atmo")

                    command = set_param(command, prof_format, "prof_format")

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

                    click.echo('\n' + command)
                    subprocess.run(shlex.split(command), shell=False)
                else:
                    click.echo("infraGA/GeoAc results already exist.  Skipping t0 = %03f..." % (time_vals[jj]))

        plt.close()
        
        # Merge output files...
        click.echo('\n' + "Applying source-time delays and merging results...")
        raypath_files = []
        arrivals_files = []
        for file in np.sort(os.listdir(tmpdirname)):
            if fnmatch.fnmatch(file, "t0_*.raypaths.dat"):
                raypath_files = raypath_files + [file]
            elif fnmatch.fnmatch(file, "t0_*.arrivals.dat"):
                arrivals_files = arrivals_files + [file]

        if write_rays:
            click.echo('\t' + "Merging ray path files...")
            raypaths = np.loadtxt(tmpdirname + "/" + raypath_files[0])
            raypaths[:, 5] = raypaths[:, 5] + float(raypath_files[0][3:11])
            for file in raypath_files[1:]:
                temp = np.loadtxt(tmpdirname + "/" + file)
                temp[:, 5] = temp[:, 5] + float(file[3:11])
                raypaths = np.vstack((raypaths, temp))
            np.savetxt(output_id + ".raypaths.dat", raypaths)

        click.echo('\t' + "Merging arrivals files...")
        # Remove files with no data (e.g., arrival files in which no rays reach the ground)
        valid_arrivals = []
        for file_name in arrivals_files:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                temp = np.atleast_2d(np.loadtxt(tmpdirname + "/" + file_name))
                if temp.shape[1] > 0:
                    valid_arrivals = valid_arrivals + [file_name]

        arrivals = np.atleast_2d(np.loadtxt(tmpdirname + "/" + valid_arrivals[0]))
        for file in valid_arrivals[1:]:
            temp = np.atleast_2d(np.loadtxt(tmpdirname + "/" + file))
            temp[:, 5] = temp[:, 5] + float(file[3:11])
            arrivals = np.vstack((arrivals, temp))
        np.savetxt(output_id + ".arrivals.dat", arrivals, header=mach_cone_arrrival_header(atmo_file, trajectory, cone_resol, traj_step, max_rng))



