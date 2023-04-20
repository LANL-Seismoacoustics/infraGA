#!which python
"""
plot_on_map.py

Methods to plot spherical coordinate
results from infraGA/GeoAc (infraga-sph)
onto a map.

usage: python plot_on_map.py output.arrivals.dat figure.png

Fill in 

Author: pblom@lanl.gov    
"""

import os 
import click
import sys
import fnmatch

import numpy as np

import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import matplotlib.ticker as mticker

from mpl_toolkits.axes_grid1 import make_axes_locatable

import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from pyproj import Geod


marker_size = 3.0
map_proj = cartopy.crs.PlateCarree()
resol = '100m'

sph_proj = Geod(ellps='sphere')

def use_offline_maps(self, pre_existing_data_dir, turn_on=True):
    # call this function to initialize the use of offline maps.  turn_on will initialize the pre_existing_data_directory
    if turn_on:
        cartopy.config['pre_existing_data_dir'] = pre_existing_data_dir
    else:
        cartopy.config['pre_existing_data_dir'] = ""


@click.command('atmo', short_help="Visualize information about an atmospheric atmo_file")
@click.option("--atmo-file", help="Atmospheric atmo_file file")
@click.option("--max-alt", help="Maximum altitude for analysis (default: 120 km)", default=None, type=float)
@click.option("--format", help="Atmospheric atmo_file format (default: 'zTuvdp')", default='zTuvdp')
def plot_atmo(atmo_file, max_alt, format):
    '''
    Visualize the sound speed, wind fields, and effective sound speed ratio ducting information for an atmospheric atmo_file


    \b
    Examples:
    \t infraga plot atmo --atmo-file examples/ToyAtmo.met
    \t infraga plot atmo --atmo-file examples/G2S_example.met

    '''

    atmo = np.loadtxt(atmo_file)
    z = atmo[:, format.find('z')]
    u = atmo[:, format.find('u')]
    v = atmo[:, format.find('v')]
    c = np.sqrt(0.14 * atmo[:, format.find('p')] / atmo[:, format.find('d')])

    grnd_ht = z[0]
    for line in open(atmo_file, 'r'):
        if "Ground Height" in line:
            grnd_ht = float(line[18:])
            break

    if max_alt is None:
        max_alt = z[-1]
    else:
        max_alt = float(max_alt)

    ht_mask = np.logical_and(grnd_ht <= z, z <= max_alt)

    f, ax = plt.subplots(1, 3, gridspec_kw={'width_ratios': [1, 1, 4]}, figsize=(12, 5))
    
    ax[0].grid(color='k', linestyle='--', linewidth=0.5)
    ax[1].grid(color='k', linestyle='--', linewidth=0.5)

    ax[0].set_ylim(grnd_ht, max_alt)
    ax[0].set_ylabel("Altitude [km]")
    ax[0].set_xlabel("Sound Speed [m/s]")

    ax[1].set_ylim(grnd_ht, max_alt)
    ax[1].set_xlabel("Wind Speed [m/s]")

    ax[1].yaxis.set_ticklabels([])

    ax[2].yaxis.set_label_position("right")
    ax[2].yaxis.tick_right()

    ax[2].set_xlim(-180.0, 180.0)
    ax[2].set_ylim(0.0, 50.0)
    ax[2].set_xticks((-180.0, -135.0, -90.0, -45.0, 0.0, 45.0, 90.0, 135.0, 180.0))
    ax[2].set_xticklabels(["S", "SW", "W", "NW", "N", "NE", "E", "SE", "S"])
    ax[2].set_xlabel("Propagation Direction")
    ax[2].set_ylabel("Inclination [deg]")

    ax[0].plot(c[ht_mask], z[ht_mask], '-k', linewidth=3.0)
    ax[1].plot(u[ht_mask], z[ht_mask], '-b', linewidth=3.0, label='Zonal')
    ax[1].plot(v[ht_mask], z[ht_mask], '-r', linewidth=3.0, label='Merid.')
    ax[1].legend(fontsize='small')

    incl_vals = np.arange(0.0, 50.0, 0.2)
    for az in np.arange(-180.0, 180.0, 1.0):
        ceff = (c + u * np.sin(np.radians(az)) + v * np.cos(np.radians(az)))[ht_mask]
        refract_ht = [z[ht_mask][np.min(np.where((ceff / ceff[0]) * np.cos(np.radians(incl)) > 1.0)[0])] if len(np.where((ceff / ceff[0]) * np.cos(np.radians(incl)) > 1.0)[0]) > 0 else max_alt for incl in incl_vals]
        sc = ax[2].scatter([az] * len(refract_ht), incl_vals, c=refract_ht, cmap=cm.jet_r, marker="s", s=5.0, alpha=0.75, edgecolor='none', vmin=grnd_ht, vmax=120.0)

    f.colorbar(sc, ax=[ax[2]], location='top', label="Estimated Refraction Altitude [km]")

    plt.show()


@click.command('azimuthal', short_help="Visualize results for a single azimuth simulation")
@click.option("--atmo-file", help="Atmospheric atmo_file file")
@click.option("--arrivals", help="Arrivals file from the simulation (optional)", default=None)
@click.option("--ray-paths", help="Ray path file from the simulation (optional)", default=None)
@click.option("--y-axis-option", help="Lower axis option (see usage info below)", default='inclination')
@click.option("--cmap-option", help="Low axis cmap option (see usage info below)", default=None)
@click.option("--reduced-tm-vel", help="Reference velocity for reduced time option", default=300.0)
@click.option("--tr-vel-ref", help="Reference velocity for trace velocity calculation", default=330.0)
@click.option("--plot-amplitudes", help="Option to plot amplitude along rays", default=True)
@click.option("--terrain-profile", help="Terrain file output from simulation", default=None)
@click.option("--figure-out", help="Name of output figure", default=None)
def plot_azimuthal(atmo_file, arrivals, ray_paths, y_axis_option, cmap_option, reduced_tm_vel, tr_vel_ref, plot_amplitudes, terrain_profile, figure_out):
    '''
    Visualize propagation results for a single azimuthal angle simulation

    \b
    Plotting Options:
    \t inclination \t\t Launch inclination angle
    \t celerity  \t\t Arrival celerity (horizontal group velocity)
    \t reduced-time \t\t Reduced arrival time (relative to --reduced-tm-vel)
    \t turning-ht \t\t Turning height
    \t trace-velocity \t Trace velocity
    \t back-azimuth \t\t Back azimuth (not available for 2d geometry)
    \t amplitude \t\t Transport equation + absorption losses

    \b
    Examples:
    \t infraga plot azimuthal --atmo-file ToyAtmo.met --y-axis-option celerity
    \t infraga plot azimuthal --atmo-file ToyAtmo.met --y-axis-option reduced-time --cmap-option trace-velocity

    '''
    print("Loading atmospheric data and simulation results for " + atmo_file)

    if arrivals is not None:
        print('\t' + "Loading specified arrivals file: " + arrivals)
        arrivals_file = arrivals
    else:
        arrivals_file = os.path.splitext(atmo_file)[0] + ".arrivals.dat"
    
    if ray_paths is not None:
        print('\t' + "Loading specified ray paths file: " + ray_paths)
        raypaths_file = ray_paths
    else:
        raypaths_file = os.path.splitext(atmo_file)[0] + ".raypaths.dat"

    if not os.path.isfile(arrivals_file):
        print('\t' + "Arrivals file (" + arrivals_file + ") not found.")
        return 0
    
    if not os.path.isfile(raypaths_file):
        print('\t' + "Ray paths file (" + raypaths_file + ") not found.")
        return 0

    # check geometry (2d, 3d, or sph)
    for line in open(arrivals_file):
        if "infraga-" in line:
            geom_info = line
        elif "azimuth:" in line:
            az_info = line
        elif "source" in line:
            src_info = line
    
    if "2d" in geom_info:
        geom = '2d'
    elif "3d" in geom_info:
        geom = '3d'
    elif "sph" in geom_info:
        geom = 'sph'

    # Extract propagation azimuth and source location
    if geom == '3d' or geom == 'sph':
        prop_az = float(az_info.split(',')[-2])
        src_loc = np.array([float(val) for val in src_info.split(":")[-1].split(",")[:2]])
    else:
        prop_az = float(az_info.split(":")[-1])
        src_loc = None

    print('\n' + "Extracted Info:")
    print('\t' + "Geometry:", geom)
    print('\t' + "Propagation azimuth:", prop_az)
    if src_loc is not None:
        print('\t' + "Source Location:", src_loc)

    atmo_data = np.loadtxt(atmo_file)
    arr_data = np.loadtxt(arrivals_file)
    ray_data = np.loadtxt(raypaths_file)
    if terrain_profile is not None:
        terrain_data = np.loadtxt(terrain_profile)

    if geom == '2d':
        ray_rngs = ray_data[:, 0]
        ray_alts = ray_data[:, 1]
        ray_amps = ray_data[:, 2] + ray_data[:, 3]

        indices = np.flatnonzero(np.gradient(ray_data[:, 4]) < 0.0)

        arr_rngs = arr_data[:, 3]
        incl_vals = arr_data[:, 0]
        tm_vals = arr_data[:, 4]
        cel_vals = arr_data[:, 5] * 1.0e3
        turn_ht_vls = arr_data[:, 6]
        arr_incl_vals = arr_data[:, 7]
        amp_vals = arr_data[:, 8] + arr_data[:, 9]

    elif geom == '3d' or geom == 'sph':
        if geom == '3d':
            ray_rngs = np.sqrt((src_loc[0] - ray_data[:, 0])**2 + (src_loc[1] - ray_data[:, 1])**2)
            arr_rngs = np.sqrt((src_loc[0] - arr_data[:, 3])**2 + (src_loc[1] - arr_data[:, 4])**2)
            if terrain_profile is not None:
                topo_rngs = np.sqrt((src_loc[0] - terrain_data[:, 0])**2 + (src_loc[1] - terrain_data[:, 1])**2)
        else:
            ray_rngs = sph_proj.inv([src_loc[1]] * len(ray_data), [src_loc[0]] * len(ray_data), ray_data[:, 1], ray_data[:, 0])[2] * 1.0e-3
            arr_rngs = sph_proj.inv([src_loc[1]] * len(arr_data), [src_loc[0]] * len(arr_data), arr_data[:, 4], arr_data[:, 3])[2] * 1.0e-3
            if terrain_profile is not None:
                topo_rngs = sph_proj.inv([src_loc[1]] * len(terrain_data), [src_loc[0]] * len(terrain_data), terrain_data[:, 1], terrain_data[:, 0])[2] * 1.0e-3

        ray_alts = ray_data[:, 2]
        ray_amps = ray_data[:, 3] + ray_data[:, 4]
        indices = np.flatnonzero(np.gradient(ray_data[:, 5]) < 0.0)

        incl_vals = arr_data[:, 0]
        tm_vals = arr_data[:, 5]
        cel_vals = arr_data[:, 6] * 1.0e3
        turn_ht_vls = arr_data[:, 7]
        arr_incl_vals = arr_data[:, 8]
        back_az_vals = arr_data[:, 9]
        amp_vals = arr_data[:, 10] + arr_data[:, 11]

    print('\n' + "Plotting...")
    fig = plt.figure(figsize=(11, 5), layout='constrained')
    spec = fig.add_gridspec(2, 7)

    ax0 = fig.add_subplot(spec[0, 0])
    ax0.set_ylabel("Altitude [km]")
    ax0.set_xlabel("Sound Speed [m/s]")
    ax0.xaxis.set_label_position('top')
    ax0.xaxis.set_ticks_position('top')

    print('\t' + "Atmospheric data...")
    ax0.plot(np.sqrt(0.14 * atmo_data[:, 5] / atmo_data[:, 4]), atmo_data[:, 0], '--k', linewidth=1.5)
    ax0.plot(np.sqrt(0.14 * atmo_data[:, 5] / atmo_data[:, 4]) + atmo_data[:, 2] * np.sin(np.radians(prop_az)) + atmo_data[:, 3] * np.cos(np.radians(prop_az)), atmo_data[:, 0], '-k', linewidth=3.0)

    ax1 = fig.add_subplot(spec[0, 1:], sharey=ax0)
    ax1.set_xlabel("Range [km]")
    ax1.xaxis.set_label_position('top')
    ax1.xaxis.set_ticks_position('top')

    ax1.set_xlim(0.0, np.max(ray_rngs))    
    ax1.set_ylim(np.min(ray_alts), np.max(ray_alts))    

    if plot_amplitudes:
        print('\t' + "Ray path info with amplitudes...")
        sc = ax1.scatter(ray_rngs, ray_alts, c=ray_amps, cmap=cm.jet, s=0.05, vmax=-20.0, vmin=-120.0)
        plt.colorbar(sc, label="Amp. (rel. 1 km) [dB]")

    else:
        print('\t' + "Ray path info without amplitudes...")
        ax1.plot(ray_rngs[:indices[0]], ray_alts[:indices[0]], '-k', linewidth=0.75)
        for n, j in enumerate(indices):
            ax1.plot(ray_rngs[indices[n - 1]:j], ray_alts[indices[n - 1]:j], '-k', linewidth=0.75)

    if terrain_profile is not None:
        ax1.fill_between(topo_rngs, terrain_data[:, 2], 0.0, color='k', alpha=1.0)

    ax2 = fig.add_subplot(spec[1, 1:], sharex=ax1)
    ax2.set_xlabel("Range [km]")

    if y_axis_option == 'inclination':
        print('\t' + "Launch inclination info...")
        ax2.set_ylabel("Launch Inclination [deg]")
        plot_vals = arr_incl_vals
    elif y_axis_option == 'celerity':
        print('\t' + "Arrival celerity info...")
        ax2.set_ylabel("Celerity [m/s]")
        plot_vals = cel_vals
    elif y_axis_option == "reduced-time":
        print('\t' + "Arrival reduced time relative to " + str(reduced_tm_vel) + " m/s...")
        ax2.set_ylabel("Reduced Time (rel. " + str(int(reduced_tm_vel)) + " m/s) [s]")
        plot_vals = tm_vals - arr_rngs / (reduced_tm_vel * 1.0e-3)
    elif y_axis_option == "turning-ht":
        print('\t' + "Turning height info...")
        ax2.set_ylabel("Turning Height [km]")
        plot_vals = turn_ht_vls
    elif y_axis_option == "trace-velocity":
        print('\t' + "Trace velocity info...")
        ax2.set_ylabel("Trace Velocity [m/s]")
        plot_vals = tr_vel_ref / np.cos(np.radians(arr_incl_vals))
    elif y_axis_option == "back-azimuth":
        if geom == '2d':
            print('\t' + "Can't plot back azimuth deviation from 2d simulation")
            return 0
        else:
            print('\t' + "Back azimuth info...")
            ax2.set_ylabel("Back Azimuth [deg]")
            plot_vals = back_az_vals
    elif y_axis_option == "amplitude":
        print('\t' + "Amplitude info...")
        ax2.set_ylabel("Amplitude (rel. 1 km) [dB]")
        plot_vals = amp_vals
    else:
        print('\t' + "WARNING!  Bad plot option provided." + '\n\t' + "Plotting launch angle inclination info...")
        ax2.set_ylabel("Launch Inclination [deg]")
        plot_vals = incl_vals

    if cmap_option == "inclination":
        cmap_label = "Launch Inclination [deg]"
        cmap_vals = incl_vals
    elif cmap_option == "celerity":
        cmap_label = "Celerity [m/s]"
        cmap_vals = cel_vals
    elif cmap_option == "reduced-time":
        cmap_label = "Reduced Time (rel. " + str(int(reduced_tm_vel)) + " m/s) [s]"
        cmap_vals = tm_vals - arr_rngs / (reduced_tm_vel * 1.0e-3)
    elif cmap_option == "turning-ht":
        cmap_label = "Turning Height [km]"
        cmap_vals = turn_ht_vls
    elif cmap_option == "trace-velocity":
        cmap_label = "Trace Velocity [m/s]"
        cmap_vals = tr_vel_ref / np.cos(np.radians(arr_incl_vals))
    elif cmap_option == "back-azimuth":
        if geom == '2d':
            print('\t' + "Can't plot back azimuth info from 2d simulation")
            cmap_option = None
        else:
            cmap_label = "Back Azimuth [deg]"
            cmap_vals = back_az_vals
    elif cmap_option == "amplitude":
        cmap_label = "Amplitude (rel. 1 km) [dB]"
        cmap_vals = amp_vals
    else:
        cmap_vals = None

    if cmap_option is not None:
        sc = ax2.scatter(arr_rngs, plot_vals, c=cmap_vals, s=7.5, cmap = cm.jet)
        plt.colorbar(sc, label=cmap_label)
    else:
        ax2.plot(arr_rngs, plot_vals, 'ok', markersize=3.0)


    if figure_out is not None:
        print('\t' + "Saving figure to " + figure_out)
        plt.savefig(figure_out, dpi=250)

    print('')
    plt.show()


@click.command("eigenray", short_help="Visualize eigenray results and predicted arrival information")
@click.option("--atmo-file", help="Atmospheric specification file")
@click.option("--arrivals", help="Arrivals file from an 'eigenray' simulation (optional)", default=None)
@click.option("--eigenrays", help="Eigenrays file from an 'eigenray' simulation (optional)", default=None)
@click.option("--y-axis-option", help="Lower axis option (see usage info below)", default='inclination')
@click.option("--tr-vel-ref", help="Reference velocity for trace velocity calculation", default=330.0)
@click.option("--figure-out", help="Name of output figure", default=None)
def plot_eigenrays(atmo_file, arrivals, eigenrays, y_axis_option, tr_vel_ref, figure_out):
    '''
    Visualize results for eigenray analysis

    \b
    Plotting Options:
    \t inclination \t\t Launch inclination angle
    \t celerity  \t\t Arrival celerity (horizontal group velocity)
    \t turning-ht \t\t Turning height
    \t trace-velocity \t Trace velocity
    \t back-azimuth \t\t Back azimuth (not available for 2d geometry)
    \t amplitude \t\t Transport equation + absorption losses

    \b
    Examples:
    \t infraga plot eigenray --atmo-file ToyAtmo.met
    \t infraga plot eigenray --atmo-file ToyAtmo.met --y-axis-option trace-velocity  
    
    '''

    # Load data for plotting
    print("Loading data...")
    atmo_data = np.loadtxt(atmo_file)

    if arrivals is not None:
        print('\t' + "Loading specified arrivals file: " + arrivals)
        arrivals_file = arrivals
    else:
        print('\t' + "Loading arrivals file using atmo file name: " + os.path.splitext(atmo_file)[0] + ".arrivals.dat")
        arrivals_file = os.path.splitext(atmo_file)[0] + ".arrivals.dat"

    if not os.path.isfile(arrivals_file):
        print('\t' + "Arrivals file (" + arrivals_file + ") not found.")
        return 0

    for line in open(arrivals_file):
        if "infraga-" in line:
            if "sph" in line:
                geom = "sph"
            else:
                geom = "3d"
        elif "source" in line:
            src_loc = np.array([float(val) for val in line.split(":")[-1].split(",")[:2]])

    arr_data = np.atleast_2d(np.loadtxt(arrivals_file))
    if len(arr_data) < 1:
        print('\t' + "Arrivals file (" + arrivals_file + ") empty.")
        return 0
    prop_az = np.average(arr_data[:, 1])


    if eigenrays is not None:
        print('\t' + "Loading specified eigenray file(s): " + eigenrays)
        if "*" in eigenrays:
            eigenray_id = eigenrays
        else:
            eigenray_id = eigenrays + "*.dat"
    else:
        print('\t' + "Loading eigenray file(s) using atmo file name: " + os.path.splitext(atmo_file)[0] + ".eigenray-*.dat")
        eigenray_id = os.path.splitext(atmo_file)[0] + ".eigenray-*.dat"

    eigenray_data = []
    if len(os.path.dirname(eigenray_id)) > 0:
        eigenray_path = os.path.dirname(eigenray_id) + "/"
    else:
        eigenray_path = ""

    if "/" in eigenray_id:
        for file in np.sort(os.listdir(os.path.dirname(eigenray_path))):
            if fnmatch.fnmatch(file, os.path.basename(eigenray_id)):
                eigenray_data = eigenray_data + [np.loadtxt(eigenray_path + file)]
    else:        
        for file in np.sort(os.listdir(".")):
            if fnmatch.fnmatch(file, os.path.basename(eigenray_id)):
                eigenray_data = eigenray_data + [np.loadtxt(file)]

    if len(eigenray_data) < 1:
        print('\t' + "Eigenray file(s) (" + eigenray_id + ") not found.")
        return 0

    # Plot data
    print('\n' + "Plotting...")
    fig = plt.figure(figsize=(9, 5), layout='constrained')
    spec = fig.add_gridspec(2, 7)

    ax0 = fig.add_subplot(spec[0, 0])
    ax0.set_ylabel("Altitude [km]")
    ax0.set_xlabel("Sound Speed [m/s]")
    ax0.xaxis.set_label_position('top')
    ax0.xaxis.set_ticks_position('top')

    print('\t' + "Atmospheric data...")
    ax0.plot(np.sqrt(0.14 * atmo_data[:, 5] / atmo_data[:, 4]), atmo_data[:, 0], '--k', linewidth=1.5)
    ax0.plot(np.sqrt(0.14 * atmo_data[:, 5] / atmo_data[:, 4]) + atmo_data[:, 2] * np.sin(np.radians(prop_az)) + atmo_data[:, 3] * np.cos(np.radians(prop_az)), atmo_data[:, 0], '-k', linewidth=3.0)

    ax1 = fig.add_subplot(spec[0, 1:], sharey=ax0)
    ax1.set_xlabel("Range [km]")
    ax1.xaxis.set_label_position('top')
    ax1.xaxis.set_ticks_position('top')  

    print('\t' + "Eigenray path(s)...")
    color_seq = ['blue', 'orange', 'green', 'red', 'cyan', 'magenta', 'black']
    for n, path in enumerate(eigenray_data):
        if geom == '3d':
            ray_rngs = np.sqrt((src_loc[0] - path[:, 0])**2 + (src_loc[1] - path[:, 1])**2)
        else:
            ray_rngs = sph_proj.inv([src_loc[1]] * len(path), [src_loc[0]] * len(path), path[:, 1], path[:, 0])[2] * 1.0e-3
        ax1.plot(ray_rngs, path[:, 2], color=color_seq[n], linewidth=2.5)  

    ax1.set_xlim(left=0.0)    
    ax1.set_ylim(bottom=0.0) 

    incl_vals = arr_data[:, 0]
    tm_vals = arr_data[:, 5]
    cel_vals = arr_data[:, 6] * 1.0e3
    turn_ht_vls = arr_data[:, 7]
    arr_incl_vals = arr_data[:, 8]
    back_az_vals = arr_data[:, 9]
    amp_vals = arr_data[:, 10] + arr_data[:, 11]

    ax2 = fig.add_subplot(spec[1, 1:])
    ax2.set_xlabel("Time (rel. origin) [s]")

    if y_axis_option == 'inclination':
        print('\t' + "Launch inclination info...")
        ax2.set_ylabel("Launch Inclination [deg]")
        plot_vals = arr_incl_vals
    elif y_axis_option == 'celerity':
        print('\t' + "Arrival celerity info...")
        ax2.set_ylabel("Celerity [m/s]")
        plot_vals = cel_vals
    elif y_axis_option == "turning-ht" or y_axis_option == "turning-height":
        print('\t' + "Turning height info...")
        ax2.set_ylabel("Turning Height [km]")
        plot_vals = turn_ht_vls
    elif y_axis_option == "trace-velocity":
        print('\t' + "Trace velocity info...")
        ax2.set_ylabel("Trace Velocity [m/s]")
        plot_vals = tr_vel_ref / np.cos(np.radians(arr_incl_vals))
    elif y_axis_option == "back-azimuth":
        if geom == '2d':
            print('\t' + "Can't plot back azimuth deviation from 2d simulation")
            return 0
        else:
            print('\t' + "Back azimuth info...")
            ax2.set_ylabel("Back Azimuth [deg]")
            plot_vals = back_az_vals
    elif y_axis_option == "amplitude":
        print('\t' + "Amplitude info...")
        ax2.set_ylabel("Amplitude (rel. 1 km) [dB]")
        plot_vals = amp_vals
    else:
        print('\t' + "WARNING!  Bad plot option provided." + '\n\t' + "Plotting launch angle inclination info...")
        ax2.set_ylabel("Launch Inclination [deg]")
        plot_vals = incl_vals

    for n, val in enumerate(plot_vals):
        ax2.plot([tm_vals[n]], [val], color=color_seq[n], markersize=7.5, marker='o')  

    if figure_out is not None:
        print('\t' + "Saving figure to " + figure_out)
        plt.savefig(figure_out, dpi=250)

    plt.show()


@click.command("eig_wvfrms", short_help="Visualize eigenrays and predicted waveform")
@click.option("--atmo-file", help="Atmospheric specification file", default=None)
@click.option("--eigenrays", help="Eigenrays file from an 'eig_wvfrm' simulation", default=None)
@click.option("--wvfrms", help="Waveforms file from an eig_wvfrm simulation", default=None)
@click.option("--tr-vel-ref", help="Reference velocity for trace velocity calculation", default=330.0)
@click.option("--y-axis-option", help="Arrival parameter to plot on y-axis", default="back-azimuth")
@click.option("--cmap-option", help="Arrival parameter to plot on colormap", default="trace-velocity")
@click.option("--figure-out", help="Name of output figure", default=None)
def plot_eig_wvfrm(atmo_file, eigenrays, wvfrms, tr_vel_ref, y_axis_option, cmap_option, figure_out):
    '''
    Visualize results for combined eigenray/waveform analysis

    \b
    Plotting Options:
    \t trace-velocity \t Trace velocity
    \t back-azimuth \t\t Back azimuth
    \t amplitude \t\t Transport equation + absorption losses

    \b
    Examples:
    \t infraga plot eig_wvfrms --eigenrays ToyAtmo.eigenrays.dat --wvfrms ToyAtmo.wvfrms.dat
        

    '''
    # Check that needed files are defined
    if atmo_file is not None:
        eigenrays = os.path.splitext(atmo_file)[0] + ".eigenrays.dat"
        wvfrms = os.path.splitext(atmo_file)[0] + ".wvfrms.dat"
    elif eigenrays is None or wvfrms is None:
        print('\t' + "Visualization requires either atmspheric specification file or eigenray and waveform output")
        return 0
    
    if not os.path.isfile(eigenrays):
        print('\t' + "Eigenray results file (" + eigenrays + ") not found.")
        return 0
    
    if not os.path.isfile(wvfrms):
        print('\t' + "Waveform results file (" + wvfrms + ") not found.")
        return 0

    # Extract eigenray arrival info:
    print("Extracting eigenray information from results...")
    file_in = open(wvfrms, 'r')
    arr_tms = []
    arr_back_az = []
    arr_trace_vel = []
    arr_amp = []

    in_arrivals = False
    for line in file_in:
        if "infraga" in line:
            if "3d" in line:
                geom = '3d'
            else:
                geom = 'sph'
        elif "source location" in line:
            src_loc = np.array([float(val) for val in line.split(":")[-1].split(",")[:2]])

        if in_arrivals:
            if len(line) < 2:
                in_arrivals = False 
                break
            line = line[2:].split(' ')
            arr_tms = arr_tms + [float(line[5])]
            arr_back_az = arr_back_az + [float(line[9])]
            arr_trace_vel = arr_trace_vel + [tr_vel_ref / np.cos(np.radians(float(line[8])))]
            arr_amp = arr_amp + [float(line[10]) + float(line[11])]

        if "# incl" in line:
            in_arrivals = True 

    # Load data for plotting
    wvfrm_data = np.loadtxt(wvfrms)
    eigenray_data = np.loadtxt(eigenrays)

    if geom == '3d':
        ray_rngs = np.sqrt(eigenray_data[:, 0]**2 + eigenray_data[:, 1]**2)
    else:
        ray_rngs = sph_proj.inv([src_loc[1]] * len(eigenray_data), [src_loc[0]] * len(eigenray_data), eigenray_data[:, 1], eigenray_data[:, 0])[2] * 1.0e-3

    print("Plotting...")
    fig = plt.figure(figsize=(7, 7), layout='constrained')
    spec = fig.add_gridspec(7)

    print('\t' + "Eigenrays...")
    ax0 = fig.add_subplot(spec[:2])
    ax0.set_ylabel("Altitude [km]")
    ax0.set_xlabel("Range [km]")

    indices = np.flatnonzero(np.gradient(eigenray_data[:, 5]) < 0.0)
    ax0.plot(ray_rngs[:indices[0]], eigenray_data[:, 2][:indices[0]], '-k', linewidth=2.5)
    for n, j in enumerate(indices):
        ax0.plot(ray_rngs[indices[n - 1]:j], eigenray_data[:, 2][indices[n - 1]:j], '-k', linewidth=2.5)
    ax0.plot(ray_rngs[indices[-1]:], eigenray_data[:, 2][indices[-1]:], '-k', linewidth=2.5)

    ax0.set_xlim(left=0.0)    
    ax0.set_ylim(bottom=0.0) 


    print('\t' + "Arrival info...")
    ax2 = fig.add_subplot(spec[5:])
    ax2.set_xlabel("Time (rel. origin time) [s]")

    if y_axis_option == "amplitude":
        ax2.set_ylabel("Amplitude (rel. 1 km) [dB]")
        y_axis_vals = arr_amp
    elif y_axis_option == "trace-velocity":
        ax2.set_ylabel("Trace Velocity [m/s]")
        y_axis_vals = arr_trace_vel
    else:
        ax2.set_ylabel("Back Azimuth [deg]")
        y_axis_vals = arr_back_az

    if cmap_option == "back-azimuth":
        cmap_label = "Back Azimuth [deg]"
        cmap_vals = arr_back_az
    elif cmap_option == "trace-velocity":
        cmap_label = "Trace Velocity [m/s]"
        cmap_vals = arr_trace_vel
    elif cmap_option == "amplitude":
        cmap_label = "Amp. (rel. 1 km) [dB]"
        cmap_vals = arr_amp
    else:
        cmap_label = None
        cmap_vals = None

    if cmap_vals is not None: 
        sc = ax2.scatter(arr_tms, y_axis_vals, c=cmap_vals, cmap=cm.jet, s=40.0)
        plt.colorbar(sc, label=cmap_label)
    else:
        ax2.plot(arr_tms, y_axis_vals, 'ok', markersize=7.5)

    print('\t' + "Waveform predictions...")
    ax1 = fig.add_subplot(spec[3:5], sharex=ax2)
    ax1.set_ylabel("Pressure [Pa]")
    ax1.set_xlabel(" ")

    ax1.set_xlim(25.0 * np.floor(min(arr_tms) / 25.0), 25.0 * np.ceil(max(arr_tms) / 25.0))
    ax1.plot(wvfrm_data[:, 0], np.sum(wvfrm_data[:, 1:], axis=1), '-k', linewidth=1.5)
    

    if figure_out is not None:
        print('\t' + "Saving figure to " + figure_out)
        plt.savefig(figure_out, dpi=250)
    
    plt.show()


@click.command('map', short_help="Visualize results on a cartopy map")
@click.option("--arrivals", help="Arrivals file from an infraga-sph simulation", default=None)
@click.option("--ray-paths", help="Ray path file from an infraga-sph simulation", default=None)
@click.option("--plot-option", help="Parameter to visualize for arrivals ('amplitude', 'turning-height', 'celerity', or 'none')", default='amplitude')
@click.option("--figure-out", help="Name of output figure", default="arrivals.png")
@click.option("--rcvrs-file", help="File containing receiver locations (optional)", default=None)
@click.option("--title", help="Title for the figure", default=None)
@click.option("--start-time", help="Propagation time [hours] for plotting sub-set of data", default=None, type=float)
@click.option("--end-time", help="Propagation time [hours] for plotting sub-set of data", default=None, type=float)
@click.option("--include-absorption", help="Include Sutherland & Bass losses", default=True)
@click.option("--offline-maps-dir", help="Use directory for offline cartopy maps", default=None)
def plot_map(arrivals, ray_paths, plot_option, figure_out, rcvrs_file, title, start_time, end_time, include_absorption, offline_maps_dir):
    '''
    Visualize arrivals or ray paths computed using infraga spherical methods on a Cartopy map

    \b
    Examples:
    \t infraga plot map --arrivals ToyAtmo.arrivals.dat --plot-option amplitude --title 'Toy Atmo arrival amplitudes' --figure-name ToyAtmo.arrivals.png
    \t infraga plot map --arrivals ToyAtmo.arrivals.dat --plot-option celerity --title 'Toy Atmo arrival celerity' --figure-name ToyAtmo.celerities.png
    \t infraga plot map --ray-paths ToyAtmo.raypaths.dat --title 'Toy Atmo ray paths' --figure-name ToyAtmo.raypaths.png

    '''
    if arrivals is not None:
        # extract source info from the header
        arrivals = open(arrivals, 'r')
        for line in arrivals: 
            if "source location" in line: 
                src_loc = [float(val) for val in line[35:-1].split(", ")[:2]] 
                break

        # load data and extract lat/lon info for the map
        arrivals = np.loadtxt(arrivals)

        lats = arrivals[:, 3]
        lons = arrivals[:, 4]

        lat_min, lat_max = np.floor(min(min(lats), src_loc[0])), np.ceil(max(max(lats), src_loc[0]))
        lon_min, lon_max = np.floor(min(min(lons), src_loc[1])), np.ceil(max(max(lons), src_loc[1]))
    else:
        src_loc = None

        # load data and extract lat/lon info for the map
        ray_paths = np.loadtxt(ray_paths)

        lats = ray_paths[:, 0]
        lons = ray_paths[:, 1]

        lat_min, lat_max = np.floor(min(lats)), np.ceil(max(lats))
        lon_min, lon_max = np.floor(min(lons)), np.ceil(max(lons))

    if offline_maps_dir is not None:
        use_offline_maps(offline_maps_dir)

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=map_proj)
    if title is not None:
        plt.title(title)

    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)

    gl = ax.gridlines(crs=map_proj, draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False

    lat_tick, lon_tick = int((lat_max - lat_min) / 5), int((lon_max - lon_min) / 5)
    gl.xlocator = mticker.FixedLocator(np.arange(lon_min - np.ceil(lon_tick / 2), lon_max + lon_tick, lon_tick))
    gl.ylocator = mticker.FixedLocator(np.arange(lat_min - np.ceil(lat_tick / 2), lat_max + lat_tick, lat_tick))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    # Add features (coast lines, borders)
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    ax.add_feature(cartopy.feature.BORDERS, linewidth=0.5)
    if (lon_max - lon_min) < 20.0:
        ax.add_feature(cartopy.feature.STATES, linewidth=0.5)
        ax.add_feature(cartopy.feature.RIVERS, edgecolor='dodgerblue', alpha=0.3)
        ax.add_feature(cartopy.feature.LAKES, facecolor='dodgerblue', edgecolor='dodgerblue', alpha=0.3)

    # Plot data
    if arrivals is not None:
        time_mask = np.ones_like(arrivals[:, 0])
        if start_time is not None:
            time_mask = np.logical_and(time_mask, start_time < arrivals[:, 5] / 3600.0)
        if end_time is not None:
            time_mask = np.logical_and(time_mask, arrivals[:, 5] / 3600.0 < end_time)
        combo_mask = np.logical_and(time_mask, arrivals[:, 7] > 80.0)

        if plot_option == "turning-height":
            print('\t' + "Generating map with turning height info....")
            ax.scatter(arrivals[:,4][combo_mask], arrivals[:,3][combo_mask], c=arrivals[:,7][combo_mask], transform=map_proj, cmap=cm.jet_r, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=0.0, vmax=130.0)

            combo_mask = np.logical_and(time_mask, np.logical_and(arrivals[:,7] > 12.0, arrivals[:, 7] < 80.0))
            ax.scatter(arrivals[:,4][combo_mask], arrivals[:,3][combo_mask], c=arrivals[:,7][combo_mask], transform=map_proj, cmap=cm.jet_r, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=0.0, vmax=130.0)

            combo_mask = np.logical_and(time_mask, np.logical_and(arrivals[:,7] > 1.8, arrivals[:, 7] < 12.0))
            sc = ax.scatter(arrivals[:,4][combo_mask], arrivals[:,3][combo_mask], c=arrivals[:, 7][combo_mask], transform=map_proj, cmap=cm.jet_r, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=0.0, vmax=130.0)

            divider = make_axes_locatable(ax)
            ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
            fig.add_axes(ax_cb)
            cbar = plt.colorbar(sc, cax=ax_cb)
            cbar.set_label('Turning Height [km]')
        elif plot_option == "amplitude":
            print('\t' + "Generating map with amplitude info....")
            if include_absorption:
                tloss = arrivals[:, 10] + arrivals[:, 11]
            else:
                tloss = arrivals[:, 10]

            ax.scatter(arrivals[:,4][combo_mask], arrivals[:,3][combo_mask], c=tloss[combo_mask], transform=map_proj, cmap=cm.jet, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=-100.0, vmax=-10.0)

            combo_mask = np.logical_and(time_mask, np.logical_and(arrivals[:,7] > 12.0, arrivals[:, 7] < 80.0))
            ax.scatter(arrivals[:,4][combo_mask], arrivals[:,3][combo_mask], c=tloss[combo_mask], transform=map_proj, cmap=cm.jet, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=-100.0, vmax=-10.0)

            combo_mask = np.logical_and(time_mask, np.logical_and(arrivals[:,7] > 1.8, arrivals[:, 7] < 12.0))
            sc = ax.scatter(arrivals[:,4][combo_mask], arrivals[:,3][combo_mask], c=tloss[combo_mask], transform=map_proj, cmap=cm.jet, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=-100.0, vmax=-10.0)

            divider = make_axes_locatable(ax)
            ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
            fig.add_axes(ax_cb)
            cbar = plt.colorbar(sc, cax=ax_cb)
            cbar.set_label('Amplitude (power rel. 1 km) [dB]')
        elif plot_option == "celerity":
            print('\t' + "Generating map with celerity info....")
            ax.scatter(arrivals[:,4][combo_mask], arrivals[:,3][combo_mask], c=arrivals[:,6][combo_mask] * 1e3, transform=map_proj, cmap=cm.jet, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=220.0, vmax=340.0)

            combo_mask = np.logical_and(time_mask, np.logical_and(arrivals[:,7] > 12.0, arrivals[:, 7] < 80.0))
            ax.scatter(arrivals[:,4][combo_mask], arrivals[:,3][combo_mask], c=arrivals[:,6][combo_mask] * 1e3, transform=map_proj, cmap=cm.jet, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=220.0, vmax=340.0)

            combo_mask = np.logical_and(time_mask, np.logical_and(arrivals[:,7] > 1.8, arrivals[:, 7] < 12.0))
            sc = ax.scatter(arrivals[:,4][combo_mask], arrivals[:,3][combo_mask], c=arrivals[:,6][combo_mask] * 1e3, transform=map_proj, cmap=cm.jet, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=220.0, vmax=340.0)

            divider = make_axes_locatable(ax)
            ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
            fig.add_axes(ax_cb)
            cbar = plt.colorbar(sc, cax=ax_cb)
            cbar.set_label('Celerity [m/s]')

        else:
            ax.plot(arrivals[:,4][time_mask], arrivals[:,3][time_mask], "b.", transform=map_proj, markersize=marker_size / 2.0)

        ax.plot([src_loc[1]], [src_loc[0]], 'r*', markersize=5.0, transform=map_proj)
    else:
        time_mask = np.ones_like(ray_paths[:, 0], dtype=bool)
        if start_time is not None:
            time_mask = np.logical_and(time_mask, start_time < ray_paths[:, 5] / 3600.0)
        if end_time is not None:
            time_mask = np.logical_and(time_mask, ray_paths[:, 5] / 3600.0 < end_time)

        sc = ax.scatter(ray_paths[:, 1][time_mask], ray_paths[:, 0][time_mask], c=(ray_paths[:,5][time_mask] / 3600.0), transform=map_proj, cmap=cm.jet_r, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=start_time, vmax=end_time)

        divider = make_axes_locatable(ax)
        ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
        fig.add_axes(ax_cb)
        cbar = plt.colorbar(sc, cax=ax_cb)
        cbar.set_label('Propagation Time [hrs]')

    if rcvrs_file:
        try:
            rcvr_locs = np.loadtxt(rcvrs_file)
            ax.plot(rcvr_locs[:, 1], rcvr_locs[:, 0], 'k^', markersize=3.0, transform=map_proj)
        except:
            print('\t\t' + "Invalid receivers file.  Omitting from plot.")

    print('\t' + "Saving map to " + figure_out)
    plt.tight_layout()
    plt.savefig(figure_out, dpi=250)
    plt.show()

