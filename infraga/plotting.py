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

import sys
import numpy as np

import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import matplotlib.ticker as mticker

from mpl_toolkits.axes_grid1 import make_axes_locatable

import cartopy.crs as crs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


marker_size = 1.0
map_proj = crs.PlateCarree()
resol = '100m'  # use data at this scale (not working at the moment)



def plot_atmo(specification, spec_format='zTuvdp', alt_max=None):

    atmo = np.loadtxt(specification)
    z = atmo[:, spec_format.find('z')]
    u = atmo[:, spec_format.find('u')]
    v = atmo[:, spec_format.find('v')]
    c = np.sqrt(0.14 * atmo[:, spec_format.find('p')] / atmo[:, spec_format.find('d')])

    grnd_ht = z[0]
    for line in open(specification, 'r'):
        if "Ground Height" in line:
            grnd_ht = float(line[18:])
            break

    if alt_max is None:
        alt_max = z[-1]
    else:
        alt_max = float(alt_max)

    ht_mask = np.logical_and(grnd_ht <= z, z <= alt_max)

    f, ax = plt.subplots(1, 3, gridspec_kw={'width_ratios': [1, 1, 4]}, figsize=(12, 5))
    
    ax[0].grid(color='k', linestyle='--', linewidth=0.5)
    ax[1].grid(color='k', linestyle='--', linewidth=0.5)

    ax[0].set_ylim(grnd_ht, alt_max)
    ax[0].set_ylabel("Altitude [km]")
    ax[0].set_xlabel("Sound Speed [m/s]")

    ax[1].set_ylim(grnd_ht, alt_max)
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
        refract_ht = [z[ht_mask][np.min(np.where((ceff / ceff[0]) * np.cos(np.radians(incl)) > 1.0)[0])] if len(np.where((ceff / ceff[0]) * np.cos(np.radians(incl)) > 1.0)[0]) > 0 else alt_max for incl in incl_vals]
        sc = ax[2].scatter([az] * len(refract_ht), incl_vals, c=refract_ht, cmap=cm.jet_r, marker="s", s=5.0, alpha=0.75, edgecolor='none', vmin=grnd_ht, vmax=120.0)

    f.colorbar(sc, ax=[ax[2]], location='top', label="Estimated Refraction Altitude [km]")

    plt.show()


def plot_2d():
    # Plot range vs. altitude results along a single azimuth

    return 0

def plot_map(arrivals_file, ray_paths_file, plot_option, file_out, rcvrs_file=None, title_text=None, time1=None, time2=None, include_absorp=True):
    if arrivals_file is not None:
        # extract source info from the header
        arrivals_file = open(arrivals_file, 'r')
        for line in arrivals_file: 
            if "source location" in line: 
                src_loc = [float(val) for val in line[35:-1].split(", ")[:2]] 
                break

        # load data and extract lat/lon info for the map
        arrivals = np.loadtxt(arrivals_file)

        lats = arrivals[:, 3]
        lons = arrivals[:, 4]

        lat_min, lat_max = np.floor(min(min(lats), src_loc[0])), np.ceil(max(max(lats), src_loc[0]))
        lon_min, lon_max = np.floor(min(min(lons), src_loc[1])), np.ceil(max(max(lons), src_loc[1]))
    else:
        src_loc = None

        # load data and extract lat/lon info for the map
        ray_paths = np.loadtxt(ray_paths_file)

        lats = ray_paths[:, 0]
        lons = ray_paths[:, 1]

        lat_min, lat_max = np.floor(min(lats)), np.ceil(max(lats))
        lon_min, lon_max = np.floor(min(lons)), np.ceil(max(lons))

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=map_proj)
    plt.title(title_text)

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
    ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    if (lon_max - lon_min) < 20.0:
        ax.add_feature(cfeature.STATES, linewidth=0.5)
        ax.add_feature(cfeature.RIVERS, edgecolor='dodgerblue', alpha=0.3)
        ax.add_feature(cfeature.LAKES, facecolor='dodgerblue', edgecolor='dodgerblue', alpha=0.3)

    # Plot data
    if arrivals_file is not None:
        time_mask = np.ones_like(arrivals[:, 0])
        if time1 is not None:
            time_mask = np.logical_and(time_mask, time1 < arrivals[:, 5] / 3600.0)
        if time2 is not None:
            time_mask = np.logical_and(time_mask, arrivals[:, 5] / 3600.0 < time2)

        if plot_option == "turning-height":
            print('\t' + "Generating map with turning height info....")
            combo_mask = np.logical_and(time_mask, arrivals[:, 7] > 80.0)
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
            combo_mask = np.logical_and(time_mask, arrivals[:, 7] > 80.0)
            if include_absorp:
                tloss = arrivals[:, 10] + arrivals[:, 11]
            else:
                tloss = arrivals[:, 10]

            ax.scatter(arrivals[:,4][combo_mask], arrivals[:,3][combo_mask], c=tloss[combo_mask], transform=map_proj, cmap=cm.jet, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=-120.0, vmax=-20.0)

            combo_mask = np.logical_and(time_mask, np.logical_and(arrivals[:,7] > 12.0, arrivals[:, 7] < 80.0))
            ax.scatter(arrivals[:,4][combo_mask], arrivals[:,3][combo_mask], c=tloss[combo_mask], transform=map_proj, cmap=cm.jet, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=-120.0, vmax=-20.0)

            combo_mask = np.logical_and(time_mask, np.logical_and(arrivals[:,7] > 1.8, arrivals[:, 7] < 12.0))
            sc = ax.scatter(arrivals[:,4][combo_mask], arrivals[:,3][combo_mask], c=tloss[combo_mask], transform=map_proj, cmap=cm.jet, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=-120.0, vmax=-20.0)

            divider = make_axes_locatable(ax)
            ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
            fig.add_axes(ax_cb)
            cbar = plt.colorbar(sc, cax=ax_cb)
            cbar.set_label('Amplitude (power rel. 1 km) [dB]')
        elif plot_option == "celerity":
            print('\t' + "Generating map with celerity info....")
            combo_mask = np.logical_and(time_mask, arrivals[:, 7] > 80.0)
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
        time_mask = np.ones_like(ray_paths[:, 0])
        if time1 is not None:
            time_mask = np.logical_and(time_mask, time1 < ray_paths[:, 5] / 3600.0)
        if time2 is not None:
            time_mask = np.logical_and(time_mask, ray_paths[:, 5] / 3600.0 < time2)

        sc = ax.scatter(ray_paths[:, 1][time_mask], ray_paths[:, 0][time_mask], c=(ray_paths[:,5][time_mask] / 3600.0), transform=map_proj, cmap=cm.jet_r, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=time1, vmax=time2)

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

    print('\t' + "Saving map to " + file_out)
    plt.savefig(file_out, dpi=250)
    plt.show()

