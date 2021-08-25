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

# Figure parameters
title_text = "infraga-sph arrival predictions"

rcvr_locs = None
# rcvr_locs = [[30, -100], [40, -110]]

marker_size = 1.0
map_proj = crs.PlateCarree()
resol = '100m'  # use data at this scale (not working at the moment)

def print_usage():
    print('\n' + "Run infraGA/GeoAc eigenray analysis for a specified source and receiver")
    print('\n' + "Usage: python -OPTION plot_on_map.py output.arrivals.dat figure.png", '\n')

    print('\t' + "Options include:")
    print('\t' + "--------------------------------------------")
    print('\t' + "-amplitude" + '\t\t' + "Plot the arrival amplitude (combined transport equation and thermo-visous losses)")
    print('\t' + "-turning-height" + '\t\t' + "Plot the turning height of arrivals (visualize tropospheric, stratospheric, and thermospheric arrivals)")
    print('\t' + "-celerity" + '\t\t' + "Plot the celerity (horizontal group velocity) of arrivals" + '\n')

if __name__ == '__main__':
    print('\n\t' + "#" * 29)
    print('\t' + "#" * 2 + "      InfraGA/GeoAc      " + "#" * 2)
    print('\t' + "#" * 2 + "     Arrival Mapping     " + "#" * 2)
    print('\t' + "#" * 29)

    if len(sys.argv) < 2:
        print_usage()
    else:
        print('\n' + "Plotting arrivals on map using cartopy...")
        # extract source info from the header
        arrivals_file = open(sys.argv[2], 'r')
        for line in arrivals_file: 
            if "source location" in line: 
                src_loc = [float(val) for val in line[35:-1].split(", ")[:2]] 
                break
        print('\t' + "Extracted source location:", src_loc)

        # load data and extract lat/lon info for the map
        arrivals = np.loadtxt(sys.argv[2])

        arrival_lats = arrivals[:, 3]
        arrival_lons = arrivals[:, 4]

        lat_min, lat_max = np.floor(min(arrival_lats)), np.ceil(max(arrival_lats))
        lon_min, lon_max = np.floor(min(arrival_lons)), np.ceil(max(arrival_lons))
        print('\t' + "Identified lat/lon bounds: [" + str(lat_min) + ", " + str(lat_max) + "], [" + str(lon_min) + ", " + str(lon_max) + "].")

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=map_proj)
        plt.title(title_text)

        ax.set_xlim(lon_min, lon_max)
        ax.set_ylim(lat_min, lat_max)

        gl = ax.gridlines(crs=map_proj, draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False

        lat_tick, lon_tick = int((lat_max - lat_min) / 5), int((lon_max - lon_min) / 5)
        gl.xlocator = mticker.FixedLocator(np.arange(lon_min - np.ceil(lon_tick / 2), lon_max + lon_tick, lon_tick))
        gl.ylocator = mticker.FixedLocator(np.arange(lat_min - np.ceil(lat_tick / 2), lat_max + lat_tick, lat_tick))
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        # Add features (coast lines, borders)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
        ax.add_feature(cfeature.STATES, linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.5)

        '''
        ax.add_feature(cfeature.COASTLINE.with_scale(resol))
        ax.add_feature(cfeature.STATES.with_scale(resol))
        ax.add_feature(cfeature.BORDERS.with_scale(resol))

        ax.add_feature(cfeature.RIVERS.with_scale(resol), edgecolor="blue")
        ax.add_feature(cfeature.LAKES.with_scale(resol), edgecolor='0.25')
        '''

        # Plot data
        if sys.argv[1] == "-turning-height":
            print('\t' + "Plotting arrival turning heights...")
            ht_mask = arrivals[:,7] > 80.0
            ax.scatter(arrivals[:,4][ht_mask], arrivals[:,3][ht_mask], c=arrivals[:,7][ht_mask], transform=map_proj, cmap=cm.jet_r, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=0.0, vmax=130.0)

            ht_mask = np.logical_and(arrivals[:,7] > 12.0, arrivals[:,7] < 80.0)
            ax.scatter(arrivals[:,4][ht_mask], arrivals[:,3][ht_mask], c=arrivals[:,7][ht_mask], transform=map_proj, cmap=cm.jet_r, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=0.0, vmax=130.0)

            ht_mask = np.logical_and(arrivals[:,7] > 1.8, arrivals[:,7] < 12.0)
            sc = ax.scatter(arrivals[:,4][ht_mask], arrivals[:,3][ht_mask], c=arrivals[:,7][ht_mask], transform=map_proj, cmap=cm.jet_r, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=0.0, vmax=130.0)

            divider = make_axes_locatable(ax)
            ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
            fig.add_axes(ax_cb)
            cbar = plt.colorbar(sc, cax=ax_cb)
            cbar.set_label('Turning Height [km]')
        elif sys.argv[1] == "-amplitude":
            print('\t' + "Plotting arrival amplitudes...")
            ht_mask = arrivals[:,7] > 80.0
            ax.scatter(arrivals[:,4][ht_mask], arrivals[:,3][ht_mask], c=(arrivals[:, 10] + arrivals[:, 11])[ht_mask], transform=map_proj, cmap=cm.jet, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=-100.0, vmax=-20.0)

            ht_mask = np.logical_and(arrivals[:,7] > 12.0, arrivals[:,7] < 80.0)
            ax.scatter(arrivals[:,4][ht_mask], arrivals[:,3][ht_mask], c=(arrivals[:, 10] + arrivals[:, 11])[ht_mask], transform=map_proj, cmap=cm.jet, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=-100.0, vmax=-20.0)

            ht_mask = np.logical_and(arrivals[:,7] > 1.8, arrivals[:,7] < 12.0)
            sc = ax.scatter(arrivals[:,4][ht_mask], arrivals[:,3][ht_mask], c=(arrivals[:, 10] + arrivals[:, 11])[ht_mask], transform=map_proj, cmap=cm.jet, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=-100.0, vmax=-20.0)

            divider = make_axes_locatable(ax)
            ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
            fig.add_axes(ax_cb)
            cbar = plt.colorbar(sc, cax=ax_cb)
            cbar.set_label('Amplitude (power rel. 1 km) [dB]')
        elif sys.argv[1] == "-celerity":
            print('\t' + "Plotting arrival celerities...")
            ht_mask = arrivals[:,7] > 80.0
            ax.scatter(arrivals[:,4][ht_mask], arrivals[:,3][ht_mask], c=arrivals[:,6][ht_mask] * 1e3, transform=map_proj, cmap=cm.jet, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=220.0, vmax=340.0)

            ht_mask = np.logical_and(arrivals[:,7] > 12.0, arrivals[:,7] < 80.0)
            ax.scatter(arrivals[:,4][ht_mask], arrivals[:,3][ht_mask], c=arrivals[:,6][ht_mask] * 1e3, transform=map_proj, cmap=cm.jet, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=220.0, vmax=340.0)

            ht_mask = np.logical_and(arrivals[:,7] > 1.8, arrivals[:,7] < 12.0)
            sc = ax.scatter(arrivals[:,4][ht_mask], arrivals[:,3][ht_mask], c=arrivals[:,6][ht_mask] * 1e3, transform=map_proj, cmap=cm.jet, marker="o", s=marker_size, alpha=0.5, edgecolor='none', vmin=220.0, vmax=340.0)

            divider = make_axes_locatable(ax)
            ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
            fig.add_axes(ax_cb)
            cbar = plt.colorbar(sc, cax=ax_cb)
            cbar.set_label('Celerity [m/s]')
        else:
            print('\t' + "Invalid -option: " + sys.argv[1] + '\n\t' + "Try: -amplitude, -turning-height, or -celerity.")



        ax.plot([src_loc[1]], [src_loc[0]], 'r*', markersize=5.0, transform=map_proj)
        if np.any(rcvr_locs):
            ax.plot(rcvr_locs[:, 1], rcvr_locs[:, 0], 'k^', markersize=3.0, transform=map_proj)

        plt.savefig(sys.argv[3], dpi=250)

