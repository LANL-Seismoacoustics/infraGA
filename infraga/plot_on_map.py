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

import cartopy 
import cartopy.crs as crs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


marker_size = 1.0
map_proj = crs.PlateCarree()
resol = '100m'  # use data at this scale (not working at the moment)

def usage():
    print('\n' + "Run infraGA/GeoAc eigenray analysis for a specified source and receiver")
    print('\n' + "Usage: python -OPTION plot_on_map.py output.arrivals.dat figure.png", '\n')

    print('\t' + "Options include:")
    print('\t' + "--------------------------------------------")
    print('\t' + "amplitude" + '\t\t' + "Plot the arrival amplitude (combined transport equation and thermo-visous losses)")
    print('\t' + "turning-height" + '\t\t' + "Plot the turning height of arrivals (visualize tropospheric, stratospheric, and thermospheric arrivals)")
    print('\t' + "celerity" + '\t\t' + "Plot the celerity (horizontal group velocity) of arrivals" + '\n')

def run(arrivals_file, ray_paths_file, plot_option, file_out, rcvrs_file=None, title_text="infraga-sph predictions", time1=None, time2=None, include_absorp=True, cartopy_data_dir=None):
    if cartopy_data_dir is not None:
        cartopy.config['pre_existing_data_dir'] = cartopy_data_dir

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

    lat_tick, lon_tick = max(0.2, int((lat_max - lat_min) / 5.0)), max(0.2, int((lon_max - lon_min) / 5.0))

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

    '''
    # these adjustable resolution calls aren't working...
    ax.add_feature(cfeature.COASTLINE.with_scale(resol))
    ax.add_feature(cfeature.STATES.with_scale(resol))
    ax.add_feature(cfeature.BORDERS.with_scale(resol))

    '''

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


if __name__ == '__main__':
    print('\n\t' + "#" * 29)
    print('\t' + "#" * 2 + "      InfraGA/GeoAc      " + "#" * 2)
    print('\t' + "#" * 2 + "     Arrival Mapping     " + "#" * 2)
    print('\t' + "#" * 29)

    if len(sys.argv) < 2:
        usage()
    else:
        print("Plotting arrival information in '" + sys.argv[2] + "' with option: '" + sys.argv[1] + "'...")
        run(sys.argv[2], sys.argv[1], sys.argv[3])


