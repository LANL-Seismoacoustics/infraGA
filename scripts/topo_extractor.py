#!which python
"""
topo_extractor.py

Methods to extract topography information
from ETOPO1 netCDF file into great circle
path (line) or a grid (Cartesian or 
latitude/longitude)

Author: pblom@lanl.gov    
"""

import os
import sys 
import numpy as np
import matplotlib.pyplot as plt

from pyproj import Geod
from scipy import interpolate
from netCDF4 import Dataset

sph_proj = Geod(ellps='sphere')

etopo_file = "ETOPO1_Ice_g_gmt4.grd"

def interp_etopo(ll_corner, ur_corner):
    """
        Loads and interpolates the ETOPO1 topography
            data within a region specified by low-left
            and upper-right corner latitude, longitudes

        Parameters
        ----------
        ll_corner : iterable
            Iterable containing the latitude and longitude
                of the lower-left corner of the region
        ur_corner : iterable
            Iterable containing the latitude and longitude
                of the upper-right corner of the region

        Returns:
        elev_interp : scipy.interpolate.interp2d
            A 2d interpolation of the elevation within 
                the specified region
    """

    # load etopo_file and extract grid information
    etopo1 = Dataset(etopo_file)
    grid_lons = etopo1.variables['x'][:]
    grid_lats = etopo1.variables['y'][:]
    grid_elev = etopo1.variables['z'][:]

    lat_mask = np.logical_and(ll_corner[0] - 1.0 <= grid_lats, grid_lats <= ur_corner[0] + 1.0).nonzero()[0]
    lon_mask = np.logical_and(ll_corner[1] - 1.0 <= grid_lons, grid_lons <= ur_corner[1] + 1.0).nonzero()[0]

    region_lat = grid_lats[lat_mask]
    region_lon = grid_lons[lon_mask]
    region_elev = grid_elev[lat_mask,:][:,lon_mask]

    # Change underwater values to sea surface
    region_elev[region_elev < 0.0] = 0.0

    return interpolate.interp2d(region_lon, region_lat, region_elev / 1000.0, kind='linear')

    
def pull_pnt2pnt(src_loc, rcvr_loc, file_out, resol=1.852):
    """
        Extract topography information along a line defined
            by source and receiver locations (latitude, 
            longitude) into file_out

        Parameters
        ----------
        src_loc : iterable
            Iterable containing the latitude and longitude
                of the source (one end of the line)
        rcvr_loc : iterable
            Iterable containing the latitude and longitude
                of the source (one end of the line)
        file_out : str
            Destination for the topography information
        resol : float
            Spatial resolution of the ETOPO1 file.  1 arc
                minute = 1.852 km, this can modified to
                increase the resolution

    """

    print('Defining planar topography from source at ' + str(src_loc[0]) + ', ' + str(src_loc[1]) + ' to receiver at ' + str(rcvr_loc[0]) + ', ' + str(rcvr_loc[1]) + '.')

    # pbuild interpolation from determine grid corners
    # note: src_loc is (lat, lon), Geod returns (lon, lat) so indices reverse
    ll_corner = [min(src_loc[0], rcvr_loc[0]) - 0.5, min(src_loc[1], rcvr_loc[1]) - 0.5]
    ur_corner = [max(src_loc[0], rcvr_loc[0]) + 0.5, max(src_loc[1], rcvr_loc[1]) + 0.5]
    elev_interp = interp_etopo(ll_corner, ur_corner)

    # Etopo1 has resolution of 1 arc minute = 1.852 kms, so write that resolution to file
    N = int((sph_proj.inv(src_loc[1], src_loc[0], rcvr_loc[1], rcvr_loc[0], radians=False)[2] / 1000.0) / resol)
    line_pnts = sph_proj.npts(src_loc[1], src_loc[0], rcvr_loc[1], rcvr_loc[0], N, radians=False)

    output = open(file_out, 'w')
    # print("# Rng [km]" + '\t' + "Elev [km]", file=output)
    # note: infraGA methods to ingest topo file can't recognize header notation yet, so no header in these files
    print(0.0, '\t', elev_interp(src_loc[1], src_loc[0])[0], file=output)
    for n in range(N):
        print(sph_proj.inv(src_loc[1], src_loc[0], line_pnts[n][0], line_pnts[n][1], radians=False)[2] / 1000.0, file=output, end='\t')
        print(elev_interp(line_pnts[n][0], line_pnts[n][1])[0], file=output)
    output.close()

    # interpolate at 2x resolution to plot
    N_highres = 2 * N
    line_pnts = sph_proj.npts(src_loc[1], src_loc[0], rcvr_loc[1], rcvr_loc[0], N_highres, radians=False)
    rng_vals, elev_vals = np.empty(N_highres), np.empty(N_highres)

    for n in range(N_highres):
        rng_vals[n] = sph_proj.inv(src_loc[1], src_loc[0], line_pnts[n][0], line_pnts[n][1], radians=False)[2] / 1000.0
        elev_vals[n] = elev_interp(line_pnts[n][0], line_pnts[n][1])[0]

    plt.fill_between(rng_vals, elev_vals, y2=0.0, color='0.25')
    plt.xlabel("Range [km]")
    plt.ylabel("Elevation [km]")
    plt.show()


def pull_line(src_loc, azimuth, rng_max, file_out, resol=1.852):
    """
        Extract topography information along a line defined
            by a great circle path from a source location 
            (latitude, longitude) along a specified azimuth
            to a maximum specified range into file_out.

        Parameters
        ----------
        src_loc : iterable
            Iterable containing the latitude and longitude
                of the source (one end of the line)
        azimuth : float
            Azimuthal direction (relative to north) along
                which to extract the topography
        rng_max : float
            Distance from the source to draw the line
                and extract topography
        file_out : str
            Destination for the topography information
        resol : float
            Spatial resolution of the ETOPO1 file.  1 arc
                minute = 1.852 km, this can modified to
                increase the resolution

    """
    print('Defining planar topography from source at ' + str(src_loc[0]) + ', ' + str(src_loc[1]) + ' along azimuth ' + str(azimuth) + ' out to a range of ' + str(rng_max))

    # project lat/lon at end of line and determine grid corners
    # note: src_loc is (lat, lon), Geod returns (lon, lat) so indices reverse
    end_loc = sph_proj.fwd(src_loc[1], src_loc[0], azimuth, rng_max * 1.0e3, radians=False)
    ll_corner = [min(src_loc[0], end_loc[1]) - 0.5, min(src_loc[1], end_loc[0]) - 0.5]
    ur_corner = [max(src_loc[0], end_loc[1]) + 0.5, max(src_loc[1], end_loc[0]) + 0.5]

    elev_interp = interp_etopo(ll_corner, ur_corner)

    # Etopo1 has resolution of 1 arc minute = 1.852 kms, so write that resolution to file
    N = int(rng_max / resol)
    line_pnts = sph_proj.npts(src_loc[1], src_loc[0], end_loc[0], end_loc[1], N, radians=False)
    rng_vals, elev_vals = np.empty(N), np.empty(N)

    output = open(file_out, 'w')
    # print("# Rng [km]" + '\t' + "Elev [km]", file=output)
    # note: infraGA methods to ingest topo file can't recognize header notation yet, so no header in these files
    print('{}\t{}\t{}'.format(src_loc[0], src_loc[1], elev_interp(src_loc[1], src_loc[0])[0]), file=output)
    for n in range(N):
        print(sph_proj.inv(src_loc[1], src_loc[0], line_pnts[n][0], line_pnts[n][1], radians=False)[2] / 1000.0, file=output, end='\t')
        print(elev_interp(line_pnts[n][0], line_pnts[n][1])[0], file=output)
    output.close()

    # interpolate at 4x resolution to plot
    N = int(rng_max / resol * 4)
    line_pnts = sph_proj.npts(src_loc[1], src_loc[0], end_loc[0], end_loc[1], N, radians=False)
    rng_vals, elev_vals = np.empty(N), np.empty(N)

    for n in range(N):
        rng_vals[n] = sph_proj.inv(src_loc[1], src_loc[0], line_pnts[n][0], line_pnts[n][1], radians=False)[2] / 1000.0
        elev_vals[n] = elev_interp(line_pnts[n][0], line_pnts[n][1])[0]
  
    plt.fill_between(rng_vals, elev_vals, y2=0.0, color='0.25')
    plt.xlabel("Range [km]")
    plt.ylabel("Elevation [km]")
    plt.show()
    
    

def pull_xy_grid(src_loc, ll_corner, ur_corner, file_out, resol=1.852):
    """
        Extract topography information across a region defined
            by the lower-left and upper-right corner latitudes
            and longitudes.  Cartesian distances are computed
            relative to a specified source location.

        Parameters
        ----------
        src_loc : iterable
            Iterable containing the latitude and longitude
                of the source (one end of the line)
        ll_corner : iterable
            Iterable containing the latitude and longitude
                of the lower-left corner of the region
        ur_corner : iterable
            Iterable containing the latitude and longitude
                of the upper-right corner of the region
        file_out : str
            Destination for the topography information
        resol : float
            Spatial resolution of the ETOPO1 file.  1 arc
                minute = 1.852 km, this can modified to
                increase the resolution

    """

    print("Extracting Cartesian grid from " + str(ll_corner[0]) + ", " + str(ll_corner[1]), end='')
    print(" to " + str(ur_corner[0]) + ", " + str(ur_corner[1]) + " with source at " + str(src_loc[0]) + ", " + str(src_loc[1]))

    # load etopo_file and interpolate within region
    elev_interp = interp_etopo(ll_corner, ur_corner)

    # define dx, dy values of corners
    az, _, dr = sph_proj.inv(src_loc[1], src_loc[0], ll_corner[1], ll_corner[0])
    ll_x, ll_y = np.round(dr * np.sin(np.radians(az)) / 1.0e3), np.round(dr * np.cos(np.radians(az)) / 1.0e3)

    az, _, dr = sph_proj.inv(src_loc[1], src_loc[0], ur_corner[1], ur_corner[0])
    ur_x, ur_y = np.round(dr * np.sin(np.radians(az)) / 1.0e3), np.round(dr * np.cos(np.radians(az)) / 1.0e3)

    x_vals = np.arange(ll_x, ur_x, resol)
    y_vals = np.arange(ll_y, ur_y, resol)
    xy_elev = np.empty((len(x_vals), len(y_vals)))
    
    output = open(file_out, 'w')
    # print("# Rng (E/W) [km]" + '\t' + "Rng (N/S) [km]"  + '\t' + "Elev [km]", file=output)
    # note: infraGA methods to ingest topo file can't recognize header notation yet, so no header in these files
    for n in range(len(x_vals)):
        for m in range(len(y_vals)):
            pnt = sph_proj.fwd(src_loc[1], src_loc[0], np.degrees(np.arctan2(x_vals[n], y_vals[m])), np.sqrt(x_vals[n]**2 + y_vals[m]**2) * 1.0e3, radians=False)
            xy_elev[n][m] = elev_interp(pnt[0], pnt[1])
            print(x_vals[n], y_vals[m], xy_elev[n][m], file=output)
    output.close()

    XX, YY = np.meshgrid(x_vals, y_vals)
    plt.pcolormesh(XX, YY, xy_elev.T, cmap=plt.cm.terrain, vmin=0.0)
    plt.xlabel("Range (E/W) [km]")
    plt.ylabel("Range (N/S) [km]")
    plt.colorbar(label="Elevation [km]")

    plt.show()


def pull_latlon_grid(ll_corner, ur_corner, file_out):
    """
        Extract topography information across a region defined
            by the lower-left and upper-right corner latitudes
            and longitudes.

        Parameters
        ----------
        ll_corner : iterable
            Iterable containing the latitude and longitude
                of the lower-left corner of the region
        ur_corner : iterable
            Iterable containing the latitude and longitude
                of the upper-right corner of the region
        file_out : str
            Destination for the topography information
    """

    print("Extracting lat/lon grid from " + str(ll_corner[0]) + ", " + str(ll_corner[1]), end='')
    print(" to " + str(ur_corner[0]) + ", " + str(ur_corner[1]))

    # load etopo_file and extract grid information
    etopo1 = Dataset(etopo_file)
    grid_lons = etopo1.variables['x'][:]
    grid_lats = etopo1.variables['y'][:]
    grid_elev = etopo1.variables['z'][:]

    lat_mask = np.logical_and(ll_corner[0] <= grid_lats, grid_lats <= ur_corner[0]).nonzero()[0]
    lon_mask = np.logical_and(ll_corner[1] <= grid_lons, grid_lons <= ur_corner[1]).nonzero()[0]

    region_lat = grid_lats[lat_mask]
    region_lon = grid_lons[lon_mask]
    region_elev = grid_elev[lat_mask,:][:,lon_mask]

    output = open(file_out, 'w')
    # print("# Latitude [deg]" + '\t' + "Longitude [deg]"  + '\t' + "Elev [km]", file=output)
    # note: infraGA methods to ingest topo file can't recognize header notation yet, so no header in these files
    for n in range(len(region_lat)):
        for m in range(len(region_lon)):
            print(region_lat[n], region_lon[m], max(region_elev[n][m], 0.0) / 1.0e3, file=output)
    output.close()

    LON, LAT = np.meshgrid(region_lon, region_lat)
    plt.pcolormesh(LON, LAT, region_elev / 1.0e3, cmap=plt.cm.terrain, vmin=0.0)
    plt.xlabel("Longitude [deg]")
    plt.ylabel("Latitude [deg]")
    plt.colorbar(label="Elevation [km]")
    plt.show()


def print_usage():
    print('\n\t' + "#" * 26)
    print('\t' + "#" * 2 + "        ETOPO1        " + "#" * 2)
    print('\t' + "#" * 2 + " Topography Extractor " + "#" * 2)
    print('\t' + "#" * 26)

    print('\n' + "Extract topography files for use in the InfraGA/GeoAc infrasound propagation")
    print('\t' + "software.  All methods require an ETOPO1 source file available at:")
    print('\t' + "https://www.ngdc.noaa.gov/mgg/global/ downloaded and placed in")
    print('\t' + "infraga/scripts/ or in a location specified in this script's header.")

    print('\n' + "Usage: python topo_extractor.py [option] [parameter values]")

    print('\n' + "Options and parameters:")
    print('\t' + "-line (extract a great circle path from source in specified direction)")
    print('\t\t' + "Parameter" + '\t\t\t' + "Units")
    print('\t\t' + "-" * 50)
    print('\t\t' + "source latitude" + '\t\t\t' + "degrees")
    print('\t\t' + "source longitude" + '\t\t' + "degrees")
    print('\t\t' + "azimuth" + '\t\t\t\t' + "degrees")
    print('\t\t' + "range" + '\t\t\t\t' + "km")
    print('\t\t' + "output file" + '\t\t\t' + "-")

    print('\n\t' + "-pnt2pnt (extract a great circle path between source and receiver locations)")
    print('\t\t' + "Parameter" + '\t\t\t' + "Units")
    print('\t\t' + "-" * 50)
    print('\t\t' + "source latitude" + '\t\t\t' + "degrees")
    print('\t\t' + "source longitude" + '\t\t' + "degrees")
    print('\t\t' + "receiver latitude" + '\t\t\t' + "degrees")
    print('\t\t' + "receiver longitude" + '\t\t' + "degrees")
    print('\t\t' + "output file" + '\t\t\t' + "-")
    
    print('\n\t' + "-xy_grid (extract a Cartesian grid over specified bounds relative to a source location)")
    print('\t\t' + "Parameter" + '\t\t\t' + "Units")
    print('\t\t' + "-" * 50)
    print('\t\t' + "source latitude" + '\t\t\t' + "degrees")
    print('\t\t' + "source longitude" + '\t\t' + "degrees")
    print('\t\t' + "lower-left corner latitude" + '\t' + "degrees")
    print('\t\t' + "lower-left corner longitude" + '\t' + "degrees")
    print('\t\t' + "upper-right corner latitude" + '\t' + "degrees")
    print('\t\t' + "upper-right corner longitude" + '\t' + "degrees")
    print('\t\t' + "output file" + '\t\t\t' + "-")

    print('\n\t' + "-latlon_grid (extract a latitude/longitude grid with specified corners)")
    print('\t\t' + "Parameter" + '\t\t\t' + "Units")
    print('\t\t' + "-" * 50)
    print('\t\t' + "lower-left corner latitude" + '\t' + "degrees")
    print('\t\t' + "lower-left corner longitude" + '\t' + "degrees")
    print('\t\t' + "upper-right corner latitude" + '\t' + "degrees")
    print('\t\t' + "upper-right corner longitude" + '\t' + "degrees")
    print('\t\t' + "output file" + '\t\t\t' + "-")

    print('\n' + "Examples:")
    print('\t' + "python topo_extractor.py -line 40.0 -102.5 -90.0 750.0 line_topo.dat")
    print('\t' + "python topo_extractor.py -pnt2pnt 40.0 -102.5 40.0 -110.0 line_topo.dat")
    print('\t' + "python topo_extractor.py -xy_grid 40.0 -102.5 35.0 -110.0 45.0 -100.0 xy_topo.dat")
    print('\t' + "python topo_extractor.py -latlon_grid 35.0 -110.0 45.0 -100.0 sph_topo.dat" + '\n')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print_usage()
    else:
        if os.path.isfile(etopo_file):
            if sys.argv[1] == "-line":
                pull_line((float(sys.argv[2]), float(sys.argv[3])), float(sys.argv[4]), float(sys.argv[5]), sys.argv[6])
            elif sys.argv[1] == "-pnt2pnt":
                pull_pnt2pnt((float(sys.argv[2]), float(sys.argv[3])), (float(sys.argv[4]), float(sys.argv[5])), sys.argv[6])
            elif sys.argv[1] == "-xy_grid":
                pull_xy_grid((float(sys.argv[2]), float(sys.argv[3])), (float(sys.argv[4]), float(sys.argv[5])), (float(sys.argv[6]), float(sys.argv[7])), sys.argv[8])
            elif sys.argv[1] == "-latlon_grid":
                pull_latlon_grid((float(sys.argv[2]), float(sys.argv[3])), (float(sys.argv[4]), float(sys.argv[5])), sys.argv[6])    
            else:
                print_usage()  
        else:
            print("Topography file not found.  Download at https://www.ngdc.noaa.gov/mgg/global/ and")
            print('\t' + "update location in script header if not placed here.")


