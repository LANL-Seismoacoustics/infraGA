#!/usr/bin/env python
"""
ecmwf_extractor.py

Methods to extract profile(s) from ECMWF netCDF format
files for use in infraGA/GeoAc analysis

Author: pblom@lanl.gov    
"""

import os 
import sys
import numpy as np
import matplotlib.pyplot as plt

from scipy.interpolate import interp1d, interp2d

from pyproj import Geod
from netCDF4 import Dataset

sph_proj = Geod(ellps='sphere')

etopo_file = "/path/to/ETOPO1_Ice_g_gmt4.grd"

##############################
##  US standard atmosphere  ##
##      polynomial fit      ##
##############################

# From: Lingevitch, J. F., Collins, M. D., & Siegmann, W. L. (1999). 
# Parabolic equations for gravity and acousto-gravity waves. 
# The Journal of the Acoustical Society of America, 105(6), 3049-3056.

gasR = 287.0
den0 = 0.001225
coeffs_A = np.array([-3.9082017e-2, -1.1526465e-3, 3.2891937e-5, -2.0494958e-7,
                        -4.7087295e-2, 1.2506387e-3, -1.5194498e-5, 6.581877e-8])
coeffs_B = np.array([-4.9244637e-3,  -1.2984142e-6, -1.5701595e-6, 1.5535974e-8,
                        -2.7221769e-2, 4.247473e-4, -3.9583181e-6, 1.7295795e-8])

def density(z):
    """
        Computes the atmospheric density according to 
            the US standard atmosphere model using a 
            polynomial fit

        Parameters
        ----------
        z : float
            Altitude above sea level [km]

        Returns:
        density : float
            Density of the atmosphere at altitude z [g/cm^3] 
    """

    poly_A, poly_B = 0.0, 1.0
    for n in range(4):
        poly_A += coeffs_A[n] * z**(n + 1)
        poly_B += coeffs_B[n] * z**(n + 1)

    return den0 * 10.0**(poly_A / poly_B)

def pressure(z, T):
    """
        Computes the atmospheric pressure according to 
            the US standard atmosphere model using a 
            polynomial fit assuming an ideal gas

        Parameters
        ----------
        z : float
            Altitude above sea level [km]

        Returns:
        pressure : float
            Pressure of the atmosphere at altitude z [mbar] 
    """
     
    return density(z) * gasR * T * 10.0


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

    return interp2d(region_lon, region_lat, region_elev / 1000.0, kind='linear')


def extract_single(ecmwf_file, lat, lon, output):
    """
        Extracts the nearest node vertical profile from
            an ECMWF netCDF4 format file at a specified
            latitude and longitude location and saves
            the profile in a specified file using the 
            zTuvdp format for infraGA/GeoAc analysis

        Parameters
        ----------
        ecmwf_file : string
            Path and name of the ecmwf netCDF4 file
        lat : float
            Latitude where the vertical profile will be extracted [deg]
        lon : float
            Longitude where the vertical profile will be extracted [deg]
        output : string
            Path and name of the zTuvdp format output file
    """

    ecmwf = Dataset(ecmwf_file)

    lat0, dlat = float(ecmwf.variables['ylat0'][:].data), float(ecmwf.variables['dy'][:].data)
    lon0, dlon = float(ecmwf.variables['xlon0'][:].data), float(ecmwf.variables['dx'][:].data)
    
    lat_vals = np.arange(lat0, lat0 + dlat * ecmwf.variables['T'].shape[1], dlat)
    lon_vals = np.arange(lon0, lon0 + dlon * ecmwf.variables['T'].shape[2], dlon)

    n_lat = np.argmin(abs(lat_vals - lat))
    n_lon = np.argmin(abs(lon_vals - lon))

    grnd_lvl_interp = interp_etopo((lat - 1.5, lon - 1.5), (lat + 1.5, lon + 1.5))
    z_gl = grnd_lvl_interp(lon, lat)[0]

    print("Extracting profile at " + str(lat) + ", " + str(lon) + " into " + output + " with ground elevation " + str(z_gl))
    
    z_vals = ecmwf.variables['height'][:].data / 1000.0 + z_gl

    atmo = z_vals
    atmo = np.vstack((atmo, ecmwf.variables['T'][:, n_lat, n_lon].data))   
    atmo = np.vstack((atmo, ecmwf.variables['U'][:, n_lat, n_lon].data))   
    atmo = np.vstack((atmo, ecmwf.variables['V'][:, n_lat, n_lon].data))   
    atmo = np.vstack((atmo, density(z_vals)))   
    atmo = np.vstack((atmo, pressure(z_vals, ecmwf.variables['T'][:, n_lat, n_lon].data)))   
    np.savetxt(output, atmo.T)

def extract_grid(ecmwf_file, lat_llc, lon_llc, lat_urc, lon_urc, output_id, grid_skip=1):
    """
        Extracts a set of vertical profiles from an 
            ECMWF netCDF4 format file onto a specified
            grid of latitude and longitude locations 
            using the zTuvdp format for infraGA/GeoAc 
            analysis.  Grid skip allows decreasing the
            resolution of the grid sampling from the
            source file.

        Parameters
        ----------
        ecmwf_file : string
            Path and name of the ecmwf netCDF4 file
        lat_llc : float
            Latitude of the lower-left corner of the grid [deg]
        lon_llc : float
            Longitude of the lower-left cornder of the grid [deg]
        lat_urc : float
            Latitude of the uppper-right corner of the grid [deg]
        lon_urc : float
            Longitude of the upper-right cornder of the grid [deg]
        output_d : string
            Path and prefix of the grid loccation files and zTuvdp format output files
        grid_skip : int
            Frequency of skips in sampling the grid (sampled numpy array as [::grid_skip])
    """

    ecmwf = Dataset(ecmwf_file)

    lat0, dlat = float(ecmwf.variables['ylat0'][:].data), float(ecmwf.variables['dy'][:].data)
    lon0, dlon = float(ecmwf.variables['xlon0'][:].data), float(ecmwf.variables['dx'][:].data)

    lat_vals = np.arange(lat0, lat0 + dlat * ecmwf.variables['T'].shape[1], dlat)
    lon_vals = np.arange(lon0, lon0 + dlon * ecmwf.variables['T'].shape[2], dlon)

    grid_lats = lat_vals[np.logical_and(lat_llc <= lat_vals, lat_vals <= lat_urc)][::grid_skip]
    grid_lons = lon_vals[np.logical_and(lon_llc <= lon_vals, lon_vals <= lon_urc)][::grid_skip]

    grid_n_lats = np.array([np.where(lat_vals==val)[0] for val in grid_lats]).flatten() 
    grid_n_lons = np.array([np.where(lon_vals==val)[0] for val in grid_lons]).flatten() 

    np.savetxt(output_id + "_lats.dat", grid_lats)
    np.savetxt(output_id + "_lons.dat", grid_lons)

    grnd_lvl_interp = interp_etopo((lat_llc, lon_llc), (lat_urc, lon_urc))

    z_min = 100.0
    for n1, n_lat in enumerate(grid_n_lats):
        for n2, n_lon in enumerate(grid_n_lons):            
            z_min = min(z_min, grnd_lvl_interp(lon_vals[n_lon], lat_vals[n_lat])[0])

    z_agl = ecmwf.variables['height'][:].data / 1000.0
    z_vals = np.linspace(z_min, z_agl[-1], len(z_agl))
    for n1, n_lat in enumerate(grid_n_lats):
        for n2, n_lon in enumerate(grid_n_lons):
            z_gl = grnd_lvl_interp(lon_vals[n_lon], lat_vals[n_lat])[0]
            print("Extracting profile at " + str(lat_vals[n_lat]) + ", " + str(lon_vals[n_lon]) + " into " + output_id + str(n1 * len(grid_n_lons) + n2) + ".met")

            T_interp = interp1d(z_gl + z_agl, ecmwf.variables['T'][:, n_lat, n_lon].data, bounds_error=False, fill_value="extrapolate")
            U_interp = interp1d(z_gl + z_agl, ecmwf.variables['U'][:, n_lat, n_lon].data, bounds_error=False, fill_value="extrapolate")
            V_interp = interp1d(z_gl + z_agl, ecmwf.variables['V'][:, n_lat, n_lon].data, bounds_error=False, fill_value="extrapolate")

            atmo = z_vals
            atmo = np.vstack((atmo, T_interp(z_vals)))
            atmo = np.vstack((atmo, U_interp(z_vals)))
            atmo = np.vstack((atmo, V_interp(z_vals)))
            atmo = np.vstack((atmo, density(z_vals)))   
            atmo = np.vstack((atmo, pressure(z_vals, T_interp(z_vals))))   

            np.savetxt(output_id + str(n1 * len(grid_n_lons) + n2) + ".met", atmo.T)

def print_usage():
    print('\n\t' + "#" * 25)
    print('\t' + "#" * 2 + "        ECMWF        " + "#" * 2)
    print('\t' + "#" * 2 + "  Profile Extractor  " + "#" * 2)
    print('\t' + "#" * 25)

    print('\n' + "Extract a single vertical profile or a grid of profiles for a specified")
    print('\t' + "atmospheric dataset from the European Centre for Medium-Range Weather")
    print('\t' + "Forecasts (ECMWF) for use in the infraGA/GeoAc ray tracing software.")

    print('\n' + "Usage: python ecmwf_extractor.py [option] ecmwf_file [parameter values]")

    print('\n' + "Options and parameters:")
    print('\t' + "-single (extract a single vertical profile)")
    print('\t\t' + "Parameter" + '\t' + "Units")
    print('\t\t' + "-" * 50)
    print('\t\t' + "latitude" + '\t' + "degrees")
    print('\t\t' + "longitude" + '\t' + "degrees")
    print('\t\t' + "output.met" + '\t' + "-")

    print('\n\t' + "-grid (extract a grid of vertical profiles)")
    print('\t\t' + "Parameter" + '\t\t' + "Units")
    print('\t\t' + "-" * 50)
    print('\t\t' + "lower-left corner latitude" + '\t' + "degrees")
    print('\t\t' + "lower-left corner longitude" + '\t' + "degrees")
    print('\t\t' + "upper-right corner latitude" + '\t' + "degrees")
    print('\t\t' + "upper-right corner longitude" + '\t' + "degrees")
    print('\t\t' + "output id" + '\t\t\t' + "-")
    print('\t\t' + "grid sampling" + '\t\t\t' + "positive integer (default 1, optional)")

    print('\n' + "Examples:")
    print('\t' + "python ecmwf_extractor.py -single ecmwf_file 30.2 -120.1 test.met")
    print('\t' + "python ecmwf_extractor.py -grid ecmwf_file 25.0 -130.0 35.0 -120.0 grid/test 2" + '\n')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print_usage()
    else:
        if os.path.isfile(etopo_file):
            if sys.argv[1] == "-single":
                extract_single(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), sys.argv[5])
            elif sys.argv[1] == "-grid":
                if len(sys.argv) < 9:
                    extract_grid(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), sys.argv[7])
                else:
                    extract_grid(sys.argv[2], float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6]), sys.argv[7], int(sys.argv[8]))
            else:
                print_usage()
        else:
            print("Topography file not found.  Download at https://www.ngdc.noaa.gov/mgg/global/ and")
            print('\t' + "update location in script header if not placed here.")
        
            





    







