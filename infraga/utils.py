#!which python
"""
topo_extractor.py

Methods to extract topography information
from ETOPO1 netCDF file into great circle
path (line) or a grid (Cartesian or 
latitude/longitude)

Author: pblom@lanl.gov    
"""

import click
import os
import sys 
import wget
import fnmatch
import subprocess

from importlib.util import find_spec, spec_from_file_location 

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1 import make_axes_locatable

import cartopy
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from pyproj import Geod
from scipy import interpolate
from netCDF4 import Dataset

sph_proj = Geod(ellps='sphere')

def use_offline_maps(self, pre_existing_data_dir, turn_on=True):
    # call this function to initialize the use of offline maps.  turn_on will initialize the pre_existing_data_directory
    if turn_on:
        cartopy.config['pre_existing_data_dir'] = pre_existing_data_dir
    else:
        cartopy.config['pre_existing_data_dir'] = ""


################################
##    Range Dependent Grid    ##
################################
@click.command('build-g2s-grid', short_help="Build grid for range dependent analysis from G2S files")
@click.option("--g2s-path", help="Path to G2S specifications", prompt="Specify directory containing G2S specifications: ")
@click.option("--output-path", help="Output dir + label", prompt="Specify output directory and label: ")
@click.option("--src-info", help="Source info (lat, lon, time) (optional)", default=None)
@click.option("--celerity-est", help="Celerity estimate [km/s] (optional)", default=0.29)
def build_g2s_grid(g2s_path, output_path, src_info=None, celerity_est=0.29):
    '''
    Construct the numbered specifications and grid files needed to run -rngdep methods.
    Assumes file format from https://g2s.ncpa.olemiss.edu/ (e.g., g2stxt_2022011506_-3.0000_-54.0000.dat)

    Inclusion of source info (location and time), an estimated celerity, and profiles at multiple reference times enables construction of a temporally varying grid where a node at a distance, r, from the source has an estimated time delay of, dt = r / cel, and uses the appropriate atmospheric information

    \b
    Examples:
    \t infraga utils build-g2s-grid --g2s-path g2s_dir/ --output-path grid/g2s_grid
    \t infraga utils build-g2s-grid --g2s-path g2s_dir/ --output-path grid/g2s_grid --src-info '[-20.56989, -175.379975, 2022-01-15T04:14:45]' --celerity-est 0.29

    '''
    # Parse profiles
    print("Parsing file list to determine grid and available datetimes...")
    dt_vals = []
    grid_lats = []
    grid_lons = []
    
    dir_files = os.listdir(g2s_path)
    for file in np.sort(dir_files):
        if fnmatch.fnmatch(file, "g2stxt_*.dat"):
            temp = file.split("_")
            dt_vals = dt_vals + [np.datetime64(temp[1][:4] + "-" + temp[1][4:6] + "-" + temp[1][6:8] + "T" + temp[1][8:10] + ":00:00")]
            grid_lats = grid_lats + [float(temp[2])]
            grid_lons = grid_lons + [float(temp[3][:-4])]

    dt_vals = np.sort(np.unique(dt_vals))
    grid_lats = np.sort(np.unique(grid_lats))
    grid_lons = np.sort(np.unique(grid_lons))

    print("File Summary:")
    print('\t' + "Unique latitudes:", grid_lats)
    print('\t' + "Unique longitudes:", grid_lons)
    print('\t' + "Unique times:", dt_vals)

    # check for longitude wrapping
    if min(grid_lons) < -179.0 and max(grid_lons) > 179.0:
        print("Detected global grid (longitudes extend to +/-180)." + '\n' + "Wrapping with extra longitude values...")
        dlon = grid_lons[1] - grid_lons[0]
        grid_lons = np.concatenate([[-180.0 - dlon], grid_lons, [180.0 + dlon]])

    # check for source info and convert
    if src_info is not None:
        for char in " ()[]":
            src_info = src_info.replace(char, "")
        src_loc = [float(val) for val in src_info.split(',')[:2]]
        src_time = np.datetime64(src_info.split(',')[2])

    print('\n' + "Mapping files to grid nodes...")
    prof_index = 0
    for lat in grid_lats:
        for lon in grid_lons:
            if lon > 180.0:
                lon = lon - 360.0
            elif lon < -180.0:
                lon = lon + 360.0

            if src_info is not None and len(dt_vals) > 1:

                prop_time = (sph_proj.inv(src_loc[1], src_loc[0], lon, lat, radians=False)[2] / 1000.0) / celerity_est
                est_arrival_time = src_time + np.timedelta64(int(prop_time), 's')
                ref_time_index = np.argmin(np.array([abs((t_ref - est_arrival_time).astype('m8[s]').astype(float)) for t_ref in dt_vals]))
                temp_datetime = dt_vals[ref_time_index].astype(object)
            else:
                temp_datetime = dt_vals[0].astype(object)


            datetime_info = "{:04d}".format(temp_datetime.year)
            datetime_info = datetime_info + "{:02d}".format(temp_datetime.month)
            datetime_info = datetime_info + "{:02d}".format(temp_datetime.day)
            datetime_info = datetime_info + "{:02d}".format(temp_datetime.hour)

            if g2s_path[-1] != "/":
                print('\t' + "Writing atmosphere " + g2s_path + "/g2stxt_" + datetime_info + "_{:.4f}".format(lat) + "_{:.4f}.dat".format(lon) + "  -->  " + output_path + "." + str(prof_index) + ".met")
                command = "cp " + g2s_path + "/g2stxt_" + datetime_info + "_{:.4f}".format(lat) + "_{:.4f}.dat".format(lon)
                command = command + " " + output_path + "." + str(prof_index) + ".met"

            else:
                print('\t' + "Writing atmosphere " + g2s_path + "g2stxt_" + datetime_info + "_{:.4f}".format(lat) + "_{:.4f}.dat".format(lon) + "  -->  " + output_path + "." + str(prof_index) + ".met")
                command = "cp " + g2s_path + "g2stxt_" + datetime_info + "_{:.4f}".format(lat) + "_{:.4f}.dat".format(lon)
                command = command + " " + output_path + "." + str(prof_index) + ".met"

            subprocess.call(command, shell=True)
            prof_index = prof_index + 1 

    file_out = open(output_path + ".lats.dat", 'w')
    for lat in grid_lats:
        print(lat, file=file_out)
    file_out.close()

    file_out = open(output_path + ".lons.dat", 'w')
    for lon in grid_lons:
        print(lon, file=file_out)
    file_out.close()

    print('\n' + "Finished grid construction." + '\n' + "Run propagation simulations using:")
    print('\t' + "infraga sph prop --atmo-prefix " + output_path + ". --grid-lats " + output_path + ".lats.dat --grid-lons " + output_path + ".lons.dat" + '\n')


############################
##    ECMWF Extraction    ##
############################

#  US standard atmosphere  #
#      polynomial fit      #
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
    etopo1 = Dataset(find_spec('infraga').submodule_search_locations[0] + "/ETOPO1_Ice_g_gmt4.grd")
    grid_lons = etopo1.variables['x'][:]
    grid_lats = etopo1.variables['y'][:]
    grid_elev = etopo1.variables['z'][:]

    lat_mask = np.logical_and(ll_corner[0] - 2.0 <= grid_lats, grid_lats <= ur_corner[0] + 2.0).nonzero()[0]
    lon_mask = np.logical_and(ll_corner[1] - 2.0 <= grid_lons, grid_lons <= ur_corner[1] + 2.0).nonzero()[0]

    region_lat = grid_lats[lat_mask]
    region_lon = grid_lons[lon_mask]
    region_elev = grid_elev[lat_mask,:][:,lon_mask]

    # Change underwater values to sea surface
    region_elev[region_elev < 0.0] = 0.0

    return interpolate.interp2d(region_lon, region_lat, region_elev / 1000.0, kind='linear')


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

    np.savetxt(output_id + "lats.dat", grid_lats)
    np.savetxt(output_id + "lons.dat", grid_lons)

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

            T_interp = interpolate.interp1d(z_gl + z_agl, ecmwf.variables['T'][:, n_lat, n_lon].data, bounds_error=False, fill_value="extrapolate")
            U_interp = interpolate.interp1d(z_gl + z_agl, ecmwf.variables['U'][:, n_lat, n_lon].data, bounds_error=False, fill_value="extrapolate")
            V_interp = interpolate.interp1d(z_gl + z_agl, ecmwf.variables['V'][:, n_lat, n_lon].data, bounds_error=False, fill_value="extrapolate")

            atmo = z_vals
            atmo = np.vstack((atmo, T_interp(z_vals)))
            atmo = np.vstack((atmo, U_interp(z_vals)))
            atmo = np.vstack((atmo, V_interp(z_vals)))
            atmo = np.vstack((atmo, density(z_vals)))   
            atmo = np.vstack((atmo, pressure(z_vals, T_interp(z_vals))))   

            np.savetxt(output_id + str(n1 * len(grid_n_lons) + n2) + ".met", atmo.T)


@click.command('extract-ecmwf', short_help="Extract atmospheric information from an ECMWF netCDF file")
@click.option("--ecmwf-file", help="ECMWF netCDF file")
@click.option("--option", help="Extraction option ('single' or 'grid')", prompt="Enter terrain option  ('single' or 'grid')")
@click.option("--lat1", help="Latitude of first point (latitude for 'single', lower-left corner for 'grid')", default=30.0)
@click.option("--lon1", help="Longitude of first point (longitude for 'single', lower-left corner for 'grid')", default=-110.0)
@click.option("--lat2", help="Latitude of second point (not used for 'single', upper-right corner for 'grid')", default=30.0)
@click.option("--lon2", help="Longitude of second point (not used for 'single', upper-right corner for grids)", default=-114.0)
@click.option("--sample_skips", help="Frequency of samples in the grid option (defaults to 1 to keep every node)", default=1)
@click.option("--output-path", help="Output file", prompt="Specify output path: ")
def extract_ecmwf(ecmwf_file, option, lat1, lon1, lat2, lon2, sample_skips, output_path):
    '''
    Construct the numbered specifications and grid files needed to run -rngdep methods.
    Assumes file format from https://g2s.ncpa.olemiss.edu/ (e.g., g2stxt_2022011506_-3.0000_-54.0000.dat)

    Inclusion of source info (location and time), an estimated celerity, and profiles at multiple reference times enables construction of a temporally varying grid where a node at a distance, r, from the source has an estimated time delay of, dt = r / cel, and uses the appropriate atmospheric information

    \b
    Examples:
    \t infraga utils extract-ecmwf --ecmwf-file EN19110100.nc --option single  --lat1 30.0 --lon1 -120.0 --output-path test.met
    \t infraga utils extract-ecmwf --ecmwf-file EN19110100.nc --option grid  --lat1 30.0 --lon1 -120.0 --lat2 40.0 --lon2 -110.0 --output-path test_grid

    '''

    """
    

    
    """

    if os.path.isfile(find_spec('infraga').submodule_search_locations[0] + "/ETOPO1_Ice_g_gmt4.grd"):
        if option == "single":
            extract_single(ecmwf_file, lat1, lon1, output_path)
        else:
            extract_grid(ecmwf_file, lat1, lon1, lat2, lon2, output_path, sample_skips)
    else:
        print("Topography file not found.  Downloading from https://www.ngdc.noaa.gov/mgg/global/")
        download_url = "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gmt4.grd.gz"
        destination = find_spec('infraga').submodule_search_locations[0] + "/ETOPO1_Ice_g_gmt4.grd.gz"
        try:
            wget.download(download_url, destination)
            print("Extracting...")
            os.system("gzip -d " + destination)
            print("ETOPO file successfully downloaded.  Running extraction...")

            if option == "single":
                extract_single(ecmwf_file, lat1, lon1, output_path)
            else:
                extract_grid(ecmwf_file, lat1, lon1, lat2, lon2, output_path, sample_skips)

        except:
            print("Download failed.")
            print("Try manual download: " + download_url)
            print("Place file in /path/to/infraGA/infraga/ (here .py files are located)")



################################
##     Terrain Extraction     ##
################################   
def pull_pnt2pnt(src_loc, rcvr_loc, file_out, resol=1.852, show_fig=True):
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

    print("# 'infraga extract-terrain --geom pnt2pnt' summary:", file=output)
    print("# src_loc:", src_loc, file=output)
    print("# rcvr_loc:", rcvr_loc, file=output)
    print('#% r, km', file=output)
    print('#% z, km', file=output)
    for n in range(N):
        print(sph_proj.inv(src_loc[1], src_loc[0], line_pnts[n][0], line_pnts[n][1], radians=False)[2] / 1000.0, file=output, end='\t')
        print(elev_interp(line_pnts[n][0], line_pnts[n][1])[0], file=output)
    output.close()

    if show_fig:
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


def pull_line(src_loc, azimuth, rng_max, file_out, resol=1.852, show_fig=True):
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
    print("# 'infraga extract-terrain --geom line' summary:", file=output)
    print("# lat: " + str(src_loc[0]), file=output)
    print("# lon: " + str(src_loc[1]), file=output)
    print("# azimuth: " + str(azimuth), file=output)
    print("# range: " + str(rng_max), file=output)
    print('#% r, km', file=output)
    print('#% z, km', file=output)

    for n in range(N):
        print(sph_proj.inv(src_loc[1], src_loc[0], line_pnts[n][0], line_pnts[n][1], radians=False)[2] / 1000.0, file=output, end='\t')
        print(elev_interp(line_pnts[n][0], line_pnts[n][1])[0], file=output)
    output.close()

    if show_fig:
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
    

def pull_xy_grid(src_loc, ll_corner, ur_corner, file_out, resol=1.852, show_fig=True):
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

    src_lat_check = np.logical_or(src_loc[0] < ll_corner[0], ur_corner[0] < src_loc[0])
    src_lon_check = np.logical_or(src_loc[1] < ll_corner[1], ur_corner[1] < src_loc[1])
    if src_lat_check or src_lon_check:
        print("Invalid source location: source must be within region.")
        return

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
    
    print("Writing terrain info into file:", file_out)
    output = open(file_out, 'w')
    # print("# Rng (E/W) [km]" + '\t' + "Rng (N/S) [km]"  + '\t' + "Elev [km]", file=output)
    # note: infraGA methods to ingest topo file can't recognize header notation yet, so no header in these files
    for n in range(len(x_vals)):
        for m in range(len(y_vals)):
            pnt = sph_proj.fwd(src_loc[1], src_loc[0], np.degrees(np.arctan2(x_vals[n], y_vals[m])), np.sqrt(x_vals[n]**2 + y_vals[m]**2) * 1.0e3, radians=False)
            xy_elev[n][m] = elev_interp(pnt[0], pnt[1])
            print(x_vals[n], y_vals[m], xy_elev[n][m], file=output)
    output.close()

    if show_fig:
        XX, YY = np.meshgrid(x_vals, y_vals)
        plt.pcolormesh(XX, YY, xy_elev.T, cmap=plt.cm.terrain, vmin=-1.4, vmax=5.0)
        plt.xlabel("Range (E/W) [km]")
        plt.ylabel("Range (N/S) [km]")
        plt.colorbar(label="Elevation [km]")

        plt.show()


def pull_latlon_grid(ll_corner, ur_corner, file_out, show_fig=True):
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
    etopo1 = Dataset(find_spec('infraga').submodule_search_locations[0] + "/ETOPO1_Ice_g_gmt4.grd")
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
    print("Writing terrain info into file:", file_out)
    for n in range(len(region_lat)):
        for m in range(len(region_lon)):
            print(region_lat[n], region_lon[m], max(region_elev[n][m], 0.0) / 1.0e3, file=output)
    output.close()

    if show_fig:
        print("Plotting terrain on map.")
        print('\t' + "Note: elevations below sea level are set to 0.0 in " + file_out + " (ray paths reflect from ocean surface)")
        map_proj = cartopy.crs.PlateCarree()
        LON, LAT = np.meshgrid(region_lon, region_lat)

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=map_proj)

        ax.set_xlim(ll_corner[1], ur_corner[1])
        ax.set_ylim(ll_corner[0], ur_corner[0])

        gl = ax.gridlines(crs=map_proj, draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
        gl.top_labels = False
        gl.right_labels = False

        lat_tick, lon_tick = int((ur_corner[0] - ll_corner[0]) / 5), int((ur_corner[1] - ll_corner[1]) / 5)
        gl.xlocator = mticker.FixedLocator(np.arange(ll_corner[0] - np.ceil(lon_tick / 2), ur_corner[0] + lon_tick, lon_tick))
        gl.ylocator = mticker.FixedLocator(np.arange(ll_corner[1] - np.ceil(lat_tick / 2), ur_corner[1] + lat_tick, lat_tick))
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER

        # Add features (coast lines, borders)
        ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
        ax.add_feature(cartopy.feature.BORDERS, linewidth=0.5)
        if (ur_corner[1] - ll_corner[0]) < 20.0:
            ax.add_feature(cartopy.feature.STATES, linewidth=0.5)

        cmesh = ax.pcolormesh(LON, LAT, region_elev / 1.0e3, cmap=plt.cm.terrain, transform=map_proj, vmin=-1.4, vmax=5.0)
        ax.set_xlabel("Longitude [deg]")
        ax.set_ylabel("Latitude [deg]")

        divider = make_axes_locatable(ax)
        ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
        fig.add_axes(ax_cb)
        cbar = plt.colorbar(cmesh, cax=ax_cb)
        cbar.set_label('Elevation [km]')

        plt.show()


@click.command('extract-terrain', short_help="Extract a line or grid of terrain information")
@click.option("--geom", help="Geometry option ('line', 'pnt2pnt', 'xy-grid' or 'latlon-grid')", prompt="Enter terrain option  ('line', 'pnt2pnt', 'xy-grid' or 'latlon-grid')")
@click.option("--lat1", help="Latitude of first point (starting point for 'pnt2pnt', lower-left corner for grids)", default=30.0)
@click.option("--lon1", help="Longitude of first point (starting point for 'pnt2pnt', lower-left corner for grids)", default=-110.0)
@click.option("--lat2", help="Latitude of second point (end point for 'pnt2pnt', upper-right corner for grids)", default=30.0)
@click.option("--lon2", help="Longitude of second point (end point for 'pnt2pnt', upper-right corner for grids)", default=-114.0)
@click.option("--ref-lat", help="Reference latitude of second point (0.0 for xy-grid option)", default=30.0)
@click.option("--ref-lon", help="Reference longitude of second point (0.0 for xy-grid option)", default=-110.0)
@click.option("--azimuth", help="Azimuth of great circle path for line option", default=-90.0)
@click.option("--range", help="Great circle distance for line option", default=1000.0)
@click.option("--output-file", help="Output file", prompt="Specify output file: ")
@click.option("--show-terrain", help="Visualize terrain results", default=True)
@click.option("--offline-maps-dir", help="Use directory for offline cartopy maps", default=None)
def extract_terrain(geom, lat1, lat2, lon1, lon2, ref_lat, ref_lon, azimuth, range, output_file, show_terrain, offline_maps_dir):
    '''
    Extract lines or grids of terrain information from an ETOPO1 file

    \b
    Examples:
    \t infraga utils extract-terrain --geom line --lat1 40.0 --lon1 -102.5 --azimuth -90.0 --range 750.0 --output-file line_topo.dat
    \t infraga utils extract-terrain --geom pnt2pnt --lat1 40.0 --lon1 -102.5 --lat2 40.0 --lon2 -110.0 --output-file line_topo.dat
    \t infraga utils extract-terrain --geom xy-grid --lat1 35.0 --lon1 -110.0 --lat2 45.0 --lon2 -100.0 --lat-ref 40.0 --lon-ref -105.0 --output-file xy_topo.dat
    \t infraga utils extract-terrain --geom latlon-grid --lat1 35.0 --lon1 -110.0 --lat2 45.0 --lon2 -100.0 --output-file sph_topo.dat

    '''

    if offline_maps_dir is not None:
        use_offline_maps(offline_maps_dir)

    if os.path.isfile(find_spec('infraga').submodule_search_locations[0] + "/ETOPO1_Ice_g_gmt4.grd"):
        if geom == "line":
            pull_line((lat1, lon1), azimuth, range, output_file, show_fig=show_terrain)
        elif geom == "pnt2pnt":
            pull_pnt2pnt((lat1, lon1), (lat2, lon2), output_file, show_fig=show_terrain)
        elif geom == "xy-grid":
            pull_xy_grid((ref_lat, ref_lon), (lat1, lon1), (lat2, lon2), output_file, show_fig=show_terrain)
        elif geom == "latlon-grid":
            pull_latlon_grid((lat1, lon1), (lat2, lon2), output_file, show_fig=show_terrain)
        else:
            print("Invalid geometry.  Options are ('line', 'pnt2pnt', 'xy-grid' or 'latlon-grid')")        
    else:
        print("Topography file not found.  Downloading from https://www.ngdc.noaa.gov/mgg/global/")
        download_url = "https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/ETOPO1_Ice_g_gmt4.grd.gz"
        destination = find_spec('infraga').submodule_search_locations[0] + "/ETOPO1_Ice_g_gmt4.grd.gz"
        try:
            wget.download(download_url, destination)
            print("Extracting...")
            os.system("gzip -d " + destination)
            print("ETOPO file successfully downloaded.  Running extraction...")

            if geom == "line":
                pull_line((lat1, lon1), azimuth, range, output_file, show_fig=show_terrain)
            elif geom == "pnt2pnt":
                pull_pnt2pnt((lat1, lon1), (lat2, lon2), output_file, show_fig=show_terrain)
            elif geom == "xy-grid":
                pull_xy_grid((ref_lat, ref_lon), (lat1, lon1), (lat2, lon2), output_file, show_fig=show_terrain)
            elif geom == "latlon-grid":
                pull_latlon_grid((lat1, lon1), (lat2, lon2), output_file, show_fig=show_terrain)  
            else:
                print("Invalid geometry.  Options are ('line', 'pnt2pnt', 'xy-grid' or 'latlon-grid')")        
        except:
            print("Download failed.")
            print("Try manual download: " + download_url)
            print("Place file in /path/to/infraGA/infraga/ (here .py files are located)")


