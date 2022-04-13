#!/usr/bin/env python

from datetime import date
import os 
import subprocess
import re 
import fnmatch

import numpy as np

from pyproj import Geod


sph_proj = Geod(ellps='sphere')


def run(profiles_path, output_path, src_info=None, celerity_est=0.29):

    # Parse profiles
    print("Parsing file list to determine grid and available datetimes...")
    dt_vals = []
    grid_lats = []
    grid_lons = []
    
    dir_files = os.listdir(profiles_path)
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

            if profiles_path[-1] != "/":
                print('\t' + "Writing atmosphere " + profiles_path + "/g2stxt_" + datetime_info + "_{:.4f}".format(lat) + "_{:.4f}.dat".format(lon) + "  -->  " + output_path + "." + str(prof_index) + ".met")
                command = "cp " + profiles_path + "/g2stxt_" + datetime_info + "_{:.4f}".format(lat) + "_{:.4f}.dat".format(lon)
                command = command + " " + output_path + "." + str(prof_index) + ".met"

            else:
                print('\t' + "Writing atmosphere " + profiles_path + "g2stxt_" + datetime_info + "_{:.4f}".format(lat) + "_{:.4f}.dat".format(lon) + "  -->  " + output_path + "." + str(prof_index) + ".met")
                command = "cp " + profiles_path + "g2stxt_" + datetime_info + "_{:.4f}".format(lat) + "_{:.4f}.dat".format(lon)
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
    print('\t' + "infraga-sph-rngdep -prop " + output_path + ". " + output_path + ".lats.dat " + output_path + ".lons.dat" + '\n')



if __name__ == '__main__':

    src_loc = [-20.56989, -175.379975]
    src_time = np.datetime64("2022-01-15T04:14:45")

    src_info = [-20.56989, -175.379975, '2022-01-15T04:14:45']

    # grid_lats = np.linspace(-90.0, 90.0, 61)
    # grid_lons = np.linspace(-180.0, 180.0, 61)

    profiles_path = "Tonga_g2s"
    output_path = "rngdep/tonga"

    celerity_est = 0.29
    use_prop_correction = True
    ref_datetime = src_time

    # Parse profiles
    print("Parsing file list to determine grid and available datetimes...")
    dt_vals = []
    grid_lats = []
    grid_lons = []
    
    dir_files = os.listdir(profiles_path)
    for file in np.sort(dir_files):
        if fnmatch.fnmatch(file, "g2stxt_*.dat"):
            temp = file.split("_")
            dt_vals = dt_vals + [np.datetime64(temp[1][:4] + "-" + temp[1][4:6] + "-" + temp[1][6:8] + "T" + temp[1][8:10] + ":00:00")]
            grid_lats = grid_lats + [float(temp[2])]
            grid_lons = grid_lons + [float(temp[3][:-4])]

    dt_vals = np.sort(np.unique(dt_vals))
    grid_lats = np.sort(np.unique(grid_lats))
    grid_lons = np.sort(np.unique(grid_lons))

    # check for longitude wrapping
    if min(grid_lons) < -179.0 and max(grid_lons) > 179.0:
        dlon = grid_lons[1] - grid_lons[0]
        grid_lons = np.concatenate([[-180.0 - dlon], grid_lons, [180.0 + dlon]])

    prof_index = 0
    for lat in grid_lats:
        for lon in grid_lons:
            print("Setting grid node at " + str(lat) + ", " + str(lon))

            if lon > 180.0:
                lon = lon - 360.0
            elif lon < -180.0:
                lon = lon + 360.0

            if use_prop_correction:
                prop_time = (sph_proj.inv(src_loc[1], src_loc[0], lon, lat, radians=False)[2] / 1000.0) / celerity_est
                est_arrival_time = src_time + np.timedelta64(int(prop_time), 's')
                ref_time_index = np.argmin(np.array([abs((t_ref - est_arrival_time).astype('m8[s]').astype(float)) for t_ref in dt_vals]))
                temp_datetime = dt_vals[ref_time_index].astype(object)

                print('\tEstimated arrival time:', est_arrival_time)
                print('\tNearest reference time:', dt_vals[ref_time_index])

            else:
                temp_datetime = ref_datetime.astype(object)

            datetime_info = "{:04d}".format(temp_datetime.year)
            datetime_info = datetime_info + "{:02d}".format(temp_datetime.month)
            datetime_info = datetime_info + "{:02d}".format(temp_datetime.day)
            datetime_info = datetime_info + "{:02d}".format(temp_datetime.hour)

            command = "cp " + profiles_path + "/g2stxt_" + datetime_info + "_{:.4f}".format(lat) + "_{:.4f}.dat".format(lon)
            command = command + " " + output_path + str(prof_index) + ".met"

            # print(command)
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
    