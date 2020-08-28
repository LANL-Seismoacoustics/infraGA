#!/usr/bin/env python
"""
multipath_wvfrm.py

Compute all eigenrays for a source-receiver pair
    and use the weakly non-linear waveform methods
    to compute each contribution to the waveform

Philip Blom (pblom@lanl.gov)

"""

import sys
import os

import numpy as np

from scipy.interpolate import interp1d

import warnings
warnings.filterwarnings('ignore', 'loadtxt: Empty input file')


# Eigenray parameters
cpu_cnt = None
incl_min, incl_max = 0.5, 45.0
incl_step_max = 0.05
iterations = 25
damping = 1.0e-3
tolerance = 0.1

# Waveform calculation parameters
wvfrm_ref = 2.0
wvfrm_ds = 0.05
wvfrm_sps = 100

wvfrm_opt = "impulse"
wvfrm_p0 = 500.0
wvfrm_t0 = 0.2
wvfrm_alpha = 4.0

# Option: Kinney & Graham overpressure and positive phase for given yield (in kg)
# If this is set to "None", then values above are used
wvfrm_yld = 50e3


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


def compute_3d_wvfrm(profile, rcvr_loc=[-400.0, 50.0, 0.0], bnc_max=1):
    """
        Computes eigenrays for a source and receiver separated
            by specified distances and then calculates waveform
            contributions for each unique eigenray and combines
            the waveforms into a single file.

        Parameters
        ----------
        profile : string
            Path and name of the atmosphere file
        rcvr_loc : iterable
            Receiver distance from source [dx, dy]
        bnc_max : int
            Maximum number of bounces to consider
    """
    
    # run eigenray analysis
    if cpu_cnt:
        command = "mpirun -np " + str(int(cpu_cnt)) + " infraga-accel-3d -eig_search " + profile 
    else :
        command = "infraga-3d -eig_search " + profile
    command = command + " rcvr_x=" + str(rcvr_loc[0]) + " rcvr_y=" + str(rcvr_loc[1])
    command = command + " z_grnd=" + str(rcvr_loc[2]) + " bnc_max=" + str(bnc_max)

    command = command + " incl_min=" + str(incl_min) + " incl_max=" + str(incl_max)
    command = command + " iterations=" + str(iterations) + " damping=" + str(damping)
    command = command + " tolerance=" + str(tolerance) + " incl_step_max=" + str(incl_step_max)

    print(command)
    os.system(command)

    # load eigenray arrivals and remove duplicates
    profile_id = os.path.splitext(profile)[0]

    eig_results = np.loadtxt(profile_id + ".arrivals.dat")
    if len(eig_results) > 0:
        os.system("rm " + profile_id + ".eigenray*")

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
        file_out = open(profile_id + ".eigenrays.dat", 'w')

        wvfrms = []
        t_lims = [np.inf, 0.0]
        for n, line in enumerate(eig_results):
            command = "infraga-3d -wnl_wvfrm " + profile
            command = command + " inclination=" + str(line[0]) + " azimuth=" + str(line[1])
            command = command + " z_grnd=" + str(rcvr_loc[2]) + " bounces=" + str(int(line[2]))
            command = command + " write_ray=true"

            if wvfrm_yld:
                p0, t0 = kg_op(wvfrm_yld, wvfrm_ref), kg_ppd(wvfrm_yld, wvfrm_ref)
                command = command + " wvfrm_p0=" + str(p0) + " wvfrm_t0="  + str(t0)
                command = command + " wvfrm_alpha=0.0 wvfrm_ref=" + str(wvfrm_ref)
            else:
                command = command + " wvfrm_p0=" + str(wvfrm_p0) + " wvfrm_t0="  + str(wvfrm_t0)
                command = command + " wvfrm_alpha=" + str(wvfrm_alpha) + " wvfrm_ref=" + str(wvfrm_ref)

            print(command)
            os.system(command)

            temp = np.loadtxt(profile_id + ".wvfrm_out.dat")
            wvfrms += [interp1d(temp[:, 0], temp[:, 1], bounds_error=False, fill_value=0.0,  kind='cubic')]
            t_lims[0] = min(t_lims[0], temp[0][0])
            t_lims[1] = max(t_lims[1], temp[-1][0])

            temp = np.loadtxt(profile_id + ".raypaths.dat")
            for line in temp:
                print(*line, file=file_out)
            print('\n', file=file_out)

            command = "rm " + profile_id + ".wvfrm_init.dat"
            command = command + " " + profile_id + ".wvfrm_out.dat"
            command = command + " " + profile_id + ".raypaths.dat"
            os.system(command)
        
        file_out.close()

        # combine waveforms and write to file
        print("Interpolating and merging waveforms...")
        t_vals = np.arange(t_lims[0], t_lims[1], 1.0 / float(wvfrm_sps))
    
        file_out = open(profile_id + ".wvfrms.dat", 'w')
        for n in range(len(t_vals)):
            print(t_vals[n], end='\t', file=file_out)
            for wvfrm in wvfrms:
                print(wvfrm(t_vals[n]), end='\t', file=file_out)
            print('', file=file_out)
        file_out.close()
    else:
        print('\n' + "No waveforms to compute.")

def compute_sph_wvfrm(profile, src_loc=[30.0, -110.0], rcvr_loc=[30.0, -114.0, 1.0], bnc_max=1):
    """
        Computes eigenrays for a source and receiver at specified
            latitude and longitude points and then calculates
            waveform contributions for each unique eigenray and
            combines the waveforms into a single file

        Parameters
        ----------
        profile : string
            Path and name of the atmosphere file
        src_loc : float
            Source location [latitude, longitude]
        rcvr_loc : iterable
            Receiver location [latitude, longitude
        bnc_max : int
            Maximum number of bounces to consider
    """
    src_lat, src_lon = src_loc[0], src_loc[1]
    rcvr_lat, rcvr_lon = rcvr_loc[0], rcvr_loc[1]
    
    # run eigenray analysis
    if cpu_cnt:
        command = "mpirun -np " + str(int(cpu_cnt)) + " infraga-accel-sph -eig_search " + profile 
    else :
        command = "infraga-sph -eig_search " + profile
    command = command + " src_lat=" + str(src_lat) + " src_lon=" + str(src_lon)
    command = command + " rcvr_lat=" + str(rcvr_lat) + " rcvr_lon=" + str(rcvr_lon)
    command = command + " z_grnd=" + str(rcvr_loc[2]) + " bnc_max=" + str(bnc_max)

    command = command + " incl_min=" + str(incl_min) + " incl_max=" + str(incl_max)
    command = command + " iterations=" + str(iterations) + " damping=" + str(damping)
    command = command + " tolerance=" + str(tolerance) + " incl_step_max=" + str(incl_step_max)
    command = command + " verbose=true"

    print(command)
    os.system(command)

    profile_id = os.path.splitext(profile)[0]
    eig_results = np.loadtxt(profile_id + ".arrivals.dat")
    if len(eig_results) > 0:
        os.system("rm " + profile_id + ".eigenray*")

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
        file_out = open(profile_id + ".eigenrays.dat", 'w')

        wvfrms = []
        t_lims = [np.inf, 0.0]
        for n, line in enumerate(eig_results):
            command = "infraga-sph -wnl_wvfrm " + profile
            command = command + " src_lat=" + str(src_lat) + " src_lon=" + str(src_lon)
            command = command + " inclination=" + str(line[0]) + " azimuth=" + str(line[1])
            command = command + " z_grnd=" + str(rcvr_loc[2]) + " bounces=" + str(int(line[2]))
            command = command + " write_ray=true"

            if wvfrm_yld:
                p0, t0 = kg_op(wvfrm_yld, wvfrm_ref), kg_ppd(wvfrm_yld, wvfrm_ref)
                command = command + " wvfrm_p0=" + str(p0) + " wvfrm_t0="  + str(t0)
                command = command + " wvfrm_alpha=0.0 wvfrm_ref=" + str(wvfrm_ref)
            else:
                command = command + " wvfrm_p0=" + str(wvfrm_p0) + " wvfrm_t0="  + str(wvfrm_t0)
                command = command + " wvfrm_alpha=" + str(wvfrm_alpha) + " wvfrm_ref=" + str(wvfrm_ref)

            print(command)
            os.system(command)

            temp = np.loadtxt(profile_id + ".wvfrm_out.dat")
            wvfrms += [interp1d(temp[:, 0], temp[:, 1], bounds_error=False, fill_value=0.0,  kind='cubic')]
            t_lims[0] = min(t_lims[0], temp[0][0])
            t_lims[1] = max(t_lims[1], temp[-1][0])

            temp = np.loadtxt(profile_id + ".raypaths.dat")
            for line in temp:
                print(*line, file=file_out)
            print('\n', file=file_out)

            command = "rm " + profile_id + ".wvfrm_init.dat"
            command = command + " " + profile_id + ".wvfrm_out.dat"
            command = command + " " + profile_id + ".raypaths.dat"
            os.system(command)
        file_out.close()
 
        # combine waveforms and write to file
        print("Interpolating and merging waveforms...")
        t_vals = np.arange(t_lims[0], t_lims[1], 1.0 / float(wvfrm_sps))
    
        file_out = open(profile_id + ".wvfrms.dat", 'w')
        for n in range(len(t_vals)):
            print(t_vals[n], end='\t', file=file_out)
            for wvfrm in wvfrms:
                print(wvfrm(t_vals[n]), end='\t', file=file_out)
            print('', file=file_out)
        file_out.close()
    else:
        print('\n' + "No waveforms to compute.")

def print_usage():
    print('\n\t' + "#" * 32)
    print('\t' + "#" * 2 + "  InfraGA/GeoAc Multi-Path  " + "#" * 2)
    print('\t' + "#" * 2 + "     Waveform Calculator    " + "#" * 2)
    print('\t' + "#" * 32)

    print('\n' + "Run infraGA/GeoAc eigenray analysis for a specified source and receiver")
    print('\t' + "pair and use a defined source waveform to compute and combine each ")
    print('\t' + "eigenray's contribution to the arrival waveform.")

    print('\n' + "Usage: python multipath_wvfrm.py [option] profile [parameter values]")

    print('\n' + "Options and parameters:")
    print('\t' + "-3d (run eigenray and waveform analysis in Cartesian geometry)")
    print('\t\t' + "Parameter" + '\t\t\t' + "Units")
    print('\t\t' + "-" * 50)
    print('\t\t' + "receiver x" + '\t\t\t' + "km")
    print('\t\t' + "receiver y" + '\t\t\t' + "km")
    print('\t\t' + "ground elevation" + '\t\t\t' + "km")
    print('\t\t' + "bnc_max" + '\t\t\t\t' + "-")

    print('\n\t' + "-sph (run eigenray and waveform analysis in global geometry)")
    print('\t\t' + "Parameter" + '\t\t\t' + "Units")
    print('\t\t' + "-" * 50)
    print('\t\t' + "source latitude" + '\t\t\t' + "degrees")
    print('\t\t' + "source longitude" + '\t\t' + "degrees")
    print('\t\t' + "receiver latitude" + '\t\t' + "degrees")
    print('\t\t' + "receiver longitude" + '\t\t' + "degrees")
    print('\t\t' + "ground elevation" + '\t\t\t' + "km")
    print('\t\t' + "bnc_max" + '\t\t\t\t' + "-")

    print('\n' + "Examples:")
    print('\t' + "python multipath_wvfrm.py -3d  ../examples/ToyAtmo.met -400.0 25.0 0.0 1")
    print('\t' + "python multipath_wvfrm.py -sph ../examples/ToyAtmo.met 30.0 -110.0 30.0 -114.0 1.0 1" + '\n')


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print_usage()
    else:
        if sys.argv[1] == "-3d":
            compute_3d_wvfrm(sys.argv[2],rcvr_loc=(float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])), bnc_max=int(sys.argv[6]))
        elif sys.argv[1] == "-sph":
            compute_sph_wvfrm(sys.argv[2], src_loc=(float(sys.argv[3]), float(sys.argv[4])), rcvr_loc=(float(sys.argv[5]), float(sys.argv[6]), float(sys.argv[7])), bnc_max=int(sys.argv[8]))
        else:
            print_usage()