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

from importlib.util import find_spec 

import numpy as np

from scipy.interpolate import interp1d

import warnings
warnings.filterwarnings('ignore', 'loadtxt: Empty input file')


# Eigenray parameters
cpu_cnt = None
incl_min, incl_max = 0.5, 45.0
incl_step_max = 0.1
rng_max = 2000.0
iterations = 25
damping = 1.0e-3
tolerance = 0.1
verbose_output = True
keep_wvfrm_init = True

# Waveform calculation parameters
wvfrm_ref = 1.0
wvfrm_ds = 0.05
wvfrm_sps = 100

wvfrm_opt = "impulse"
wvfrm_p0 = 250.0
wvfrm_t0 = 0.1
wvfrm_alpha = 2.0

# Option: Kinney & Graham overpressure and positive phase for given yield (in kg)
# If this is set to "None", then values above are used
wvfrm_yld = None

# Option: Load a waveform from file for use.
# If this is set to "None", then values above are used
wvfrm_file = None


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


def compute_3d_wvfrm(profile, src_alt=0.0, rcvr_loc=[-400.0, 50.0, 0.0], bnc_max=1, incl_min=incl_min, incl_max=incl_max, incl_step_max=incl_step_max, rng_max=rng_max, verbose_output=verbose_output, wvfrm_ref=wvfrm_ref, wvfrm_yld=wvfrm_yld, wvfrm_file=wvfrm_file, cpu_cnt=None):
    """
        Computes eigenrays for a source and receiver separated
            by specified distances and then calculates waveform
            contributions for each unique eigenray and combines
            the waveforms into a single file.

        Parameters
        ----------
        profile : string
            Path and name of the atmosphere file
        src_alt : float
            Altitude of the source
        rcvr_loc : iterable
            Receiver distance from source and ground elevation [dx, dy, z_gr]
        bnc_max : int
            Maximum number of bounces to consider
    """
    
    # run eigenray analysis
    if cpu_cnt:
        command = "mpirun -np " + str(int(cpu_cnt)) + " " + find_spec('infraga').submodule_search_locations[0][:-8] + "/bin/infraga-accel-3d -eig_search " + profile 
    else :
        command =  find_spec('infraga').submodule_search_locations[0][:-8] + "/bin/infraga-3d -eig_search " + profile

    command = command + " src_alt=" + str(src_alt)
    command = command + " rcvr_x=" + str(rcvr_loc[0]) + " rcvr_y=" + str(rcvr_loc[1])
    command = command + " z_grnd=" + str(rcvr_loc[2]) + " bnc_max=" + str(bnc_max) + " rng_max=" + str(rng_max)

    command = command + " incl_min=" + str(incl_min) + " incl_max=" + str(incl_max)
    command = command + " iterations=" + str(iterations) + " damping=" + str(damping)
    command = command + " tolerance=" + str(tolerance) + " incl_step_max=" + str(incl_step_max)
    if verbose_output:
        command = command + " verbose=true"

    print(command)
    os.system(command)

    # load eigenray arrivals and remove duplicates
    profile_id = os.path.splitext(profile)[0]

    eig_results = np.loadtxt(profile_id + ".arrivals.dat")
    if len(eig_results) > 0:
        eig_results = np.atleast_2d(eig_results)
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
        print("# multipath_wvfrm.py -3d eigenrays", '\n#', file=file_out)
        print("# Source Alt:", src_alt, file=file_out)
        print("# Receiver:", rcvr_loc, file=file_out)
        print("# Bounce max:", bnc_max, '\n#', file=file_out)

        print("# Inclination range: (" + str(incl_min) + ", " + str(incl_max) + ")", file=file_out)
        print("# Inclination step max:", incl_step_max, file=file_out)
        print("# Range max:", rng_max, file=file_out)
        print("# Iterations:", iterations, file=file_out)
        print("# Damping:", damping, file=file_out)
        print("# Tolerance:", tolerance, '\n', file=file_out)

        wvfrms = []
        t_lims = [np.inf, 0.0]
        for n, line in enumerate(eig_results):
            command = "infraga-3d -wnl_wvfrm " + profile
            command = command + " inclination=" + str(line[0]) + " azimuth=" + str(line[1]) + " src_alt=" + str(src_alt)
            command = command + " z_grnd=" + str(rcvr_loc[2]) + " bounces=" + str(int(line[2]))
            command = command + " write_ray=true"

            if wvfrm_file:
                command = command + " wvfrm_file=" + wvfrm_file + " wvfrm_ref=" + str(wvfrm_ref)
            else:
                wvfrm_ref = 0.035 * wvfrm_yld**(1.0 / 3.0)
                p0, t0 = kg_op(wvfrm_yld, wvfrm_ref), kg_ppd(wvfrm_yld, wvfrm_ref)
                command = command + " wvfrm_p0=" + str(p0) + " wvfrm_t0="  + str(t0)
                command = command + " wvfrm_alpha=0.0 wvfrm_ref=" + str(wvfrm_ref)

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

            command = "rm " + profile_id + ".wvfrm_out.dat"
            command = command + " " + profile_id + ".raypaths.dat"
            command = command + " " + profile_id + ".wvfrm_init.dat"
            command = command + " " + profile_id + ".arrivals.dat"
            os.system(command)
        
        file_out.close()

        # combine waveforms and write to file
        print("Interpolating and merging waveforms...")
        t_vals = np.arange(t_lims[0], t_lims[1], 1.0 / float(wvfrm_sps))
    
        file_out = open(profile_id + ".wvfrms.dat", 'w')
        print("# multipath_wvfrm.py -3d waveform prediction", '\n#', file=file_out)
        print("# Source Alt:", src_alt, file=file_out)
        print("# Receiver:", rcvr_loc, file=file_out)
        print("# Bounce max:", bnc_max, '\n#', file=file_out)

        print("# Inclination range: (" + str(incl_min) + ", " + str(incl_max) + ")", file=file_out)
        print("# Inclination step max:", incl_step_max, file=file_out)
        print("# Range max:", rng_max, file=file_out)
        print("# Iterations:", iterations, file=file_out)
        print("# Damping:", damping, file=file_out)
        print("# Tolerance:", tolerance, '\n#', file=file_out)

        # Waveform calculation parameters
        print("# Waveform reference distance:", wvfrm_ref, file=file_out)
        print("# Waveform ds:", wvfrm_ds, file=file_out)
        print("# Waveform sps:", wvfrm_sps, file=file_out)

        if wvfrm_yld:
            print("# Source yield:", wvfrm_yld, file=file_out)
        else:
            print("# Waveform option:", wvfrm_opt, file=file_out)
            print("# Waveform Pk OP:", wvfrm_p0, file=file_out)
            print("# Waveform Time Sc.:", wvfrm_t0, file=file_out)
            print("# Waveform Shaping Param.:", wvfrm_alpha, file=file_out)

        print("#", file=file_out)
        print("# Eigenray arrivals:", file=file_out)
        print("# incl [deg]	az [deg]	n_b	x_0 [deg]	y_0 [deg]	time [s]	cel [km/s]	turning ht [km]	inclination [deg]	back azimuth [deg]	trans. coeff. [dB]	absorption [dB]", file=file_out)
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

def compute_sph_wvfrm(profile, src_loc=[30.0, -110.0, 0.0], rcvr_loc=[30.0, -114.0, 0.0], bnc_max=1, incl_min=incl_min, incl_max=incl_max, incl_step_max=incl_step_max, rng_max=rng_max, verbose_output=verbose_output, wvfrm_ref=wvfrm_ref, wvfrm_yld=wvfrm_yld, wvfrm_file=wvfrm_file, cpu_cnt=None):
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
            Source location [latitude, longitude, altitude]
        rcvr_loc : iterable
            Receiver location and ground elevation [latitude, longitude, z_gr]
        bnc_max : int
            Maximum number of bounces to consider
    """
    src_lat, src_lon, src_alt = src_loc[0], src_loc[1], src_loc[2]
    rcvr_lat, rcvr_lon = rcvr_loc[0], rcvr_loc[1]
    
    # run eigenray analysis
    if cpu_cnt:
        command = "mpirun -np " + str(int(cpu_cnt)) +  " " + find_spec('infraga').submodule_search_locations[0][:-8] + "/bin/infraga-accel-sph -eig_search " + profile 
    else :
        command =  find_spec('infraga').submodule_search_locations[0][:-8] + "/bin/infraga-sph -eig_search " + profile

    command = command + " src_lat=" + str(src_lat) + " src_lon=" + str(src_lon) + " src_alt=" + str(src_alt)
    command = command + " rcvr_lat=" + str(rcvr_lat) + " rcvr_lon=" + str(rcvr_lon) + " z_grnd=" + str(rcvr_loc[2])
    command = command + " bnc_max=" + str(bnc_max) + " rng_max=" + str(rng_max)

    command = command + " incl_min=" + str(incl_min) + " incl_max=" + str(incl_max)
    command = command + " iterations=" + str(iterations) + " damping=" + str(damping)
    command = command + " tolerance=" + str(tolerance) + " incl_step_max=" + str(incl_step_max)
    if verbose_output:
        command = command + " verbose=true"

    print(command)
    os.system(command)

    profile_id = os.path.splitext(profile)[0]
    eig_results = np.loadtxt(profile_id + ".arrivals.dat")
    if len(eig_results) > 0:
        eig_results = np.atleast_2d(eig_results)

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
        print("# multipath_wvfrm.py -sph eigenrays", '\n#', file=file_out)
        print("# Source:", src_loc, file=file_out)
        print("# Receiver:", rcvr_loc, file=file_out)
        print("# Bounce max:", bnc_max, '\n#', file=file_out)

        print("# Inclination range: (" + str(incl_min) + ", " + str(incl_max) + ")", file=file_out)
        print("# Inclination step max:", incl_step_max, file=file_out)
        print("# Range max:", rng_max, file=file_out)
        print("# Iterations:", iterations, file=file_out)
        print("# Damping:", damping, file=file_out)
        print("# Tolerance:", tolerance, '\n', file=file_out)

        wvfrms = []
        t_lims = [np.inf, 0.0]
        for n, line in enumerate(eig_results):
            command = "infraga-sph -wnl_wvfrm " + profile
            command = command + " src_lat=" + str(src_lat) + " src_lon=" + str(src_lon) + " src_alt=" + str(src_alt)
            command = command + " inclination=" + str(line[0]) + " azimuth=" + str(line[1]) + " z_grnd=" + str(rcvr_loc[2])
            command = command + " bounces=" + str(int(line[2])) + " write_ray=true"

            if wvfrm_file:
                command = command + " wvfrm_file=" + wvfrm_file + " wvfrm_ref=" + str(wvfrm_ref)
            else:
                wvfrm_ref = 0.035 * wvfrm_yld**(1.0 / 3.0)
                p0, t0 = kg_op(wvfrm_yld, wvfrm_ref), kg_ppd(wvfrm_yld, wvfrm_ref)
                command = command + " wvfrm_p0=" + str(p0) + " wvfrm_t0="  + str(t0)
                command = command + " wvfrm_alpha=0.1 wvfrm_ref=" + str(wvfrm_ref)

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

            command = "rm " + profile_id + ".wvfrm_out.dat"
            command = command + " " + profile_id + ".raypaths.dat"
            command = command + " " + profile_id + ".wvfrm_init.dat"
            command = command + " " + profile_id + ".arrivals.dat"
            os.system(command)
        file_out.close()
 
        # combine waveforms and write to file
        print("Interpolating and merging waveforms...")
        print('\t' + "Eigenrays written into " + profile_id + ".eigenrays.dat")
        print('\t' + "Arrival waveform written into " + profile_id + ".wvfrms.dat")


        t_vals = np.arange(t_lims[0], t_lims[1], 1.0 / float(wvfrm_sps))
    
        file_out = open(profile_id + ".wvfrms.dat", 'w')
        print("# multipath_wvfrm.py -sph waveform prediction", '\n#', file=file_out)
        print("# Source:", src_loc, file=file_out)
        print("# Receiver:", rcvr_loc, file=file_out)
        print("# Bounce max:", bnc_max, '\n#', file=file_out)

        print("# Inclination range: (" + str(incl_min) + ", " + str(incl_max) + ")", file=file_out)
        print("# Inclination step max:", incl_step_max, file=file_out)
        print("# Range max:", rng_max, file=file_out)
        print("# Iterations:", iterations, file=file_out)
        print("# Damping:", damping, file=file_out)
        print("# Tolerance:", tolerance, '\n#', file=file_out)

        # Waveform calculation parameters
        print("# Waveform reference distance:", wvfrm_ref, file=file_out)
        print("# Waveform ds:", wvfrm_ds, file=file_out)
        print("# Waveform sps:", wvfrm_sps, file=file_out)

        if wvfrm_yld:
            print("# Source yield:", wvfrm_yld, file=file_out)
        else:
            print("# Waveform option:", wvfrm_opt, file=file_out)
            print("# Waveform Pk OP:", wvfrm_p0, file=file_out)
            print("# Waveform Time Sc.:", wvfrm_t0, file=file_out)
            print("# Waveform Shaping Param.:", wvfrm_alpha, file=file_out)

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
        print('\n' + "No waveforms to compute.")

def usage():
    print('\n\t' + "#" * 32)
    print('\t' + "#" * 2 + "  InfraGA/GeoAc Multi-Path  " + "#" * 2)
    print('\t' + "#" * 2 + "     Waveform Calculator    " + "#" * 2)
    print('\t' + "#" * 32)

    print('\n' + "Run infraGA/GeoAc eigenray analysis for a specified source and receiver")
    print('\t' + "pair and use a defined source waveform to compute and combine each ")
    print('\t' + "eigenray's contribution to the arrival waveform.")

    print('\n' + "Usage: python multipath_wvfrm.py [option] profile [parameter values]")

    print('\n' + "Options and parameters (all parameters required):")
    print('\t' + "-3d (run eigenray and waveform analysis in Cartesian geometry)")
    print('\t\t' + "Parameter" + '\t\t\t' + "Units")
    print('\t\t' + "-" * 40)
    print('\t\t' + "source altitude" + '\t\t\t' + "km")
    print('\t\t' + "receiver dx" + '\t\t\t' + "km")
    print('\t\t' + "receiver dy" + '\t\t\t' + "km")
    print('\t\t' + "ground elevation" + '\t\t' + "km")
    print('\t\t' + "bnc_max" + '\t\t\t\t' + "-")

    print('\n\t' + "-sph (run eigenray and waveform analysis in global geometry)")
    print('\t\t' + "Parameter" + '\t\t\t' + "Units")
    print('\t\t' + "-" * 40)
    print('\t\t' + "source latitude" + '\t\t\t' + "degrees")
    print('\t\t' + "source longitude" + '\t\t' + "degrees")
    print('\t\t' + "source altitude" + '\t\t\t' + "km")

    print('\t\t' + "receiver latitude" + '\t\t' + "degrees")
    print('\t\t' + "receiver longitude" + '\t\t' + "degrees")
    print('\t\t' + "ground elevation" + '\t\t' + "km")
    print('\t\t' + "bnc_max" + '\t\t\t\t' + "-")

    print('\n' + "Examples:")
    print('\t' + "python multipath_wvfrm.py -3d  ../examples/ToyAtmo.met 0.0 -400.0 25.0 0.0 1")
    print('\t' + "python multipath_wvfrm.py -sph ../examples/ToyAtmo.met 30.0 -110.0 0.0 30.0 -114.0 1.0 1" + '\n')




def run(specification, option, src_lat, src_lon, src_alt, rcvr_x, rcvr_y, rcvr_lat, rcvr_lon, z_grnd, bnc_max, incl_min=0.5, incl_max=45.0, incl_step_max=0.1, rng_max=2000.0, verbose_output=True, wvfrm_ref=1.0, wvfrm_yld=10.0e3, wvfrm_file=None, cpu_cnt=None):
    if option == "3d":
        compute_3d_wvfrm(specification, src_alt=src_alt, rcvr_loc=(rcvr_x, rcvr_y, z_grnd), bnc_max=bnc_max, incl_min=incl_min, incl_max=incl_max, incl_step_max=incl_step_max, rng_max=rng_max, verbose_output=verbose_output, wvfrm_ref=wvfrm_ref, wvfrm_yld=wvfrm_yld, wvfrm_file=wvfrm_file, cpu_cnt=cpu_cnt)
    elif option == "sph":
        compute_sph_wvfrm(specification, src_loc=(src_lat, src_lon, src_alt), rcvr_loc=(rcvr_lat, rcvr_lon, z_grnd), bnc_max=bnc_max, incl_min=incl_min, incl_max=incl_max, incl_step_max=incl_step_max, rng_max=rng_max, verbose_output=verbose_output, wvfrm_ref=wvfrm_ref, wvfrm_yld=wvfrm_yld, wvfrm_file=wvfrm_file, cpu_cnt=cpu_cnt)




if __name__ == '__main__':
    if len(sys.argv) < 2:
        usage()
    else:
        if sys.argv[1] == "-3d":
            compute_3d_wvfrm(sys.argv[2], src_alt=float(sys.argv[3]), rcvr_loc=(float(sys.argv[4]), float(sys.argv[5]), float(sys.argv[6])), bnc_max=int(sys.argv[7]))
        elif sys.argv[1] == "-sph":
            compute_sph_wvfrm(sys.argv[2], src_loc=(float(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])), rcvr_loc=(float(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8])), bnc_max=int(sys.argv[9]))
        else:
            usage()