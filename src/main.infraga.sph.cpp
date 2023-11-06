#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <sstream>
#include <math.h>

#include "atmo/atmo_state.h"
#include "atmo/atmo_io.sph.strat.h"

#include "geoac/geoac.params.h"
#include "geoac/geoac.eqset.h"
#include "geoac/geoac.interface.h"
#include "geoac/geoac.eigenray.h"

#include "util/fileIO.h"
#include "util/rk4solver.h"
#include "util/interpolation.h"
#include "util/globe.h"
#include "util/waveforms.h"

using namespace std;

void version(){
    cout << '\n' << '\t' << "infraga v 1.0.3" << '\n';
    cout << '\t' << "Copyright (c) 2014, Triad National Security, LLC.  All rights reserved." << '\n';
    cout << '\t' << "This software was produced under U.S. Government contract 89233218CNA000001 "<< '\n';
    cout << '\t' << "for Los Alamos National Laboratory (LANL), which is operated by Triad" << '\n';
    cout << '\t' << "National Security, LLC, for the U.S. Department of Energy/National Nuclear" << '\n';
    cout << '\t' << "Security Administration." << '\n' << '\n';
    
    cout << '\t' << "All rights in the program are reserved by Triad National Security, LLC, and " << '\n';
    cout << '\t' << "the U.S. Department of Energy/National Nuclear Security Administration. The" << '\n';
    cout << '\t' << "Government is granted for itself and others acting on its behalf a nonexclusive," << '\n';
    cout << '\t' << "paid-up, irrevocable worldwide license in this material to reproduce, prepare " << '\n';
    cout << '\t' << "derivative works, distribute copies to the public, perform publicly and display" << '\n';
    cout << '\t' << "publicly, and to permit others to do so." << '\n' << '\n';
    
    cout << '\t' << "License: Open Source MIT <http://opensource.org/licenses/MIT>" << '\n';
    cout << '\t' << "This is free software: you are free to change and redistribute it." << '\n';
    cout << '\t' << "The software is provided 'AS IS', without warranty of any kind, expressed or implied." << '\n';
    cout << '\t' << "See manual for full licence information." << '\n' << '\n';
}

void usage(){
    cout << '\n';
    cout << '\t' << "#####################################################" << '\n';
    cout << '\t' << "####                 infraga-sph                 ####" << '\n';
    cout << '\t' << "####   Ray Tracing Using Spherical Coordinates   ####" << '\n';
    cout << '\t' << "####       and a Stratified, Moving Medium       ####" << '\n';
    cout << '\t' << "#####################################################" << '\n' << '\n';
    
    cout << '\n';
    cout << "Usage: infraga-sph [option] profile.met [parameters]" << '\n';
    cout << '\t' << '\t' << "Enter only 1 option." << '\n';
    cout << '\t' << '\t' << "Profile.met is expected to contain columns describing {z[km]  T[K]  u (zonal wind) [m/s]  v (meridional wind) [m/s]  density[g/cm^3]  p[mbar]} " << '\n';
    cout << '\t' << '\t' << '\t' << "Profile format can be modified, see manual for details." << '\n';
    cout << '\t' << '\t' << "Parameter calls are expected using the format: parameter_name=value." << '\n' << '\n';
    
    cout << "Options and parameters are:" << '\n';
    cout << '\t' << "-prop (generate ray paths for propagations at multiple azimuth and inclination angles)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units" << '\t' << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "incl_min"          << '\t' << "degrees"            << '\t' << '\t' << "0.5" << '\n';
    cout << '\t' << '\t' << "incl_max"          << '\t' << "degrees"            << '\t' << '\t' << "45.0" << '\n';
    cout << '\t' << '\t' << "incl_step"         << '\t' << "degrees"            << '\t' << '\t' << "0.5"  << '\n';
    cout << '\t' << '\t' << "inclination"       << '\t' << "see manual"         << '\t' << "0.5" << '\n';
    cout << '\t' << '\t' << "az_min"            << '\t' << '\t' << "degrees"    << '\t' << '\t' << "-90.0" << '\n';
    cout << '\t' << '\t' << "az_max"            << '\t' << '\t' << "degrees"    << '\t' << '\t' << "-90.0" << '\n';
    cout << '\t' << '\t' << "az_step"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "1.0"  << '\n';
    cout << '\t' << '\t' << "azimuth"           << '\t' << '\t' << "see manual" << '\t' << "-90.0" << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "10" << '\n';

    cout << '\t' << '\t' << "src_lat"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "30.0" << '\n';
    cout << '\t' << '\t' << "src_lon"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "src_alt"           << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';

    cout << '\t' << '\t' << "turn_ht_min"       << '\t' << "km"                 << '\t' << '\t' << "0.2" << '\n';
    cout << '\t' << '\t' << "write_rays"        << '\t' << "true/false"         << '\t' << "true" << '\n';
    cout << '\t' << '\t' << "write_topo"        << '\t' << "true/false"         << '\t' << "false" << '\n' << '\n';
    /*
    cout << '\t' << "-back_proj (back project from a receiver towards a potential source)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "inclination"       << '\t' << "degrees"            << '\t' << '\t' << "15.0"  << '\n';
    cout << '\t' << '\t' << "azimuth"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "-90.0"  << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "0"  << '\n';
    cout << '\t' << '\t' << "rcvr_lat"          << '\t' << '\t' << "degrees"    << '\t' << '\t' << "30.0" << '\n';
    cout << '\t' << '\t' << "rcvr_lon"          << '\t' << '\t' << "degrees"    << '\t' << '\t' << "0.0" << '\n' << '\n';
    */
    cout << '\t' << "-eig_search (Search for all eigenrays connecting a source at (src_lat, src_lon, src_alt) to a receiver " << '\n';
    cout << '\t' << '\t' << '\t' << "at (rcvr_lat, rcvr_lon, z_grnd) which have inclinations and ground reflections within specified limits)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units" << '\t' << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "incl_min"          << '\t' << "degrees"            << '\t' << '\t' << "1.0" << '\n';
    cout << '\t' << '\t' << "incl_max"          << '\t' << "degrees"            << '\t' << '\t' << "45.0" << '\n';
    cout << '\t' << '\t' << "bnc_min"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "0" << '\n';
    cout << '\t' << '\t' << "bnc_max"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "0" << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "see manual" << '\t' << "0" << '\n';
    cout << '\t' << '\t' << "src_lat"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "30.0" << '\n';
    cout << '\t' << '\t' << "src_lon"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "src_alt"           << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "rcvr_lat"          << '\t' << "degrees"            << '\t' << '\t' << "30" << '\n';
    cout << '\t' << '\t' << "rcvr_lon"          << '\t' << "degrees"            << '\t' << '\t' << "-2.5" << '\n';
    cout << '\t' << '\t' << "verbose"           << '\t' << '\t' << "true/false" << '\t' << "false" << '\n';
    cout << '\t' << '\t' << "iterations"        << '\t' << "integer"            << '\t' << '\t' << "25" << '\n';
    cout << '\t' << '\t' << "damping"           << '\t' << '\t' << "scalar"     << '\t' << '\t' << "1.0e-3" << '\n';
    cout << '\t' << '\t' << "tolerance"         << '\t' << "km"                 << '\t' << '\t' << "0.1" << '\n';
    cout << '\t' << '\t' << "az_dev_lim"        << '\t' << "degrees"            << '\t' << '\t' << "2.0" << '\n';
    cout << '\t' << '\t' << "incl_step_min"     << '\t' << "degrees"            << '\t' << '\t' << "0.001" << '\n';
    cout << '\t' << '\t' << "incl_step_max"     << '\t' << "degrees"            << '\t' << '\t' << "0.1" << '\n' << '\n';
    
    cout << '\t' << "-eig_direct (Search for a single eigenray connecting a source at (src_lat, src_lon, src_alt) to a receiver " << '\n';
    cout << '\t' << '\t' << '\t' << "at (rcvr_lat, rcvr_lon, z_grnd) near an estiamted azimuth and inclination assuming a specific number of ground reflections)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units" << '\t' << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "incl_est"          << '\t' << "degrees"            << '\t' << '\t' << "15.0" << '\n';
    cout << '\t' << '\t' << "az_est"            << '\t' << '\t' << "degrees"    << '\t' << '\t' << "great circle bearing" << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "2" << '\n';
    cout << '\t' << '\t' << "src_lat"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "30.0" << '\n';
    cout << '\t' << '\t' << "src_lon"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "src_alt"           << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "rcvr_lat"          << '\t' << "degrees"            << '\t' << '\t' << "30" << '\n';
    cout << '\t' << '\t' << "rcvr_lon"          << '\t' << "degrees"            << '\t' << '\t' << "10.0" << '\n';
    cout << '\t' << '\t' << "verbose"           << '\t' << '\t' << "true/false" << '\t' << "false" << '\n';
    cout << '\t' << '\t' << "iterations"        << '\t' << "integer"            << '\t' << '\t' << "25" << '\n';
    cout << '\t' << '\t' << "damping"           << '\t' << '\t' << "scalar"     << '\t' << '\t' << "1.0e-3" << '\n';
    cout << '\t' << '\t' << "tolerance"         << '\t' << "km"                 << '\t' << '\t' << "0.1" << '\n' << '\n';

    cout << '\t' << "-wnl_wvfrm (compute the weakly non-linear waveform along a specific ray path)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "inclination"       << '\t' << "degrees"            << '\t' << '\t' << "15.0"  << '\n';
    cout << '\t' << '\t' << "azimuth"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "90.0"  << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "0"  << '\n';
    cout << '\t' << '\t' << "src_lat"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "src_lon"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "src_alt"           << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "write_ray"         << '\t' << "true/false"         << '\t' << "false" << '\n';

    cout << '\t' << '\t' << "wvfrm_file"        << '\t' << "see manual"         << '\t' << "none" << '\n';
    cout << '\t' << '\t' << "wvfrm_opt"         << '\t' << "see manual"         << '\t' << "impulse" << '\n';
    cout << '\t' << '\t' << "wvfrm_p0"          << '\t' << "Pa"                 << '\t' << '\t' << "10.0" << '\n';
    cout << '\t' << '\t' << "wvfrm_t0"          << '\t' << "seconds"            << '\t' << '\t' << "1.0" << '\n';
    cout << '\t' << '\t' << "wvfrm_alpha"       << '\t' << "-"                  << '\t' << '\t' << "1.0" << '\n';
    cout << '\t' << '\t' << "wvfrm_ref"         << '\t' << "km"                 << '\t' << '\t' << "1.0" << '\n';
    cout << '\t' << '\t' << "wvfrm_out_step"    << '\t' << "km"                 << '\t' << '\t' << "none" << '\n' << '\n';
    
    cout << '\t' << '\t' << "wvfrm_ds"          << '\t' << "see manual"         << '\t' << "0.5" << '\n';
    cout << '\t' << '\t' << "wvfrm_N_rise"      << '\t' << "see manual"         << '\t' << "10" << '\n';
    cout << '\t' << '\t' << "wvfrm_N_period"    << '\t' << "see manual"         << '\t' << "40" << '\n';
    cout << '\t' << '\t' << "wvfrm_len"         << '\t' << "see manual"         << '\t' << "pow(2, 13)" << '\n';
    cout << '\t' << '\t' << "wvfrm_filt1_g"     << '\t' << "see manual"         << '\t' << "0.5" << '\n';
    cout << '\t' << '\t' << "wvfrm_filt1_n1"    << '\t' << "see manual"         << '\t' << "4.0" << '\n';
    cout << '\t' << '\t' << "wvfrm_filt1_n2"    << '\t' << "see manual"         << '\t' << "0.5" << '\n';
    cout << '\t' << '\t' << "wvfrm_filt2_g"     << '\t' << "see manual"         << '\t' << "0.75" << '\n';
    cout << '\t' << '\t' << "wvfrm_filt2_n1"    << '\t' << "see manual"         << '\t' << "2.0" << '\n';
    cout << '\t' << '\t' << "wvfrm_filt2_n2"    << '\t' << "see manual"         << '\t' << "1.0" << '\n' << '\n';

    cout << '\t' << "Additional Parameters"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << "freq"              << '\t' << '\t' << '\t' << "Hz"         << '\t' << '\t' << "0.1" << '\n';
    cout << '\t' << "abs_coeff"         << '\t' << '\t' << "scalar"             << '\t' << '\t' << "1.0" << '\n';
    cout << '\t' << "z_grnd"            << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << "write_atmo"        << '\t' << '\t' << "true/false"         << '\t' << "false" << '\n';
    cout << '\t' << "prof_format"       << '\t' << '\t' << "see manual"         << '\t' << "zTuvdp" << '\n';
    cout << '\t' << "reverse_winds"     << '\t' << '\t' << "true/false"         << '\t' << "false" << '\n';
    cout << '\t' << "output_id"         << '\t' << '\t' << "see manual"         << '\t' << "from profile.met" << '\n';
    cout << '\t' << "write_caustics"    << '\t' << '\t' << "true/false"         << '\t' << "false" << '\n';
    cout << '\t' << "calc_amp"          << '\t' << '\t' << "true/false"         << '\t' << "true" << '\n';
    
    cout << '\t' << "max_alt"           << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "interpolation max" << '\n';
    cout << '\t' << "max_rng"           << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "1000.0" << '\n';
    cout << '\t' << "max_time"          << '\t' << '\t' << '\t' << "hr"         << '\t' << '\t' << "10.0" << '\n';
    cout << '\t' << "min_lat"           << '\t' << '\t' << '\t' << "degrees"    << '\t' << '\t' << "-90.0" << '\n';
    cout << '\t' << "max_lat"           << '\t' << '\t' << '\t' << "degrees"    << '\t' << '\t' << "+90.0" << '\n';
    cout << '\t' << "min_lon"           << '\t' << '\t' << '\t' << "degrees"    << '\t' << '\t' << "-180.0" << '\n';
    cout << '\t' << "max_lon"           << '\t' << '\t' << '\t' << "degrees"    << '\t' << '\t' << "+180.0" << '\n';
    cout << '\t' << "min_ds"            << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "0.001" << '\n';
    cout << '\t' << "max_ds"            << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "0.05" << '\n';
    cout << '\t' << "max_s"             << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "1000.0" << '\n';
    cout << '\t' << "topo_file"         << '\t' << '\t' << "see manual"         << '\t' << "none" << '\n';
    cout << '\t' << "topo_use_BLw"      << '\t' << '\t' << "see manual"         << '\t' << "false" << '\n' << '\n';

    cout << "Output (see output files or manual for units):" << '\n';
    cout << '\t' << "atmo.dat -> z[km] : c [m/s]  : u (zonal winds) [m/s] : v (meridional winds) [m/s] : density[g/cm^3] : ceff [km/s]" << '\n';
    cout << '\t' << "{...}.raypaths.dat -> lat : lon : z : trans. coeff. : absorption : time " << '\n';
    cout << '\t' << "{...}.arrivals.dat -> incl : az : n_bnc : lat : lon : time : cel : z_max : arrival incl : back az : trans. coeff. : absorption" << '\n' << '\n';
    // cout << '\t' << "{...}.projection.dat -> lat : lon : z : t : LAT_incl : LON_incl : Z_incl : T_incl : LAT_az : LON_az : Z_az : T_az" << '\n' << '\n';

    cout << "Examples:" << '\n';
    cout << '\t' << "./bin/infraga-sph -prop examples/ToyAtmo.met src_lat=30.0 src_lon=-100.0 azimuth=-90.0 max_rng=750.0" << '\n';
    cout << '\t' << "./bin/infraga-sph -eig_search examples/ToyAtmo.met src_lat=30.0 src_lon=-100.0 rcvr_lat=30.25 rcvr_lon=-104.25 bnc_max=1 verbose=true" << '\n';
    cout << '\t' << "./bin/infraga-sph -eig_direct examples/ToyAtmo.met src_lat=30.0 src_lon=-100 rcvr_lat=30.25 rcvr_lon=-104.25 bounces=1 incl_est=5.0 verbose=true" << '\n';
    //cout << '\t' << "./bin/infraga-sph -back_proj examples/ToyAtmo.met rcvr_lat=30.25 rcvr_lon=-104.25 azimuth=92.896723 inclination=4.1689992" << '\n';
    cout << '\t' << "./bin/infraga-sph -wnl_wvfrm examples/ToyAtmo.met src_lat=30.0 src_lon=-100.0 inclination=4.1314876 azimuth=-84.969455 bounces=1 wvfrm_opt=impulse wvfrm_p0=500.0" << '\n' << '\n';
}

void run_prop(char* inputs[], int count){
    cout << '\n';
    cout << '\t' << "###########################################" << '\n';
    cout << '\t' << "####        Running infraga-sph        ####" << '\n';
    cout << '\t' << "####            Propagation            ####" << '\n';
    cout << '\t' << "###########################################" << '\n' << '\n';
    
    double theta_min = 0.5, theta_max = 45.0, theta_step = 0.5;
    double phi_min = -90.0, phi_max = -90.0, phi_step = 1.0;
    int bounces=10, file_check;
    double lat_src = 30.0, lon_src = 0.0, z_src = 0.0;
    bool write_atmo=false, write_rays=true, write_caustics=false, write_topo=false, custom_output_id=false, print_resid=false, reverse_winds=false;
    double freq = 0.1, turn_ht_min = 0.2;
    char* prof_format = "zTuvdp";
    char* topo_file = "None";
    char* output_id;
    char input_check;
    
    topo::z0 = 0.0;
    atmo::tweak_abs = 1.0;

    geoac::calc_amp = true;
    geoac::is_topo = false;
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "prof_format=",12) == 0){            prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){            topo::z0 = atof(inputs[i] + 7);}
        else if (strncmp(inputs[i], "reverse_winds=", 14) == 0){    reverse_winds = string2bool(inputs[i] + 14);}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){        topo_file = inputs[i] + 10; geoac::is_topo=true;}
    }

    if (geoac::is_topo){    file_check = set_region(inputs[2], topo_file, prof_format, reverse_winds);}
    else{                   file_check = set_region(inputs[2], prof_format, reverse_winds);}
    if(!file_check){
        return;
    }

    if (geoac::is_topo){
        lat_src = (geoac::lat_max + geoac::lat_min) / 2.0 * (180.0 / Pi);
        lon_src = (geoac::lon_max + geoac::lon_min) / 2.0 * (180.0 / Pi);
    }
    
    for(int i = 3; i < count; i++){
        if ((strncmp(inputs[i], "incl_min=", 9) == 0) || (strncmp(inputs[i], "min_incl=", 9) == 0)){        theta_min = atof(inputs[i] + 9);}
        else if ((strncmp(inputs[i], "incl_max=", 9) == 0) || (strncmp(inputs[i], "max_incl=", 9) == 0)){   theta_max = atof(inputs[i] + 9);}
        else if ((strncmp(inputs[i], "az_min=", 7) == 0) || (strncmp(inputs[i], "min_az=", 7) == 0)){       phi_min = atof(inputs[i] + 7);}
        else if ((strncmp(inputs[i], "az_max=", 7) == 0) || (strncmp(inputs[i], "max_az=", 7) == 0)){       phi_max = atof(inputs[i] + 7);}
        
        else if (strncmp(inputs[i], "incl_step=", 10) == 0){                                                theta_step = atof(inputs[i] + 10);}
        else if (strncmp(inputs[i], "az_step=", 8) == 0){                                                   phi_step = atof(inputs[i] + 8);}
        
        else if (strncmp(inputs[i], "inclination=", 12) == 0){                                              theta_min = atof(inputs[i] + 12);
            theta_max = atof(inputs[i] + 12);
            theta_step = 1.0;}
        else if (strncmp(inputs[i], "azimuth=", 8) == 0){                                                   phi_min = atof(inputs[i] + 8);
            phi_max = atof(inputs[i] + 8);
            phi_step = 1.0;}
        
        else if (strncmp(inputs[i], "bounces=", 8) == 0){                                                   bounces = atoi(inputs[i] + 8);}
        
        else if ((strncmp(inputs[i], "src_alt=", 8) == 0) || (strncmp(inputs[i], "alt_src=", 8) == 0)){     z_src = atof(inputs[i] + 8);}
        else if ((strncmp(inputs[i], "src_lat=", 8) == 0) || (strncmp(inputs[i], "lat_src=", 8) == 0)){     lat_src = atof(inputs[i] + 8);}
        else if ((strncmp(inputs[i], "src_lon=", 8) == 0) || (strncmp(inputs[i], "lon_src=", 8) == 0)){     lon_src = atof(inputs[i] + 8);}

        else if (strncmp(inputs[i], "turn_ht_min=", 12) == 0){                                              turn_ht_min = atof(inputs[i] + 12);}
        else if (strncmp(inputs[i], "write_rays=", 11) == 0){                                               write_rays = string2bool(inputs[i] + 11);}
                 
        else if (strncmp(inputs[i], "freq=", 5) == 0){                                                      freq = atof(inputs[i] + 5);}
        else if (strncmp(inputs[i], "abs_coeff=", 10) == 0){                                                atmo::tweak_abs = max(0.0, atof(inputs[i] + 10));}
        else if (strncmp(inputs[i], "write_caustics=", 15) == 0){                                           write_caustics = string2bool(inputs[i] + 15);}
        else if (strncmp(inputs[i], "calc_amp=", 9) == 0){                                                  geoac::calc_amp = string2bool(inputs[i] + 9);}
        
        else if ((strncmp(inputs[i], "max_alt=", 8) == 0) || (strncmp(inputs[i], "alt_max=", 8) == 0)){     geoac::alt_max = atof(inputs[i] + 8);   cout << '\t' << "User selected altitude maximum = " << geoac::alt_max << '\n';}
        else if ((strncmp(inputs[i], "max_rng=", 8) == 0) || (strncmp(inputs[i], "rng_max=", 8) == 0)){     geoac::rng_max = atof(inputs[i] + 8);   cout << '\t' << "User selected range maximum = " << geoac::rng_max << '\n';}
        else if ((strncmp(inputs[i], "max_time=", 9) == 0) || (strncmp(inputs[i], "time_max=", 9) == 0)){   geoac::time_max = atof(inputs[i] + 9);   cout << '\t' << "User selected propagation time maximum = " << geoac::time_max << '\n';}

        else if ((strncmp(inputs[i], "min_lat=", 8) == 0) || (strncmp(inputs[i], "lat_min=", 8) == 0)){     geoac::lat_min = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected latitude minimum = " << geoac::lat_min * (180.0 / Pi) << '\n';}
        else if ((strncmp(inputs[i], "max_lat=", 8) == 0) || (strncmp(inputs[i], "lat_max=", 8) == 0)){     geoac::lat_max = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected latitude maximum = " << geoac::lat_max * (180.0 / Pi) << '\n';}
        else if ((strncmp(inputs[i], "min_lon=", 8) == 0) || (strncmp(inputs[i], "lon_min=", 8) == 0)){     geoac::lon_min = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected longitude minimum = " << geoac::lon_min * (180.0 / Pi) << '\n';}
        else if ((strncmp(inputs[i], "max_lon=", 8) == 0) || (strncmp(inputs[i], "lon_max=", 8) == 0)){     geoac::lon_max = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected longitude maximum = " << geoac::lon_max * (180.0 / Pi) << '\n';}
                 
        else if ((strncmp(inputs[i], "min_ds=", 7) == 0) || (strncmp(inputs[i], "ds_min=", 7) == 0)){       geoac::ds_min = atof(inputs[i] + 7);    cout << '\t' << "User selected ds minimum = " << geoac::ds_min << '\n';}
        else if ((strncmp(inputs[i], "max_ds=", 7) == 0) || (strncmp(inputs[i], "ds_max=", 7) == 0)){       geoac::ds_max = atof(inputs[i] + 7);    cout << '\t' << "User selected ds maximum = " << geoac::ds_max << '\n';}
        else if ((strncmp(inputs[i], "max_s=", 6) == 0) || (strncmp(inputs[i], "s_max=", 6) == 0)){         geoac::s_max = atof(inputs[i] + 6);     cout << '\t' << "User selected s maximum (path length between reflections) = " << geoac::s_max << '\n';}

        else if (strncmp(inputs[i], "write_atmo=", 11) == 0){                                               write_atmo = string2bool(inputs[i] + 11);}
        else if (strncmp(inputs[i], "prof_format=", 12) == 0){                                              prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "reverse_winds=", 14) == 0){                                            reverse_winds = string2bool(inputs[i] + 14);}

        else if (strncmp(inputs[i], "output_id=", 10) == 0){                                                custom_output_id = true; 
                                                                                                            output_id = inputs[i] + 10;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){                                                    if(!geoac::is_topo){
                                                                                                                topo::z0 = atof(inputs[i] + 7);
                                                                                                                topo::z_max = topo::z0;
                                                                                                                topo::z_bndlyr = topo::z0 + 2.0;
                                                                                                            } else {
                                                                                                                cout << '\t' << "Note: cannot adjust ground elevation with topography." << '\n';
                                                                                                            }}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){                                                topo_file = inputs[i] + 10; geoac::is_topo=true;}
        else if (strncmp(inputs[i], "topo_use_BLw=", 13) == 0){                                             topo::use_BLw = string2bool(inputs[i] + 13);}
        else if (strncmp(inputs[i], "write_topo=", 11) == 0){                                               write_topo = string2bool(inputs[i] + 11);}

        else if (strncmp(inputs[i], "print_resid=", 12) == 0){                                              print_resid = string2bool(inputs[i] + 12);}

        else{
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }

    lat_src *= Pi / 180.0;
    lon_src *= Pi / 180.0;
    z_src = max(z_src, topo::z(lat_src, lon_src) - globe::r0);
    if(write_atmo)      geoac::write_prof("atmo.dat", lat_src, lon_src, (90.0 - phi_min) * Pi / 180.0);
    if(write_caustics)  geoac::calc_amp=true;
    geoac::configure();

    cout << '\n' << "Parameter summary:" << '\n';
    cout << '\t' << "inclination: " << theta_min << ", " << theta_max << ", " << theta_step << '\n';
    cout << '\t' << "azimuth: " << phi_min << ", " << phi_max << ", " << phi_step << '\n';
    cout << '\t' << "bounces: " << bounces << '\n';
    cout << '\t' << "source location (lat, lon, alt): " << lat_src * (180.0 / Pi) << ", " << lon_src * (180.0 / Pi) << ", " << z_src << '\n';
    if(!geoac::is_topo){
        cout << '\t' << "ground elevation: " << topo::z0 << '\n';
    }
    cout << '\t' << "turn_ht_min: " << turn_ht_min << '\n';
    cout << '\t' << "frequency: " << freq << '\n';
    cout << '\t' << "S&B atten coeff: " << atmo::tweak_abs << '\n';
    
    if (write_atmo){        cout << '\t' << "write_atmo: true" << '\n';}        else { cout << '\t' << "write_atmo: false" << '\n';}
    if (write_rays){        cout << '\t' << "write_rays: true" << '\n';}        else { cout << '\t' << "write_rays: false" << '\n';}
    if (write_topo){        cout << '\t' << "write_topo: true" << '\n';}        else { cout << '\t' << "write_topo: false" << '\n';}
    if (write_caustics){    cout << '\t' << "write_caustics: true" << '\n';}    else { cout << '\t' << "write_caustics: false" << '\n';}
    if (geoac::calc_amp){   cout << '\t' << "calc_amp: true" << '\n';}          else { cout << '\t' << "calc_amp: false" << '\n';}
    cout << '\n';
    
    // Extract the file name from the input and use it to distinguish the resulting output
    char output_buffer [512];
    if(!custom_output_id){
        output_id = inputs[2];
        for(int m = strlen(output_id); m >= 0; m--){
            if(output_id[m] == '.'){
                output_id[m] = '\0';
                break;
            }
        }
    }

    // Define variables used for analysis
	double D, D_prev, travel_time_sum, attenuation, r_max, inclination, back_az;
	int k, length = int(geoac::s_max / geoac::ds_min);
	bool break_check;

	ofstream results, raypath, caustics, topo_out;

    sprintf(output_buffer, "%s.arrivals.dat", output_id);
    results.open(output_buffer);

    results << "# infraga-sph -prop summary:" << '\n';
    results << "#" << '\t' << "profile: " << inputs[2] << '\n';
    results << "#" << '\t' << "inclination: " << theta_min << ", " << theta_max << ", " << theta_step << '\n';
    results << "#" << '\t' << "azimuth: " << phi_min << ", " << phi_max << ", " << phi_step << '\n';
    results << "#" << '\t' << "bounces: " << bounces << '\n';
    results << "#" << '\t' << "source location (lat, lon, alt): " << lat_src * (180.0 / Pi) << ", " << lon_src * (180.0 / Pi) << ", " << z_src << '\n';
    if(!geoac::is_topo){
        results << "#" << '\t' << "ground elevation: " << topo::z0 << '\n';
    } else {
        results << "#" << '\t' << "topo file:" << topo_file << '\n';
    }
    results << "#" << '\t' << "turn_ht_min: " << turn_ht_min << '\n';
    results << "#" << '\t' << "frequency: " << freq << '\n';
    results << "#" << '\t' << "S&B atten coeff: " << atmo::tweak_abs << '\n' << '\n';

    results << "# incl [deg]";
    results << '\t' << "az [deg]";
    results << '\t' << "n_b";
    results << '\t' << "lat_0 [deg]";
    results << '\t' << "lon_0 [deg]";
    results << '\t' << "time [s]";
    results << '\t' << "cel [km/s]";
    results << '\t' << "turning ht [km]";
    results << '\t' << "inclination [deg]";
    results << '\t' << "back azimuth [deg]";
    results << '\t' << "trans. coeff. [dB]";
    results << '\t' << "absorption [dB]";
    results << '\n';
    
	if(write_rays){
        sprintf(output_buffer, "%s.raypaths.dat", output_id);
        raypath.open(output_buffer);
        
        raypath << "# lat [deg]";
        raypath << '\t' << "lon [deg]";
        raypath << '\t' << "z [km]";
        raypath << '\t' << "trans. coeff. [dB]";
        raypath << '\t' << "absorption [dB]";
        raypath << '\t' << "time [s]";
        raypath << '\n';
        
    }
    if(write_caustics){
        for (int bnc = 0; bnc <= bounces; bnc++){
            sprintf(output_buffer, "%s.caustics-%i.dat", output_id, bnc);
            caustics.open(output_buffer);
            
            caustics << "# lat [deg]";
            caustics << '\t' << "lon [deg]";
            caustics << '\t' << "z [km]";
            caustics << '\t' << "time [s]";
            caustics << '\n';
            caustics.close();
        }
    }

	double** solution;
	geoac::build_solution(solution, length);
    
    for(double phi = phi_min; phi <= phi_max; phi+=phi_step){
        geoac::phi = Pi / 2.0 - phi * (Pi / 180.0);

        double dzgdx, dzgdy, theta_grnd;
        if (z_src - (topo::z(lat_src, lon_src) - globe::r0) < 0.01){
            if (geoac::is_topo){
                dzgdx = 1.0 / (globe::r0 * cos(lat_src)) * topo::dz(lat_src, lon_src, 1);
                dzgdy = 1.0 / globe::r0 * topo::dz(lat_src, lon_src, 0);
                theta_grnd = atan(dzgdx * sin(phi * Pi / 180.0) + dzgdy * cos(phi * Pi / 180.0)) * (180.0 / Pi) + 0.1;
            } else {
                theta_grnd = 0.1;
            }
        } else {
            theta_grnd=-89.9;
        }

        for(double theta = max(theta_min, theta_grnd); theta <= max(theta_max, theta_grnd); theta+=theta_step){
	    	cout << "Calculating ray path: " << theta << " degrees inclination, " << phi << " degrees azimuth." << '\n';
            geoac::theta = theta * (Pi / 180.0);
        
            geoac::set_initial(solution, z_src, lat_src, lon_src);
            travel_time_sum = 0.0;
            attenuation = 0.0;
            r_max = 0.0;

            if((fabs(theta - max(theta_min, theta_grnd)) < theta_step) && write_topo){
                sprintf(output_buffer, "%s.terrain.dat", output_id);
                topo_out.open(output_buffer);
            }

            for(int bnc_cnt = 0; bnc_cnt <= bounces; bnc_cnt++){
                k = geoac::prop_rk4(solution, break_check, length);
                if(print_resid){
                    cout << '\t' << "Residuals at Ray Termination (region boundary or reflection):" << '\n';
                    cout << '\t' << '\t' << "Eikonal residual: " << geoac::eval_eikonal(solution, k - 1) << '\n';
                    cout << '\t' << '\t' << "Aux. param. residual: " << geoac::eval_eikonal_deriv(solution, k - 1) << '\n';
                }

                if(write_rays || write_caustics){
                    if(write_caustics){
                        sprintf(output_buffer, "%s.caustics-%i.dat", output_id, bnc_cnt);
                        caustics.open(output_buffer,fstream::app);
                    
                        D_prev = geoac::jacobian(solution,1);
                    }

                    for(int m = 1; m < k ; m++){
                        if(write_caustics){
                            D = geoac::jacobian(solution, m);
                        }
                        geoac::travel_time(travel_time_sum, solution, m - 1, m);
                        geoac::atten(attenuation, solution, m - 1, m, freq);
                        
                        if(write_rays && (m == 1 || m % 15 == 0)){
                            raypath << setprecision(8) << solution[m][1] * (180.0 / Pi);
                            raypath << '\t' << setprecision(8) << solution[m][2] * (180.0 / Pi);
                            raypath << '\t' << solution[m][0] - globe::r0;
                            if(geoac::calc_amp){
                                raypath << '\t' << 20.0 * log10(geoac::amp(solution, m));
                            } else{
                                raypath << '\t' << 0.0;
                            }
                            raypath << '\t' << -2.0 * attenuation;
                            raypath << '\t' << travel_time_sum << '\n';
                        }
                        
                        if(write_caustics && D*D_prev < 0.0){
                            caustics << setprecision(8) << solution[m][1] * (180.0 / Pi);
                            caustics << '\t' << setprecision(8) << solution[m][2] * (180.0 / Pi);
                            caustics << '\t' << solution[m][0] - globe::r0;
                            caustics << '\t' << travel_time_sum << '\n';
                        }
                        if(write_caustics){
                            D_prev = D;
                        }
                    }
                    if(write_caustics){
                        caustics.close();
                    }
                } else {
                    travel_time_sum+= geoac::travel_time(solution, k);
                    attenuation+= geoac::atten(solution,k,freq);
                }

                if((fabs(theta - max(theta_min, theta_grnd)) < theta_step) && write_topo){
                    for(int m = 1; m < k ; m+=10){                        
                        topo_out << setprecision(8) << solution[m][1] * (180.0 / Pi);
                        topo_out << '\t' << setprecision(8) << solution[m][2] * (180.0 / Pi);
                        topo_out << '\t' << topo::z(solution[m][1], solution[m][2]) - globe::r0 << '\n';
                    }
                }

                for(int m = 0; m < k ; m++){
                    r_max = max (r_max, solution[m][0] - globe::r0);
                }

                
                if(!geoac::is_topo && r_max < (topo::z0 + turn_ht_min)){
                    break;
                }

                if(break_check || k < 2 || (travel_time_sum / 3600.0) > geoac::time_max){
                    break;
                }
            
                inclination = - asin(atmo::c(solution[k][0], solution[k][1], solution[k][2]) / atmo::c(globe::r0 + z_src, lat_src, lon_src) * solution[k][3]) * (180.0 / Pi);
                back_az = 90.0 - atan2(-solution[k][4], -solution[k][5]) * (180.0 / Pi);
                if(back_az < -180.0){
                    back_az +=360.0;
                } else if(back_az >  180.0){
                    back_az -=360.0;
                }
                
                results << theta;
                results << '\t' << phi;
                results << '\t' << bnc_cnt;
                results << '\t' << setprecision(8) << solution[k][1] * (180.0 / Pi);
                results << '\t' << setprecision(8) << solution[k][2] * (180.0 / Pi);
                results << '\t' << travel_time_sum;
                results << '\t' << globe::gc_dist(solution[k][1], solution[k][2], lat_src, lon_src) / travel_time_sum;
                results << '\t' << r_max;
                results << '\t' << inclination;
                results << '\t' << back_az;
                if(geoac::calc_amp){    results << '\t' << 20.0 * log10(geoac::amp(solution,k));}
                else{                   results << '\t' << 0.0;}
                results << '\t' << -2.0 * attenuation;
                results << '\n';

            
                geoac::set_refl(solution,k);
            }
            if(write_rays){raypath << '\n';}
            geoac::clear_solution(solution,k);

            if((fabs(theta - max(theta_min, theta_grnd)) < theta_step) && write_topo){
                topo_out.close();
            }

        }
    results << '\n';
    }
	
    raypath.close();
    results.close();
    geoac::delete_solution(solution, length);
    clear_region();
}


void run_back_proj(char* inputs[], int count){
    cout << '\n';
    cout << '\t' << "##########################################" << '\n';
    cout << '\t' << "####        Running infraga-sph       ####" << '\n';
    cout << '\t' << "####          Back Projection         ####" << '\n';
    cout << '\t' << "##########################################" << '\n' << '\n';
    
    double lat_rcvr = 30.0, lon_rcvr = 0.0, z_rcvr = 0.0;
    double freq = 0.1, D, D_prev;
    int bounces = 0, file_check;
    bool write_atmo = false, custom_output_id=false;
    char* prof_format = "zTuvdp";
    char* topo_file = "None";
    char* output_id;
    char input_check;
    
    topo::z0 = 0.0;
    atmo::tweak_abs = 1.0;
    
    geoac::calc_amp = true;
    geoac::is_topo = false;
    
    geoac::theta = 15.0 * (Pi / 180.0);
    geoac::phi = Pi / 2.0 - 90.0 * (Pi / 180.0);
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "prof_format=", 12) == 0){       prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){        topo::z0 = atof(inputs[i] + 7);}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){    topo_file = inputs[i] + 10; geoac::is_topo=true;}
    }
    
    if (geoac::is_topo){    file_check = set_region(inputs[2], topo_file, prof_format, true);}
    else{                   file_check = set_region(inputs[2], prof_format, true);}
    if(!file_check){
        return;
    }

    if (geoac::is_topo){
        lat_rcvr = (geoac::lat_max + geoac::lat_min) / 2.0 * (180.0 / Pi),
        lon_rcvr = (geoac::lon_max + geoac::lon_min) / 2.0 * (180.0 / Pi),
        z_rcvr = 0.0;
    }
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "inclination=", 12) == 0){                                                   geoac::theta = atof(inputs[i] + 12) * (Pi / 180.0);}
        else if (strncmp(inputs[i], "azimuth=", 8) == 0){                                                   geoac::phi = Pi / 2.0 - atof(inputs[i] + 8) * (Pi / 180.0);}
        else if (strncmp(inputs[i], "bounces=", 8) == 0){                                                   bounces = max(0, atoi(inputs[i] + 8));}
        
        else if ((strncmp(inputs[i], "rcvr_lat=", 9) == 0) || (strncmp(inputs[i], "lat_rcvr=", 9) == 0)){   lat_rcvr = atof(inputs[i] + 9);}
        else if ((strncmp(inputs[i], "rcvr_lon=", 9) == 0) || (strncmp(inputs[i], "lon_rcvr=", 9) == 0)){   lon_rcvr = atof(inputs[i] + 9);}
        
        else if (strncmp(inputs[i], "freq=", 5) == 0){                                                      freq = atof(inputs[i] + 5);}
        else if (strncmp(inputs[i], "abs_coeff=", 10) == 0){                                                atmo::tweak_abs = max(0.0, atof(inputs[i]+10));}
        
        else if ((strncmp(inputs[i], "max_alt=", 8) == 0) || (strncmp(inputs[i], "alt_max=", 8) == 0)){     geoac::alt_max = atof(inputs[i] + 8);   cout << '\t' << "User selected altitude maximum = " << geoac::alt_max << '\n';}
        else if ((strncmp(inputs[i], "max_rng=", 8) == 0) || (strncmp(inputs[i], "rng_max=", 8) == 0)){     geoac::rng_max = atof(inputs[i] + 8);   cout << '\t' << "User selected range maximum = " << geoac::rng_max << '\n';}
        else if ((strncmp(inputs[i], "max_time=", 9) == 0) || (strncmp(inputs[i], "time_max=", 9) == 0)){   geoac::time_max = atof(inputs[i] + 9);   cout << '\t' << "User selected propagation time maximum = " << geoac::time_max << '\n';}

        else if ((strncmp(inputs[i], "min_lat=", 8) == 0) || (strncmp(inputs[i], "lat_min=", 8) == 0)){     geoac::lat_min = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected latitude minimum = " << geoac::lat_min * (180.0 / Pi) << '\n';}
        else if ((strncmp(inputs[i], "max_lat=", 8) == 0) || (strncmp(inputs[i], "lat_max=", 8) == 0)){     geoac::lat_max = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected latitude maximum = " << geoac::lat_max * (180.0 / Pi) << '\n';}
        else if ((strncmp(inputs[i], "min_lon=", 8) == 0) || (strncmp(inputs[i], "lon_min=", 8) == 0)){     geoac::lon_min = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected longitude minimum = " << geoac::lon_min * (180.0 / Pi) << '\n';}
        else if ((strncmp(inputs[i], "max_lon=", 8) == 0) || (strncmp(inputs[i], "lon_max=", 8) == 0)){     geoac::lon_max = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected longitude maximum = " << geoac::lon_max * (180.0 / Pi) << '\n';}
        
        else if ((strncmp(inputs[i], "min_ds=", 7) == 0) || (strncmp(inputs[i], "ds_min=", 7) == 0)){       geoac::ds_min = atof(inputs[i] + 7);    cout << '\t' << "User selected ds minimum = " << geoac::ds_min << '\n';}
        else if ((strncmp(inputs[i], "max_ds=", 7) == 0) || (strncmp(inputs[i], "ds_max=", 7) == 0)){       geoac::ds_max = atof(inputs[i] + 7);    cout << '\t' << "User selected ds maximum = " << geoac::ds_max << '\n';}
        else if ((strncmp(inputs[i], "max_s=", 6) == 0) || (strncmp(inputs[i], "s_max=", 6) == 0)){         geoac::s_max = atof(inputs[i] + 6);     cout << '\t' << "User selected s maximum (path length between reflections) = " << geoac::s_max << '\n';}
        
        else if (strncmp(inputs[i], "write_atmo=", 11) == 0){                                               write_atmo = string2bool(inputs[i] + 11);}
        else if (strncmp(inputs[i], "prof_format=", 12) == 0){                                              prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "output_id=", 10) == 0){                                                custom_output_id = true; 
                                                                                                            output_id = inputs[i] + 10;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){                                                    if(!geoac::is_topo){
                                                                                                                topo::z0 = atof(inputs[i] + 7);
                                                                                                                topo::z_max = topo::z0;
                                                                                                                topo::z_bndlyr = topo::z0 + 2.0;
                                                                                                            } else {
                                                                                                                cout << '\t' << "Note: cannot adjust ground elevation with topography." << '\n';
                                                                                                            }}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){                                                topo_file = inputs[i] + 10; geoac::is_topo=true;}
        else if (strncmp(inputs[i], "topo_use_BLw=", 13) == 0){                                             topo::use_BLw = string2bool(inputs[i] + 13);}       
        else{
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }
    lat_rcvr *= Pi / 180.0;
    lon_rcvr *= Pi / 180.0;
    z_rcvr = max(z_rcvr, topo::z(lat_rcvr, lon_rcvr) - globe::r0);
    if(write_atmo)      geoac::write_prof("atmo.dat", lat_rcvr, lon_rcvr, geoac::phi);
    geoac::configure();
    
    cout << '\n' << "Parameter summary:" << '\n';
    cout << '\t' << "inclination: " << geoac::theta * (180.0 / Pi) << '\n';
    cout << '\t' << "azimuth: " << 90.0 - geoac::phi * (180.0 / Pi) << '\n';
    cout << '\t' << "bounces: " << bounces << '\n';
    cout << '\t' << "receiver location (lat, lon, alt): " << lat_rcvr * (180.0 / Pi) << ", " << lon_rcvr * (180.0 / Pi) << ", " << z_rcvr << '\n';
    if(!geoac::is_topo){
        cout << '\t' << "ground elevation: " << topo::z0 << '\n';
    }
    cout << '\t' << "frequency: " << freq << '\n';
    cout << '\t' << "S&B atten coeff: " << atmo::tweak_abs << '\n' << '\n';
    
    // Define variables used for analysis
    double travel_time_sum, travel_time_var_incl, travel_time_var_az, attenuation, r_max;
    double back_az, inclination;
    
    int k, length = int(geoac::s_max / geoac::ds_min);
    bool break_check;
    ofstream projection;
    
    // Extract the file name from the input and use it to distinguish the resulting output
    char output_buffer [512];
    if(!custom_output_id){
        output_id = inputs[2];
        for(int m = strlen(output_id); m >= 0; m--){
            if(output_id[m] == '.'){
                output_id[m] = '\0';
                break;
            }
        }
    }
    
    sprintf(output_buffer, "%s.projection.dat", output_id);
    projection.open(output_buffer);
    projection << "# lat [deg]";
    projection << '\t' << "lon [deg]";
    projection << '\t' << "z [km]";
    projection << '\t' << "time [s]";
    projection << '\t' << "LAT^{(incl)} [-]";
    projection << '\t' << "LON^{(incl)} [-]";
    projection << '\t' << "Z^{(incl)} [km/deg]";
    projection << '\t' << "T^{(incl)} [s/deg]";
    projection << '\t' << "LAT^{(az)} [-]";
    projection << '\t' << "LON^{(az)} [-]";
    projection << '\t' << "Z^{(az)} [km/deg]";
    projection << '\t' << "T^{(az)} [s/deg]";
    projection << '\n';
    
    double** solution;
    geoac::build_solution(solution,length);
    geoac::set_initial(solution, z_rcvr, lat_rcvr, lon_rcvr);
    
    travel_time_sum = 0.0;
    attenuation = 0.0;
    r_max = 0.0;
    
    for(int bnc_cnt = 0; bnc_cnt <= bounces; bnc_cnt++){
        k = geoac::prop_rk4(solution, break_check, length);
        
        for(int m = 1; m < k ; m++){     // write profiles to data files and vector arrays
            geoac::travel_time_var(travel_time_sum, travel_time_var_incl, travel_time_var_az, solution, m - 1, m);
            geoac::atten(attenuation, solution, m - 1, m, freq);
            r_max = max(r_max, solution[m][0] - globe::r0);
            
            if(m == 1 || m % 5 == 0){
                projection << setprecision(8) << solution[m][1] * (180.0 / Pi);
                projection << '\t' << setprecision(8) << solution[m][2] * (180.0 / Pi);
                projection << '\t' << solution[m][0] - globe::r0;
                projection << '\t' << travel_time_sum;
                
                projection << '\t' << setprecision(8) << solution[m][7];
                projection << '\t' << setprecision(8) << solution[m][8];
                projection << '\t' << solution[m][6]*  (Pi / 180.0);
                projection << '\t' << travel_time_var_incl * (Pi / 180.0);

                projection << '\t' << setprecision(8) << solution[m][13];
                projection << '\t' << setprecision(8) << solution[m][14];
                projection << '\t' << solution[m][12] * (Pi / 180.0);
                projection << '\t' << travel_time_var_az  * (Pi / 180.0);

                projection << '\n';
            }
        }
        if(break_check) break;
        geoac::set_refl(solution,k);
    }
    
    inclination = - asin(atmo::c(solution[k][0], solution[k][1], solution[k][2]) / atmo::c(globe::r0 + z_rcvr, lat_rcvr, lon_rcvr) * solution[k][3]) * (180.0 / Pi);
    back_az = 90.0 - atan2(-solution[k][4], -solution[k][5]) * (180.0 / Pi);
    while(back_az < -180.0) back_az += 360.0;
    while(back_az >  180.0) back_az -= 360.0;
    
    if(!break_check){
        cout << '\n';
        cout << '\t' << "Back projected arrival:" << '\n';
        cout << '\t' << '\t' << "latitude [deg] = " << setprecision(8) << solution[k][1] * (180.0 / Pi) << '\n';
        cout << '\t' << '\t' << "longitude [deg] = " << setprecision(8) << solution[k][2] * (180.0 / Pi) << '\n';
        cout << '\t' << '\t' << "time [s] = " << travel_time_sum << '\n';
        cout << '\t' << '\t' << "celerity [km/s] = " << globe::gc_dist(solution[k][1], solution[k][2], lat_rcvr, lon_rcvr) / travel_time_sum << '\n';
        cout << '\t' << '\t' << "turning height [km] = " << r_max << '\n';
        cout << '\t' << '\t' << "arrival inclination [deg] = " << inclination << '\n';
        cout << '\t' << '\t' << "back azimuth [deg] = " << back_az << '\n';
        cout << '\t' << '\t' << "trans. coeff. [dB] = " << 20.0 * log10(geoac::amp(solution,k)) << '\n';
        cout << '\t' << '\t' << "absorption [dB] = " << -2.0 * attenuation << '\n' << '\n';
    } else {
        cout << '\n';
        cout << '\t' << '\t' << "Ray path does not return to the ground." << '\n' << '\n';
    }
    
    projection.close();
    
    geoac::delete_solution(solution, length);
    clear_region();
}


void run_eig_search(char* inputs[], int count){
    cout << '\n';
    cout << '\t' << "###########################################" << '\n';
    cout << '\t' << "####        Running infraga-sph        ####" << '\n';
    cout << '\t' << "####          Eigenray Search          ####" << '\n';
    cout << '\t' << "###########################################" << '\n' << '\n';

    double src [3]   = {0.0, 30.0, 0.0};
    double rcvr [2] = {30.0, -2.5};
    double theta_min = 0.5, theta_max = 45.0;
    int bnc_min = 0, bnc_max = 0;
    int iterations=25, file_check;
    double az_err_lim=2.0;
    double freq = 0.1;
    char* prof_format = "zTuvdp";
    char* topo_file = "None";
    char* output_id;
    char input_check;
    bool write_atmo = false, custom_output_id=false, reverse_winds=false;
    
    topo::z0 = 0.0;
    atmo::tweak_abs = 1.0;

    geoac::is_topo = false;
    geoac::verbose = false;
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "prof_format=", 12) == 0){       prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){        topo::z0 = atof(inputs[i] + 7);}
        else if (strncmp(inputs[i], "reverse_winds=", 14) == 0){    reverse_winds = string2bool(inputs[i] + 14);}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){    topo_file = inputs[i] + 10; geoac::is_topo=true;}
    }
    
    if (geoac::is_topo){    file_check = set_region(inputs[2], topo_file, prof_format, reverse_winds);}
    else{                   file_check = set_region(inputs[2], prof_format, reverse_winds);}
    if(!file_check){
        return;
    }

    if (geoac::is_topo){
        src[1] = (geoac::lat_max + geoac::lat_min) / 2.0 * (180.0 / Pi);    rcvr[0] = src[1] + 0.0;
        src[2] = (geoac::lon_max + geoac::lon_min) / 2.0 * (180.0 / Pi);    rcvr[1] = src[2] + 2.5;
    }
    
    for(int i = 3; i < count; i++){
        if ((strncmp(inputs[i], "incl_min=", 9) == 0) || (strncmp(inputs[i], "min_incl=", 9) == 0)){                    theta_min = atof(inputs[i] + 9);}
        else if ((strncmp(inputs[i], "incl_max=", 9) == 0) || (strncmp(inputs[i], "max_incl=", 9) == 0)){               theta_max = atof(inputs[i] + 9);}
        else if ((strncmp(inputs[i], "bnc_min=", 8) == 0) || (strncmp(inputs[i], "bnc_min=", 8) == 0)){                 bnc_min = atoi(inputs[i] + 8);}
        else if ((strncmp(inputs[i], "bnc_max=", 8) == 0) || (strncmp(inputs[i], "bnc_max=", 8) == 0)){                 bnc_max = atoi(inputs[i] + 8);}
        else if (strncmp(inputs[i], "bounces=", 8) == 0){                                                               bnc_min = atoi(inputs[i] + 8);
                                                                                                                        bnc_max = atoi(inputs[i] + 8);}
        
        else if ((strncmp(inputs[i], "src_alt=", 8) == 0) || (strncmp(inputs[i], "alt_src=", 8) == 0)){                 src[0] = atof(inputs[i] + 8);}
        else if ((strncmp(inputs[i], "src_lat=", 8) == 0) || (strncmp(inputs[i], "lat_src=", 8) == 0)){                 src[1] = atof(inputs[i] + 8);}
        else if ((strncmp(inputs[i], "src_lon=", 8) == 0) || (strncmp(inputs[i], "lon_src=", 8) == 0)){                 src[2] = atof(inputs[i] + 8);}
        else if ((strncmp(inputs[i], "rcvr_lat=", 9) == 0) || (strncmp(inputs[i], "lat_rcvr=", 9) == 0)){               rcvr[0] = atof(inputs[i] + 9);}
        else if ((strncmp(inputs[i], "rcvr_lon=", 9) == 0) || (strncmp(inputs[i], "lon_rcvr=", 9) == 0)){               rcvr[1] = atof(inputs[i] + 9);}
        
        else if (strncmp(inputs[i], "verbose=", 8) == 0){                                                               geoac::verbose = string2bool(inputs[i] + 8);}
        else if (strncmp(inputs[i], "iterations=", 11) == 0){                                                           iterations=atof(inputs[i] + 11);}
        else if (strncmp(inputs[i], "damping=", 8) == 0){                                                               geoac::damping = atof(inputs[i] + 8);}
        else if (strncmp(inputs[i], "tolerance=", 10) == 0){                                                            geoac::tolerance = atof(inputs[i] + 10);}
        else if (strncmp(inputs[i], "az_dev_lim=", 11) == 0){                                                           az_err_lim = atof(inputs[i] + 11);}
        else if ((strncmp(inputs[i], "incl_step_min=", 14) == 0) || (strncmp(inputs[i], "min_incl_step=", 14) == 0)){   geoac::dth_sml = atof(inputs[i] + 14) * (Pi / 180.0);}
        else if ((strncmp(inputs[i], "incl_step_max=", 14) == 0) || (strncmp(inputs[i], "max_incl_step=", 14) == 0)){   geoac::dth_big = atof(inputs[i] + 14) * (Pi / 180.0);}

        else if (strncmp(inputs[i], "freq=", 5) == 0){                                                                  freq = atof(inputs[i] + 5);}
        else if (strncmp(inputs[i], "abs_coeff=", 10) == 0){                                                            atmo::tweak_abs = max(0.0, atof(inputs[i] + 10));}
        
        else if ((strncmp(inputs[i], "max_alt=", 8) == 0) || (strncmp(inputs[i], "alt_max=", 8) == 0)){                 geoac::alt_max = atof(inputs[i] + 8);   cout << '\t' << "User selected altitude maximum = " << geoac::alt_max << '\n';}
        else if ((strncmp(inputs[i], "max_rng=", 8) == 0) || (strncmp(inputs[i], "rng_max=", 8) == 0)){                 geoac::rng_max = atof(inputs[i] + 8);   cout << '\t' << "User selected range maximum = " << geoac::rng_max << '\n';}
        else if ((strncmp(inputs[i], "max_time=", 9) == 0) || (strncmp(inputs[i], "time_max=", 9) == 0)){               geoac::time_max = atof(inputs[i] + 9);   cout << '\t' << "User selected propagation time maximum = " << geoac::time_max << '\n';}

        else if ((strncmp(inputs[i], "min_lat=", 8) == 0) || (strncmp(inputs[i], "lat_min=", 8) == 0)){                 geoac::lat_min = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected latitude minimum = " << geoac::lat_min * (180.0 / Pi) << '\n';}
        else if ((strncmp(inputs[i], "max_lat=", 8) == 0) || (strncmp(inputs[i], "lat_max=", 8) == 0)){                 geoac::lat_max = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected latitude maximum = " << geoac::lat_max * (180.0 / Pi) << '\n';}
        else if ((strncmp(inputs[i], "min_lon=", 8) == 0) || (strncmp(inputs[i], "lon_min=", 8) == 0)){                 geoac::lon_min = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected longitude minimum = " << geoac::lon_min * (180.0 / Pi) << '\n';}
        else if ((strncmp(inputs[i], "max_lon=", 8) == 0) || (strncmp(inputs[i], "lon_max=", 8) == 0)){                 geoac::lon_max = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected longitude maximum = " << geoac::lon_max * (180.0 / Pi) << '\n';}
        
        else if ((strncmp(inputs[i], "min_ds=", 7) == 0) || (strncmp(inputs[i], "ds_min=", 7) == 0)){                   geoac::ds_min = atof(inputs[i] + 7);    cout << '\t' << "User selected ds minimum = " << geoac::ds_min << '\n';}
        else if ((strncmp(inputs[i], "max_ds=", 7) == 0) || (strncmp(inputs[i], "ds_max=", 7) == 0)){                   geoac::ds_max = atof(inputs[i] + 7);    cout << '\t' << "User selected ds maximum = " << geoac::ds_max << '\n';}
        else if ((strncmp(inputs[i], "max_s=", 6) == 0) || (strncmp(inputs[i], "s_max=", 6) == 0)){                     geoac::s_max = atof(inputs[i] + 6);     cout << '\t' << "User selected s maximum (path length between reflections) = " << geoac::s_max << '\n';}

        else if (strncmp(inputs[i], "write_atmo=", 11) == 0){                                                           write_atmo = string2bool(inputs[i] + 11);}
        else if (strncmp(inputs[i], "prof_format=",12) == 0){                                                           prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "reverse_winds=", 14) == 0){                                                        reverse_winds = string2bool(inputs[i] + 14);}

        else if (strncmp(inputs[i], "output_id=", 10) == 0){                                                            custom_output_id = true; 
                                                                                                                        output_id = inputs[i] + 10;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){                                                                if(!geoac::is_topo){
                                                                                                                            topo::z0 = atof(inputs[i] + 7);
                                                                                                                            topo::z_max = topo::z0;
                                                                                                                            topo::z_bndlyr = topo::z0 + 2.0;
                                                                                                                        } else {
                                                                                                                            cout << '\t' << "Note: cannot adjust ground elevation with topography." << '\n';
                                                                                                                        }}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){                                                            topo_file = inputs[i] + 10; geoac::is_topo=true;}       
        else if (strncmp(inputs[i], "topo_use_BLw=", 13) == 0){                                                         topo::use_BLw = string2bool(inputs[i] + 13);}
        else{
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }
    src[1] *= (Pi / 180.0);   rcvr[0] *= (Pi / 180.0);
    src[2] *= (Pi / 180.0);   rcvr[1] *= (Pi / 180.0);
    src[0] = max(src[0], topo::z(src[1], src[2]) - globe::r0);
    if(write_atmo)      geoac::write_prof("atmo.dat", src[1], src[2], globe::bearing(src[1], src[2], rcvr[0], rcvr[1]));

    theta_min *= (Pi / 180.0);
    theta_max *= (Pi / 180.0);

    cout << '\n' << "Parameter summary:" << '\n';
    cout << '\t' << "source location (lat, lon, alt): " << src[1] * (180.0 / Pi) << ", " << src[2] * (180.0 / Pi) << ", " << src[0] << '\n';
    cout << '\t' << "receiver location (lat, lon, alt): " << rcvr[0] * (180.0 / Pi) << ", " << rcvr[1] * (180.0 / Pi) << ", " << topo::z(rcvr[0], rcvr[1]) - globe::r0 <<'\n';
    cout << '\t' << "inclination range: " << theta_min * (180.0 / Pi) << ", " << theta_max * (180.0 / Pi) << '\n';
    cout << '\t' << "bounces: " << bnc_min << ", " << bnc_max << '\n';
    if(!geoac::is_topo){
        cout << '\t' << "ground elevation: " << topo::z0 << '\n';
    }
    cout << '\t' << "damping:" << geoac::damping << '\n';
    cout << '\t' << "frequency: " << freq << '\n';
    cout << '\t' << "S&B atten coeff: " << atmo::tweak_abs << '\n';
    if (geoac::verbose){
        cout << '\t' << "verbose: true" << '\n';
    } else {
        cout << '\t' << "verbose: false" << '\n';
    }
    cout << '\n';

    // Extract the file name from the input and use it to distinguish the resulting output
    char output_buffer [512];
    if(!custom_output_id){
        output_id = inputs[2];
        for(int m = strlen(output_id); m >= 0; m--){
            if(output_id[m] == '.'){
                output_id[m] = '\0';
                break;
            }
        }
    }
    
    double theta_start, theta_next, theta_est, phi_est;
    bool estimate_success;
    
    sprintf(output_buffer, "%s.arrivals.dat", output_id);
    geoac::eig_results.open(output_buffer);
    
    geoac::eig_results << "# infraga-sph -eig_search summary:" << '\n';
    geoac::eig_results << "# " << '\t' << "profile: " << inputs[2] << '\n';
    geoac::eig_results << "# " << '\t' << "source location (lat, lon, alt): " << src[1] * (180.0 / Pi) << ", " << src[2] * (180.0 / Pi) << ", " << src[0] << '\n';
    geoac::eig_results << "# " << '\t' << "receiver location (lat, lon, alt): " << rcvr[0] * (180.0 / Pi) << ", " << rcvr[1] * (180.0 / Pi) << ", " << topo::z(rcvr[0], rcvr[1]) - globe::r0 <<'\n';
    geoac::eig_results << "# " << '\t' << "inclination range: " << theta_min * (180.0 / Pi) << ", " << theta_max * (180.0 / Pi) << '\n';
    geoac::eig_results << "# " << '\t' << "bounces: " << bnc_min << ", " << bnc_max << '\n';
    if(!geoac::is_topo){
        geoac::eig_results << "#" << '\t' << "ground elevation: " << topo::z0 << '\n';
    } else {
        geoac::eig_results << "#" << '\t' << "topo file:" << topo_file << '\n';
    }
    geoac::eig_results << "# " << '\t' << "damping: " << geoac::damping << '\n';
    geoac::eig_results << "# " << '\t' << "frequency: " << freq << '\n';
    geoac::eig_results << "# " << '\t' << "S&B atten coeff: " << atmo::tweak_abs << '\n';
    
    geoac::eig_results << "# incl [deg]";
    geoac::eig_results << '\t' << "az [deg]";
    geoac::eig_results << '\t' << "n_b";
    geoac::eig_results << '\t' << "lat_0 [deg]";
    geoac::eig_results << '\t' << "lon_0 [deg]";
    geoac::eig_results << '\t' << "time [s]";
    geoac::eig_results << '\t' << "cel [km/s]";
    geoac::eig_results << '\t' << "turning ht [km]";
    geoac::eig_results << '\t' << "inclination [deg]";
    geoac::eig_results << '\t' << "back azimuth [deg]";
    geoac::eig_results << '\t' << "trans. coeff. [dB]";
    geoac::eig_results << '\t' << "absorption [dB]";
    geoac::eig_results << '\n';
    
    for(int n_bnc = bnc_min; n_bnc <= bnc_max; n_bnc++){
        cout << "Searching for " << n_bnc << " bounce eigenrays." << '\n';
        
        theta_start = theta_min;
        while(theta_start < theta_max){
            estimate_success = geoac::est_eigenray(src, rcvr, theta_start, theta_max, theta_est, phi_est, theta_next, n_bnc, az_err_lim);
            if(estimate_success){
                geoac::find_eigenray(src, rcvr, theta_est, phi_est, freq, n_bnc, iterations, output_id);
            }
            theta_start = theta_next;
        }
    }

    geoac::eig_results.close();
    cout << "Identified " << geoac::eigenray_cnt << " eigenray(s)." << '\n';
    clear_region();
}

void run_eig_direct(char* inputs[], int count){
    cout << '\n';
    cout << '\t' << "###########################################" << '\n';
    cout << '\t' << "####        Running infraga-sph        ####" << '\n';
    cout << '\t' << "####          Eigenray Direct          ####" << '\n';
    cout << '\t' << "###########################################" << '\n' << '\n';
 
    double src [3]   = {0.0, 30.0, 0.0};
    double rcvr [2] = {30.0, -2.5};
    double theta_est = 10.0 * (Pi / 180.0), phi_est = 0.0;
    int bounces = 0, iterations = 25, file_check;
    double freq = 0.1;
    char* prof_format = "zTuvdp";
    char* topo_file = "None";
    char* output_id;
    char input_check;
    bool write_atmo = false, custom_output_id=false, reverse_winds=false;
    
    topo::z0 = 0.0;
    atmo::tweak_abs = 1.0;
    
    geoac::calc_amp = true;
    geoac::is_topo = false;
    geoac::verbose = false;
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "prof_format=", 12) == 0){           prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){            topo::z0 = atof(inputs[i] + 7);}
        else if (strncmp(inputs[i], "reverse_winds=", 14) == 0){    reverse_winds = string2bool(inputs[i] + 14);}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){        topo_file = inputs[i] + 10; geoac::is_topo=true;}
    }
    
    if (geoac::is_topo){    file_check = set_region(inputs[2], topo_file, prof_format, reverse_winds);}
    else{                   file_check = set_region(inputs[2], prof_format, reverse_winds);}
    if(!file_check){
        return;
    }
    
    if (geoac::is_topo){
        src[1] = (geoac::lat_max + geoac::lat_min) / 2.0 * (180.0 / Pi);    rcvr[0] = src[1] + 0.0;
        src[2] = (geoac::lon_max + geoac::lon_min) / 2.0 * (180.0 / Pi);    rcvr[1] = src[2] + 2.5;
    }
    
    for(int i = 3; i < count; i++){
        if ((strncmp(inputs[i], "incl_est=", 9) == 0) || (strncmp(inputs[i], "est_incl=", 9) == 0)){        theta_est = atof(inputs[i] + 9) * (Pi / 180.0);}
        else if ((strncmp(inputs[i], "az_est=", 7) == 0) || (strncmp(inputs[i], "est_az=", 7) == 0)){       phi_est = Pi / 2.0 - atof(inputs[i] + 7)  * (Pi / 180.0);}
        else if (strncmp(inputs[i], "bounces=", 8) == 0){                                                   bounces = atoi(inputs[i] + 8);}
        
        else if ((strncmp(inputs[i], "src_alt=", 8) == 0) || (strncmp(inputs[i], "alt_src=", 8) == 0)){     src[0] = atof(inputs[i] + 8);}
        else if ((strncmp(inputs[i], "src_lat=", 8) == 0) || (strncmp(inputs[i], "lat_src=", 8) == 0)){     src[1] = atof(inputs[i] + 8);}
        else if ((strncmp(inputs[i], "src_lon=", 8) == 0) || (strncmp(inputs[i], "lon_src=", 8) == 0)){     src[2] = atof(inputs[i] + 8);}
        else if ((strncmp(inputs[i], "rcvr_lat=", 9) == 0) || (strncmp(inputs[i], "lat_rcvr=", 9) == 0)){   rcvr[0] = atof(inputs[i] + 9);}
        else if ((strncmp(inputs[i], "rcvr_lon=", 9) == 0) || (strncmp(inputs[i], "lon_rcvr=", 9) == 0)){   rcvr[1] = atof(inputs[i] + 9);}
        
        else if (strncmp(inputs[i], "verbose=", 8) == 0){                                                   geoac::verbose = string2bool(inputs[i]+8);}
        else if (strncmp(inputs[i], "iterations=", 11) == 0){                                               iterations = atof(inputs[i]+11);}
        else if (strncmp(inputs[i], "damping=", 8) == 0){                                                   geoac::damping = atof(inputs[i] + 8);}
        else if (strncmp(inputs[i], "tolerance=", 10) == 0){                                                geoac::tolerance = atof(inputs[i] + 10);}
        
        else if (strncmp(inputs[i], "freq=", 5) == 0){                                                      freq = atof(inputs[i] + 5);}
        else if (strncmp(inputs[i], "abs_coeff=", 10) == 0){                                                atmo::tweak_abs = max(0.0, atof(inputs[i] + 10));}
        
        else if ((strncmp(inputs[i], "max_alt=", 8) == 0) || (strncmp(inputs[i], "alt_max=", 8) == 0)){     geoac::alt_max = atof(inputs[i] + 8);   cout << '\t' << "User selected altitude maximum = " << geoac::alt_max << '\n';}
        else if ((strncmp(inputs[i], "max_rng=", 8) == 0) || (strncmp(inputs[i], "rng_max=", 8) == 0)){     geoac::rng_max = atof(inputs[i] + 8);   cout << '\t' << "User selected range maximum = " << geoac::rng_max << '\n';}
        else if ((strncmp(inputs[i], "max_time=", 9) == 0) || (strncmp(inputs[i], "time_max=", 9) == 0)){   geoac::time_max = atof(inputs[i] + 9);   cout << '\t' << "User selected propagation time maximum = " << geoac::time_max << '\n';}

        else if ((strncmp(inputs[i], "min_lat=", 8) == 0) || (strncmp(inputs[i], "lat_min=", 8) == 0)){     geoac::lat_min = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected latitude minimum = " << geoac::lat_min * (180.0 / Pi) << '\n';}
        else if ((strncmp(inputs[i], "max_lat=", 8) == 0) || (strncmp(inputs[i], "lat_max=", 8) == 0)){     geoac::lat_max = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected latitude maximum = " << geoac::lat_max * (180.0 / Pi) << '\n';}
        else if ((strncmp(inputs[i], "min_lon=", 8) == 0) || (strncmp(inputs[i], "lon_min=", 8) == 0)){     geoac::lon_min = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected longitude minimum = " << geoac::lon_min * (180.0 / Pi) << '\n';}
        else if ((strncmp(inputs[i], "max_lon=", 8) == 0) || (strncmp(inputs[i], "lon_max=", 8) == 0)){     geoac::lon_max = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected longitude maximum = " << geoac::lon_max * (180.0 / Pi) << '\n';}
        
        else if ((strncmp(inputs[i], "min_ds=", 7) == 0) || (strncmp(inputs[i], "ds_min=", 7) == 0)){       geoac::ds_min = atof(inputs[i] + 7);    cout << '\t' << "User selected ds minimum = " << geoac::ds_min << '\n';}
        else if ((strncmp(inputs[i], "max_ds=", 7) == 0) || (strncmp(inputs[i], "ds_max=", 7) == 0)){       geoac::ds_max = atof(inputs[i] + 7);    cout << '\t' << "User selected ds maximum = " << geoac::ds_max << '\n';}
        else if ((strncmp(inputs[i], "max_s=", 6) == 0) || (strncmp(inputs[i], "s_max=", 6) == 0)){         geoac::s_max = atof(inputs[i] + 6);     cout << '\t' << "User selected s maximum (path length between reflections) = " << geoac::s_max << '\n';}

        else if (strncmp(inputs[i], "write_atmo=", 11) == 0){                                               write_atmo = string2bool(inputs[i] + 11);}
        else if (strncmp(inputs[i], "prof_format=", 12) == 0){                                              prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "reverse_winds=", 14) == 0){                                            reverse_winds = string2bool(inputs[i] + 14);}

        else if (strncmp(inputs[i], "output_id=", 10) == 0){                                                custom_output_id = true; 
                                                                                                            output_id = inputs[i] + 10;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){                                                    if(!geoac::is_topo){
                                                                                                                topo::z0 = atof(inputs[i] + 7);
                                                                                                                topo::z_max = topo::z0;
                                                                                                                topo::z_bndlyr = topo::z0 + 2.0;
                                                                                                            } else {
                                                                                                                cout << '\t' << "Note: cannot adjust ground elevation with topography." << '\n';
                                                                                                            }}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){                                                topo_file = inputs[i] + 10; geoac::is_topo=true;}
        else if (strncmp(inputs[i], "topo_use_BLw=", 13) == 0){                                             topo::use_BLw = string2bool(inputs[i] + 13);}
        else{
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }
    src[1] *= (Pi / 180.0);   rcvr[0] *= (Pi / 180.0);
    src[2] *= (Pi / 180.0);   rcvr[1] *= (Pi / 180.0);
    src[0] = max(src[0], topo::z(src[1], src[2]) - globe::r0);
    if(write_atmo)  geoac::write_prof("atmo.dat", src[1], src[2], globe::bearing(src[1], src[2], rcvr[0], rcvr[1]));

    // Set the phi_start angle by the great circle bearing to the receiver and then change it if it's been input
    phi_est = Pi / 2.0 - globe::bearing(src[1], src[2], rcvr[0], rcvr[1]);
    for(int i = 3; i < count; i++){
        if((strncmp(inputs[i], "az_est=", 7) == 0) || (strncmp(inputs[i], "est_az=", 7) == 0)){
            phi_est = Pi / 2.0 - atof(inputs[i] + 7) * (Pi / 180.0);
        }
    }
    
    cout << '\n' << "Parameter summary:" << '\n';
    cout << '\t' << "source location (lat, lon, alt): " << src[1] * (180.0 / Pi) << ", " << src[2] * (180.0 / Pi) << ", " << src[0] << '\n';
    cout << '\t' << "receiver location (lat, lon, alt): " << rcvr[0] * (180.0 / Pi) << ", " << rcvr[1] * (180.0 / Pi) << ", " << topo::z(rcvr[0], rcvr[1]) - globe::r0 <<'\n';
    cout << '\t' << "inclination estimate: " << theta_est * (180.0 / Pi) << '\n';
    cout << '\t' << "azimuth estimate: " << 90.0 - phi_est * (180.0 / Pi) << '\n';
    cout << '\t' << "bounces: " << bounces << '\n';
    if(!geoac::is_topo){
        cout << '\t' << "ground elevation: " << topo::z0 << '\n';
    }
    cout << '\t' << "damping:" << geoac::damping << '\n';
    cout << '\t' << "frequency: " << freq << '\n';
    cout << '\t' << "S&B atten coeff: " << atmo::tweak_abs << '\n';
    if (geoac::verbose){
        cout << '\t' << "verbose: true" << '\n';
    } else {
        cout << '\t' << "verbose: false" << '\n';
    }
    cout << '\n';
    
    // Extract the file name from the input and use it to distinguish the resulting output
    if(!custom_output_id){
        output_id = inputs[2];
        for(int m = strlen(output_id); m >= 0; m--){
            if(output_id[m] == '.'){
                output_id[m] = '\0';
                break;
            }
        }
    }

    geoac::find_eigenray(src, rcvr, theta_est, phi_est, freq, bounces, iterations, output_id);
    clear_region();
}


void run_wnl_wvfrm(char* inputs[], int count){
    cout << '\n';
    cout << '\t' << "###########################################" << '\n';
    cout << '\t' << "####        Running infraga-sph        ####" << '\n';
    cout << '\t' << "####     Weakly Non-Linear Waveform    ####" << '\n';
    cout << '\t' << "###########################################" << '\n' << '\n';
    
    double lat_src = 30.0, lon_src = 0.0, z_src = 0.0;
    double freq = 0.1, D, D_prev;
    int bounces = 0, file_check;
    bool write_atmo = false, write_rays=false, custom_output_id=false, start_wvfrm=false, reverse_winds=false;
    char* prof_format = "zTuvdp";
    char* topo_file = "None";
    char* output_id;
    char input_check;

    double wvfrm_ref=1.0, wvfrm_out_step=1.0e10;
    double ray_length, c0, nu0, cg0r, cg0th, cg0ph, cg0, rho0, D0, p0;
    double** wvfrm_array;
    int wvfrm_ref_k = 0;
    char* wvfrm_file = "none";
    char* wvfrm_opt = "impulse";
    ofstream wvfrm_out;

    topo::z0 = 0.0;
    atmo::tweak_abs = 1.0;
    
    geoac::calc_amp = true;
    geoac::is_topo = false;

    geoac::theta = 15.0 * (Pi / 180.0);
    geoac::phi = Pi / 2.0 - 90.0 * (Pi / 180.0);
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "prof_format=", 12) == 0){           prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){            topo::z0 = atof(inputs[i] + 7);}
        else if (strncmp(inputs[i], "reverse_winds=", 14) == 0){    reverse_winds = string2bool(inputs[i] + 14);}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){        topo_file = inputs[i] + 10; geoac::is_topo=true;}
    }
    
    if (geoac::is_topo){    file_check = set_region(inputs[2], topo_file, prof_format, reverse_winds);}
    else{                   file_check = set_region(inputs[2], prof_format, reverse_winds);}
    if(!file_check){
        return;
    }

    if (geoac::is_topo){
        lat_src = (geoac::lat_max + geoac::lat_min) / 2.0 * (180.0 / Pi),
        lon_src = (geoac::lon_max + geoac::lon_min) / 2.0 * (180.0 / Pi),
        z_src = 0.0;
    }
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "inclination=", 12) == 0){                                                   geoac::theta = atof(inputs[i] + 12) * (Pi / 180.0);}
        else if (strncmp(inputs[i], "azimuth=", 8) == 0){                                                   geoac::phi = Pi / 2.0 - atof(inputs[i] + 8) * (Pi / 180.0);}
        else if (strncmp(inputs[i], "bounces=", 8) == 0){                                                   bounces = max(0, atoi(inputs[i] + 8));}

        else if ((strncmp(inputs[i], "src_alt=", 8) == 0) || (strncmp(inputs[i], "alt_src=", 8) == 0)){     z_src = atof(inputs[i] + 8);}
        else if ((strncmp(inputs[i], "src_lat=", 8) == 0) || (strncmp(inputs[i], "lat_src=", 8) == 0)){     lat_src = atof(inputs[i] + 8);}
        else if ((strncmp(inputs[i], "src_lon=", 8) == 0) || (strncmp(inputs[i], "lon_src=", 8) == 0)){     lon_src = atof(inputs[i] + 8);}
        else if (strncmp(inputs[i], "write_ray=", 10) == 0){                                                write_rays=string2bool(inputs[i] + 10);}

        else if (strncmp(inputs[i], "wvfrm_file=", 11) == 0){                                               wvfrm_file = inputs[i] + 11;}
        else if (strncmp(inputs[i], "wvfrm_opt=", 10) == 0){                                                wvfrm_opt = inputs[i] + 10;}
        else if (strncmp(inputs[i], "wvfrm_p0=", 9) == 0){                                                  wvfrm::p0 = atof(inputs[i] + 9);}
        else if (strncmp(inputs[i], "wvfrm_t0=", 9) == 0){                                                  wvfrm::t0 = atof(inputs[i] + 9);}
        else if (strncmp(inputs[i], "wvfrm_alpha=", 12) == 0){                                              wvfrm::alpha = atof(inputs[i] + 12);}
        else if (strncmp(inputs[i], "wvfrm_ref=", 10) == 0){                                                wvfrm_ref = atof(inputs[i] + 10);}
        else if (strncmp(inputs[i], "wvfrm_out_step=", 15) == 0){                                           wvfrm_out_step = atof(inputs[i] + 15);}
        
        else if (strncmp(inputs[i], "wvfrm_ds=", 9) == 0){                                                  geoac::ds_wvfrm = atof(inputs[i] + 9);}
        else if (strncmp(inputs[i], "wvfrm_N_rise=", 13) == 0){                                             wvfrm::N_rise = atof(inputs[i] + 15);}
        else if (strncmp(inputs[i], "wvfrm_N_period=", 15) == 0){                                           wvfrm::N_period = atof(inputs[i] + 15);}
        else if (strncmp(inputs[i], "wvfrm_len=", 10) == 0){                                                wvfrm::len = atof(inputs[i] + 10);}
        else if (strncmp(inputs[i], "wvfrm_filt1_g=", 14) == 0){                                            wvfrm::filt1_g = atof(inputs[i] + 14);}
        else if (strncmp(inputs[i], "wvfrm_filt1_n1=", 15) == 0){                                           wvfrm::filt1_n1 = atof(inputs[i] + 15);}
        else if (strncmp(inputs[i], "wvfrm_filt1_n2=", 15) == 0){                                           wvfrm::filt1_n2 = atof(inputs[i] + 15);}
        else if (strncmp(inputs[i], "wvfrm_filt2_g=", 14) == 0){                                            wvfrm::filt2_g = atof(inputs[i] + 14);}
        else if (strncmp(inputs[i], "wvfrm_filt2_n1=", 15) == 0){                                           wvfrm::filt2_n1 = atof(inputs[i] + 15);}
        else if (strncmp(inputs[i], "wvfrm_filt2_n2=", 15) == 0){                                           wvfrm::filt2_n2 = atof(inputs[i] + 15);}

        else if (strncmp(inputs[i], "freq=", 5) == 0){                                                      freq = atof(inputs[i] + 5);}
        else if (strncmp(inputs[i], "abs_coeff=", 10) == 0){                                                atmo::tweak_abs = max(0.0, atof(inputs[i]+10));}
        
        else if ((strncmp(inputs[i], "max_alt=", 8) == 0) || (strncmp(inputs[i], "alt_max=", 8) == 0)){     geoac::alt_max = atof(inputs[i] + 8);   cout << '\t' << "User selected altitude maximum = " << geoac::alt_max << '\n';}
        else if ((strncmp(inputs[i], "max_rng=", 8) == 0) || (strncmp(inputs[i], "rng_max=", 8) == 0)){     geoac::rng_max = atof(inputs[i] + 8);   cout << '\t' << "User selected range maximum = " << geoac::rng_max << '\n';}
        else if ((strncmp(inputs[i], "max_time=", 9) == 0) || (strncmp(inputs[i], "time_max=", 9) == 0)){   geoac::time_max = atof(inputs[i] + 9);   cout << '\t' << "User selected propagation time maximum = " << geoac::time_max << '\n';}

        else if ((strncmp(inputs[i], "min_lat=", 8) == 0) || (strncmp(inputs[i], "lat_min=", 8) == 0)){     geoac::lat_min = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected latitude minimum = " << geoac::lat_min * (180.0 / Pi) << '\n';}
        else if ((strncmp(inputs[i], "max_lat=", 8) == 0) || (strncmp(inputs[i], "lat_max=", 8) == 0)){     geoac::lat_max = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected latitude maximum = " << geoac::lat_max * (180.0 / Pi) << '\n';}
        else if ((strncmp(inputs[i], "min_lon=", 8) == 0) || (strncmp(inputs[i], "lon_min=", 8) == 0)){     geoac::lon_min = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected longitude minimum = " << geoac::lon_min * (180.0 / Pi) << '\n';}
        else if ((strncmp(inputs[i], "max_lon=", 8) == 0) || (strncmp(inputs[i], "lon_max=", 8) == 0)){     geoac::lon_max = atof(inputs[i] + 8) * (Pi / 180.0);   cout << '\t' << "User selected longitude maximum = " << geoac::lon_max * (180.0 / Pi) << '\n';}
        
        else if ((strncmp(inputs[i], "min_ds=", 7) == 0) || (strncmp(inputs[i], "ds_min=", 7) == 0)){       geoac::ds_min = atof(inputs[i] + 7);    cout << '\t' << "User selected ds minimum = " << geoac::ds_min << '\n';}
        else if ((strncmp(inputs[i], "max_ds=", 7) == 0) || (strncmp(inputs[i], "ds_max=", 7) == 0)){       geoac::ds_max = atof(inputs[i] + 7);    cout << '\t' << "User selected ds maximum = " << geoac::ds_max << '\n';}
        else if ((strncmp(inputs[i], "max_s=", 6) == 0) || (strncmp(inputs[i], "s_max=", 6) == 0)){         geoac::s_max = atof(inputs[i] + 6);     cout << '\t' << "User selected s maximum (path length between reflections) = " << geoac::s_max << '\n';}

        else if (strncmp(inputs[i], "write_atmo=", 11) == 0){                                               write_atmo = string2bool(inputs[i] + 11);}
        else if (strncmp(inputs[i], "prof_format=", 12) == 0){                                              prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "reverse_winds=", 14) == 0){                                            reverse_winds = string2bool(inputs[i] + 14);}

        else if (strncmp(inputs[i], "output_id=", 10) == 0){                                                custom_output_id = true; 
                                                                                                            output_id = inputs[i] + 10;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){                                                    if(!geoac::is_topo){
                                                                                                                topo::z0 = atof(inputs[i] + 7);
                                                                                                                topo::z_max = topo::z0;
                                                                                                                topo::z_bndlyr = topo::z0 + 2.0;
                                                                                                            } else {
                                                                                                                cout << '\t' << "Note: cannot adjust ground elevation with topography." << '\n';
                                                                                                            }}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){                                                topo_file = inputs[i] + 10; geoac::is_topo=true;}
        else if (strncmp(inputs[i], "topo_use_BLw=", 13) == 0){                                             topo::use_BLw = string2bool(inputs[i] + 13);}        
        else{
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }

    cout << "s_max --> " << geoac::s_max << '\n';

    lat_src *= Pi / 180.0;
    lon_src *= Pi / 180.0;
    z_src = max(z_src, topo::z(lat_src, lon_src) - globe::r0);
    if(write_atmo)  geoac::write_prof("atmo.dat", lat_src, lon_src, geoac::phi);
    geoac::configure();
    
    cout << '\n' << "Parameter summary:" << '\n';
    cout << '\t' << "inclination: " << geoac::theta * (180.0 / Pi) << '\n';
    cout << '\t' << "azimuth: " << 90.0 - geoac::phi * (180.0 / Pi) << '\n';
    cout << '\t' << "bounces: " << bounces << '\n';
    cout << '\t' << "source location (lat, lon, alt): " << lat_src * (180.0 / Pi) << ", " << lon_src * (180.0 / Pi) << ", " << z_src << '\n';
    if(!geoac::is_topo){
        cout << '\t' << "ground elevation: " << topo::z0 << '\n';
    }
    cout << '\t' << "frequency: " << freq << '\n';
    cout << '\t' << "S&B atten coeff: " << atmo::tweak_abs << '\n';
    if (strncmp(wvfrm_file, "none", 4) != 0){
        cout << '\t' << "waveform source file: " << wvfrm_file << '\n';
    } else {
        cout << '\t' << "waveform option: " << wvfrm_opt << '\n';
        cout << '\t' << "waveform p0: " << wvfrm::p0 << '\n';
        cout << '\t' << "waveform t0: " << wvfrm::t0 << '\n';
        cout << '\t' << "waveform alpha: " << wvfrm::alpha << '\n';
    }
    cout << '\t' << "waveform reference location: " << wvfrm_ref << '\n' << '\n'; 
    
    // Extract the file name from the input and use it to distinguish the resulting output
    char output_buffer [512];
    if(!custom_output_id){
        output_id = inputs[2];
        for(int m = strlen(output_id); m >= 0; m--){
            if(output_id[m] == '.'){
                output_id[m] = '\0';
                break;
            }
        }
    }

    // Define variables used for analysis
    double travel_time_sum, attenuation, r_max;
    double back_az, inclination;
    
	int k, length = int(geoac::s_max / geoac::ds_min);
    bool break_check;
    
    ofstream raypath;
    if(write_rays){
        sprintf(output_buffer, "%s.raypaths.dat", output_id);
        raypath.open(output_buffer);
        raypath << "# z [km]";
        raypath << '\t' << "lat [deg]";
        raypath << '\t' << "lon [deg]";
        raypath << '\t' << "trans. coeff. [dB]";
        raypath << '\t' << "absorption [dB]";
        raypath << '\t' << "time [s]";
        raypath << '\n';
    }

    cout << "Calculating ray path geometry and weakly non-linear waveform evolution..." << '\n';

    double** solution;
    geoac::build_solution(solution,length);
    geoac::set_initial(solution, z_src, lat_src, lon_src);

    if (strncmp(wvfrm_file, "none", 4) != 0){
        wvfrm::load_wvfrm(wvfrm_array, wvfrm_file);
    } else {
        wvfrm::build_wvfrm(wvfrm_array, wvfrm_opt);
    }

    travel_time_sum = 0.0;
    attenuation = 0.0;
    r_max = 0.0;

    ray_length = 0.0;

    for(int bnc_cnt = 0; bnc_cnt <= bounces; bnc_cnt++){
        k = geoac::prop_rk4(solution, break_check, length);            

        // write profiles to data files and vector arrays
        for(int m = 1; m < k ; m++){
            geoac::travel_time(travel_time_sum, solution, m - 1, m);
            geoac::atten(attenuation, solution, m - 1, m, freq);
            r_max = max(r_max, solution[m][0] - globe::r0);

            if(write_rays && (m == 1 || m % 15 == 0)){
                raypath << setprecision(8) << solution[m][1] * (180.0 / Pi);
                raypath << '\t' << setprecision(8) << solution[m][2] * (180.0 / Pi);
                raypath << '\t' << solution[m][0] - globe::r0;
                raypath << '\t' << 20.0 * log10(geoac::amp(solution, m));
                raypath << '\t' << -2.0 * attenuation;
                raypath << '\t' << travel_time_sum;
                raypath << '\n';
            }
        }

        if(!start_wvfrm){
            for(wvfrm_ref_k = 0; wvfrm_ref_k < k; wvfrm_ref_k++){
                double dr, dt, dp, ds, r, t;
            
                dr = solution[wvfrm_ref_k + 1][0] - solution[wvfrm_ref_k][0];   r = solution[wvfrm_ref_k][0] + dr / 2.0;
                dt = solution[wvfrm_ref_k + 1][1] - solution[wvfrm_ref_k][1];   t = solution[wvfrm_ref_k][1] + dt / 2.0;
                dp = solution[wvfrm_ref_k + 1][2] - solution[wvfrm_ref_k][2];
            
                ray_length += sqrt(pow(dr, 2) + pow(r * dt, 2) + pow(r * cos(t) * dp, 2));
                if (ray_length >= wvfrm_ref){

                    c0 = atmo::c(solution[wvfrm_ref_k + 1][0], solution[wvfrm_ref_k + 1][1], solution[wvfrm_ref_k + 1][2]);
                    rho0 = atmo::rho(solution[wvfrm_ref_k + 1][0], solution[wvfrm_ref_k + 1][1], solution[wvfrm_ref_k + 1][2]);

                    nu0 = sqrt(pow(solution[wvfrm_ref_k + 1][3], 2) + pow(solution[wvfrm_ref_k + 1][4], 2) + pow(solution[wvfrm_ref_k + 1][5], 2));
                    cg0r  = c0 * solution[wvfrm_ref_k + 1][3] / nu0 + atmo::w(solution[wvfrm_ref_k + 1][0], solution[wvfrm_ref_k + 1][1], solution[wvfrm_ref_k + 1][2]);
                    cg0th = c0 * solution[wvfrm_ref_k + 1][4] / nu0 + atmo::v(solution[wvfrm_ref_k + 1][0], solution[wvfrm_ref_k + 1][1], solution[wvfrm_ref_k + 1][2]);
                    cg0ph = c0 * solution[wvfrm_ref_k + 1][5] / nu0 + atmo::u(solution[wvfrm_ref_k + 1][0], solution[wvfrm_ref_k + 1][1], solution[wvfrm_ref_k + 1][2]);
                    cg0 = sqrt(pow(cg0r, 2) + pow(cg0th, 2) + pow(cg0ph, 2));

                    D0 = geoac::jacobian(solution, wvfrm_ref_k + 1);

                    sprintf(output_buffer, "%s.wvfrm_init.dat", output_id);
                    wvfrm_out.open(output_buffer);
                    wvfrm_out << "# t [sec]" << '\t' << "p(t) [Pa]" << '\n';
                    for (int n = 0; n < wvfrm::len; n++){
                        wvfrm_out << setprecision(8) << geoac::travel_time(solution, wvfrm_ref_k + 1) + wvfrm_array[n][0] << '\t' << wvfrm_array[n][1] << '\n';
                    }
                    wvfrm_out.close();
            
                    p0 = 0.0;
                    for (int n = 0; n < wvfrm::len; n++){    p0 = max(p0, fabs(wvfrm_array[n][1]));}
                    for (int n = 0; n < wvfrm::len; n++){    wvfrm_array[n][1] /= p0;}
                    start_wvfrm = true;

                    ray_length = geoac::wnl_wvfrm(solution, wvfrm_array, wvfrm_ref_k + 1, k, wvfrm_ref, c0, cg0, nu0, rho0, D0, p0, wvfrm_out_step);                   
                    break;
                }
            }
        } else {
            ray_length += geoac::wnl_wvfrm(solution, wvfrm_array, 0, k, ray_length, c0, cg0, nu0, rho0, D0, p0, wvfrm_out_step);
        }

        if(break_check) break;
        geoac::set_refl(solution, k);
    }
        
    inclination = - asin(atmo::c(solution[k][0], solution[k][1], solution[k][2]) / atmo::c(globe::r0 + z_src, lat_src, lon_src) * solution[k][3]) * (180.0 / Pi);
    back_az = 90.0 - atan2(-solution[k][4], -solution[k][5]) * (180.0 / Pi);
    while(back_az < -180.0) back_az += 360.0;
    while(back_az >  180.0) back_az -= 360.0;
        
    if(!break_check){
        cout << '\n';
        cout << '\t' << "Arrival summary:" << '\n';
        cout << '\t' << '\t' << "latitude [deg] = " << setprecision(8) << solution[k][1] * (180.0 / Pi) << '\n';
        cout << '\t' << '\t' << "longitude [deg] = " << setprecision(8) << solution[k][2] * (180.0 / Pi) << '\n';
        cout << '\t' << '\t' << "time [s] = " << travel_time_sum << '\n';
        cout << '\t' << '\t' << "celerity [km/s] = " << globe::gc_dist(solution[k][1], solution[k][2], lat_src, lon_src) / travel_time_sum << '\n';
        cout << '\t' << '\t' << "turning height [km] = " << r_max << '\n';
        cout << '\t' << '\t' << "arrival inclination [deg] = " << inclination << '\n';
        cout << '\t' << '\t' << "back azimuth [deg] = " << back_az << '\n';
        cout << '\t' << '\t' << "trans. coeff. [dB] = " << 20.0 * log10(geoac::amp(solution,k)) << '\n';
        cout << '\t' << '\t' << "absorption [dB] = " << -2.0 * attenuation << '\n' << '\n';
    } else {
        cout << '\n';
        cout << '\t' << '\t' << "Ray path does not return to the ground." << '\n' << '\n';
    }

    double c = atmo::c(solution[k][0], solution[k][1], solution[k][2]);
    double rho = atmo::rho(solution[k][0], solution[k][1], solution[k][2]);
    double nu = sqrt(pow(solution[k][3], 2) + pow(solution[k][4], 2) + pow(solution[k][5], 2));

    double cgr  = c * solution[k][3] / nu + atmo::w(solution[k][0], solution[k][1], solution[k][2]);
    double cgth = c * solution[k][4] / nu + atmo::v(solution[k][0], solution[k][1], solution[k][2]);
    double cgph = c * solution[k][5] / nu + atmo::u(solution[k][0], solution[k][1], solution[k][2]);
    double cg = sqrt(pow(cgr, 2) + pow(cgth, 2) + pow(cgph, 2));

    sprintf(output_buffer, "%s.wvfrm_out.dat", output_id);
    wvfrm_out.open(output_buffer);

    wvfrm_out << "# infraga-sph -wnl_wvfrm summary:" << '\n';
    wvfrm_out << "#" << '\t' << "profile: " << inputs[2] << '\n';
    wvfrm_out << "#" << '\t' << "inclination: " << geoac::theta * (180.0 / Pi)  << '\n';
    wvfrm_out << "#" << '\t' << "azimuth: " << geoac::phi * (180.0 / Pi)  << '\n';
    wvfrm_out << "#" << '\t' << "bounces: " << bounces << '\n';
    wvfrm_out << "#" << '\t' << "source (lat, lon, alt): " << lat_src * (180.0 / Pi) << ", " << lon_src * (180.0 / Pi) << ", " << z_src << '\n';
    if(!geoac::is_topo){
        wvfrm_out << "#" << '\t' << "ground elevation: " << topo::z0 << '\n';
    } else {
        wvfrm_out << "#" << '\t' << "topo file:" << topo_file << '\n';
    }    
    wvfrm_out << "#" << '\t' << "S&B atten coeff: " << atmo::tweak_abs << '\n';
    if (strncmp(wvfrm_file, "none", 4) != 0){
        wvfrm_out << "#" << '\t' << "waveform source file: " << wvfrm_file << '\n';
    } else {
        wvfrm_out << "#" << '\t' << "waveform option: " << wvfrm_opt << '\n';
        wvfrm_out << "#" << '\t' << "waveform p0: " << wvfrm::p0 << '\n';
        wvfrm_out << "#" << '\t' << "waveform t0: " << wvfrm::t0 << '\n';
        wvfrm_out << "#" << '\t' << "waveform alpha: " << wvfrm::alpha << '\n';
    }
    wvfrm_out << "#" << '\t' << "waveform reference location: " << wvfrm_ref << '\n' << '\n';

    wvfrm_out << "# t [sec]" << '\t' << "p(t) [Pa]" << '\n';
    for (int n = 0; n < wvfrm::len; n++){
        wvfrm_out << setprecision(8) << travel_time_sum + wvfrm_array[n][0];
        wvfrm_out << '\t' << wvfrm_array[n][1] * p0 * sqrt(fabs( (rho * pow(c, 3) * nu) / (rho0 * pow(c0, 3) * nu0) * (cg0 * D0) / (cg * geoac::jacobian(solution, k)) ));
        wvfrm_out << '\n';
    }
    wvfrm_out.close();

    if(write_rays){
        raypath.close();
    }
    
    wvfrm_out.close();
    wvfrm::delete_wvfrm(wvfrm_array);

    geoac::delete_solution(solution, length);
    clear_region();
}

void run_region_test(char* inputs[], int count){
    cout << '\n';
    cout << '\t' << "###########################################" << '\n';
    cout << '\t' << "####        Running infraga-sph        ####" << '\n';
    cout << '\t' << "####            Region Test            ####" << '\n';
    cout << '\t' << "###########################################" << '\n' << '\n';
    
    ofstream file_out;
    double alt_max = 120.0, lat_src=0.0, lon_src=0.0, eps=1.0e-6;
    int file_check;
    
    char* prof_format = "zTuvdp";
    char* topo_file = "None";
    
    topo::z0 = 0.0;
    atmo::tweak_abs = 1.0;
    geoac::calc_amp = true;
    geoac::is_topo = false;
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "prof_format=", 12) == 0){       prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){        topo::z0 = atof(inputs[i] + 7);}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){    topo_file = inputs[i] + 10; geoac::is_topo=true;}
        else if (strncmp(inputs[i], "eps=", 4) == 0){           eps = atof(inputs[i] + 4);}  
        else if (strncmp(inputs[i], "alt_max=", 8) == 0){       alt_max = atof(inputs[i] + 8);}
        else if (strncmp(inputs[i], "topo_use_BLw=", 13) == 0){ topo::use_BLw = string2bool(inputs[i] + 13);}
    }
    
    if (geoac::is_topo){    file_check = set_region(inputs[2], topo_file, prof_format, false);}
    else{                   file_check = set_region(inputs[2], prof_format, false);}
    if(!file_check){
        return;
    }

    if (geoac::is_topo){
        lat_src = (geoac::lat_max + geoac::lat_min) / 2.0 * (180.0 / Pi);
        lon_src = (geoac::lon_max + geoac::lon_min) / 2.0 * (180.0 / Pi);
    }
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "lat_src=", 8) == 0){       lat_src = atof(inputs[i] + 8);}
        else if (strncmp(inputs[i], "lon_src=", 8) == 0){   lon_src = atof(inputs[i] + 8);}
    }
    
    lat_src *= Pi / 180.0;
    lon_src *= Pi / 180.0;

    // Write the sound speed and its derivatives (analytic and defined)
    cout << '\n' << "Running region test comparing defined derivatives to finite differences (make sure you've created 'region_tests/sph/' directory)..." << '\n';
    cout << '\t' << "Evaluating sound speed..." << '\n';
    file_out.open("region_tests/sph/c.analytic.dat");
    for(int m = 0; m < int(alt_max * 100); m++){        
        double r0 = topo::z(lat_src, lon_src) + m / 100.0;

        if(r0 - globe::r0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;
        
        file_out << r0 - globe::r0 << '\t';
        file_out << atmo::c(r0, lat_src, lon_src) << '\t';

        file_out << atmo::dc(r0, lat_src, lon_src, 0) << '\t';
        file_out << atmo::dc(r0, lat_src, lon_src, 1) << '\t';
        file_out << atmo::dc(r0, lat_src, lon_src, 2) << '\t';
        
        file_out << atmo::ddc(r0, lat_src, lon_src, 0, 0) << '\t';
        file_out << atmo::ddc(r0, lat_src, lon_src, 1, 1) << '\t';
        file_out << atmo::ddc(r0, lat_src, lon_src, 2, 2) << '\t';
        
        file_out << atmo::ddc(r0, lat_src, lon_src, 0, 1) << '\t';
        file_out << atmo::ddc(r0, lat_src, lon_src, 0, 2) << '\t';
        file_out << atmo::ddc(r0, lat_src, lon_src, 1, 2) << '\t';

        file_out << '\n';
    }
    file_out.close();

    file_out.open("region_tests/sph/c.finite.dat");
    for(int m = 1; m < int(alt_max * 100); m++){
        double r0 = topo::z(lat_src, lon_src) + m / 100.0;
        if(r0 - globe::r0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;
        
        file_out << r0 - globe::r0 << '\t';
        file_out << atmo::c(r0, lat_src, lon_src) << '\t';

        file_out << (atmo::c(r0 + eps, lat_src, lon_src) - atmo::c(r0 - eps, lat_src, lon_src)) / (2.0 * eps) << '\t';
        file_out << (atmo::c(r0, lat_src + eps, lon_src) - atmo::c(r0, lat_src - eps, lon_src)) / (2.0 * eps) << '\t';
        file_out << (atmo::c(r0, lat_src, lon_src + eps) - atmo::c(r0, lat_src, lon_src - eps)) / (2.0 * eps) << '\t';
        
        file_out << (atmo::c(r0 + eps, lat_src, lon_src) + atmo::c(r0 - eps, lat_src, lon_src) - 2.0 * atmo::c(r0, lat_src, lon_src)) / pow(eps, 2.0) << '\t';
        file_out << (atmo::c(r0, lat_src + eps, lon_src) + atmo::c(r0, lat_src - eps, lon_src) - 2.0 * atmo::c(r0, lat_src, lon_src)) / pow(eps, 2.0) << '\t';
        file_out << (atmo::c(r0, lat_src, lon_src + eps) + atmo::c(r0, lat_src, lon_src - eps) - 2.0 * atmo::c(r0, lat_src, lon_src)) / pow(eps, 2.0) << '\t';
        
        file_out << (atmo::c(r0 + eps, lat_src + eps, lon_src) + atmo::c(r0 - eps, lat_src - eps, lon_src) - atmo::c(r0 + eps, lat_src - eps, lon_src) - atmo::c(r0 - eps, lat_src + eps, lon_src)) / pow(2.0 * eps, 2.0) << '\t';
        file_out << (atmo::c(r0 + eps, lat_src, lon_src + eps) + atmo::c(r0 - eps, lat_src, lon_src - eps) - atmo::c(r0 + eps, lat_src, lon_src - eps) - atmo::c(r0 - eps, lat_src, lon_src + eps)) / pow(2.0 * eps, 2.0) << '\t';
        file_out << (atmo::c(r0, lat_src + eps, lon_src + eps) + atmo::c(r0, lat_src - eps, lon_src - eps) - atmo::c(r0, lat_src + eps, lon_src - eps) - atmo::c(r0, lat_src - eps, lon_src + eps)) / pow(2.0 * eps, 2.0) << '\t';

        file_out << '\n';
    }
    file_out.close();

    file_out.open("region_tests/sph/c.calc_all.dat");
    for(int m = 1; m < int(alt_max * 100); m++){
        double r0 = topo::z(lat_src, lon_src) + m / 100.0;
        if(r0 - globe::r0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;
        
        double c, dc, ddc;
        interp::eval_all(r0 - globe::r0, atmo::c_spline, c, dc, ddc);

        file_out << r0 - globe::r0 << '\t';
        file_out << c << '\t';

        file_out << dc << '\t';
        file_out << 0.0 << '\t';
        file_out << 0.0 << '\t';
        
        file_out << ddc << '\t';
        file_out << 0.0 << '\t';
        file_out << 0.0 << '\t';
        
        file_out << 0.0 << '\t';
        file_out << 0.0 << '\t';
        file_out << 0.0 << '\t';

        file_out << '\n';
    }
    file_out.close();

    // Write the (u) wind component and its derivatives (analytic and finite diff)
    cout << '\t' << "Evaluating winds (u)..." << '\n';
    file_out.open("region_tests/sph/u.analytic.dat");
    for(int m = 0; m < int(alt_max * 100); m++){
        double r0 = topo::z(lat_src, lon_src) + m/100.0;
        if(r0 - globe::r0 > atmo::u_spline.x_vals[atmo::u_spline.length - 1]) break;
        
        file_out << r0 - globe::r0 << '\t';
        file_out << atmo::u(r0, lat_src, lon_src) << '\t';
        
        file_out << atmo::du(r0, lat_src, lon_src, 0) << '\t';
        file_out << atmo::du(r0, lat_src, lon_src, 1) << '\t';
        file_out << atmo::du(r0, lat_src, lon_src, 2) << '\t';
        
        file_out << atmo::ddu(r0, lat_src, lon_src, 0, 0) << '\t';
        file_out << atmo::ddu(r0, lat_src, lon_src, 1, 1) << '\t';
        file_out << atmo::ddu(r0, lat_src, lon_src, 2, 2) << '\t';
        
        file_out << atmo::ddu(r0, lat_src, lon_src, 0, 1) << '\t';
        file_out << atmo::ddu(r0, lat_src, lon_src, 0, 2) << '\t';
        file_out << atmo::ddu(r0, lat_src, lon_src, 1, 2) << '\t';
        file_out << '\n';
    }
    file_out.close();
    
    file_out.open("region_tests/sph/u.finite.dat");
    for(int m = 1; m < int(alt_max * 100); m++){
        double r0 = topo::z(lat_src, lon_src) + m/100.0;
        if(r0 - globe::r0 > atmo::u_spline.x_vals[atmo::u_spline.length - 1]) break;
        
        file_out << r0 - globe::r0 << '\t';
        file_out << atmo::u(r0, lat_src, lon_src) << '\t';
        
        file_out << (atmo::u(r0 + eps, lat_src, lon_src) - atmo::u(r0 - eps, lat_src, lon_src)) / (2.0 * eps) << '\t';
        file_out << (atmo::u(r0, lat_src + eps, lon_src) - atmo::u(r0, lat_src - eps, lon_src)) / (2.0 * eps) << '\t';
        file_out << (atmo::u(r0, lat_src, lon_src + eps) - atmo::u(r0, lat_src, lon_src - eps)) / (2.0 * eps) << '\t';
        
        file_out << (atmo::u(r0 + eps, lat_src, lon_src) + atmo::u(r0 - eps, lat_src, lon_src) - 2.0 * atmo::u(r0, lat_src, lon_src)) / pow(eps, 2.0) << '\t';
        file_out << (atmo::u(r0, lat_src + eps, lon_src) + atmo::u(r0, lat_src - eps, lon_src) - 2.0 * atmo::u(r0, lat_src, lon_src)) / pow(eps, 2.0) << '\t';
        file_out << (atmo::u(r0, lat_src, lon_src + eps) + atmo::u(r0, lat_src, lon_src - eps) - 2.0 * atmo::u(r0, lat_src, lon_src)) / pow(eps, 2.0) << '\t';
        
        file_out << (atmo::u(r0 + eps, lat_src + eps, lon_src) + atmo::u(r0 - eps, lat_src - eps, lon_src) - atmo::u(r0 + eps, lat_src - eps, lon_src) - atmo::u(r0 - eps, lat_src + eps, lon_src)) / pow(2.0 * eps, 2.0) << '\t';
        file_out << (atmo::u(r0 + eps, lat_src, lon_src + eps) + atmo::u(r0 - eps, lat_src, lon_src - eps) - atmo::u(r0 + eps, lat_src, lon_src - eps) - atmo::u(r0 - eps, lat_src, lon_src + eps)) / pow(2.0 * eps, 2.0) << '\t';
        file_out << (atmo::u(r0, lat_src + eps, lon_src + eps) + atmo::u(r0, lat_src - eps, lon_src - eps) - atmo::u(r0, lat_src + eps, lon_src - eps) - atmo::u(r0, lat_src - eps, lon_src + eps)) / pow(2.0 * eps, 2.0) << '\t';
        
        file_out << '\n';
    }
    file_out.close();

    file_out.open("region_tests/sph/u.calc_uvw.dat");
    for(int m = 0; m < int(alt_max * 100); m++){
        double r0 = topo::z(lat_src, lon_src) + m/100.0;
        if(r0 - globe::r0 > atmo::u_spline.x_vals[atmo::u_spline.length - 1]) break;
        
        double u, v, w, du[3], dv[3], dw[3], ddu[6], ddv[6], ddw[6];
        atmo::calc_uvw(r0, lat_src, lon_src, u, v, w, du, dv, dw, ddu, ddv, ddw);
        atmo::calc_uvw(r0, lat_src, lon_src, u, v, w, du, dv, dw);
            
        file_out << r0 - globe::r0 << '\t';
        file_out << u << '\t';
        
        file_out << du[0] << '\t';
        file_out << du[1] << '\t';
        file_out << du[2] << '\t';
        
        file_out << ddu[0] << '\t';
        file_out << ddu[1] << '\t';
        file_out << ddu[2] << '\t';
        
        file_out << ddu[3] << '\t';
        file_out << ddu[4] << '\t';
        file_out << ddu[5] << '\t';
        file_out << '\n';
    }
    file_out.close();    


    // Write the meridional wind component and its derivatives
    cout << '\t' << "Evaluating winds (v)..." << '\n';
    file_out.open("region_tests/sph/v.analytic.dat");
    for(int m = 0; m < int(alt_max * 100); m++){
        double r0 = topo::z(lat_src, lon_src) + m/100.0;
        if(r0 - globe::r0 > atmo::v_spline.x_vals[atmo::v_spline.length - 1]) break;
        
        file_out << r0 - globe::r0 << '\t';
        file_out << atmo::v(r0, lat_src, lon_src) << '\t';
        
        file_out << atmo::dv(r0, lat_src, lon_src, 0) << '\t';
        file_out << atmo::dv(r0, lat_src, lon_src, 1) << '\t';
        file_out << atmo::dv(r0, lat_src, lon_src, 2) << '\t';
        
        file_out << atmo::ddv(r0, lat_src, lon_src, 0, 0) << '\t';
        file_out << atmo::ddv(r0, lat_src, lon_src, 1, 1) << '\t';
        file_out << atmo::ddv(r0, lat_src, lon_src, 2, 2) << '\t';
        
        file_out << atmo::ddv(r0, lat_src, lon_src, 0, 1) << '\t';
        file_out << atmo::ddv(r0, lat_src, lon_src, 0, 2) << '\t';
        file_out << atmo::ddv(r0, lat_src, lon_src, 1, 2) << '\t';
        file_out << '\n';
    }
    file_out.close();
    
    file_out.open("region_tests/sph/v.finite.dat");
    for(int m = 1; m < int(alt_max * 100); m++){
        double r0 = topo::z(lat_src, lon_src) + m/100.0;
        if(r0 - globe::r0 > atmo::v_spline.x_vals[atmo::v_spline.length - 1]) break;
        
        file_out << r0 - globe::r0 << '\t';
        file_out << atmo::v(r0, lat_src, lon_src) << '\t';
        
        file_out << (atmo::v(r0 + eps, lat_src, lon_src) - atmo::v(r0 - eps, lat_src, lon_src)) / (2.0 * eps) << '\t';
        file_out << (atmo::v(r0, lat_src + eps, lon_src) - atmo::v(r0, lat_src - eps, lon_src)) / (2.0 * eps) << '\t';
        file_out << (atmo::v(r0, lat_src, lon_src + eps) - atmo::v(r0, lat_src, lon_src - eps)) / (2.0 * eps) << '\t';
        
        file_out << (atmo::v(r0 + eps, lat_src, lon_src) + atmo::v(r0 - eps, lat_src, lon_src) - 2.0 * atmo::v(r0, lat_src, lon_src)) / pow(eps, 2.0) << '\t';
        file_out << (atmo::v(r0, lat_src + eps, lon_src) + atmo::v(r0, lat_src - eps, lon_src) - 2.0 * atmo::v(r0, lat_src, lon_src)) / pow(eps, 2.0) << '\t';
        file_out << (atmo::v(r0, lat_src, lon_src + eps) + atmo::v(r0, lat_src, lon_src - eps) - 2.0 * atmo::v(r0, lat_src, lon_src)) / pow(eps, 2.0) << '\t';
        
        file_out << (atmo::v(r0 + eps, lat_src + eps, lon_src) + atmo::v(r0 - eps, lat_src - eps, lon_src) - atmo::v(r0 + eps, lat_src - eps, lon_src) - atmo::v(r0 - eps, lat_src + eps, lon_src)) / pow(2.0 * eps, 2.0) << '\t';
        file_out << (atmo::v(r0 + eps, lat_src, lon_src + eps) + atmo::v(r0 - eps, lat_src, lon_src - eps) - atmo::v(r0 + eps, lat_src, lon_src - eps) - atmo::v(r0 - eps, lat_src, lon_src + eps)) / pow(2.0 * eps, 2.0) << '\t';
        file_out << (atmo::v(r0, lat_src + eps, lon_src + eps) + atmo::v(r0, lat_src - eps, lon_src - eps) - atmo::v(r0, lat_src + eps, lon_src - eps) - atmo::v(r0, lat_src - eps, lon_src + eps)) / pow(2.0 * eps, 2.0) << '\t';
        
        file_out << '\n';
    }
    file_out.close();
    
    file_out.open("region_tests/sph/v.calc_uvw.dat");
    for(int m = 0; m < int(alt_max * 100); m++){
        double r0 = topo::z(lat_src, lon_src) + m/100.0;
        if(r0 - globe::r0 > atmo::u_spline.x_vals[atmo::u_spline.length - 1]) break;
        
        double u, v, w, du[3], dv[3], dw[3], ddu[6], ddv[6], ddw[6];
        atmo::calc_uvw(r0, lat_src, lon_src, u, v, w, du, dv, dw, ddu, ddv, ddw);
        atmo::calc_uvw(r0, lat_src, lon_src, u, v, w, du, dv, dw);
            
        file_out << r0 - globe::r0 << '\t';
        file_out << v << '\t';
        
        file_out << dv[0] << '\t';
        file_out << dv[1] << '\t';
        file_out << dv[2] << '\t';
        
        file_out << ddv[0] << '\t';
        file_out << ddv[1] << '\t';
        file_out << ddv[2] << '\t';
        
        file_out << ddv[3] << '\t';
        file_out << ddv[4] << '\t';
        file_out << ddv[5] << '\t';
        file_out << '\n';
    }
    file_out.close(); 

    // Write the vertical wind component and its derivatives
    cout << '\t' << "Evaluating winds (w)..." << '\n';
    file_out.open("region_tests/sph/w.analytic.dat");
    for(int m = 0; m < int(alt_max * 100); m++){
        double r0 = topo::z(lat_src, lon_src) + m/100.0;
        if(r0 - globe::r0 > atmo::v_spline.x_vals[atmo::v_spline.length - 1]) break;
        
        file_out << r0 - globe::r0 << '\t';
        file_out << atmo::w(r0, lat_src, lon_src) << '\t';
        
        file_out << atmo::dw(r0, lat_src, lon_src, 0) << '\t';
        file_out << atmo::dw(r0, lat_src, lon_src, 1) << '\t';
        file_out << atmo::dw(r0, lat_src, lon_src, 2) << '\t';
        
        file_out << atmo::ddw(r0, lat_src, lon_src, 0, 0) << '\t';
        file_out << atmo::ddw(r0, lat_src, lon_src, 1, 1) << '\t';
        file_out << atmo::ddw(r0, lat_src, lon_src, 2, 2) << '\t';
        
        file_out << atmo::ddw(r0, lat_src, lon_src, 0, 1) << '\t';
        file_out << atmo::ddw(r0, lat_src, lon_src, 0, 2) << '\t';
        file_out << atmo::ddw(r0, lat_src, lon_src, 1, 2) << '\t';
        file_out << '\n';
    }
    file_out.close();
    
    file_out.open("region_tests/sph/w.finite.dat");
    for(int m = 1; m < int(alt_max * 100); m++){
        double r0 = topo::z(lat_src, lon_src) + m/100.0;
        if(r0 - globe::r0 > atmo::u_spline.x_vals[atmo::u_spline.length - 1]) break;
        
        file_out << r0 - globe::r0 << '\t';
        file_out << atmo::w(r0, lat_src, lon_src) << '\t';
        
        file_out << (atmo::w(r0 + eps, lat_src, lon_src) - atmo::w(r0 - eps, lat_src, lon_src)) / (2.0 * eps) << '\t';
        file_out << (atmo::w(r0, lat_src + eps, lon_src) - atmo::w(r0, lat_src - eps, lon_src)) / (2.0 * eps) << '\t';
        file_out << (atmo::w(r0, lat_src, lon_src + eps) - atmo::w(r0, lat_src, lon_src - eps)) / (2.0 * eps) << '\t';
        
        file_out << (atmo::w(r0 + eps, lat_src, lon_src) + atmo::w(r0 - eps, lat_src, lon_src) - 2.0 * atmo::w(r0, lat_src, lon_src)) / pow(eps, 2.0) << '\t';
        file_out << (atmo::w(r0, lat_src + eps, lon_src) + atmo::w(r0, lat_src - eps, lon_src) - 2.0 * atmo::w(r0, lat_src, lon_src)) / pow(eps, 2.0) << '\t';
        file_out << (atmo::w(r0, lat_src, lon_src + eps) + atmo::w(r0, lat_src, lon_src - eps) - 2.0 * atmo::w(r0, lat_src, lon_src)) / pow(eps, 2.0) << '\t';
        
        file_out << (atmo::w(r0 + eps, lat_src + eps, lon_src) + atmo::w(r0 - eps, lat_src - eps, lon_src) - atmo::w(r0 + eps, lat_src - eps, lon_src) - atmo::w(r0 - eps, lat_src + eps, lon_src)) / pow(2.0 * eps, 2.0) << '\t';
        file_out << (atmo::w(r0 + eps, lat_src, lon_src + eps) + atmo::w(r0 - eps, lat_src, lon_src - eps) - atmo::w(r0 + eps, lat_src, lon_src - eps) - atmo::w(r0 - eps, lat_src, lon_src + eps)) / pow(2.0 * eps, 2.0) << '\t';
        file_out << (atmo::w(r0, lat_src + eps, lon_src + eps) + atmo::w(r0, lat_src - eps, lon_src - eps) - atmo::w(r0, lat_src + eps, lon_src - eps) - atmo::w(r0, lat_src - eps, lon_src + eps)) / pow(2.0 * eps, 2.0) << '\t';
        
        file_out << '\n';
    }
    file_out.close();
    
    file_out.open("region_tests/sph/w.calc_uvw.dat");
    for(int m = 0; m < int(alt_max * 100); m++){
        double r0 = topo::z(lat_src, lon_src) + m/100.0;
        if(r0 - globe::r0 > atmo::u_spline.x_vals[atmo::u_spline.length - 1]) break;
        
        double u, v, w, du[3], dv[3], dw[3], ddu[6], ddv[6], ddw[6];
        atmo::calc_uvw(r0, lat_src, lon_src, u, v, w, du, dv, dw, ddu, ddv, ddw);
        // atmo::calc_uvw(r0, lat_src, lon_src, u, v, w, du, dv, dw);
            
        file_out << r0 - globe::r0 << '\t';
        file_out << w << '\t';
        
        file_out << dw[0] << '\t';
        file_out << dw[1] << '\t';
        file_out << dw[2] << '\t';
        
        file_out << ddw[0] << '\t';
        file_out << ddw[1] << '\t';
        file_out << ddw[2] << '\t';
        
        file_out << ddw[3] << '\t';
        file_out << ddw[4] << '\t';
        file_out << ddw[5] << '\t';
        file_out << '\n';
    }
    file_out.close(); 
   
    // Write the topography and its derivatives (analytic and finite diff)
    if(geoac::is_topo){
        cout << '\t' << "Evaluating topography..." << '\n';
        eps=1.0e-3 * Pi / 180.0;
        
        file_out.open("region_tests/sph/topo.dat");
        for(int m1 = 0; m1 < 2 * topo::spline.length_x; m1++){
            double lat0 = topo::spline.x_vals[0] + m1 / (2.0 * topo::spline.length_x) * (topo::spline.x_vals[topo::spline.length_x - 1] - topo::spline.x_vals[0]);
            for(int m2 = 0; m2 < 2 * topo::spline.length_y; m2++){
                double lon0 = topo::spline.y_vals[0] + m2 / (2.0 * topo::spline.length_y) * (topo::spline.y_vals[topo::spline.length_y - 1] - topo::spline.y_vals[0]);

                double zg, dzg[2], ddzg[3];
                interp::eval_all(lat0, lon0, topo::spline, zg, dzg, ddzg);

                file_out << setprecision(8) << lat0 * 180.0 / Pi << '\t';
                file_out << setprecision(8) << lon0 * 180.0 / Pi << '\t';

                file_out << zg << '\t';
                file_out << dzg[0] << '\t';
                file_out << dzg[1] << '\t';
            
                file_out << ddzg[0] << '\t';
                file_out << ddzg[1] << '\t';
                file_out << ddzg[2] << '\t';

                file_out << '\n';
            }
        }
        file_out.close();

        file_out.open("region_tests/sph/topo.lat.analytic.dat");
        for(int m = 0; m < 2 * topo::spline.length_x; m++){
            double r0 = topo::z(lat_src, lon_src);
            double lat0 = topo::spline.x_vals[0] + m / (2.0 * topo::spline.length_x) * (topo::spline.x_vals[topo::spline.length_x - 1] - topo::spline.x_vals[0]);
            double lon0 = lon_src;
            
            file_out << setprecision(8) << lat0 * 180.0 / Pi << '\t';
            file_out << topo::z(lat0, lon0) - globe::r0 << '\t';
            file_out << topo::dz(lat0, lon0, 0) << '\t';
            file_out << topo::dz(lat0, lon0, 1) << '\t';
            
            file_out << topo::ddz(lat0, lon0, 0, 0) << '\t';
            file_out << topo::ddz(lat0, lon0, 0, 1) << '\t';
            file_out << topo::ddz(lat0, lon0, 1, 1) << '\t';

            file_out << topo::dddz(lat0, lon0, 0, 0, 0) << '\t';
            file_out << topo::dddz(lat0, lon0, 0, 0, 1) << '\t';
            file_out << topo::dddz(lat0, lon0, 0, 1, 1) << '\t';
            file_out << topo::dddz(lat0, lon0, 1, 1, 1) << '\t';

            file_out << '\n';
        }
        file_out.close();
        file_out.open("region_tests/sph/topo.lon.analytic.dat");
        for(int m = 0; m < 2 * topo::spline.length_y; m++){
            double r0 = topo::z(lat_src, lon_src);
            double lat0 = lat_src;
            double lon0 = topo::spline.y_vals[0] + m / (2.0 * topo::spline.length_y) * (topo::spline.y_vals[topo::spline.length_y - 1] - topo::spline.y_vals[0]);
            
            file_out << setprecision(8) << lon0 * 180.0 / Pi << '\t';
            file_out << topo::z(lat0, lon0) - globe::r0 << '\t';
            file_out << topo::dz(lat0, lon0, 0) << '\t';
            file_out << topo::dz(lat0, lon0, 1) << '\t';
            
            file_out << topo::ddz(lat0, lon0, 0, 0) << '\t';
            file_out << topo::ddz(lat0, lon0, 0, 1) << '\t';
            file_out << topo::ddz(lat0, lon0, 1, 1) << '\t';

            file_out << topo::dddz(lat0, lon0, 0, 0, 0) << '\t';
            file_out << topo::dddz(lat0, lon0, 0, 0, 1) << '\t';
            file_out << topo::dddz(lat0, lon0, 0, 1, 1) << '\t';
            file_out << topo::dddz(lat0, lon0, 1, 1, 1) << '\t';

            file_out << '\n';
        }
        file_out.close();
                
        file_out.open("region_tests/sph/topo.lat.finite.dat");
        for(int m = 1; m < 2 * topo::spline.length_x - 1; m++){
            double r0 = topo::z(lat_src, lon_src);
            double lat0 = topo::spline.x_vals[0] + m / (2.0 * topo::spline.length_x) * (topo::spline.x_vals[topo::spline.length_x - 1] - topo::spline.x_vals[0]);
            double lon0 = lon_src;
            
            file_out << setprecision(8) << lat0 * 180.0 / Pi << '\t';
            file_out << topo::z(lat0, lon0) - globe::r0 << '\t';
            file_out << (topo::z(lat0 + eps, lon0) - topo::z(lat0 - eps, lon0)) / (2.0 * eps) << '\t';
            file_out << (topo::z(lat0, lon0 + eps) - topo::z(lat0, lon0 - eps)) / (2.0 * eps) << '\t';
            
            file_out << (topo::z(lat0 + eps, lon0) + topo::z(lat0 - eps, lon0) - 2.0 * topo::z(lat0, lon0)) / pow(eps, 2.0) << '\t';
            file_out << (topo::z(lat0 + eps, lon0 + eps) + topo::z(lat0 - eps, lon0 - eps) - topo::z(lat0 + eps, lon0 - eps) - topo::z(lat0 - eps, lon0 + eps)) / pow(2.0 * eps, 2.0) << '\t';
            file_out << (topo::z(lat0, lon0 + eps) + topo::z(lat0, lon0 - eps) - 2.0 * topo::z(lat0, lon0)) / pow(eps, 2.0) << '\t';
            
            file_out << '\n';
        }
        file_out.close();
        file_out.open("region_tests/sph/topo.lon.finite.dat");
        for(int m = 0; m < 2 * topo::spline.length_y; m++){
            double r0 = topo::z(lat_src, lon_src);
            double lat0 = lat_src;
            double lon0 = topo::spline.y_vals[0] + m / (2.0 * topo::spline.length_y) * (topo::spline.y_vals[topo::spline.length_y - 1] - topo::spline.y_vals[0]);
            
            file_out << setprecision(8) << lon0 * 180.0 / Pi << '\t';
            file_out << topo::z(lat0, lon0) - globe::r0 << '\t';
            file_out << (topo::z(lat0 + eps, lon0) - topo::z(lat0 - eps, lon0)) / (2.0 * eps) << '\t';
            file_out << (topo::z(lat0, lon0 + eps) - topo::z(lat0, lon0 - eps)) / (2.0 * eps) << '\t';
            
            file_out << (topo::z(lat0 + eps, lon0) + topo::z(lat0 - eps, lon0) - 2.0 * topo::z(lat0, lon0)) / pow(eps, 2.0) << '\t';
            file_out << (topo::z(lat0 + eps, lon0 + eps) + topo::z(lat0 - eps, lon0 - eps) - topo::z(lat0 + eps, lon0 - eps) - topo::z(lat0 - eps, lon0 + eps)) / pow(2.0 * eps, 2.0) << '\t';
            file_out << (topo::z(lat0, lon0 + eps) + topo::z(lat0, lon0 - eps) - 2.0 * topo::z(lat0, lon0)) / pow(eps, 2.0) << '\t';
            
            file_out << '\n';
        }
        file_out.close();

        file_out.open("region_tests/sph/topo.lat.calc_all.dat");
        for(int m = 0; m < 2 * topo::spline.length_x - 1; m++){
            double r0 = topo::z(lat_src, lon_src);
            double lat0 = topo::spline.x_vals[0] + m / (2.0 * topo::spline.length_x) * (topo::spline.x_vals[topo::spline.length_x - 1] - topo::spline.x_vals[0]);
            double lon0 = lon_src;
            
            double zg, dzg[2], ddzg[3], dddzg[4];
            interp::eval_all(lat0, lon0, topo::spline, zg, dzg, ddzg, dddzg);

            file_out << setprecision(8) << lat0 * 180.0 / Pi << '\t';
            file_out << zg << '\t';
            file_out << dzg[0] << '\t';
            file_out << dzg[1] << '\t';
            
            file_out << ddzg[0] << '\t';
            file_out << ddzg[1] << '\t';
            file_out << ddzg[2] << '\t';
            
            file_out << dddzg[0] << '\t';
            file_out << dddzg[1] << '\t';
            file_out << dddzg[2] << '\t';
            file_out << dddzg[3] << '\t';

            file_out << '\n';
        }
        file_out.close();
        file_out.open("region_tests/sph/topo.lon.calc_all.dat");
        for(int m = 0; m < 2 * topo::spline.length_y; m++){
            double r0 = topo::z(lat_src, lon_src);
            double lat0 = lat_src;
            double lon0 = topo::spline.y_vals[0] + m / (2.0 * topo::spline.length_y) * (topo::spline.y_vals[topo::spline.length_y - 1] - topo::spline.y_vals[0]);
            
            double zg, dzg[2], ddzg[3], dddzg[4];
            interp::eval_all(lat0, lon0, topo::spline, zg, dzg, ddzg, dddzg);

            file_out << setprecision(8) << lon0 * 180.0 / Pi << '\t';
            file_out << zg << '\t';
            file_out << dzg[0] << '\t';
            file_out << dzg[1] << '\t';
            
            file_out << ddzg[0] << '\t';
            file_out << ddzg[1] << '\t';
            file_out << ddzg[2] << '\t';

            file_out << dddzg[0] << '\t';
            file_out << dddzg[1] << '\t';
            file_out << dddzg[2] << '\t';
            file_out << dddzg[3] << '\t';

            
            file_out << '\n';
        }
        file_out.close();


    }

    clear_region();
}


int main(int argc, char* argv[]){
    if (argc < 3){
        usage();
        return 0;
    } else {
        if ((strncmp(argv[1], "-version", 8) == 0) || (strncmp(argv[1], "-v", 2) == 0)){            version();}
        else if ((strncmp(argv[1], "-usage", 6) == 0) || (strncmp(argv[1], "-u", 2) == 0)){         usage();}
        else if ((strncmp(argv[1], "-prop", 5) == 0) || (strncmp(argv[1], "-p", 2) == 0)){          run_prop(argv, argc);}
        else if ((strncmp(argv[1], "-back_proj", 10) == 0) || (strncmp(argv[1], "-b", 2) == 0)){    run_back_proj(argv, argc);}
        else if ((strncmp(argv[1], "-eig_search", 11) == 0) || (strncmp(argv[1], "-s", 2) == 0)){   run_eig_search(argv, argc);}
        else if ((strncmp(argv[1], "-eig_direct", 11) == 0) || (strncmp(argv[1], "-d", 2) == 0)){   run_eig_direct(argv, argc);}
        else if ((strncmp(argv[1], "-wnl_wvfrm", 10) == 0) || (strncmp(argv[1], "-w", 2) == 0)){    run_wnl_wvfrm(argv, argc);}
        else if (strncmp(argv[1], "-region_test", 12) == 0){                                        run_region_test(argv, argc);} // Development use only to evaluate interpolation
        else {                                                                                      cout << "Unrecognized option." << '\n';}
    }
    return 0;
}
