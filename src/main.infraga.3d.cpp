#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <sstream>
#include <math.h>

#include "atmo/atmo_state.h"
#include "atmo/atmo_io.3d.strat.h"

#include "geoac/geoac.params.h"
#include "geoac/geoac.eqset.h"
#include "geoac/geoac.interface.h"
#include "geoac/geoac.eigenray.h"

#include "util/fileIO.h"
#include "util/rk4solver.h"
#include "util/interpolation.h"
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
    cout << '\t' << "#############################################" << '\n';
    cout << '\t' << "####              infraga-3d             ####" << '\n';
    cout << '\t' << "####    Three-Dimensional Ray Tracing    ####" << '\n';
    cout << '\t' << "#### Through a Stratified, Moving Medium ####" << '\n';
    cout << '\t' << "#############################################" << '\n' << '\n';
    
    cout << '\n';
    cout << "Usage: infraga-3d [option] profile.met [parameters]" << '\n';
    cout << '\t' << '\t' << "Enter only 1 option." << '\n';
    cout << '\t' << '\t' << "Profile.met is expected to contain columns describing {z[km]  T[K]  u (zonal wind) [m/s]  v (meridional wind) [m/s]  density[g/cm^3]  p[mbar]} " << '\n';
    cout << '\t' << '\t' << '\t' << "Profile format can be modified, see manual document for details." << '\n';
    cout << '\t' << '\t' << "Parameter calls are expected using the format: parameter_name=value." << '\n' << '\n';
    
    cout << "Options and parameters are:" << '\n';
    cout << '\t' << "-prop (generate ray paths for propagations at multiple azimuth and inclination angles)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "incl_min"          << '\t' << "degrees"            << '\t' << '\t' << "0.5" << '\n';
    cout << '\t' << '\t' << "incl_max"          << '\t' << "degrees"            << '\t' << '\t' << "45.0" << '\n';
    cout << '\t' << '\t' << "incl_step"         << '\t' << "degrees"            << '\t' << '\t' << "0.5"  << '\n';
    cout << '\t' << '\t' << "inclination"       << '\t' << "see manual"         << '\t' << "0.5" << '\n';

    cout << '\t' << '\t' << "az_min"            << '\t' << '\t' << "degrees"    << '\t' << '\t' << "-90.0" << '\n';
    cout << '\t' << '\t' << "az_max"            << '\t' << '\t' << "degrees"    << '\t' << '\t' << "-90.0" << '\n';
    cout << '\t' << '\t' << "az_step"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "1.0"  << '\n';
    cout << '\t' << '\t' << "azimuth"           << '\t' << '\t' << "see manual" << '\t' << "-90.0" << '\n';
    
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "2" << '\n';
    
    cout << '\t' << '\t' << "src_x"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "src_y"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
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
    cout << '\t' << '\t' << "rcvr_x"            << '\t' << '\t' << "km"         << '\t' << '\t' << "250.0" << '\n';
    cout << '\t' << '\t' << "rcvr_y"            << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n' << '\n';
    */
    cout << '\t' << "-eig_search (Search for all eigenrays connecting a source at (src_x, src_y, src_alt) to a receiver " << '\n';
    cout << '\t' << '\t' << '\t' << "at (rcvr_x, rcvr_y, z_grnd) which have inclinations and ground reflections within specified limits)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "incl_min"         << '\t' << "degrees"            << '\t' << '\t' << "0.5" << '\n';
    cout << '\t' << '\t' << "incl_max"         << '\t' << "degrees"            << '\t' << '\t' << "45.0" << '\n';
    cout << '\t' << '\t' << "bnc_min"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "0" << '\n';
    cout << '\t' << '\t' << "bnc_max"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "0" << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "see manual" << '\t' << "0" << '\n';
    cout << '\t' << '\t' << "src_x"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "src_y"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "src_alt"           << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "rcvr_x"            << '\t' << '\t' << "km"         << '\t' << '\t' << "250.0" << '\n';
    cout << '\t' << '\t' << "rcvr_y"            << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "verbose"           << '\t' << '\t' << "true/false" << '\t' << "false" << '\n';
    cout << '\t' << '\t' << "iterations"        << '\t' << "integer"            << '\t' << '\t' << "25" << '\n';
    cout << '\t' << '\t' << "damping"           << '\t' << '\t' << "scalar"     << '\t' << '\t' << "1.0e-3" << '\n';
    cout << '\t' << '\t' << "tolerance"         << '\t' << "km"                 << '\t' << '\t' << "0.1" << '\n';
    cout << '\t' << '\t' << "az_dev_lim"        << '\t' << "degrees"            << '\t' << '\t' << "2.0" << '\n';
    cout << '\t' << '\t' << "incl_step_min"     << '\t' << "degrees"            << '\t' << '\t' << "0.001" << '\n';
    cout << '\t' << '\t' << "incl_step_max"     << '\t' << "degrees"            << '\t' << '\t' << "0.1" << '\n' << '\n';
  
    cout << '\t' << "-eig_direct (Search for a single eigenray connecting a source at (src_x, src_y, src_alt) to a receiver at (rcvr_x, rcvr_y, z_grnd)" << '\n';
    cout << '\t' << '\t' << '\t' << "near an estimated azimuth and inclination pair assuming a specific number of ground reflections)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "incl_est"          << '\t' << "degrees"            << '\t' << '\t' << "15.0" << '\n';
    cout << '\t' << '\t' << "az_est"            << '\t' <<'\t' << "degrees"     << '\t' << '\t' << "atan2(rvcr_y, rcvr_x)" << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "0" << '\n';
    cout << '\t' << '\t' << "src_x"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "src_y"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "src_alt"           << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "rcvr_x"            << '\t' << '\t' << "km"         << '\t' << '\t' << "250.0" << '\n';
    cout << '\t' << '\t' << "rcvr_y"            << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "verbose"           << '\t' << '\t' << "true/false" << '\t' << "true" << '\n';
    cout << '\t' << '\t' << "iterations"        << '\t' << "integer"            << '\t' << '\t' << "25" << '\n';
    cout << '\t' << '\t' << "damping"           << '\t' << '\t' << "scalar"     << '\t' << '\t' << "1.0e-3" << '\n';
    cout << '\t' << '\t' << "tolerance"         << '\t' << "km"                 << '\t' << '\t' << "0.1" << '\n' << '\n';

    cout << '\t' << "-wnl_wvfrm (compute the weakly non-linear waveform along a specific ray path)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "inclination"       << '\t' << "degrees"            << '\t' << '\t' << "15.0"  << '\n';
    cout << '\t' << '\t' << "azimuth"           << '\t' << '\t' << "degrees"    << '\t' << '\t' << "90.0"  << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "integer"    << '\t' << '\t' << "0"  << '\n';
    cout << '\t' << '\t' << "src_x"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
    cout << '\t' << '\t' << "src_y"             << '\t' << '\t' << "km"         << '\t' << '\t' << "0.0" << '\n';
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
    cout << '\t' << "output_id"         << '\t' << '\t' << "see manual"         << '\t' << "from profile.met" << '\n';
    cout << '\t' << "write_caustics"    << '\t' << '\t' << "true/false"         << '\t' << "false" << '\n';
    cout << '\t' << "calc_amp"          << '\t' << '\t' << "true/false"         << '\t' << "true" << '\n';
    cout << '\t' << "max_alt"           << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "interpolation max" << '\n';
    cout << '\t' << "max_rng"           << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "2000.0" << '\n';
    cout << '\t' << "min_x"             << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "-2000.0" << '\n';
    cout << '\t' << "min_y"             << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "-2000.0" << '\n';
    cout << '\t' << "max_x"             << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "2000.0" << '\n';
    cout << '\t' << "max_y"             << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "2000.0" << '\n';
    cout << '\t' << "min_ds"            << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "0.001" << '\n';
    cout << '\t' << "max_ds"            << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "0.05" << '\n';
    cout << '\t' << "max_s"             << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "2500.0" << '\n';
    cout << '\t' << "topo_file"         << '\t' << '\t' << "see manual"         << '\t' << "none" << '\n';
    cout << '\t' << "topo_use_BLw"      << '\t' << '\t' << "see manual"         << '\t' << "false" << '\n' << '\n';

    cout << "Output (see output files or manual for units):" << '\n';
    cout << '\t' << "atmo.dat -> z[km] : c [m/s]  : u (zonal winds) [m/s] : v (meridional winds) [m/s] : density[g/cm^3] : ceff [km/s]" << '\n';
    cout << '\t' << "{...}.raypaths.dat -> x : y : z : trans. coeff. : absorption : time " << '\n';
    cout << '\t' << "{...}.arrivals.dat -> incl : az : n_bnc : x : y : time : cel : z_max : arrival incl : back az : trans. coeff. : absorption" << '\n' << '\n';
    // cout << '\t' << "{...}.projection.dat -> x : y : z : t : X_incl : Y_incl : Z_incl : T_incl : X_az : Y_az : Z_az : T_az" << '\n' << '\n';

    cout << "Examples:" << '\n';
    cout << '\t' << "./bin/infraga-3d -prop examples/ToyAtmo.met incl_step=1.0 bounces=2 azimuth=-75.0 max_rng=500.0" << '\n';
    cout << '\t' << "./bin/infraga-3d -eig_search examples/ToyAtmo.met rcvr_x=-175.0 rcvr_y=75.0 verbose=true incl_step_max=0.2" << '\n';
    cout << '\t' << "./bin/infraga-3d -eig_direct examples/ToyAtmo.met rcvr_x=-175.0 rcvr_y=75.0 incl_est=8.0" << '\n';
    // cout << '\t' << "./bin/infraga-3d -back_proj examples/ToyAtmo.met rcvr_x=-175.0 rcvr_y=75.0 azimuth=113.87504 inclination=11.164744" << '\n';
    cout << '\t' << "./bin/infraga-3d -wnl_wvfrm examples/ToyAtmo.met azimuth=-66.124958 inclination=11.164744 wvfrm_opt=impulse wvfrm_p0=500.0" << '\n' << '\n';
}

void run_prop(char* inputs[], int count){
    cout << '\n';
    cout << '\t' << "##########################################" << '\n';
    cout << '\t' << "####        Running infraga-3d        ####" << '\n';
    cout << '\t' << "####            Propagation           ####" << '\n';
    cout << '\t' << "##########################################" << '\n' << '\n';
    
    double theta_min = 0.5, theta_max=45.0, theta_step=0.5;
    double phi_min=-90.0, phi_max=-90.0, phi_step=1.0;
    int bounces=2;
    double x_src=0.0, y_src=0.0, z_src=0.0;
    bool  write_atmo=false, write_rays=true, write_caustics=false, write_topo=false, custom_output_id=false;
    double freq=0.1, turn_ht_min=0.2;
    char* prof_format = "zTuvdp";
    char* topo_file = "None";
    char* output_id;
    char input_check;

    topo::z0 = 0.0;
    atmo::tweak_abs = 1.0;
    geoac::calc_amp = true;
    geoac::is_topo = false;
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "prof_format=",12) == 0){        prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){        topo::z0 = atof(inputs[i] + 7);}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){    topo_file = inputs[i] + 10; geoac::is_topo=true;}
    }
    
    if (count < 3){         cout << "You have to specify an atmosphere file..." << '\n'; return;}
    if (geoac::is_topo){    set_region(inputs[2], topo_file, prof_format, false);}
    else{                   set_region(inputs[2], prof_format, false);}
    
    if (geoac::is_topo){
        x_src = (geoac::x_max + geoac::x_min) / 2.0;
        y_src = (geoac::y_max + geoac::y_min) / 2.0;
    }
    
    for(int i = 3; i < count; i++){
        if ((strncmp(inputs[i], "incl_min=", 9) == 0) || (strncmp(inputs[i], "min_incl=", 9) == 0)){        theta_min = atof(inputs[i] + 9);}
        else if ((strncmp(inputs[i], "incl_max=", 9) == 0) || (strncmp(inputs[i], "max_incl=", 9) == 0)){   theta_max = atof(inputs[i] + 9);}
        else if ((strncmp(inputs[i], "az_min=", 7) == 0) || (strncmp(inputs[i], "max_az=", 7) == 0)){       phi_min = atof(inputs[i] + 7);}
        else if ((strncmp(inputs[i], "az_max=", 7) == 0) || (strncmp(inputs[i], "min_az=", 7) == 0)){       phi_max = atof(inputs[i] + 7);}

        else if (strncmp(inputs[i], "incl_step=", 10) == 0){                                                theta_step = atof(inputs[i] + 10);}
        else if (strncmp(inputs[i], "az_step=", 8) == 0){                                                   phi_step = atof(inputs[i] + 8);}

        else if (strncmp(inputs[i], "inclination=", 12) == 0){                                              theta_min = atof(inputs[i] + 12);
                                                                                                            theta_max = atof(inputs[i] + 12);
                                                                                                            theta_step = 1.0;}
        else if (strncmp(inputs[i], "azimuth=", 8) == 0){                                                   phi_min = atof(inputs[i] + 8);
                                                                                                            phi_max = atof(inputs[i] + 8);
                                                                                                            phi_step = 1.0;}

        else if (strncmp(inputs[i], "bounces=", 8) == 0){                                                   bounces = atoi(inputs[i] + 8);}
        
        else if ((strncmp(inputs[i], "src_x=", 6) == 0) || (strncmp(inputs[i], "x_src=", 6) == 0)){         x_src = atof(inputs[i] + 6);}
        else if ((strncmp(inputs[i], "src_y=", 6) == 0) || (strncmp(inputs[i], "y_src=", 6) == 0)){         y_src = atof(inputs[i] + 6);}
        else if ((strncmp(inputs[i], "src_alt=", 8) == 0) || (strncmp(inputs[i], "alt_src=", 6) == 0)){     z_src = atof(inputs[i] + 8);}
        
        else if (strncmp(inputs[i], "turn_ht_min=", 12) == 0){                                              turn_ht_min = atof(inputs[i] + 12);}
        else if (strncmp(inputs[i], "write_rays=", 11) == 0){                                               write_rays = string2bool(inputs[i] + 11);}
        
        else if (strncmp(inputs[i], "freq=",5) == 0){                                                       freq = atof(inputs[i] + 5);}
        else if (strncmp(inputs[i], "abs_coeff=", 10) == 0){                                                atmo::tweak_abs = max(0.0, atof(inputs[i] + 10));}
        else if (strncmp(inputs[i], "write_caustics=", 15) == 0){                                           write_caustics = string2bool(inputs[i] + 15);}
        else if (strncmp(inputs[i], "calc_amp=", 9) == 0){                                                  geoac::calc_amp = string2bool(inputs[i] + 9);}
        
        else if ((strncmp(inputs[i], "max_alt=", 8) == 0) || (strncmp(inputs[i], "alt_max=", 8) == 0)){     geoac::alt_max = atof(inputs[i] + 8);   cout << '\t' << "User selected altitude maximum = " << geoac::alt_max << '\n';}
        else if ((strncmp(inputs[i], "max_rng=", 8) == 0) || (strncmp(inputs[i], "rng_max=", 8) == 0)){     geoac::rng_max = atof(inputs[i] + 8);   cout << '\t' << "User selected range maximum = " << geoac::rng_max << '\n';}
        else if ((strncmp(inputs[i], "min_x=", 6) == 0) || (strncmp(inputs[i], "x_min=", 6) == 0)){         geoac::x_min = atof(inputs[i] + 6);     cout << '\t' << "User selected x minimum = " << geoac::x_min << '\n';}
        else if ((strncmp(inputs[i], "max_x=", 6) == 0) || (strncmp(inputs[i], "x_max=", 6) == 0)){         geoac::x_max = atof(inputs[i] + 6);     cout << '\t' << "User selected x maximum = " << geoac::x_max << '\n';}
        else if ((strncmp(inputs[i], "min_y=", 6) == 0) || (strncmp(inputs[i], "y_min=", 6) == 0)){         geoac::y_min = atof(inputs[i] + 6);     cout << '\t' << "User selected y minimum = " << geoac::y_min << '\n';}
        else if ((strncmp(inputs[i], "max_y=", 6) == 0) || (strncmp(inputs[i], "y_max=", 6) == 0)){         geoac::y_max = atof(inputs[i] + 6);     cout << '\t' << "User selected y maximum = " << geoac::y_max << '\n';}

        else if ((strncmp(inputs[i], "min_ds=", 7) == 0) || (strncmp(inputs[i], "ds_min=", 7) == 0)){       geoac::ds_min = atof(inputs[i] + 7);   cout << '\t' << "User selected ds minimum = " << geoac::ds_min << '\n';}
        else if ((strncmp(inputs[i], "max_ds=", 7) == 0) || (strncmp(inputs[i], "ds_max=", 7) == 0)){       geoac::ds_max = atof(inputs[i] + 7);   cout << '\t' << "User selected ds maximum = " << geoac::ds_max << '\n';}
        else if ((strncmp(inputs[i], "max_s=", 6) == 0) || (strncmp(inputs[i], "s_max=", 6) == 0)){         geoac::s_max = atof(inputs[i] + 6);    cout << '\t' << "User selected s maximum (path length between reflections) = " << geoac::s_max << '\n';}

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
        else if (strncmp(inputs[i], "write_topo=", 11) == 0){                                               write_topo = string2bool(inputs[i] + 11);}
        else{
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }
    
    z_src = max(z_src,  topo::z(x_src,y_src));
    if(write_atmo)      geoac::write_prof("atmo.dat", x_src, y_src, (90.0 - phi_min) * Pi / 180.0);
    if(write_caustics)  geoac::calc_amp = true;
    geoac::configure();

    cout << '\n' << "Parameter summary:" << '\n';
    cout << '\t' << "inclination: " << theta_min << ", " << theta_max << ", " << theta_step << '\n';
    cout << '\t' << "azimuth: " << phi_min << ", " << phi_max << ", " << phi_step << '\n';
    cout << '\t' << "bounces: " << bounces << '\n';
    cout << '\t' << "source location: " << x_src << ", " << y_src << ", " << z_src << '\n';
    if(!geoac::is_topo){
        cout << '\t' << "ground elevation: " << topo::z0 << '\n';
    }
    cout << '\t' << "frequency: " << freq << '\n';
    cout << '\t' << "S&B atten coeff: " << atmo::tweak_abs << '\n' << '\n';
    
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
	double D, D_prev, travel_time_sum, attenuation, z_max, inclination, back_az;
	int k, length = int(geoac::s_max / geoac::ds_min);
	bool break_check;
    
	ofstream results, raypath, caustics, topo_out;
    
    sprintf(output_buffer, "%s.arrivals.dat", output_id);
    results.open(output_buffer);
    results << "# infraga-3d -prop summary:" << '\n';
    results << "#" << '\t' << "profile: " << inputs[2] << '\n';
    results << "#" << '\t' << "inclination: " << theta_min << ", " << theta_max << ", " << theta_step << '\n';
    results << "#" << '\t' << "azimuth: " << phi_min << ", " << phi_max << ", " << phi_step << '\n';
    results << "#" << '\t' << "bounces: " << bounces << '\n';
    results << "#" << '\t' << "source location: " << x_src << ", " << y_src << ", " << z_src << '\n';
    if(!geoac::is_topo){
        results << "#" << '\t' << "ground elevation: " << topo::z0 << '\n';
    } else {
        results << "#" << '\t' << "topo file:" << topo_file << '\n';
    }
    results << "#" << '\t' << "frequency: " << freq << '\n';
    results << "#" << '\t' << "S&B atten coeff: " << atmo::tweak_abs << '\n' << '\n';

    results << "# incl [deg]";
    results << '\t' << "az [deg]";
    results << '\t' << "n_b";
    results << '\t' << "x_0 [km]";
    results << '\t' << "y_0 [km]";
    results << '\t' << "time [s]";
    results << '\t' << "cel [km/s]";
    results << '\t' << "turning ht [km]";
    results << '\t' << "arrival incl [deg]";
    results << '\t' << "back az [deg]";
    results << '\t' << "trans. coeff. [dB]";
    results << '\t' << "absorption [dB]";
    results << '\n';
    
	if(write_rays){
        sprintf(output_buffer, "%s.raypaths.dat", output_id);
        raypath.open(output_buffer);
        raypath << "# x [km]";
        raypath << '\t' << "y [km]";
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
            
            caustics << "# x [km]";
            caustics << '\t' << "y [km]";
            caustics << '\t' << "z [km]";
            caustics << '\t' << "time [s]";
            caustics << '\n';
            
            caustics.close();
        }
    }
    
	double** solution;
    geoac::build_solution(solution,length);
    
    for(double phi = phi_min; phi <= phi_max; phi+=phi_step){
        geoac::phi = Pi / 2.0 - phi * Pi / 180.0;

        double dzgdx, dzgdy, theta_grnd;
        if (z_src - topo::z(x_src, y_src) < 0.01){
            if (geoac::is_topo){
                theta_grnd = atan(topo::dz(x_src, y_src, 0) * sin(phi * Pi / 180.0) + topo::dz(x_src, y_src, 1) * cos(phi * Pi / 180.0)) * (180.0 / Pi) + 0.1;
            } else {
                theta_grnd = 0.1;
            }
        } else {
            theta_grnd = -89.9;
        }

        for(double theta = max(theta_min, theta_grnd); theta <= max(theta_max, theta_grnd); theta+=theta_step){
            geoac::theta = theta * Pi / 180.0;

            geoac::set_initial(solution, x_src, y_src, z_src);
            travel_time_sum = 0.0;
            attenuation = 0.0;
            z_max = topo::z0;

            if((fabs(theta - max(theta_min, theta_grnd)) < theta_step) && write_topo){
                topo_out.open("topography.dat");
            }

    		cout << "Calculating ray path: " << theta << " degrees inclination, " << phi << " degrees azimuth." << '\n';
            for(int bnc_cnt = 0; bnc_cnt <= bounces; bnc_cnt++){
                k = geoac::prop_rk4(solution, break_check, length);
                
                if(write_rays || write_caustics){
                    if(write_caustics){
                        sprintf(output_buffer, "%s.caustics-%i.dat", output_id, bnc_cnt);
                        caustics.open(output_buffer,fstream::app);

                        D_prev = geoac::jacobian(solution,1);
                    }
                
                    for(int m = 1; m < k ; m++){
                        geoac::travel_time(travel_time_sum, solution, m - 1, m);
                        geoac::atten(attenuation, solution, m - 1, m, freq);
                        if(write_caustics){
                            D = geoac::jacobian(solution, m);
                        }
                    
                        if(write_rays && (m == 1 || m % 15 == 0)){
                            raypath << solution[m][0];
                            raypath << '\t' << solution[m][1];
                            raypath << '\t' << max(solution[m][2], topo::z(solution[m][0], solution[m][1]));
                            if(geoac::calc_amp){
                                raypath << '\t' << 20.0 * log10(geoac::amp(solution, m));
                            } else{
                                raypath << '\t' << 0.0;
                            }
                            raypath << '\t' << -2.0 * attenuation;
                            raypath << '\t' << travel_time_sum;
                            raypath << '\n';
                        }

                        if(write_caustics && D * D_prev < 0.0){
                            caustics << solution[m][0];
                            caustics << '\t' << solution[m][1];
                            caustics << '\t' << max(solution[m][2], topo::z(solution[m][0], solution[m][1]));
                            caustics <<'\t' << travel_time_sum << '\n';
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
                    attenuation+= geoac::atten(solution, k, freq);
                }

                if((fabs(theta - max(theta_min, theta_grnd)) < theta_step) && write_topo){
                    for(int m = 1; m < k ; m+=10){                        
                        topo_out << solution[m][0];
                        topo_out << '\t' << solution[m][1];
                        topo_out << '\t' << topo::z(solution[m][0], solution[m][1]) << '\n';
                    }
                }


                for(int m = 0; m < k ; m++){
                    z_max = max(z_max, solution[m][2]);
                }
                
                if(!geoac::is_topo && z_max < (topo::z0 + turn_ht_min)){
                    break;
                }

                if(break_check || k < 2){
                    break;
                }
                
                inclination = - asin(atmo::c(solution[k][0], solution[k][1], topo::z(solution[k][0], solution[k][1])) / atmo::c(x_src, y_src, z_src) * solution[k][5]) * 180.0 / Pi;
                back_az = 90.0 - atan2(-solution[k][4], -solution[k][3]) * 180.0 / Pi;
                if(back_az < -180.0){
                    back_az +=360.0;
                } else if(back_az >  180.0){
                    back_az -=360.0;
                }
                
                results << theta;
                results << '\t' << phi;
                results << '\t' << bnc_cnt;
                results << '\t' << solution[k][0];
                results << '\t' << solution[k][1];
                results << '\t' << travel_time_sum;
                results << '\t' << sqrt(pow(solution[k][0]- x_src, 2) + pow(solution[k][1] - y_src, 2)) / travel_time_sum;
                results << '\t' << z_max;
                results << '\t' << inclination;
                results << '\t' << back_az;
                if(geoac::calc_amp){
                    results << '\t' << 20.0 * log10(geoac::amp(solution,k));
                } else{
                    results << '\t' << 0.0;
                }
                results << '\t' << -2.0 * attenuation;
                results << '\n';           

                geoac::set_refl(solution,k);
            }
            if(write_rays){
                raypath << '\n';
            }
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
    cout << '\t' << "####        Running infraga-3d        ####" << '\n';
    cout << '\t' << "####          Back Projection         ####" << '\n';
    cout << '\t' << "##########################################" << '\n' << '\n';

    double x_rcvr = 0.0, y_rcvr = 0.0, z_rcvr = 0.0;
    double freq = 0.1, D, D_prev;
    int bounces = 0;
    bool write_atmo = false, custom_output_id=false;
    char* prof_format = "zTuvdp";
    char* topo_file = "None";
    char* output_id;
    char input_check;
    
    geoac::theta = 15.0 * (Pi / 180.0);
    geoac::phi = Pi / 2.0 - 90.0 * (Pi / 180.0);
    
    topo::z0 = 0.0;
    atmo::tweak_abs = 1.0;
    geoac::calc_amp = true;
    geoac::is_topo = false;
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "prof_format=", 12) == 0){       prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){        topo::z0 = atof(inputs[i] + 7);}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){    topo_file = inputs[i] + 10; geoac::is_topo=true;}
    }
    
    if (count < 3){         cout << "You have to specify an atmosphere file..." << '\n'; return;}
    if (geoac::is_topo){    set_region(inputs[2], topo_file, prof_format, true);}
    else{                   set_region(inputs[2], prof_format, true);}
    
    if (geoac::is_topo){
        x_rcvr = (geoac::x_max + geoac::x_min) / 2.0;
        y_rcvr = (geoac::y_max + geoac::y_min) / 2.0;
    }
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "inclination=", 12) == 0){                                                   geoac::theta = atof(inputs[i] + 12) * (Pi / 180.0);}
        else if (strncmp(inputs[i], "azimuth=", 8) == 0){                                                   geoac::phi = Pi / 2.0 - atof(inputs[i] + 8) * (Pi / 180.0);}
        else if (strncmp(inputs[i], "bounces=", 8) == 0){                                                   bounces = max(0, atoi(inputs[i] + 8));}
        
        else if ((strncmp(inputs[i], "rcvr_x=", 7) == 0) || (strncmp(inputs[i], "x_rcvr=", 7) == 0)){       x_rcvr = atof(inputs[i] + 7);}
        else if ((strncmp(inputs[i], "rcvr_y=", 7) == 0) || (strncmp(inputs[i], "y_rcvr=", 7) == 0)){       y_rcvr = atof(inputs[i] + 7);}
        
        else if (strncmp(inputs[i], "freq=", 5) == 0){                                                      freq = atof(inputs[i] + 5);}
        else if (strncmp(inputs[i], "abs_coeff=", 10) == 0){                                                atmo::tweak_abs = max(0.0, atof(inputs[i] + 10));}
        
        else if ((strncmp(inputs[i], "max_alt=", 8) == 0) || (strncmp(inputs[i], "alt_max=", 8) == 0)){     geoac::alt_max = atof(inputs[i] + 8);   cout << '\t' << "User selected altitude maximum = " << geoac::alt_max << '\n';}
        else if ((strncmp(inputs[i], "max_rng=", 8) == 0) || (strncmp(inputs[i], "rng_max=", 8) == 0)){     geoac::rng_max = atof(inputs[i] + 8);   cout << '\t' << "User selected range maximum = " << geoac::rng_max << '\n';}
        else if ((strncmp(inputs[i], "min_x=", 6) == 0) || (strncmp(inputs[i], "x_min=", 6) == 0)){         geoac::x_min = atof(inputs[i] + 6);     cout << '\t' << "User selected x minimum = " << geoac::x_min << '\n';}
        else if ((strncmp(inputs[i], "max_x=", 6) == 0) || (strncmp(inputs[i], "x_max=", 6) == 0)){         geoac::x_max = atof(inputs[i] + 6);     cout << '\t' << "User selected x maximum = " << geoac::x_max << '\n';}
        else if ((strncmp(inputs[i], "min_y=", 6) == 0) || (strncmp(inputs[i], "y_min=", 6) == 0)){         geoac::y_min = atof(inputs[i] + 6);     cout << '\t' << "User selected y minimum = " << geoac::y_min << '\n';}
        else if ((strncmp(inputs[i], "max_y=", 6) == 0) || (strncmp(inputs[i], "y_max=", 6) == 0)){         geoac::y_max = atof(inputs[i] + 6);     cout << '\t' << "User selected y maximum = " << geoac::y_max << '\n';}
        
        else if ((strncmp(inputs[i], "min_ds=", 7) == 0) || (strncmp(inputs[i], "ds_min=", 7) == 0)){       geoac::ds_min = atof(inputs[i] + 7);   cout << '\t' << "User selected ds minimum = " << geoac::ds_min << '\n';}
        else if ((strncmp(inputs[i], "max_ds=", 7) == 0) || (strncmp(inputs[i], "ds_max=", 7) == 0)){       geoac::ds_max = atof(inputs[i] + 7);   cout << '\t' << "User selected ds maximum = " << geoac::ds_max << '\n';}
        else if ((strncmp(inputs[i], "max_s=", 6) == 0) || (strncmp(inputs[i], "s_max=", 6) == 0)){         geoac::s_max = atof(inputs[i] + 6);    cout << '\t' << "User selected s maximum (path length between reflections) = " << geoac::s_max << '\n';}
        
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
        else {
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }
    z_rcvr = max(topo::z(x_rcvr, y_rcvr), z_rcvr);
    if(write_atmo)      geoac::write_prof("atmo.dat", x_rcvr, y_rcvr, geoac::phi);
    geoac::configure();
    
    cout << '\n' << "Parameter summary:" << '\n';
    cout << '\t' << "inclination: " << geoac::theta * (180.0 / Pi) << '\n';
    cout << '\t' << "azimuth: " << 90.0 - geoac::phi * (180.0 / Pi) << '\n';
    cout << '\t' << "bounces: " << bounces << '\n';
    cout << '\t' << "receiver location: " << x_rcvr << ", " << y_rcvr << ", " << z_rcvr << '\n';
    if(!geoac::is_topo){
        cout << '\t' << "ground elevation: " << topo::z0 << '\n';
    }
    cout << '\t' << "frequency: " << freq << '\n';
    cout << '\t' << "S&B atten coeff: " << atmo::tweak_abs << '\n' << '\n';
    
    // Define variables used for analysis
    double travel_time_sum, travel_time_var_incl, travel_time_var_az, attenuation, z_max, inclination, back_az;
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
    projection << "# x [km]";
    projection << '\t' << "y [km]";
    projection << '\t' << "z [km]";
    projection << '\t' << "time [s]";
    projection << '\t' << "X^{(incl)} [km/deg]";
    projection << '\t' << "Y^{(incl)} [km/deg]";
    projection << '\t' << "Z^{(incl)} [km/deg]";
    projection << '\t' << "T^{(incl)} [s/deg]";
    projection << '\t' << "X^{(az)} [km/deg]";
    projection << '\t' << "Y^{(az)} [km/deg]";
    projection << '\t' << "Z^{(az)} [km/deg]";
    projection << '\t' << "T^{(az)} [s/deg]";
    projection << '\n';
    
    double** solution;
    geoac::build_solution(solution, length);

    geoac::set_initial(solution, x_rcvr, y_rcvr, z_rcvr);
    travel_time_sum = 0.0;
    travel_time_var_incl = 0.0;
    travel_time_var_az = 0.0;
    
    attenuation = 0.0;
    z_max = 0.0;
    
    for(int bnc_cnt = 0; bnc_cnt <= bounces; bnc_cnt++){
        k = geoac::prop_rk4(solution, break_check, length);
        
        for(int m = 1; m < k ; m++){
            geoac::travel_time_var(travel_time_sum, travel_time_var_incl, travel_time_var_az, solution, m - 1, m);
            geoac::atten(attenuation, solution, m - 1, m, freq);
            z_max = max (z_max, solution[m][2]);
            
            if(m == 1 || m % 5 == 0){
                projection << solution[m][0];
                projection << '\t' << solution[m][1];
                projection << '\t' << max(solution[m][2], topo::z(solution[m][0], solution[m][1]));
                projection << '\t' << travel_time_sum;
                
                projection << '\t' << solution[m][6] * (Pi / 180.0);
                projection << '\t' << solution[m][7] * (Pi / 180.0);
                projection << '\t' << solution[m][8] * (Pi / 180.0);
                projection << '\t' << travel_time_var_incl  * (Pi / 180.0);

                projection << '\t' << solution[m][12] * (Pi / 180.0);
                projection << '\t' << solution[m][13] * (Pi / 180.0);
                projection << '\t' << solution[m][14] * (Pi / 180.0);
                projection << '\t' << travel_time_var_az * (Pi / 180.0);

                projection << '\n';
            }
        }
        
        if(break_check) break;
        geoac::set_refl(solution,k);
    }
    
    inclination = - asin(atmo::c(solution[k][0], solution[k][1], topo::z(solution[k][0], solution[k][1])) / atmo::c(x_rcvr, y_rcvr, z_rcvr) * solution[k][5]) * 180.0 / Pi;
    back_az = 90.0 - atan2(-solution[k][4], -solution[k][3]) * 180.0 / Pi;
    while(back_az < -180.0) back_az +=360.0;
    while(back_az >  180.0) back_az -=360.0;
    
    if(!break_check){
        cout << '\t' << "Back projected arrival:" << '\n';
        cout << '\t' << '\t' << "x [km E-W] = " << solution[k][0] << '\n';
        cout << '\t' << '\t' << "y [km N-S] = " << solution[k][1] << '\n';
        cout << '\t' << '\t' << "time [s] = " << travel_time_sum << '\n';
        cout << '\t' << '\t' << "celerity [km/s] = " << sqrt(pow(solution[k][0] - x_rcvr, 2) + pow(solution[k][1] - y_rcvr, 2)) / travel_time_sum << '\n';
        cout << '\t' << '\t' << "turning height [km] = " << z_max << '\n';
        cout << '\t' << '\t' << "arrival inclination [deg] = " << inclination << '\n';
        cout << '\t' << '\t' << "back azimuth = " << back_az  << '\n';
        cout << '\t' << '\t' << "trans. coeff. [dB] = " << 20.0 * log10(geoac::amp(solution, k)) << '\n';
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
    cout << '\t' << "##########################################" << '\n';
    cout << '\t' << "####        Running infraga-3d        ####" << '\n';
    cout << '\t' << "####          Eigenray Search         ####" << '\n';
    cout << '\t' << "##########################################" << '\n' << '\n';

    double src [3]   = {0.0, 0.0, 0.0};
    double rcvr [2] = {-250.0, 0.0};
    double theta_min = 0.5, theta_max = 45.0;
    bool write_atmo = false, custom_output_id=false;
    int bnc_min = 0, bnc_max = 0;
    int iterations = 25;
    double az_err_lim = 2.0;
    double freq = 0.1;
    char* prof_format = "zTuvdp";
    char* topo_file = "None";
    char* output_id;
    char input_check;

    topo::z0 = 0.0;
    atmo::tweak_abs = 1.0;
    geoac::calc_amp = true;
    geoac::is_topo = false;
    geoac::verbose = false;
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "prof_format=", 12) == 0){       prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){        topo::z0 = atof(inputs[i] + 7);}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){    topo_file = inputs[i] + 10; geoac::is_topo=true;}
    }
    
    if (count < 3){         cout << "You have to specify an atmosphere file..." << '\n'; return;}
    if (geoac::is_topo){    set_region(inputs[2], topo_file, prof_format, false);}
    else{                   set_region(inputs[2], prof_format, false);}
    
    if (geoac::is_topo){
        src[0] = (geoac::x_max + geoac::x_min) / 2.0;    rcvr[0] = src[0] - 250.0;
        src[1] = (geoac::y_max + geoac::y_min) / 2.0;    rcvr[1] = src[1];
    }
    
    for(int i = 3; i < count; i++){
        if ((strncmp(inputs[i], "incl_min=", 9) == 0) || (strncmp(inputs[i], "min_incl=", 9) == 0)){                    theta_min = atof(inputs[i] + 9);}
        else if ((strncmp(inputs[i], "incl_max=", 9) == 0) || (strncmp(inputs[i], "max_incl=", 9) == 0)){               theta_max = atof(inputs[i] + 9);}
        else if ((strncmp(inputs[i], "bnc_min=", 8) == 0) || (strncmp(inputs[i], "min_bnc=", 8) == 0)){                 bnc_min = atoi(inputs[i] + 8);}
        else if ((strncmp(inputs[i], "bnc_max=", 8) == 0) || (strncmp(inputs[i], "max_bnc=", 8) == 0)){                 bnc_max = atoi(inputs[i] + 8);}
        else if (strncmp(inputs[i], "bounces=", 8) == 0){                                                               bnc_min = atoi(inputs[i] + 8);
                                                                                                                        bnc_max = atoi(inputs[i] + 8);}
                 
        else if ((strncmp(inputs[i], "src_x=", 6) == 0) || (strncmp(inputs[i], "x_src=", 6) == 0)){                     src[0] = atof(inputs[i] + 6);}
        else if ((strncmp(inputs[i], "src_y=", 6) == 0) || (strncmp(inputs[i], "y_src_y=", 6) == 0)){                   src[1] = atof(inputs[i] + 6);}
        else if ((strncmp(inputs[i], "src_alt=", 8) == 0) || (strncmp(inputs[i], "alt_src=", 8) == 0)){                 src[2] = atof(inputs[i] + 8);}
        else if ((strncmp(inputs[i], "rcvr_x=", 7) == 0) || (strncmp(inputs[i], "x_rcvr=", 7) == 0)){                   rcvr[0] = atof(inputs[i] + 7);}
        else if ((strncmp(inputs[i], "rcvr_y=", 7) == 0) || (strncmp(inputs[i], "y_rcvr=", 7) == 0)){                   rcvr[1] = atof(inputs[i] + 7);}
        
        else if (strncmp(inputs[i], "verbose=", 8) == 0){                                                               geoac::verbose = string2bool(inputs[i] + 8);}
        else if (strncmp(inputs[i], "iterations=", 11) == 0){                                                           iterations = atof(inputs[i] + 11);}
        else if (strncmp(inputs[i], "damping=", 8) == 0){                                                               geoac::damping = atof(inputs[i] + 8);}
        else if (strncmp(inputs[i], "tolerance=", 10) == 0){                                                            geoac::tolerance = atof(inputs[i] + 10);}
        else if (strncmp(inputs[i], "az_dev_lim=", 11) == 0){                                                           az_err_lim = atof(inputs[i] + 11);}
        else if ((strncmp(inputs[i], "incl_step_min=", 14) == 0) || (strncmp(inputs[i], "min_incl_step=", 14) == 0)){   geoac::dth_sml = atof(inputs[i] + 14) * (Pi / 180.0);}
        else if ((strncmp(inputs[i], "incl_step_max=", 14) == 0) || (strncmp(inputs[i], "max_incl_step=", 14) == 0)){   geoac::dth_big = atof(inputs[i] + 14) * (Pi / 180.0);}
        
        else if (strncmp(inputs[i], "freq=", 5) == 0){                                                                  freq = atof(inputs[i] + 5);}
        else if (strncmp(inputs[i], "abs_coeff=", 10) == 0){                                                            atmo::tweak_abs = max(0.0, atof(inputs[i] + 10));}
        
        else if ((strncmp(inputs[i], "max_alt=", 8) == 0) || (strncmp(inputs[i], "alt_max=", 8) == 0)){                 geoac::alt_max = atof(inputs[i] + 8);   cout << '\t' << "User selected altitude maximum = " << geoac::alt_max << '\n';}
        else if ((strncmp(inputs[i], "max_rng=", 8) == 0) || (strncmp(inputs[i], "rng_max=", 8) == 0)){                 geoac::rng_max = atof(inputs[i] + 8);   cout << '\t' << "User selected range maximum = " << geoac::rng_max << '\n';}
        else if ((strncmp(inputs[i], "min_x=", 6) == 0) || (strncmp(inputs[i], "x_min=", 6) == 0)){                     geoac::x_min = atof(inputs[i] + 6);     cout << '\t' << "User selected x minimum = " << geoac::x_min << '\n';}
        else if ((strncmp(inputs[i], "max_x=", 6) == 0) || (strncmp(inputs[i], "x_max=", 6) == 0)){                     geoac::x_max = atof(inputs[i] + 6);     cout << '\t' << "User selected x maximum = " << geoac::x_max << '\n';}
        else if ((strncmp(inputs[i], "min_y=", 6) == 0) || (strncmp(inputs[i], "y_min=", 6) == 0)){                     geoac::y_min = atof(inputs[i] + 6);     cout << '\t' << "User selected y minimum = " << geoac::y_min << '\n';}
        else if ((strncmp(inputs[i], "max_y=", 6) == 0) || (strncmp(inputs[i], "y_max=", 6) == 0)){                     geoac::y_max = atof(inputs[i] + 6);     cout << '\t' << "User selected y maximum = " << geoac::y_max << '\n';}
        
        else if ((strncmp(inputs[i], "min_ds=", 7) == 0) || (strncmp(inputs[i], "ds_min=", 7) == 0)){                   geoac::ds_min = atof(inputs[i] + 7);   cout << '\t' << "User selected ds minimum = " << geoac::ds_min << '\n';}
        else if ((strncmp(inputs[i], "max_ds=", 7) == 0) || (strncmp(inputs[i], "ds_max=", 7) == 0)){                   geoac::ds_max = atof(inputs[i] + 7);   cout << '\t' << "User selected ds maximum = " << geoac::ds_max << '\n';}
        else if ((strncmp(inputs[i], "max_s=", 6) == 0) || (strncmp(inputs[i], "s_max=", 6) == 0)){                     geoac::s_max = atof(inputs[i] + 6);    cout << '\t' << "User selected s maximum (path length between reflections) = " << geoac::s_max << '\n';}

        else if (strncmp(inputs[i], "write_atmo=", 11) == 0){                                                           write_atmo = string2bool(inputs[i] + 11);}
        else if (strncmp(inputs[i], "prof_format=", 12) == 0){                                                          prof_format = inputs[i]+15;}
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
    src[2] = max(topo::z(src[0], src[1]), src[2]);
    if(write_atmo)  geoac::write_prof("atmo.dat", src[0], src[1], atan2(rcvr[1] - src[1], rcvr[0] - src[0]));

    theta_min *= Pi / 180.0;
    theta_max *= Pi / 180.0;
    
    cout << '\n' << "Parameter summary:" << '\n';
    cout << '\t' << "source location: " << src[0] << ", " << src[1] << ", " << src[2] << '\n';
    cout << '\t' << "receiver location: " << rcvr[0] << ", " << rcvr[1] << ", " << topo::z(rcvr[0], rcvr[1]) << '\n';
    cout << '\t' << "inclination range: " << theta_min * (180.0 / Pi) << ", " << theta_max * (180.0 / Pi) << '\n';
    cout << '\t' << "bounces: " << bnc_min << ", " << bnc_max << '\n';
    if(!geoac::is_topo){
        cout << '\t' << "ground elevation: " << topo::z0 << '\n';
    }
    cout << '\t' << "damping:" << geoac::damping << '\n';
    cout << '\t' << "frequency: " << freq << '\n';
    cout << '\t' << "S&B atten coeff: " << atmo::tweak_abs << '\n' << '\n';
    
    // Extract the file name from the input and use it to distinguish the output
	char output_buffer[512];
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
    
    geoac::eig_results << "# infraga-3d -eig_search summary:" << '\n';
    geoac::eig_results << "# " << '\t' << "profile: " << inputs[2] << '\n';
    geoac::eig_results << "# " << '\t' << "source location: " << src[0] << ", " << src[1] << ", " << src[2] << '\n';
    geoac::eig_results << "# " << '\t' << "receiver location: " << rcvr[0] << ", " << rcvr[1] << ", " << topo::z(rcvr[0], rcvr[1]) << '\n';
    geoac::eig_results << "# " << '\t' << "inclination range: " << theta_min * (180.0 / Pi) << ", " << theta_max * (180.0 / Pi) << '\n';
    geoac::eig_results << "# " << '\t' << "bounces: " << bnc_min << ", " << bnc_max << '\n';
    if(!geoac::is_topo){
        geoac::eig_results << "#" << '\t' << "ground elevation: " << topo::z0 << '\n';
    } else {
        geoac::eig_results << "#" << '\t' << "topo file:" << topo_file << '\n';
    }
    geoac::eig_results << "# " << '\t' << "damping: " << geoac::damping << '\n';
    geoac::eig_results << "# " << '\t' << "frequency: " << freq << '\n';
    geoac::eig_results << "# " << '\t' << "S&B atten coeff: " << atmo::tweak_abs << '\n' << '\n';
    
    geoac::eig_results << "# incl [deg]";
    geoac::eig_results << '\t' << "az [deg]";
    geoac::eig_results << '\t' << "n_b";
    geoac::eig_results << '\t' << "x_0 [km]";
    geoac::eig_results << '\t' << "y_0 [km]";
    geoac::eig_results << '\t' << "time [s]";
    geoac::eig_results << '\t' << "cel [km/s]";
    geoac::eig_results << '\t' << "turning ht [km]";
    geoac::eig_results << '\t' << "arrival incl [deg]";
    geoac::eig_results << '\t' << "back az [deg]";
    geoac::eig_results << '\t' << "trans. coeff. [dB]";
    geoac::eig_results << '\t' << "absorption [dB]";
    geoac::eig_results << '\n';
    
    for(int n_bnc = bnc_min; n_bnc <= bnc_max; n_bnc++){
        cout << "Searching for " << n_bnc << " bounce eigenray(s) with inclination between " << theta_min * 180.0 / Pi << " and " << theta_max * 180.0 / Pi << "." << '\n';
        
        theta_start = theta_min;
        while(theta_start < theta_max){
            estimate_success = geoac::est_eigenray(src, rcvr, theta_start, theta_max, theta_est, phi_est, theta_next, n_bnc, az_err_lim);
            if(estimate_success) geoac::find_eigenray(src, rcvr, theta_est, phi_est, freq, n_bnc, iterations, output_id);
            
            theta_start = theta_next;
        }
    }
    geoac::eig_results.close();
    cout << '\t' << "Identified " << geoac::eigenray_cnt << " eigenray(s)." << '\n';
    clear_region();
}

void run_eig_direct(char* inputs[], int count){
    cout << '\n';
    cout << '\t' << "##########################################" << '\n';
    cout << '\t' << "####        Running infraga-3d        ####" << '\n';
    cout << '\t' << "####          Eigenray Direct         ####" << '\n';
    cout << '\t' << "##########################################" << '\n' << '\n';
 
    double src [3]  = {0.0, 0.0, 0.0};
    double rcvr [2] = {-250.0, 0.0};
    double theta_est = 10.0 * (Pi / 180.0), phi_est = 45.0;
    int bounces = 0, iterations = 25;
    double freq = 0.1;
    char* prof_format = "zTuvdp";
    char* topo_file = "None";
    char* output_id;
    char input_check;
    bool write_atmo = false, custom_output_id=false;

    topo::z0 = 0.0;
    atmo::tweak_abs = 1.0;
    geoac::calc_amp = true;
    geoac::is_topo = false;
    geoac::verbose = true;
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "prof_format=", 12) == 0){       prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){        topo::z0 = atof(inputs[i] + 7);}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){    topo_file = inputs[i] + 10; geoac::is_topo=true;}
    }
    
    if (count < 3){         cout << "You have to specify an atmosphere file..." << '\n'; return;}
    if (geoac::is_topo){    set_region(inputs[2], topo_file, prof_format, false);}
    else{                   set_region(inputs[2], prof_format, false);}
    
    if (geoac::is_topo){
        src[0] = (geoac::x_max + geoac::x_min) / 2.0;    rcvr[0] = src[0] - 250.0;
        src[1] = (geoac::y_max + geoac::y_min) / 2.0;    rcvr[1] = src[1];
    }
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "incl_est=", 9) == 0){                                                       theta_est = atof(inputs[i] + 9) * (Pi / 180.0);}
        else if (strncmp(inputs[i], "az_est=", 7) == 0){                                                    phi_est = atof(inputs[i] + 7) * (Pi / 180.0);}
        else if (strncmp(inputs[i], "bounces=", 8) == 0){                                                   bounces = atoi(inputs[i] + 8);}

        else if ((strncmp(inputs[i], "src_x=", 6) == 0) || (strncmp(inputs[i], "x_src=", 6) == 0)){         src[0] = atof(inputs[i] + 6);}
        else if ((strncmp(inputs[i], "src_y=", 6) == 0) || (strncmp(inputs[i], "y_src_y=", 6) == 0)){       src[1] = atof(inputs[i] + 6);}
        else if ((strncmp(inputs[i], "src_alt=", 8) == 0) || (strncmp(inputs[i], "alt_src=", 8) == 0)){     src[2] = atof(inputs[i] + 8);}
        else if ((strncmp(inputs[i], "rcvr_x=", 7) == 0) || (strncmp(inputs[i], "x_rcvr=", 7) == 0)){       rcvr[0] = atof(inputs[i] + 7);}
        else if ((strncmp(inputs[i], "rcvr_y=", 7) == 0) || (strncmp(inputs[i], "y_rcvr=", 7) == 0)){       rcvr[1] = atof(inputs[i] + 7);}
        
        else if (strncmp(inputs[i], "verbose=", 8) == 0){                                                   geoac::verbose = string2bool(inputs[i] + 8);}
        else if (strncmp(inputs[i], "iterations=", 11) == 0){                                               iterations = atof(inputs[i] + 11);}
        else if (strncmp(inputs[i], "damping=", 8) == 0){                                                   geoac::damping = atof(inputs[i] + 8);}
        else if (strncmp(inputs[i], "tolerance=", 10) == 0){                                                geoac::tolerance = atof(inputs[i] + 10);}
                 
        else if (strncmp(inputs[i], "freq=", 5) == 0){                                                      freq = atof(inputs[i] + 5);}
        else if (strncmp(inputs[i], "abs_coeff=", 10) == 0){                                                atmo::tweak_abs = max(0.0, atof(inputs[i] + 10));}
        
        else if ((strncmp(inputs[i], "max_alt=", 8) == 0) || (strncmp(inputs[i], "alt_max=", 8) == 0)){     geoac::alt_max = atof(inputs[i] + 8);   cout << '\t' << "User selected altitude maximum = " << geoac::alt_max << '\n';}
        else if ((strncmp(inputs[i], "max_rng=", 8) == 0) || (strncmp(inputs[i], "rng_max=", 8) == 0)){     geoac::rng_max = atof(inputs[i] + 8);   cout << '\t' << "User selected range maximum = " << geoac::rng_max << '\n';}
        else if ((strncmp(inputs[i], "min_x=", 6) == 0) || (strncmp(inputs[i], "x_min=", 6) == 0)){         geoac::x_min = atof(inputs[i] + 6);     cout << '\t' << "User selected x minimum = " << geoac::x_min << '\n';}
        else if ((strncmp(inputs[i], "max_x=", 6) == 0) || (strncmp(inputs[i], "x_max=", 6) == 0)){         geoac::x_max = atof(inputs[i] + 6);     cout << '\t' << "User selected x maximum = " << geoac::x_max << '\n';}
        else if ((strncmp(inputs[i], "min_y=", 6) == 0) || (strncmp(inputs[i], "y_min=", 6) == 0)){         geoac::y_min = atof(inputs[i] + 6);     cout << '\t' << "User selected y minimum = " << geoac::y_min << '\n';}
        else if ((strncmp(inputs[i], "max_y=", 6) == 0) || (strncmp(inputs[i], "y_max=", 6) == 0)){         geoac::y_max = atof(inputs[i] + 6);     cout << '\t' << "User selected y maximum = " << geoac::y_max << '\n';}
        
        else if ((strncmp(inputs[i], "min_ds=", 7) == 0) || (strncmp(inputs[i], "ds_min=", 7) == 0)){       geoac::ds_min = atof(inputs[i] + 7);   cout << '\t' << "User selected ds minimum = " << geoac::ds_min << '\n';}
        else if ((strncmp(inputs[i], "max_ds=", 7) == 0) || (strncmp(inputs[i], "ds_max=", 7) == 0)){       geoac::ds_max = atof(inputs[i] + 7);   cout << '\t' << "User selected ds maximum = " << geoac::ds_max << '\n';}
        else if ((strncmp(inputs[i], "max_s=", 6) == 0) || (strncmp(inputs[i], "s_max=", 6) == 0)){         geoac::s_max = atof(inputs[i] + 6);    cout << '\t' << "User selected s maximum (path length between reflections) = " << geoac::s_max << '\n';}

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
    src[2] = max(topo::z(src[0], src[1]), src[2]);

    phi_est = atan2(rcvr[1] - src[1], rcvr[0] - src[0]);
    for(int i = 3; i < count; i++){
        if((strncmp(inputs[i], "az_est=", 7) == 0) || (strncmp(inputs[i], "est_az=", 7) == 0)){
            phi_est = Pi / 2.0 - atof(inputs[i] + 7) * (Pi / 180.0);
        }
    }
    if(write_atmo)  geoac::write_prof("atmo.dat", src[0], src[1], phi_est);

    cout << '\n' << "Parameter summary:" << '\n';
    cout << '\t' << "source location: " << src[0] << ", " << src[1] << ", " << src[2] << '\n';
    cout << '\t' << "receiver location: " << rcvr[0] << ", " << rcvr[1] << ", " << topo::z(rcvr[0], rcvr[1]) << '\n';
    cout << '\t' << "inclination estimate: " << theta_est * (180.0 / Pi) << '\n';
    cout << '\t' << "azimuth estimate: " << 90.0 - phi_est * 180.0 / Pi << '\n';
    cout << '\t' << "bounces: " << bounces << '\n';
    if(!geoac::is_topo){
        cout << '\t' << "ground elevation: " << topo::z0 << '\n';
    }
    cout << '\t' << "damping:" << geoac::damping << '\n';
    cout << '\t' << "frequency: " << freq << '\n';
    cout << '\t' << "S&B atten coeff: " << atmo::tweak_abs << '\n' << '\n';
    
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
    cout << '\t' << "##########################################" << '\n';
    cout << '\t' << "####        Running infraga-3d        ####" << '\n';
    cout << '\t' << "####    Weakly Non-Linear Waveform    ####" << '\n';
    cout << '\t' << "##########################################" << '\n' << '\n';
    
    double x_src = 0.0, y_src = 0.0, z_src = 0.0;
    double freq = 0.1, D, D_prev;
    int bounces = 0;
    bool write_atmo=false, write_rays=false, custom_output_id=false;
    char* prof_format = "zTuvdp";
    char* topo_file = "None";
    char* output_id;
    char input_check;

    double wvfrm_ref=1.0, wvfrm_out_step=1.0e10;
    double ray_length, c0, nu0, cg0x, cg0y, cg0z, cg0, rho0, D0, p0;
    double** wvfrm_array;
    int wvfrm_ref_k = 0;
    char* wvfrm_file = "none";
    char* wvfrm_opt = "impulse";
    ofstream wvfrm_out;

    geoac::theta = 15.0 * (Pi / 180.0);
    geoac::phi = Pi / 2.0 - 90.0 * (Pi / 180.0);

    topo::z0 = 0.0;
    atmo::tweak_abs = 1.0;
    geoac::calc_amp = true;
    geoac::is_topo = false;
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "prof_format=", 12) == 0){       prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){        topo::z0 = atof(inputs[i] + 7);}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){    topo_file = inputs[i] + 10; geoac::is_topo=true;}
    }
    
    if (count < 3){         cout << "You have to specify an atmosphere file..." << '\n'; return;}
    if (geoac::is_topo){    set_region(inputs[2], topo_file, prof_format, false);}
    else{                   set_region(inputs[2], prof_format, false);}

    if (geoac::is_topo){
        x_src = (geoac::x_max + geoac::x_min) / 2.0;
        y_src = (geoac::y_max + geoac::y_min) / 2.0;
    }
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "inclination=", 12) == 0){                                                   geoac::theta = atof(inputs[i] + 12) * (Pi / 180.0);}
        else if (strncmp(inputs[i], "azimuth=", 8) == 0){                                                   geoac::phi = Pi / 2.0 - atof(inputs[i] + 8) * (Pi / 180.0);}
        else if (strncmp(inputs[i], "bounces=", 8) == 0){                                                   bounces = max(0, atoi(inputs[i] + 8));}

        else if ((strncmp(inputs[i], "src_x=", 6) == 0) || (strncmp(inputs[i], "x_src=", 6) == 0)){         x_src = atof(inputs[i] + 6);}
        else if ((strncmp(inputs[i], "src_y=", 6) == 0) || (strncmp(inputs[i], "y_src_y=", 6) == 0)){       y_src = atof(inputs[i] + 6);}
        else if ((strncmp(inputs[i], "src_alt=", 8) == 0) || (strncmp(inputs[i], "alt_src=", 8) == 0)){     z_src = atof(inputs[i] + 8);}
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
        else if (strncmp(inputs[i], "abs_coeff=", 10) == 0){                                                atmo::tweak_abs = max(0.0, atof(inputs[i] + 10));}

        else if ((strncmp(inputs[i], "max_alt=", 8) == 0) || (strncmp(inputs[i], "alt_max=", 8) == 0)){     geoac::alt_max = atof(inputs[i] + 8);   cout << '\t' << "User selected altitude maximum = " << geoac::alt_max << '\n';}
        else if ((strncmp(inputs[i], "max_rng=", 8) == 0) || (strncmp(inputs[i], "rng_max=", 8) == 0)){     geoac::rng_max = atof(inputs[i] + 8);   cout << '\t' << "User selected range maximum = " << geoac::rng_max << '\n';}
        else if ((strncmp(inputs[i], "min_x=", 6) == 0) || (strncmp(inputs[i], "x_min=", 6) == 0)){         geoac::x_min = atof(inputs[i] + 6);     cout << '\t' << "User selected x minimum = " << geoac::x_min << '\n';}
        else if ((strncmp(inputs[i], "max_x=", 6) == 0) || (strncmp(inputs[i], "x_max=", 6) == 0)){         geoac::x_max = atof(inputs[i] + 6);     cout << '\t' << "User selected x maximum = " << geoac::x_max << '\n';}
        else if ((strncmp(inputs[i], "min_y=", 6) == 0) || (strncmp(inputs[i], "y_min=", 6) == 0)){         geoac::y_min = atof(inputs[i] + 6);     cout << '\t' << "User selected y minimum = " << geoac::y_min << '\n';}
        else if ((strncmp(inputs[i], "max_y=", 6) == 0) || (strncmp(inputs[i], "y_max=", 6) == 0)){         geoac::y_max = atof(inputs[i] + 6);     cout << '\t' << "User selected y maximum = " << geoac::y_max << '\n';}
        
        else if ((strncmp(inputs[i], "min_ds=", 7) == 0) || (strncmp(inputs[i], "ds_min=", 7) == 0)){       geoac::ds_min = atof(inputs[i] + 7);   cout << '\t' << "User selected ds minimum = " << geoac::ds_min << '\n';}
        else if ((strncmp(inputs[i], "max_ds=", 7) == 0) || (strncmp(inputs[i], "ds_max=", 7) == 0)){       geoac::ds_max = atof(inputs[i] + 7);   cout << '\t' << "User selected ds maximum = " << geoac::ds_max << '\n';}
        else if ((strncmp(inputs[i], "max_s=", 6) == 0) || (strncmp(inputs[i], "s_max=", 6) == 0)){         geoac::s_max = atof(inputs[i] + 6);    cout << '\t' << "User selected s maximum (path length between reflections) = " << geoac::s_max << '\n';}

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
        else {
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }
    z_src = max(topo::z(x_src, y_src), z_src);
    if(write_atmo)  geoac::write_prof("atmo.dat", x_src, y_src, geoac::phi);
    geoac::configure();
    
    cout << '\n' << "Parameter summary:" << '\n';
    cout << '\t' << "inclination: " << geoac::theta * (180.0 / Pi) << '\n';
    cout << '\t' << "azimuth: " << 90.0 - geoac::phi * (180.0 / Pi) << '\n';
    cout << '\t' << "bounces: " << bounces << '\n';
    cout << '\t' << "source location: " << x_src << ", " << y_src << ", " << z_src << '\n';
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
    double travel_time_sum, attenuation, z_max, inclination, back_az;
	int k, length = int(geoac::s_max / geoac::ds_min);
    bool break_check;
    
    double** solution;
    geoac::build_solution(solution,length);

    ofstream raypath;
    if(write_rays){
        sprintf(output_buffer, "%s.raypaths.dat", output_id);
        raypath.open(output_buffer);
        raypath << "# x [km]";
        raypath << '\t' << "y [km]";
        raypath << '\t' << "z [km]";
        raypath << '\t' << "trans. coeff. [dB]";
        raypath << '\t' << "absorption [dB]";
        raypath << '\t' << "time [s]";
        // raypath << '\t' << "density [g/cm^3]";
        // raypath << '\t' << "snd spd [km/s]";
        // raypath << '\t' << "u [km/s]";
        // raypath << '\t' << "v [km/s]";
        // raypath << '\t' << "w [km/s]";
        // raypath << '\t' << "s_x [s/km]";
        // raypath << '\t' << "s_y [s/km]";
        // raypath << '\t' << "s_z [s/km]";
        // raypath << '\t' << "D [km^2/rad^2]";
        raypath << '\n';
    }

    cout << "Calculating ray path geometry and weakly non-linear waveform evolution..." << '\n';
    if (strncmp(wvfrm_file, "none", 4) != 0){   wvfrm::load_wvfrm(wvfrm_array, wvfrm_file);}
    else {                                      wvfrm::build_wvfrm(wvfrm_array, wvfrm_opt);}
        
    geoac::set_initial(solution, x_src, y_src, z_src);
    travel_time_sum = 0.0;
    attenuation = 0.0;
    z_max = 0.0;
		
    for(int bnc_cnt = 0; bnc_cnt <= bounces; bnc_cnt++){
        k = geoac::prop_rk4(solution, break_check, length);

        for(int m = 1; m < k ; m++){
            geoac::travel_time(travel_time_sum, solution, m - 1,m);
            geoac::atten(attenuation, solution, m - 1, m, freq);
            z_max = max (z_max, solution[m][2]);

            if(write_rays && (m == 1 || m % 15 == 0)){
                raypath << solution[m][0];
                raypath << '\t' << solution[m][1];
                raypath << '\t' << max(solution[m][2], topo::z(solution[m][0], solution[m][1]));
                raypath << '\t' << 20.0 * log10(geoac::amp(solution, m));
                raypath << '\t' << -2.0 * attenuation;
                raypath << '\t' << travel_time_sum;

                // Uncomment these to output the necessary items to compare with NCPAprop wnl_wvfrm methods
                // raypath << '\t' << atmo::rho(solution[m][0], solution[m][1], solution[m][2]);
                // raypath << '\t' << atmo::c(solution[m][0], solution[m][1], solution[m][2]);
                // raypath << '\t' << atmo::u(solution[m][0], solution[m][1], solution[m][2]);
                // raypath << '\t' << atmo::v(solution[m][0], solution[m][1], solution[m][2]);
                // raypath << '\t' << atmo::w(solution[m][0], solution[m][1], solution[m][2]);
                // raypath << '\t' << solution[m][3] / atmo::c(solution[0][0], solution[0][1], solution[0][2]);
                // raypath << '\t' << solution[m][4] / atmo::c(solution[0][0], solution[0][1], solution[0][2]);
                // raypath << '\t' << solution[m][5] / atmo::c(solution[0][0], solution[0][1], solution[0][2]);
                // raypath << '\t' << geoac::jacobian(solution, m);
                raypath << '\n';
            }
        }
        if(bnc_cnt == 0){
            wvfrm_ref_k = 0;
            ray_length = 0.0;
            for(wvfrm_ref_k = 0; wvfrm_ref_k < k; wvfrm_ref_k++){
                ray_length += sqrt(pow(solution[wvfrm_ref_k + 1][0] - solution[wvfrm_ref_k][0], 2) + pow(solution[wvfrm_ref_k + 1][1] - solution[wvfrm_ref_k][1], 2) + pow(solution[wvfrm_ref_k + 1][2] - solution[wvfrm_ref_k][2], 2));
                if (ray_length >= wvfrm_ref) break;
            }

            c0 = atmo::c(solution[wvfrm_ref_k + 1][0], solution[wvfrm_ref_k + 1][1], solution[wvfrm_ref_k + 1][2]);
            rho0 = atmo::rho(solution[wvfrm_ref_k + 1][0], solution[wvfrm_ref_k + 1][1], solution[wvfrm_ref_k + 1][2]);

            nu0 = sqrt(pow(solution[wvfrm_ref_k + 1][3], 2) + pow(solution[wvfrm_ref_k + 1][4], 2) + pow(solution[wvfrm_ref_k + 1][5], 2));
            cg0x = c0 * solution[wvfrm_ref_k + 1][3] / nu0 + atmo::u(solution[wvfrm_ref_k + 1][0], solution[wvfrm_ref_k + 1][1], solution[wvfrm_ref_k + 1][2]);
            cg0y = c0 * solution[wvfrm_ref_k + 1][4] / nu0 + atmo::v(solution[wvfrm_ref_k + 1][0], solution[wvfrm_ref_k + 1][1], solution[wvfrm_ref_k + 1][2]);
            cg0z = c0 * solution[wvfrm_ref_k + 1][5] / nu0 + atmo::w(solution[wvfrm_ref_k + 1][0], solution[wvfrm_ref_k + 1][1], solution[wvfrm_ref_k + 1][2]);
            cg0 = sqrt(pow(cg0x, 2) + pow(cg0y, 2) + pow(cg0z, 2));

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

            ray_length = geoac::wnl_wvfrm(solution, wvfrm_array, wvfrm_ref_k + 1, k, 0.0, c0, cg0, nu0, rho0, D0, p0, wvfrm_out_step);
        } else {
            ray_length += geoac::wnl_wvfrm(solution, wvfrm_array, 0, k, ray_length, c0, cg0, nu0, rho0, D0, p0, wvfrm_out_step);
        }

        if(break_check) break;
        geoac::set_refl(solution,k);
    }
        
    inclination = - asin(atmo::c(solution[k][0], solution[k][1], topo::z(solution[k][0], solution[k][1])) / atmo::c(x_src, y_src, z_src) * solution[k][5]) * 180.0 / Pi;
    back_az = 90.0 - atan2(-solution[k][4], -solution[k][3]) * 180.0 / Pi;
    while(back_az < -180.0) back_az +=360.0;
    while(back_az >  180.0) back_az -=360.0;
        
    if(!break_check){
        cout << '\n';
        cout << '\t' << "Arrival summary:" << '\n';
        cout << '\t' << '\t' << "x [km E-W] = " << solution[k][0] << '\n';
        cout << '\t' << '\t' << "y [km N-S] = " << solution[k][1] << '\n';
        cout << '\t' << '\t' << "time [s] = " << travel_time_sum << '\n';
        cout << '\t' << '\t' << "celerity [km/s] = " << sqrt(pow(solution[k][0] - x_src, 2) + pow(solution[k][1] - y_src, 2)) / travel_time_sum << '\n';
        cout << '\t' << '\t' << "turning height [km] = " << z_max << '\n';
        cout << '\t' << '\t' << "arrival inclination [deg] = " << inclination << '\n';
        cout << '\t' << '\t' << "back azimuth = " << back_az  << '\n';
        cout << '\t' << '\t' << "attenuation (geometric) [dB] = " << 20.0 * log10(geoac::amp(solution, k)) << '\n';
        cout << '\t' << '\t' << "absorption [dB] = " << -2.0 * attenuation << '\n' << '\n';
    } else {
        cout << '\n';
        cout << '\t' << '\t' << "Ray path does not return to the ground." << '\n' << '\n';
    }

    double c = atmo::c(solution[k][0], solution[k][1], solution[k][2]);
    double rho = atmo::rho(solution[k][0], solution[k][1], solution[k][2]);
    double nu = sqrt(pow(solution[k][3], 2) + pow(solution[k][4], 2) + pow(solution[k][5], 2));

    double cgx = c * solution[k][3] / nu + atmo::u(solution[k][0], solution[k][1], solution[k][2]);
    double cgy = c * solution[k][4] / nu + atmo::v(solution[k][0], solution[k][1], solution[k][2]);
    double cgz = c * solution[k][5] / nu + atmo::w(solution[k][0], solution[k][1], solution[k][2]);
    double cg = sqrt(pow(cgx, 2) + pow(cgy, 2) + pow(cgz, 2));

    sprintf(output_buffer, "%s.wvfrm_out.dat", output_id);
    wvfrm_out.open(output_buffer);

    wvfrm_out << "# infraga-3d -wnl_wvfrm summary:" << '\n';
    wvfrm_out << "#" << '\t' << "profile: " << inputs[2] << '\n';
    wvfrm_out << "#" << '\t' << "inclination: " << geoac::theta * (180.0 / Pi)  << '\n';
    wvfrm_out << "#" << '\t' << "azimuth: " << geoac::phi * (180.0 / Pi)  << '\n';
    wvfrm_out << "#" << '\t' << "bounces: " << bounces << '\n';
    wvfrm_out << "#" << '\t' << "source location: " << x_src << ", " << y_src << ", " << z_src << '\n';
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
        wvfrm_out << '\t' << wvfrm_array[n][1] * p0 * sqrt(fabs((rho * pow(c, 3) * nu) / (rho0 * pow(c0, 3) * nu0) * (cg0 * D0) / (cg * geoac::jacobian(solution, k))));
        wvfrm_out << '\n';
    }
    wvfrm_out.close();

    if(write_rays){
        raypath.close();
    }
    wvfrm::delete_wvfrm(wvfrm_array);

    geoac::delete_solution(solution, length);
    clear_region();
}


void run_region_test(char* inputs[], int count){
    cout << '\n';
    cout << '\t' << "##########################################" << '\n';
    cout << '\t' << "####        Running infraga-3d        ####" << '\n';
    cout << '\t' << "####            Region Test           ####" << '\n';
    cout << '\t' << "##########################################" << '\n' << '\n';
    
    ofstream file_out;
    double x_src=0.0, y_src=0.0, alt_max=129.0, eps=1.0e-3;
    
    char* prof_format = "zTuvdp";
    char* topo_file = "None";
    
    topo::z0 = 0.0;
    atmo::tweak_abs = 1.0;
    geoac::calc_amp = true;
    geoac::is_topo = false;
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "prof_format=",12) == 0){        prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){        topo::z0 = atof(inputs[i] + 7);}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){    topo_file = inputs[i] + 10; geoac::is_topo=true;}
        else if (strncmp(inputs[i], "eps=", 4) == 0){           eps = atof(inputs[i] + 4);}  
        else if (strncmp(inputs[i], "alt_max=", 8) == 0){       alt_max = atof(inputs[i] + 8);}
        else if (strncmp(inputs[i], "topo_use_BLw=", 13) == 0){ topo::use_BLw = string2bool(inputs[i] + 13);}
    }
    
    if (count < 3){         cout << "You have to specify an atmosphere file..." << '\n'; return;}
    if (geoac::is_topo){    set_region(inputs[2], topo_file, prof_format, false);}
    else{                   set_region(inputs[2], prof_format, false);}
    
    if (geoac::is_topo){
        x_src = (geoac::x_max + geoac::x_min) / 2.0;
        y_src = (geoac::y_max + geoac::y_min) / 2.0;
    }
    
    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "x_src=", 6) == 0){              x_src = atof(inputs[i] + 6);}
        else if (strncmp(inputs[i], "y_src=", 6) == 0){         y_src = atof(inputs[i] + 6);}
    }   

    // Write the sound speed and its derivatives (analytic and defined)
    cout << '\n' << "Running region test comparing defined derivatives to finite differences (make sure you've created 'region_tests/3D/' directory)..." << '\n';
    cout << '\t' << "Evaluating sound speed..." << '\n';
    file_out.open("region_tests/3D/c.analytic.dat");
    for(int m = 0; m < int(alt_max * 100); m++){
        double z0 = topo::z(x_src, y_src) + m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;
        
        file_out << z0 << '\t';
        file_out << atmo::c(x_src, y_src, z0) << '\t';
        file_out << atmo::dc(x_src, y_src, z0, 0) << '\t';
        file_out << atmo::dc(x_src, y_src, z0, 1) << '\t';
        file_out << atmo::dc(x_src, y_src, z0, 2) << '\t';
        file_out << atmo::ddc(x_src, y_src, z0, 0, 0) << '\t';
        file_out << atmo::ddc(x_src, y_src, z0, 1, 1) << '\t';
        file_out << atmo::ddc(x_src, y_src, z0, 2, 2) << '\t';
        file_out << atmo::ddc(x_src, y_src, z0, 0, 1) << '\t';
        file_out << atmo::ddc(x_src, y_src, z0, 0, 2) << '\t';
        file_out << atmo::ddc(x_src, y_src, z0, 1, 2) << '\t';
        file_out << '\n';
    }
    file_out.close();
    
    file_out.open("region_tests/3D/c.finite.dat");
     for(int m = 1; m < int(alt_max * 100); m++){
       double z0 = topo::z(x_src, y_src) + m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;
        
        file_out << z0 << '\t';
        file_out << atmo::c(x_src, y_src, z0) << '\t';
        file_out << (atmo::c(x_src + eps, y_src, z0) - atmo::c(x_src - eps, y_src, z0)) / (2.0 * eps) << '\t';
        file_out << (atmo::c(x_src, y_src + eps, z0) - atmo::c(x_src, y_src - eps, z0)) / (2.0 * eps) << '\t';
        file_out << (atmo::c(x_src, y_src, z0 + eps) - atmo::c(x_src, y_src, z0 - eps)) / (2.0 * eps) << '\t';
        
        file_out << (atmo::c(x_src + eps, y_src, z0) + atmo::c(x_src - eps, y_src, z0) - 2.0 * atmo::c(x_src, y_src, z0)) / pow(eps, 2.0) << '\t';
        file_out << (atmo::c(x_src, y_src + eps, z0) + atmo::c(x_src, y_src - eps, z0) - 2.0 * atmo::c(x_src, y_src, z0)) / pow(eps, 2.0) << '\t';
        file_out << (atmo::c(x_src, y_src, z0 + eps) + atmo::c(x_src, y_src, z0 - eps) - 2.0 * atmo::c(x_src, y_src, z0)) / pow(eps, 2.0) << '\t';

        file_out << (atmo::c(x_src + eps, y_src + eps, z0) + atmo::c(x_src - eps, y_src - eps, z0) - atmo::c(x_src + eps, y_src - eps, z0) - atmo::c(x_src - eps, y_src + eps, z0)) / pow(2.0 * eps, 2.0) << '\t';
        file_out << (atmo::c(x_src + eps, y_src, z0 + eps) + atmo::c(x_src - eps, y_src, z0 - eps) - atmo::c(x_src + eps, y_src, z0 - eps) - atmo::c(x_src - eps, y_src, z0 + eps)) / pow(2.0 * eps, 2.0) << '\t';
        file_out << (atmo::c(x_src, y_src + eps, z0 + eps) + atmo::c(x_src, y_src - eps, z0 - eps) - atmo::c(x_src, y_src + eps, z0 - eps) - atmo::c(x_src, y_src - eps, z0 + eps)) / pow(2.0 * eps, 2.0) << '\t';
        
        file_out << '\n';
    }
    file_out.close();

    // Write the zonal wind component and its derivatives
    cout << '\t' << "Evaluating winds (u)..." << '\n';
    file_out.open("region_tests/3D/u.analytic.dat");
    for(int m = 0; m < int(alt_max * 100); m++){
        double z0 = topo::z(x_src, y_src) + m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;
        
        file_out << z0 << '\t';
        file_out << atmo::u(x_src, y_src, z0) << '\t';
        file_out << atmo::du(x_src, y_src, z0, 0) << '\t';
        file_out << atmo::du(x_src, y_src, z0, 1) << '\t';
        file_out << atmo::du(x_src, y_src, z0, 2) << '\t';
        file_out << atmo::ddu(x_src, y_src, z0, 0, 0) << '\t';
        file_out << atmo::ddu(x_src, y_src, z0, 1, 1) << '\t';
        file_out << atmo::ddu(x_src, y_src, z0, 2, 2) << '\t';
        file_out << atmo::ddu(x_src, y_src, z0, 0, 1) << '\t';
        file_out << atmo::ddu(x_src, y_src, z0, 0, 2) << '\t';
        file_out << atmo::ddu(x_src, y_src, z0, 1, 2) << '\t';
        file_out << '\n';
    }
    file_out.close();
    
    file_out.open("region_tests/3D/u.finite.dat");
    for(int m = 1; m < int(alt_max * 100); m++){
        double z0 = topo::z(x_src, y_src) + m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;
        
        file_out << z0 << '\t';
        file_out << atmo::u(x_src, y_src, z0) << '\t';
        file_out << (atmo::u(x_src + eps, y_src, z0) - atmo::u(x_src - eps, y_src, z0)) / (2.0 * eps) << '\t';
        file_out << (atmo::u(x_src, y_src + eps, z0) - atmo::u(x_src, y_src - eps, z0)) / (2.0 * eps) << '\t';
        file_out << (atmo::u(x_src, y_src, z0 + eps) - atmo::u(x_src, y_src, z0 - eps)) / (2.0 * eps) << '\t';
        
        file_out << (atmo::u(x_src + eps, y_src, z0) + atmo::u(x_src - eps, y_src, z0) - 2.0 * atmo::u(x_src, y_src, z0)) / pow(eps, 2.0) << '\t';
        file_out << (atmo::u(x_src, y_src + eps, z0) + atmo::u(x_src, y_src - eps, z0) - 2.0 * atmo::u(x_src, y_src, z0)) / pow(eps, 2.0) << '\t';
        file_out << (atmo::u(x_src, y_src, z0 + eps) + atmo::u(x_src, y_src, z0 - eps) - 2.0 * atmo::u(x_src, y_src, z0)) / pow(eps, 2.0) << '\t';
        
        file_out << (atmo::u(x_src + eps, y_src + eps, z0) + atmo::u(x_src - eps, y_src - eps, z0) - atmo::u(x_src + eps, y_src - eps, z0) - atmo::u(x_src - eps, y_src + eps, z0)) / pow(2.0 * eps, 2.0) << '\t';
        file_out << (atmo::u(x_src + eps, y_src, z0 + eps) + atmo::u(x_src - eps, y_src, z0 - eps) - atmo::u(x_src + eps, y_src, z0 - eps) - atmo::u(x_src - eps, y_src, z0 + eps)) / pow(2.0 * eps, 2.0) << '\t';
        file_out << (atmo::u(x_src, y_src + eps, z0 + eps) + atmo::u(x_src, y_src - eps, z0 - eps) - atmo::u(x_src, y_src + eps, z0 - eps) - atmo::u(x_src, y_src - eps, z0 + eps)) / pow(2.0 * eps, 2.0) << '\t';
        
        file_out << '\n';
    }
    file_out.close();

    file_out.open("region_tests/3D/u.calc_uvw.dat");
    for(int m = 0; m < int(alt_max * 100); m++){
        double z0 = topo::z(x_src, y_src) + m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;
        
        double u, v, w, du[3], dv[3], dw[3], ddu[6], ddv[6], ddw[6];
        atmo::calc_uvw(x_src, y_src, z0, u, v, w, du, dv, dw, ddu, ddv, ddw);    
        atmo::calc_uvw(x_src, y_src, z0, u, v, w, du, dv, dw);
            
        file_out << z0 << '\t';
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
    file_out.open("region_tests/3D/v.analytic.dat");
    for(int m = 0; m < int(alt_max * 100); m++){
        double z0 = topo::z(x_src, y_src) + m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;
        
        file_out << z0 << '\t';
        file_out << atmo::v(x_src, y_src, z0) << '\t';
        file_out << atmo::dv(x_src, y_src, z0, 0) << '\t';
        file_out << atmo::dv(x_src, y_src, z0, 1) << '\t';
        file_out << atmo::dv(x_src, y_src, z0, 2) << '\t';
        file_out << atmo::ddv(x_src, y_src, z0, 0, 0) << '\t';
        file_out << atmo::ddv(x_src, y_src, z0, 1, 1) << '\t';
        file_out << atmo::ddv(x_src, y_src, z0, 2, 2) << '\t';
        file_out << atmo::ddv(x_src, y_src, z0, 0, 1) << '\t';
        file_out << atmo::ddv(x_src, y_src, z0, 0, 2) << '\t';
        file_out << atmo::ddv(x_src, y_src, z0, 1, 2) << '\t';
        file_out << '\n';
    }
    file_out.close();
    
    file_out.open("region_tests/3D/v.finite.dat");
    for(int m = 1; m < int(alt_max * 100); m++){
        double z0 = topo::z(x_src, y_src) + m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;
        
        file_out << z0 << '\t';
        file_out << atmo::v(x_src, y_src, z0) << '\t';
        file_out << (atmo::v(x_src + eps, y_src, z0) - atmo::v(x_src - eps, y_src, z0)) / (2.0 * eps) << '\t';
        file_out << (atmo::v(x_src, y_src + eps, z0) - atmo::v(x_src, y_src - eps, z0)) / (2.0 * eps) << '\t';
        file_out << (atmo::v(x_src, y_src, z0 + eps) - atmo::v(x_src, y_src, z0 - eps)) / (2.0 * eps) << '\t';
        
        file_out << (atmo::v(x_src + eps, y_src, z0) + atmo::v(x_src - eps, y_src, z0) - 2.0 * atmo::v(x_src, y_src, z0)) / pow(eps, 2.0) << '\t';
        file_out << (atmo::v(x_src, y_src + eps, z0) + atmo::v(x_src, y_src - eps, z0) - 2.0 * atmo::v(x_src, y_src, z0)) / pow(eps, 2.0) << '\t';
        file_out << (atmo::v(x_src, y_src, z0 + eps) + atmo::v(x_src, y_src, z0 - eps) - 2.0 * atmo::v(x_src, y_src, z0)) / pow(eps, 2.0) << '\t';
        
        file_out << (atmo::v(x_src + eps, y_src + eps, z0) + atmo::v(x_src - eps, y_src - eps, z0) - atmo::v(x_src + eps, y_src - eps, z0) - atmo::v(x_src - eps, y_src + eps, z0)) / pow(2.0 * eps, 2.0) << '\t';
        file_out << (atmo::v(x_src + eps, y_src, z0 + eps) + atmo::v(x_src - eps, y_src, z0 - eps) - atmo::v(x_src + eps, y_src, z0 - eps) - atmo::v(x_src - eps, y_src, z0 + eps)) / pow(2.0 * eps, 2.0) << '\t';
        file_out << (atmo::v(x_src, y_src + eps, z0 + eps) + atmo::v(x_src, y_src - eps, z0 - eps) - atmo::v(x_src, y_src + eps, z0 - eps) - atmo::v(x_src, y_src - eps, z0 + eps)) / pow(2.0 * eps, 2.0) << '\t';
        
        file_out << '\n';
    }
    file_out.close();

    file_out.open("region_tests/3D/v.calc_uvw.dat");
    for(int m = 0; m < int(alt_max * 100); m++){
        double z0 = topo::z(x_src, y_src) + m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;
        
        double u, v, w, du[3], dv[3], dw[3], ddu[6], ddv[6], ddw[6];
        atmo::calc_uvw(x_src, y_src, z0, u, v, w, du, dv, dw, ddu, ddv, ddw);    
        atmo::calc_uvw(x_src, y_src, z0, u, v, w, du, dv, dw);
            
        file_out << z0 << '\t';
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
    file_out.open("region_tests/3D/w.analytic.dat");
    for(int m = 0; m < int(alt_max * 100); m++){
        double z0 = topo::z(x_src, y_src) + m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;

        file_out << z0 << '\t';
        file_out << atmo::w(x_src, y_src, z0) << '\t';
        file_out << atmo::dw(x_src, y_src, z0, 0) << '\t';
        file_out << atmo::dw(x_src, y_src, z0, 1) << '\t';
        file_out << atmo::dw(x_src, y_src, z0, 2) << '\t';
        file_out << atmo::ddw(x_src, y_src, z0, 0, 0) << '\t';
        file_out << atmo::ddw(x_src, y_src, z0, 1, 1) << '\t';
        file_out << atmo::ddw(x_src, y_src, z0, 2, 2) << '\t';
        file_out << atmo::ddw(x_src, y_src, z0, 0, 1) << '\t';
        file_out << atmo::ddw(x_src, y_src, z0, 0, 2) << '\t';
        file_out << atmo::ddw(x_src, y_src, z0, 1, 2) << '\t';
        file_out << '\n';
    }
    file_out.close();

    file_out.open("region_tests/3D/w.finite.dat");
    for(int m = 1; m < int(alt_max * 100); m++){
        double z0 = topo::z(x_src, y_src) + m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;

        file_out << z0 << '\t';
        file_out << atmo::w(x_src, y_src, z0) << '\t';
        file_out << (atmo::w(x_src + eps, y_src, z0) - atmo::w(x_src - eps, y_src, z0)) / (2.0 * eps) << '\t';
        file_out << (atmo::w(x_src, y_src + eps, z0) - atmo::w(x_src, y_src - eps, z0)) / (2.0 * eps) << '\t';
        file_out << (atmo::w(x_src, y_src, z0 + eps) - atmo::w(x_src, y_src, z0 - eps)) / (2.0 * eps) << '\t';

        file_out << (atmo::w(x_src + eps, y_src, z0) + atmo::w(x_src - eps, y_src, z0) - 2.0 * atmo::w(x_src, y_src, z0)) / pow(eps, 2.0) << '\t';
        file_out << (atmo::w(x_src, y_src + eps, z0) + atmo::w(x_src, y_src - eps, z0) - 2.0 * atmo::w(x_src, y_src, z0)) / pow(eps, 2.0) << '\t';
        file_out << (atmo::w(x_src, y_src, z0 + eps) + atmo::w(x_src, y_src, z0 - eps) - 2.0 * atmo::w(x_src, y_src, z0)) / pow(eps, 2.0) << '\t';

        file_out << (atmo::w(x_src + eps, y_src + eps, z0) + atmo::w(x_src - eps, y_src - eps, z0) - atmo::w(x_src + eps, y_src - eps, z0) - atmo::w(x_src - eps, y_src + eps, z0)) / pow(2.0 * eps, 2.0) << '\t';
        file_out << (atmo::w(x_src + eps, y_src, z0 + eps) + atmo::w(x_src - eps, y_src, z0 - eps) - atmo::w(x_src + eps, y_src, z0 - eps) - atmo::w(x_src - eps, y_src, z0 + eps)) / pow(2.0 * eps, 2.0) << '\t';
        file_out << (atmo::w(x_src, y_src + eps, z0 + eps) + atmo::w(x_src, y_src - eps, z0 - eps) - atmo::w(x_src, y_src + eps, z0 - eps) - atmo::w(x_src, y_src - eps, z0 + eps)) / pow(2.0 * eps, 2.0) << '\t';

        file_out << '\n';
    }
    file_out.close();
  
    file_out.open("region_tests/3D/w.calc_uvw.dat");
    for(int m = 0; m < int(alt_max * 100); m++){
        double z0 = topo::z(x_src, y_src) + m / 100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;
        
        double u, v, w, du[3], dv[3], dw[3], ddu[6], ddv[6], ddw[6];
        atmo::calc_uvw(x_src, y_src, z0, u, v, w, du, dv, dw, ddu, ddv, ddw);    
        atmo::calc_uvw(x_src, y_src, z0, u, v, w, du, dv, dw);
            
        file_out << z0 << '\t';
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
        eps=1.0e-1;

        cout << '\t' << '\t' << "field..." << '\n';
        file_out.open("region_tests/3D/topo.dat");
        for(int m1 = 0; m1 < topo::spline.length_x / 4; m1++){
            double x0 = topo::spline.x_vals[0] + m1 / (0.25 * topo::spline.length_x) * (topo::spline.x_vals[topo::spline.length_x - 1] - topo::spline.x_vals[0]);
            for(int m2 = 0; m2 < topo::spline.length_y / 4.0; m2++){
                double y0 = topo::spline.y_vals[0] + m2 / (0.25 * topo::spline.length_y) * (topo::spline.y_vals[topo::spline.length_y - 1] - topo::spline.y_vals[0]);

                double zg, dzg[2], ddzg[3];
                interp::eval_all(x0, y0, topo::spline, zg, dzg, ddzg);

                file_out << x0 << '\t';
                file_out << y0 << '\t';

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

        cout << '\t' << '\t' << "analytic..." << '\n';
        file_out.open("region_tests/3D/topo.x.analytic.dat");
        for(int m = 0; m < 2 * topo::spline.length_x; m++){
            double x0 = topo::spline.x_vals[0] + m / (2.0 * topo::spline.length_x) * (topo::spline.x_vals[topo::spline.length_x - 1] - topo::spline.x_vals[0]);
            double y0 = y_src;
        
            file_out << x0 << '\t';
            file_out << topo::z(x0, y0) << '\t';
            file_out << topo::dz(x0, y0, 0) << '\t';
            file_out << topo::dz(x0, y0, 1) << '\t';

            file_out << topo::ddz(x0, y0, 0, 0) << '\t';
            file_out << topo::ddz(x0, y0, 0, 1) << '\t';
            file_out << topo::ddz(x0, y0, 1, 1) << '\t';
            file_out << '\n';
        }
        file_out.close();

        file_out.open("region_tests/3D/topo.y.analytic.dat");
        for(int m = 0; m < 2 * topo::spline.length_y; m++){
            double x0 = x_src;
            double y0 = topo::spline.y_vals[0] + m / (2.0 * topo::spline.length_y) * (topo::spline.y_vals[topo::spline.length_y - 1] - topo::spline.y_vals[0]);
            
            file_out << y0 << '\t';
            file_out << topo::z(x0, y0) << '\t';
            file_out << topo::dz(x0, y0, 0) << '\t';
            file_out << topo::dz(x0, y0, 1) << '\t';
            
            file_out << topo::ddz(x0, y0, 0, 0) << '\t';
            file_out << topo::ddz(x0, y0, 0, 1) << '\t';
            file_out << topo::ddz(x0, y0, 1, 1) << '\t';
            file_out << '\n';
        }
        file_out.close();

        cout << '\t' << '\t' << "finite diff..." << '\n';        
        file_out.open("region_tests/3D/topo.x.finite.dat");
        for(int m = 0; m < 2 * topo::spline.length_x; m++){
            double x0 = topo::spline.x_vals[0] + m / (2.0 * topo::spline.length_x) * (topo::spline.x_vals[topo::spline.length_x - 1] - topo::spline.x_vals[0]);
            double y0 = y_src;
        
            file_out << x0 << '\t';
            file_out << topo::z(x0, y0) << '\t';
            file_out << (topo::z(x0 + eps, y0) - topo::z(x0 - eps, y0)) / (2.0 * eps) << '\t';
            file_out << (topo::z(x0, y0 + eps) - topo::z(x0, y0 - eps)) / (2.0 * eps) << '\t';
        
            file_out << (topo::z(x0 + eps, y0) + topo::z(x0 - eps, y0) - 2.0 * topo::z(x0, y0)) / pow(eps, 2.0) << '\t';
            file_out << (topo::z(x0 + eps, y0 + eps) + topo::z(x0 - eps, y0 - eps) - topo::z(x0 + eps, y0 - eps) - topo::z(x0 - eps, y0 + eps)) / pow(2.0 * eps, 2.0) << '\t';
            file_out << (topo::z(x0, y0 + eps) + topo::z(x0, y0 - eps) - 2.0 * topo::z(x0, y0)) / pow(eps, 2.0) << '\t';
        
            file_out << '\n';
        }
        file_out.close();

        file_out.open("region_tests/3D/topo.y.finite.dat");
        for(int m = 0; m < 2 * topo::spline.length_y; m++){
            double x0 = x_src;
            double y0 = topo::spline.y_vals[0] + m / (2.0 * topo::spline.length_y) * (topo::spline.y_vals[topo::spline.length_y - 1] - topo::spline.y_vals[0]);
            
            file_out << y0 << '\t';
            file_out << topo::z(x0, y0) << '\t';
            file_out << (topo::z(x0 + eps, y0) - topo::z(x0 - eps, y0)) / (2.0 * eps) << '\t';
            file_out << (topo::z(x0, y0 + eps) - topo::z(x0, y0 - eps)) / (2.0 * eps) << '\t';
            
            file_out << (topo::z(x0 + eps, y0) + topo::z(x0 - eps, y0) - 2.0 * topo::z(x0, y0)) / pow(eps, 2.0) << '\t';
            file_out << (topo::z(x0 + eps, y0 + eps) + topo::z(x0 - eps, y0 - eps) - topo::z(x0 + eps, y0 - eps) - topo::z(x0 - eps, y0 + eps)) / pow(2.0 * eps, 2.0) << '\t';
            file_out << (topo::z(x0, y0 + eps) + topo::z(x0, y0 - eps) - 2.0 * topo::z(x0, y0)) / pow(eps, 2.0) << '\t';
            
            file_out << '\n';
        }
        file_out.close();

        cout << '\t' << '\t' << "calc_all..." << '\n';
        file_out.open("region_tests/3D/topo.x.calc_all.dat");
        for(int m = 0; m < 2 * topo::spline.length_x; m++){
            double x0 = topo::spline.x_vals[0] + m / (2.0 * topo::spline.length_x) * (topo::spline.x_vals[topo::spline.length_x - 1] - topo::spline.x_vals[0]);
            double y0 = y_src;

            double zg, dzg[2], ddzg[3], dddzg[4];
            interp::eval_all(x0, y0, topo::spline, zg, dzg, ddzg, dddzg);

            file_out << x0 << '\t';
            file_out << zg << '\t';
            file_out << dzg[0] << '\t';
            file_out << dzg[1] << '\t';
            
            file_out << ddzg[0] << '\t';
            file_out << ddzg[1] << '\t';
            file_out << ddzg[2] << '\t';

            file_out << '\n';
        }
        file_out.close();

        file_out.open("region_tests/3D/topo.y.calc_all.dat");
        for(int m = 0; m < 2 * topo::spline.length_y; m++){
            double x0 = x_src;
            double y0 = topo::spline.y_vals[0] + m / (2.0 * topo::spline.length_y) * (topo::spline.y_vals[topo::spline.length_y - 1] - topo::spline.y_vals[0]);

            double zg, dzg[2], ddzg[3], dddzg[4];
            interp::eval_all(x0, y0, topo::spline, zg, dzg, ddzg, dddzg);

            file_out << y0 << '\t';
            file_out << zg << '\t';
            file_out << dzg[0] << '\t';
            file_out << dzg[1] << '\t';
            
            file_out << ddzg[0] << '\t';
            file_out << ddzg[1] << '\t';
            file_out << ddzg[2] << '\t';

            file_out << '\n';
        }
        file_out.close();


    }
    clear_region();
}

int main(int argc, char* argv[]){
    if (argc < 2){ usage();}
    else {
        if ((strncmp(argv[1], "--version", 9) == 0) || (strncmp(argv[1], "-v", 2) == 0)){       version();  return 0;}
        else if ((strncmp(argv[1], "--usage", 7) == 0) || (strncmp(argv[1], "-u", 2) == 0)){    usage();    return 0;}
        
        else if (strncmp(argv[1], "-prop", 5) == 0){                                            run_prop(argv, argc);           return 0;}
        else if (strncmp(argv[1], "-back_proj", 10) == 0){                                      run_back_proj(argv, argc);      return 0;}
        else if (strncmp(argv[1], "-eig_search", 11) == 0){                                     run_eig_search(argv, argc);     return 0;}
        else if (strncmp(argv[1], "-eig_direct", 11) == 0){                                     run_eig_direct(argv, argc);     return 0;}
        else if (strncmp(argv[1], "-wnl_wvfrm", 10) == 0){                                      run_wnl_wvfrm(argv, argc);      return 0;}

        // Test methods...development use only.
        else if (strncmp(argv[1], "-region_test", 12) == 0){                                    run_region_test(argv, argc);    return 0;}

        else {                                                                                  cout << "Unrecognized option." << '\n';  return 0;}
    }
    return 0;
}
