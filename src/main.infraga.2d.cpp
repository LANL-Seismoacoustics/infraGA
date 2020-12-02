#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <sstream>
#include <math.h>

#include "atmo/atmo_state.h"
#include "atmo/atmo_io.2d.h"

#include "geoac/geoac.params.h"
#include "geoac/geoac.eqset.h"
#include "geoac/geoac.interface.h"

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
    cout << '\t' << "########################################################" << '\n';
    cout << '\t' << "####                   infraga-2d                   ####" << '\n';
    cout << '\t' << "####       Two-Dimensional Planar Ray Tracing       ####" << '\n';
    cout << '\t' << "#### Using the Effective Sound Speed Approximation  ####" << '\n';
    cout << '\t' << "########################################################" << '\n' << '\n';
    
    cout << "Usage: infraga-2d [option] profile.met [parameters]" << '\n';
    cout << '\t' << '\t' << "Enter only 1 option." << '\n';
    cout << '\t' << '\t' << "Profile.met is expected to contain columns describing {z[km]  T[K]  u (zonal winds) [m/s]  v (meridional winds) [m/s]  density[g/cm^3]  p[mbar]} " << '\n';
    cout << '\t' << '\t' << '\t' << "Profile format can be modified, see manual document for details." << '\n';
    cout << '\t' << '\t' << "Parameter calls are expected using the format: parameter_name=value." << '\n' << '\n';
    
    cout << "Options and parameters are:" << '\n';
    cout << '\t' << "-prop (generate ray paths for propagations at fixed azimuth using multiple inclination angles)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "incl_min"          << '\t' << "degrees"            << '\t'  << '\t' << "0.5" << '\n';
    cout << '\t' << '\t' << "incl_max"          << '\t' << "degrees"            << '\t'  << '\t' << "45.0" << '\n';
    cout << '\t' << '\t' << "incl_step"         << '\t' << "degrees"            << '\t'  << '\t' << "0.5"  << '\n';
    cout << '\t' << '\t' << "inclination"       << '\t' << "see manual"         << '\t' << "0.5" << '\n';
    cout << '\t' << '\t' << "azimuth"           << '\t' << '\t' << "degrees"    << '\t'  << '\t' << "-90.0"  << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "integer"    << '\t'  << '\t' << "2" << '\n';
    cout << '\t' << '\t' << "src_alt"           << '\t' << '\t' << "km"         << '\t'  << '\t' << "0.0" << '\n' << '\n';

    cout << '\t' << "-wnl_wvfrm (compute the weakly non-linear waveform along a specific ray path)" << '\n';
    cout << '\t' << '\t' << "Parameter"  << '\t' << "Units/Options" << '\t' << "Default Value" << '\n';
    cout << '\t' << '\t' << "---------------------------------------------" << '\n';
    cout << '\t' << '\t' << "inclination"       << '\t' << "degrees"            << '\t'  << '\t' << "15.0"  << '\n';
    cout << '\t' << '\t' << "azimuth"           << '\t' << '\t' << "degrees"    << '\t'  << '\t' << "90.0"  << '\n';
    cout << '\t' << '\t' << "bounces"           << '\t' << '\t' << "integer"    << '\t'  << '\t' << "0"  << '\n';
    cout << '\t' << '\t' << "src_alt"           << '\t' << '\t' << "km"         << '\t'  << '\t' << "0.0" << '\n';
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
    cout << '\t' << '\t' << "wvfrm_len"         << '\t' << "see manual"         << '\t' << "2^{13}" << '\n';
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
    cout << '\t' << "write_atmo"        << '\t' << '\t' << "true/false"         << '\t' << "false" << '\n' ;
    cout << '\t' << "prof_format"       << '\t' << '\t' << "see manual"         << '\t' << "zTuvdp" << '\n';
    cout << '\t' << "write_caustics"    << '\t' << '\t' << "true/false"         << '\t' << "false" << '\n';
    cout << '\t' << "calc_amp"          << '\t' << '\t' << "true/false"         << '\t' << "true" << '\n';
    cout << '\t' << "max_alt"           << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "interpolation max" << '\n';
    cout << '\t' << "max_rng"           << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "2000.0" << '\n';
    cout << '\t' << "min_ds"            << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "0.001" << '\n';
    cout << '\t' << "max_ds"            << '\t' << '\t' << '\t' << "km"         << '\t' << '\t' << "0.05" << '\n';
    cout << '\t' << "topo_file"         << '\t' << '\t' << "see manual" << '\t' << "none" << '\n' << '\n';

    cout << "Output (see output files or manual for units):" << '\n';
    cout << '\t' << "atmo.dat -> z[km] : c [m/s]  : u (zonal winds) [m/s] : v (meridional winds) [m/s] : density[g/cm^3] : ceff [km/s]" << '\n';
    cout << '\t' << "{...}.raypaths.dat -> r : z : geo atten : absorption : time " << '\n';
    cout << '\t' << "{...}.arrivals.dat -> incl : az : n_bnc : r : time : cel : z_max : arrival incl : geo atten : absorption" << '\n' << '\n';

    cout << "Examples:" << '\n';
    cout << '\t' << "./bin/infraga-2d -prop examples/ToyAtmo.met incl_step=1.0 bounces=2 max_rng=500.0" << '\n';
    cout << '\t' << "./bin/infraga-2d -wnl_wvfrm examples/ToyAtmo.met azimuth=-90.0 inclination=12.0 wvfrm_opt=impulse wvfrm_p0=500.0" << '\n' << '\n';
}

void run_prop(char* inputs[], int count){
    cout << '\n';
    cout << '\t' << "##########################################" << '\n';
    cout << '\t' << "####        Running infraga-2d        ####" << '\n';
    cout << '\t' << "####            Propagation           ####" << '\n';
    cout << '\t' << "##########################################" << '\n' << '\n';
    
    double theta_min = 0.5, theta_max=45.0, theta_step=0.5, azimuth=-90.0;
    double r_src = 0.0, z_src = 0.0, freq = 0.1;
    int bounces=2;
    bool write_atmo=false, write_caustics=false;
    char* prof_format = "zTuvdp";
    char* topo_file = "None";
    char input_check;
    
    topo::z0 = 0.0;
    atmo::tweak_abs = 1.0;
    geoac::calc_amp = true;
    geoac::is_topo = false;

    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "prof_format=", 12) == 0){       prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){    topo_file = inputs[i] + 10; geoac::is_topo=true;}
    }

    if (count < 3){         cout << "You have to specify an atmosphere file..." << '\n'; return;}
    if (geoac::is_topo){    set_region(inputs[2], topo_file, prof_format, false);}
    else{                   set_region(inputs[2], prof_format, false);}
    
    for(int i = 3; i < count; i++){
        if ((strncmp(inputs[i], "incl_min=", 9) == 0) || (strncmp(inputs[i], "min_incl=", 9) == 0)){        theta_min = atof(inputs[i] + 9);}
        else if ((strncmp(inputs[i], "incl_max=", 9) == 0) || (strncmp(inputs[i], "max_incl=", 9) == 0)){   theta_max = atof(inputs[i] + 9);}
        else if (strncmp(inputs[i], "incl_step=", 10) == 0){                                                theta_step = atof(inputs[i] + 10);}
        else if (strncmp(inputs[i], "inclination=", 12) == 0){                                              theta_min = atof(inputs[i] + 12);
                                                                                                            theta_max = atof(inputs[i] + 12);
                                                                                                            theta_step = 1.0;}

        else if (strncmp(inputs[i], "azimuth=", 8) == 0){                                                   azimuth = atof(inputs[i] + 8);}
        else if (strncmp(inputs[i], "bounces=", 8) == 0){                                                   bounces = atoi(inputs[i] + 8);}
        
        else if (strncmp(inputs[i], "src_alt=", 8) == 0){                                                   z_src = atof(inputs[i] + 8);}
        
        else if (strncmp(inputs[i], "freq=", 5) == 0){                                                      freq = atof(inputs[i] + 5);}
        else if (strncmp(inputs[i], "abs_coeff=", 10) == 0){                                                atmo::tweak_abs = max(0.0, atof(inputs[i] + 10));}
        else if (strncmp(inputs[i], "write_caustics=", 15) == 0){                                           write_caustics = string2bool(inputs[i] + 15);}
        else if (strncmp(inputs[i], "calc_amp=", 9) == 0){                                                  geoac::calc_amp = string2bool(inputs[i] + 9);}
        
        else if ((strncmp(inputs[i], "max_rng=", 8) == 0) || (strncmp(inputs[i], "rng_max=", 8) == 0)){     geoac::rng_max = atof(inputs[i] + 8); cout << '\t' << "User selected range maximum = " << geoac::rng_max << '\n';}
        else if ((strncmp(inputs[i], "max_alt=", 8) == 0) || (strncmp(inputs[i], "alt_max=", 8) == 0)){     geoac::alt_max = atof(inputs[i] + 8); cout << '\t' << "User selected altitude maximum = " << geoac::alt_max << '\n';}

        else if ((strncmp(inputs[i], "max_ds=", 7) == 0) || (strncmp(inputs[i], "ds_max=", 7) == 0)){       geoac::ds_max = atof(inputs[i] + 7);   cout << '\t' << "User selected ds maximum = " << geoac::ds_max << '\n';}
        else if ((strncmp(inputs[i], "min_ds=", 7) == 0) || (strncmp(inputs[i], "ds_min=", 7) == 0)){       geoac::ds_min = atof(inputs[i] + 7);   cout << '\t' << "User selected ds minimum = " << geoac::ds_min << '\n';}

        else if (strncmp(inputs[i], "write_atmo=", 11) == 0){                                               write_atmo = string2bool(inputs[i] + 11);}
        else if (strncmp(inputs[i], "prof_format=", 12) == 0){                                              prof_format = inputs[i] + 12;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){                                                    if(!geoac::is_topo){
                                                                                                                topo::z0 = atof(inputs[i] + 7);
                                                                                                                topo::z_max = topo::z0;
                                                                                                                topo::z_bndlyr = topo::z0 + 2.0;
                                                                                                            } else {
                                                                                                                cout << '\t' << "Note: cannot adjust ground elevation with topography." << '\n';
                                                                                                            }}
        else if (strncmp(inputs[i], "topo_file=",10) == 0){                                                 topo_file = inputs[i] + 10; geoac::is_topo=true;}
        else{
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }
    z_src = max(z_src, topo::z(r_src));
    if(write_atmo)     geoac::write_prof("atmo.dat", (90.0 - azimuth) * Pi / 180.0);

    if(write_caustics) geoac::calc_amp=true;
    geoac::configure();
    
    cout << '\n' << "Parameter summary:" << '\n';
    cout << '\t' << "inclination: " << theta_min << ", " << theta_max << ", " << theta_step << '\n';
    cout << '\t' << "azimuth: " << azimuth << '\n';
    cout << '\t' << "bounces: " << bounces << '\n';
    cout << '\t' << "source location: " << r_src << ", " << z_src << '\n';
    if(!geoac::is_topo){
        cout << '\t' << "ground elevation: " << topo::z0 << '\n';
    }
    cout << '\t' << "frequency: " << freq << '\n';
    cout << '\t' << "S&B atten coeff: " << atmo::tweak_abs << '\n' << '\n';
    
    if (write_caustics){    cout << '\t' << "write_caustics: true" << '\n';} else { cout << '\t' << "write_caustics: false" << '\n';}
    if (geoac::calc_amp){   cout << '\t' << "calc_amp: true" << '\n';} else {       cout << '\t' << "calc_amp: false" << '\n';}
    cout << '\n';
    
    // Extract the file name from the input and use it to ID the output
    char output_buffer [512];
    char* file_id = inputs[2];
    for(int m = strlen(file_id); m >= 0; m--){
        if(file_id[m]=='.'){
            file_id[m] = '\0'; break;
        }
    }
    
    // Define variables used for analysis
	double D, D_prev, travel_time_sum, attenuation, z_max;
	int k, length = geoac::s_max * int(1.0 / (geoac::ds_min * 10));
	bool break_check;

	ofstream results, raypath, caustics;    
    sprintf(output_buffer, "%s.arrivals.dat", file_id);
    results.open(output_buffer);

    results << "# infraga-2d -prop summary:" << '\n';
    results << "#" << '\t' << "profile: " << inputs[2] << ".met" << '\n';
    results << "#" << '\t' << "inclination: " << theta_min << ", " << theta_max << ", " << theta_step << '\n';
    results << "#" << '\t' << "azimuth: " << azimuth << '\n';
    results << "#" << '\t' << "bounces: " << bounces << '\n';
    results << "#" << '\t' << "source location: " << r_src << ", " << z_src << '\n';
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
    results << '\t' << "r_0 [km]";
    results << '\t' << "time [s]";
    results << '\t' << "celerity [km/s]";
    results << '\t' << "turning ht [km]";
    results << '\t' << "arrival incl. [deg]";
    results << '\t' << "geo. atten. [dB]";
    results << '\t' << "absorption [dB]";
    results << '\t' << "perp. dist. [km]";
    results << '\n';
    
    sprintf(output_buffer, "%s.raypaths.dat", file_id);
    raypath.open(output_buffer);
    raypath << "# r [km]";
    raypath << '\t' << "z [km]";
    raypath << '\t' << "geo. atten. [dB]";
    raypath << '\t' << "absorption [dB]";
    raypath << '\t' << "time [s]";
    raypath << '\n';
    
    if(write_caustics){
        for (int bnc = 0; bnc <= bounces; bnc++){
            sprintf(output_buffer, "%s.caustics-%i.dat", file_id, bnc);
            caustics.open(output_buffer);
            caustics << "# r [km]";
            caustics << '\t' << "z [km]";
            caustics << '\t' << "time [s]";
            caustics << '\n';
            caustics.close();
        }
    }
    
	double** solution;
    geoac::build_solution(solution, length);

    double theta_grnd;
    if (z_src - topo::z(0.0) < 0.01){
        if (geoac::is_topo){
            theta_grnd = atan(topo::dz(r_src)) * (180.0 / Pi) + 0.1;
        } else {
            theta_grnd = 0.1;
        }
    } else {
        theta_grnd = -89.9;
    }

	for(double theta = max(theta_min, theta_grnd); theta <= max(theta_max, theta_grnd); theta+=theta_step){
		cout << "Calculating ray path: " << theta << " degrees inclination, " << azimuth << " degrees azimuth." << '\n';
        geoac::theta = theta * Pi / 180.0;
        geoac::phi = Pi / 2.0 - azimuth * Pi / 180.0;
        
        geoac::set_initial(solution, r_src, z_src);
        travel_time_sum = 0.0;
        attenuation = 0.0;
        z_max = 0.0;
		
		for(int bnc_cnt = 0; bnc_cnt <= bounces; bnc_cnt++){
            
			k = geoac::prop_rk4(solution, break_check);

            if(write_caustics) {
                sprintf(output_buffer, "%s.caustics-%i.dat", file_id, bnc_cnt);
                caustics.open(output_buffer,fstream::app);
                D_prev = geoac::jacobian(solution, 1);
            }
            
            
            for(int m = 1; m < k ; m++){
                geoac::travel_time(travel_time_sum, solution, m - 1,m);
                geoac::atten(attenuation, solution, m - 1, m, freq);
                if(write_caustics) D = geoac::jacobian(solution, m);
                
   				if(m == 1 || m % 15 == 0){
                    raypath << solution[m][0];
                    raypath << '\t' << max(solution[m][1], topo::z(solution[m][1]));
                    if(geoac::calc_amp){    raypath << '\t' << 10.0 * log10(geoac::amp(solution, m));}
                    else{                   raypath << '\t' << 0.0;}
                    raypath << '\t' << -attenuation;
                    raypath << '\t' << travel_time_sum;
                    raypath << '\n';
                }
                if((write_caustics) && (D * D_prev < 0.0)){
                    caustics << solution[m][0];
                    caustics << '\t' << solution[m][1];
                    caustics << '\t' << travel_time_sum << '\n';
                    caustics << '\n';
                }
                if(write_caustics) D_prev = D;
            }
            if(write_caustics) caustics.close();
            
    		if(break_check || k < 2){
                break;
            }

            for(int m = 0; m < k ; m++){ z_max = max (z_max, solution[m][1]);}
            
            results << theta;
			results << '\t' << azimuth;
            results << '\t' << bnc_cnt;
			results << '\t' << solution[k][0];
			results << '\t' << travel_time_sum;
            results << '\t' << solution[k][0] / travel_time_sum;
            results << '\t' << z_max;
            results << '\t' <<  - asin(atmo::c(0.0, 0.0, topo::z(solution[k][0])) / atmo::c(0.0, 0.0, z_src) * solution[k][3]) * (180.0 / Pi);
            if(geoac::calc_amp){    results << '\t' << 10.0 * log10(geoac::amp(solution, k));}
            else{                   results << '\t' << 0.0;}
            results << '\t' << -attenuation;
            results << '\t' << geoac::est_dev(solution,k);
			results << '\n';
            
            geoac::set_refl(solution,k);
		}
        geoac::clear_solution(solution,k);
        raypath << '\n';
	}
	
	raypath.close(); results.close();
    geoac::delete_solution(solution, length);
    clear_region();
}


void run_wnl_wvfrm(char* inputs[], int count){
    cout << '\t' << "##########################################" << '\n';
    cout << '\t' << "####        Running infraga-2d        ####" << '\n';
    cout << '\t' << "####    Weakly Non-Linear Waveform    ####" << '\n';
    cout << '\t' << "##########################################" << '\n' << '\n';
    
    double r_src = 0.0, z_src = 0.0, freq = 0.1, D, D_prev;
    int bounces = 0;
    bool break_check, write_atmo=false, write_rays=false;
    char* prof_format = "zTuvdp";
    char* topo_file = "none";
    char input_check;

    double wvfrm_ref=1.0, wvfrm_out_step=1.0e10;
    double ray_length, c0, rho0, D0, p0;
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
        else if (strncmp(inputs[i], "topo_file=", 10) == 0){    topo_file = inputs[i] + 10; geoac::is_topo=true;}
    }

    if (count < 3){         cout << "You have to specify an atmosphere file..." << '\n'; return;}
    if (geoac::is_topo){    set_region(inputs[2], topo_file, prof_format, false);}
    else{                   set_region(inputs[2], prof_format, false);}

    for(int i = 3; i < count; i++){
        if (strncmp(inputs[i], "inclination=", 12) == 0){                                                   geoac::theta = atof(inputs[i] + 12) * (Pi / 180.0);}
        else if (strncmp(inputs[i], "azimuth=", 8) == 0){                                                   geoac::phi = Pi / 2.0 - atof(inputs[i] + 8) * (Pi / 180.0);}
        else if (strncmp(inputs[i], "bounces=", 8) == 0){                                                   bounces = max(0, atoi(inputs[i] + 8));}
        else if (strncmp(inputs[i], "src_alt=", 8) == 0){                                                   z_src = atof(inputs[i] + 8);}
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

        else if ((strncmp(inputs[i], "max_rng=", 8) == 0) || (strncmp(inputs[i], "rng_max=", 8) == 0)){     geoac::rng_max = atof(inputs[i] + 8); cout << '\t' << "User selected range maximum = " << geoac::rng_max << '\n';}
        else if ((strncmp(inputs[i], "max_alt=", 8) == 0) || (strncmp(inputs[i], "alt_max=", 8) == 0)){     geoac::alt_max = atof(inputs[i] + 8); cout << '\t' << "User selected altitude maximum = " << geoac::alt_max << '\n';}
        
        else if ((strncmp(inputs[i], "max_ds=", 7) == 0) || (strncmp(inputs[i], "ds_max=", 7) == 0)){       geoac::ds_max = atof(inputs[i] + 7);   cout << '\t' << "User selected ds maximum = " << geoac::ds_max << '\n';}
        else if ((strncmp(inputs[i], "min_ds=", 7) == 0) || (strncmp(inputs[i], "ds_min=", 7) == 0)){       geoac::ds_min = atof(inputs[i] + 7);   cout << '\t' << "User selected ds minimum = " << geoac::ds_min << '\n';}

        else if (strncmp(inputs[i], "write_atmo=", 11) == 0){                                               write_atmo = string2bool(inputs[i] + 11);}
        else if (strncmp(inputs[i], "prof_format=", 12) == 0){                                              prof_format = inputs[i] + 15;}
        else if (strncmp(inputs[i], "z_grnd=", 7) == 0){                                                    if(!geoac::is_topo){
                                                                                                                topo::z0 = atof(inputs[i] + 7);
                                                                                                                topo::z_max = topo::z0;
                                                                                                                topo::z_bndlyr = topo::z0 + 2.0;
                                                                                                            } else {
                                                                                                                cout << '\t' << "Note: cannot adjust ground elevation with topography." << '\n';
                                                                                                            }}
        else if (strncmp(inputs[i], "topo_file=",10) == 0){                                                 topo_file = inputs[i] + 10; geoac::is_topo=true;}
        else {
            cout << "***WARNING*** Unrecognized parameter entry: " << inputs[i] << '\n';
            cout << "Continue? (y/n):"; cin >> input_check;
            if(input_check!='y' && input_check!='Y') return;
        }
    }
    z_src = max(z_src, topo::z(r_src));
    if(write_atmo)     geoac::write_prof("atmo.dat", geoac::phi);

    geoac::configure();

    cout << '\n' << "Parameter summary:" << '\n';
    cout << '\t' << "inclination: " << geoac::theta * (180.0 / Pi) << '\n';
    cout << '\t' << "azimuth: " << 90.0 - geoac::phi * (180.0 / Pi) << '\n';
    cout << '\t' << "bounces: " << bounces << '\n';
    cout << '\t' << "source location: " << r_src << ", " << z_src << '\n';
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
    char* file_id = inputs[2];
    for(int m = strlen(file_id); m >= 0; m--){
        if(file_id[m]=='.'){
            file_id[m] = '\0';
            break;
        }
    }
    
    // Define variables used for analysis
    double travel_time_sum, attenuation, z_max;
	int k, length = geoac::s_max * int(1.0 / (geoac::ds_min * 10));

    double** solution;
    geoac::build_solution(solution,length);
    
    ofstream raypath;
    if(write_rays){
        sprintf(output_buffer, "%s.raypaths.dat", file_id);
        raypath.open(output_buffer);
        raypath << "# r [km]";
        raypath << '\t' << "z [km]";
        raypath << '\t' << "geo. atten. [dB]";
        raypath << '\t' << "absorption [dB]";
        raypath << '\t' << "time [s]";
        raypath << '\n';
    }

    cout << "Calculating ray path geometry and weakly non-linear waveform evolution..." << '\n';
    if (strncmp(wvfrm_file, "none", 4) != 0){   wvfrm::load_wvfrm(wvfrm_array, wvfrm_file);}
    else {                                      wvfrm::build_wvfrm(wvfrm_array, wvfrm_opt);}

    geoac::set_initial(solution, r_src, z_src);
    travel_time_sum = 0.0;
    attenuation = 0.0;
    z_max = 0.0;

    for(int bnc_cnt = 0; bnc_cnt <= bounces; bnc_cnt++){
        k = geoac::prop_rk4(solution, break_check);

        for(int m = 1; m < k ; m++){
            geoac::travel_time(travel_time_sum, solution, m - 1, m);
            geoac::atten(attenuation, solution, m - 1, m, freq);
            z_max = max (z_max, solution[m][1]);

            if(write_rays && (m == 1 || m % 15 == 0)){
                raypath << solution[m][0];
                raypath << '\t' << max(solution[m][1], topo::z(solution[m][1]));
                raypath << '\t' << 10.0 * log10(geoac::amp(solution, m));
                raypath << '\t' << -attenuation;
                raypath << '\t' << travel_time_sum;
                raypath << '\t' << geoac::est_dev(solution, k);
                raypath << '\n';
            }
        }
        if(bnc_cnt == 0){
            wvfrm_ref_k = 0;
            ray_length = 0.0;
            for(wvfrm_ref_k = 0; wvfrm_ref_k < k; wvfrm_ref_k++){
                ray_length += sqrt(pow(solution[wvfrm_ref_k + 1][0] - solution[wvfrm_ref_k][0], 2) + pow(solution[wvfrm_ref_k + 1][1] - solution[wvfrm_ref_k][1], 2));
                if (ray_length >= wvfrm_ref) break;
            }

            c0 = atmo::c(0.0, 0.0, solution[wvfrm_ref_k + 1][1]) + atmo::u(0.0, 0.0, solution[wvfrm_ref_k + 1][1]) * cos(geoac::phi) + atmo::v(0.0, 0.0, solution[wvfrm_ref_k + 1][1]) * sin(geoac::phi);
            rho0 = atmo::rho(0.0, 0.0, solution[wvfrm_ref_k + 1][1]);
            D0 = geoac::jacobian(solution, wvfrm_ref_k + 1);

            sprintf(output_buffer, "%s.wvfrm_init.dat", file_id);
            wvfrm_out.open(output_buffer);
            wvfrm_out << "# t [sec]" << '\t' << "p(t) [Pa]" << '\n';
            for (int n = 0; n < wvfrm::len; n++){
                wvfrm_out << setprecision(8) << geoac::travel_time(solution, wvfrm_ref_k + 1) + wvfrm_array[n][0] << '\t' << wvfrm_array[n][1] << '\n';
            }
            wvfrm_out.close();

            p0 = 0.0;
            for (int n = 0; n < wvfrm::len; n++){    p0 = max(p0, fabs(wvfrm_array[n][1]));}
            for (int n = 0; n < wvfrm::len; n++){    wvfrm_array[n][1] /= p0;}

            ray_length = geoac::wnl_wvfrm(solution, wvfrm_array, wvfrm_ref_k + 1, k, 0.0, c0, rho0, D0, p0, wvfrm_out_step);
        } else {
            ray_length += geoac::wnl_wvfrm(solution, wvfrm_array, 0, k, ray_length, c0, rho0, D0, p0, wvfrm_out_step);
        }

        if(break_check) break;
        geoac::set_refl(solution,k);
    }

    if(!break_check){
        cout << '\n';
        cout << '\t' << "Arrival summary:" << '\n';
        cout << '\t' << '\t' << "range [km] = " << solution[k][0] << '\n';
        cout << '\t' << '\t' << "time [s] = " << travel_time_sum << '\n';
        cout << '\t' << '\t' << "celerity [km/s] = " << solution[k][0] / travel_time_sum << '\n';
        cout << '\t' << '\t' << "turning height [km] = " << z_max << '\n';
        cout << '\t' << '\t' << "inclination [deg] = " << - asin(atmo::c(0.0, 0.0, topo::z(solution[k][0])) / atmo::c(0.0, 0.0, z_src) * solution[k][3]) * (180.0 / Pi) << '\n';
        cout << '\t' << '\t' << "attenuation (geometric) [dB] = " << 10.0 * log10(geoac::amp(solution, k)) << '\n';
        cout << '\t' << '\t' << "absorption [dB] = " << -attenuation << '\n' << '\n';
    } else {
        cout << '\t' << '\t' << "Ray path does not return to the ground." << '\n' << '\n';
    }

    sprintf(output_buffer, "%s.wvfrm_out.dat", file_id);
    wvfrm_out.open(output_buffer);
    wvfrm_out << "# infraga-2d -wnl_wvfrm summary:" << '\n';
    wvfrm_out << "#" << '\t' << "profile: " << inputs[2] << ".met" << '\n';
    wvfrm_out << "#" << '\t' << "inclination: " << geoac::theta * (180.0 / Pi)  << '\n';
    wvfrm_out << "#" << '\t' << "azimuth: " << geoac::phi * (180.0 / Pi)  << '\n';
    wvfrm_out << "#" << '\t' << "bounces: " << bounces << '\n';
    wvfrm_out << "#" << '\t' << "source location: " << r_src << ", " << z_src << '\n';
    if(!geoac::is_topo){
        wvfrm_out << "#" << '\t' << "ground elevation: " << topo::z0 << '\n';
    } else {
        wvfrm_out << "#" << '\t' << "topo file:" << topo_file << '\n';
    }

    wvfrm_out << "#" << '\t' << "ground evevation: " << topo::z0 << '\n';
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
        wvfrm_out << '\t' << wvfrm_array[n][1] * p0 * sqrt(fabs((atmo::rho(0.0, 0.0, solution[k][1]) * (atmo::c(0.0, 0.0, solution[k][1]) + atmo::u(0.0, 0.0, solution[k][1]) * cos(geoac::phi) + atmo::v(0.0, 0.0, solution[k][1]) * sin(geoac::phi))) / (rho0 * c0) * D0 / geoac::jacobian(solution, k)));
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
    cout << '\t' << "##########################################" << '\n';
    cout << '\t' << "####        Running infraga-2d        ####" << '\n';
    cout << '\t' << "####            Region Test           ####" << '\n';
    cout << '\t' << "##########################################" << '\n' << '\n';
    
    ofstream file_out;
    double alt_max=120.0, eps=1.0e-3;
    
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
    }
    
    if (count < 3){         cout << "You have to specify an atmosphere file..." << '\n'; return;}
    if (geoac::is_topo){    set_region(inputs[2], topo_file, prof_format, false);}
    else{                   set_region(inputs[2], prof_format, false);}

    // Write the sound speed and its derivatives (analytic and defined)
    cout << '\n' << "Running region test comparing defined derivatives to finite differences (make sure you've created 'region_tests/2D/' directory)..." << '\n';
    cout << '\t' << "Evaluating sound speed..." << '\n';
    file_out.open("region_tests/2D/c.analytic.dat");
    for(int m = 0; m < int(alt_max * 100); m++){
        double z0 = m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;
        
        file_out << z0 << '\t';
        file_out << atmo::c(0.0, 0.0, z0) << '\t';
        file_out << atmo::dc(0.0, 0.0, z0, 2) << '\t';
        file_out << atmo::ddc(0.0, 0.0, z0, 2, 2) << '\n';
    }
    file_out.close();
    
    file_out.open("region_tests/2D/c.finite.dat");
    for(int m = 1; m < int(alt_max * 100); m++){
        double z0 = m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;
        
        file_out << z0 << '\t';
        file_out << atmo::c(0.0, 0.0, z0) << '\t';
        file_out << (atmo::c(0.0, 0.0, z0 + eps) - atmo::c(0.0, 0.0, z0 - eps)) / (2.0 * eps) << '\t';
        file_out << (atmo::c(0.0, 0.0, z0 + eps) + atmo::c(0.0, 0.0, z0 - eps) - 2.0 * atmo::c(0.0, 0.0, z0)) / pow(eps, 2.0) << '\n';
    }
    file_out.close();
    
    file_out.open("region_tests/2D/c.eval_all.dat");
    for(int m = 1; m < int(alt_max * 100); m++){
        double z0 = m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;

        double c, dc, ddc;
        interp::eval_all(z0, atmo::c_spline, c, dc, ddc);
        
        file_out << z0 << '\t';
        file_out << c << '\t';
        file_out << dc << '\t';
        file_out << ddc << '\n';
    }
    file_out.close();

    // Write the (u) wind component and its derivatives (analytic and finite diff)
    cout << '\t' << "Evaluating winds (u)..." << '\n';
    file_out.open("region_tests/2D/u.analytic.dat");
    for(int m = 0; m < int(alt_max * 100); m++){
        double z0 = m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;
        
        file_out << z0 << '\t';
        file_out << atmo::u(0.0, 0.0, z0) << '\t';
        file_out << atmo::du(0.0, 0.0, z0, 2) << '\t';
        file_out << atmo::ddu(0.0, 0.0, z0, 2, 2) << '\n';
    }
    file_out.close();
    
    file_out.open("region_tests/2D/u.finite.dat");
    for(int m = 1; m < int(alt_max * 100); m++){
        double z0 = m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;

        file_out << z0 << '\t';
        file_out << atmo::u(0.0, 0.0, z0) << '\t';
        file_out << (atmo::u(0.0, 0.0, z0 + eps) - atmo::u(0.0, 0.0, z0 - eps)) / (2.0 * eps) << '\t';
        file_out << (atmo::u(0.0, 0.0, z0 + eps) + atmo::u(0.0, 0.0, z0 - eps) - 2.0 * atmo::u(0.0, 0.0, z0)) / pow(eps, 2.0) << '\n';
    }
    file_out.close();

    file_out.open("region_tests/2D/u.eval_all.dat");
    for(int m = 1; m < int(alt_max * 100); m++){
        double z0 = m / 100.0;
        if(z0 > atmo::u_spline.x_vals[atmo::u_spline.length - 1]) break;

        double u, du, ddu;
        interp::eval_all(z0, atmo::u_spline, u, du, ddu);
        
        file_out << z0 << '\t';
        file_out << u << '\t';
        file_out << du << '\t';
        file_out << ddu << '\n';
    }
    file_out.close();

    // Write the (v) wind component and its derivatives (analytic and finite diff)
    cout << '\t' << "Evaluating winds (v)..." << '\n';
    file_out.open("region_tests/2D/v.analytic.dat");
    for(int m = 0; m < int(alt_max * 100); m++){
        double z0 = m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;

        file_out << z0 << '\t';
        file_out << atmo::v(0.0, 0.0, z0) << '\t';
        file_out << atmo::dv(0.0, 0.0, z0, 2) << '\t';
        file_out << atmo::ddv(0.0, 0.0, z0, 2, 2) << '\n';
    }
    file_out.close();
    
    file_out.open("region_tests/2D/v.finite.dat");
    for(int m = 1; m < int(alt_max * 100); m++){
        double z0 = m/100.0;
        if(z0 > atmo::c_spline.x_vals[atmo::c_spline.length - 1]) break;

        file_out << z0 << '\t';
        file_out << atmo::v(0.0, 0.0, z0) << '\t';
        file_out << (atmo::v(0.0, 0.0, z0 + eps) - atmo::v(0.0, 0.0, z0 - eps)) / (2.0 * eps) << '\t';
        file_out << (atmo::v(0.0, 0.0, z0 + eps) + atmo::v(0.0, 0.0, z0 - eps) - 2.0 * atmo::v(0.0, 0.0, z0)) / pow(eps, 2.0) << '\n';
    }
    file_out.close();
        
    file_out.open("region_tests/2D/v.eval_all.dat");
    for(int m = 1; m < int(alt_max * 100); m++){
        double z0 = m / 100.0;
        if(z0 > atmo::v_spline.x_vals[atmo::v_spline.length - 1]) break;

        double v, dv, ddv;
        interp::eval_all(z0, atmo::u_spline, v, dv, ddv);
        
        file_out << z0 << '\t';
        file_out << v << '\t';
        file_out << dv << '\t';
        file_out << ddv << '\n';
    }
    file_out.close();

    
    // Write the topography and its derivatives to file (analytic and finite diff)
    if(geoac::is_topo){
        cout << '\t' << "Evaluating topography..." << '\n';
        eps=1.0e-1;
        file_out.open("region_tests/2D/topo.analytic.dat");
        for(int m = 0; m < 5000; m++){
            double r0 = topo::spline.x_vals[0] + m / 5000.0 * (topo::spline.x_vals[topo::spline.length - 1] - topo::spline.x_vals[0]);
            
            file_out << r0 << '\t';
            file_out << topo::z(r0) << '\t';
            file_out << topo::dz(r0) << '\t';
            file_out << topo::ddz(r0) << '\n';
        }
        file_out.close();
        
        
        file_out.open("region_tests/2D/topo.finite.dat");
        for(int m = 1; m < 5000; m++){
            double r0 = topo::spline.x_vals[0] + m / 5000.0 * (topo::spline.x_vals[topo::spline.length - 1] - topo::spline.x_vals[0]);
            
            file_out << r0 << '\t';
            file_out << topo::z(r0) << '\t';
            file_out << (topo::z(r0 + eps) - topo::z(r0 - eps)) / (2.0 * eps) << '\t';
            file_out << (topo::z(r0 + eps) + topo::z(r0 - eps) - 2.0 * topo::z(r0)) / pow(eps, 2.0) << '\n';
        }
        file_out.close();

        file_out.open("region_tests/2D/topo.eval_all.dat");
        for(int m = 1; m < 5000; m++){
            double r0 = topo::spline.x_vals[0] + m / 5000.0 * (topo::spline.x_vals[topo::spline.length - 1] - topo::spline.x_vals[0]);
            
            double zg, dzg, ddzg;
            interp::eval_all(r0, topo::spline, zg, dzg, ddzg);
        
            file_out << r0 << '\t';
            file_out << zg << '\t';
            file_out << dzg << '\t';
            file_out << ddzg << '\n';
        }
        file_out.close();

    }
    clear_region();
}


int main(int argc, char* argv[]){
    if (argc < 2){ usage();}
    else {
        if ((strncmp(argv[1], "--version", 9) == 0) || (strncmp(argv[1], "-v", 2) == 0)){       version();}
        else if ((strncmp(argv[1], "--usage", 7) == 0) || (strncmp(argv[1], "-u", 2) == 0)){    usage();}
        
        else if (strncmp(argv[1], "-prop",5) == 0){                                             run_prop(argv, argc);}
        else if (strncmp(argv[1], "-wnl_wvfrm",10) == 0){                                       run_wnl_wvfrm(argv, argc);}

        // Test methods...development use only.
        else if (strncmp(argv[1], "-region_test",12) == 0){                                     run_region_test(argv, argc);}
        else {                                                                                  cout << "Unrecognized option." << '\n';}
    }
    return 0;
}
