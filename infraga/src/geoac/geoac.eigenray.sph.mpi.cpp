#ifndef _GEOAC_EIGENRAY_SPH_MPI_CPP_
#define _GEOAC_EIGENRAY_SPH_MPI_CPP_

#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>

#include "geoac.params.h"
#include "geoac.eqset.h"
#include "geoac.interface.h"
#include "geoac.eigenray.mpi.h"

#include "../atmo/atmo_state.h"

#include "../util/rk4solver.h"
#include "../util/globe.h"

// NOTE: All angles are assumed to be radians (convert deg --> rad in {}_main.cpp file)
// Also, the source array contains (altitude, latitude, longitude) while the reciever
// array contains (latitude, longitude).  Latitudes and longitudes are converted to
// radians in the main.cpp file and altitude is relative to the wgs84 ellipsoid.

using namespace std;

bool geoac::verbose = false;
int geoac::verbose_opt = 0;

int geoac::eigenray_cnt = 0;

double geoac::damping = 1.0e-3;
double geoac::tolerance = 0.1;

double geoac::dth_big = 0.1 * (Pi / 180.0);
double geoac::dth_sml = 0.001 * (Pi / 180.0);

ofstream geoac::eig_results;

double geoac::mod_dth(double dr, double dr_dtheta){
    return dth_big - (dth_big - dth_sml) * exp(- 1.0 / 2.0 * pow(dr / dr_dtheta, 2));
}

bool geoac::est_eigenray(double src[3], double rcvr[2], double th_min, double th_max, double & th_est, double & ph_est, double & th_next, int bncs, double az_err_lim, MPI_Comm comm, int rank, int size, int rcvr_id){
    double rcvr_rng = globe::gc_dist(src[1], src[2], rcvr[0], rcvr[1]);
    double rcvr_az = globe::bearing(src[1], src[2], rcvr[0], rcvr[1]) * (180.0 / Pi);
    
    if(verbose && verbose_opt == rcvr_id && rank == 0){
        cout << '\t' << "Estimating eigenray angles for source-receiver at great circle distance " << rcvr_rng << " km and azimuth ";
        cout << rcvr_az << " degrees.  Inclination limits: [" << th_min * (180.0 / Pi) << ", " << th_max * (180.0 / Pi) << "]." << '\n';
    }
    
    int	iterations = 0, k, length = int(s_max / ds_min), success[2];
    double dth = dth_big, dph = 100.0, arrival_rng, arrival_az_dev, results_buffer[3], rngs[size], az_devs[size], prev_rng = rng_max, prev_az_dev = 1000.0;
    bool break_check, th_max_reached;

    calc_amp = false;
    configure();
    
    double** solution;
    build_solution(solution, length);
    
    phi = Pi / 2.0 - rcvr_az * (Pi / 180.0);

    th_max_reached = false;
    while(fabs(dph) > az_err_lim && iterations < 5){
        success[0] = 0;
        success[1] = 0;

        for(double th = th_min; th <= th_max; th += dth * size){
            if(th + dth * size >= th_max){
                th_max_reached = true;
            }
            
            theta =	th + dth * rank;
            set_initial(solution, src[0], src[1], src[2]);
            if(verbose && verbose_opt == rcvr_id && rank == 0){
                cout << '\t' << '\t' << "Computing ray paths with inclinations " << theta * (180.0 / Pi);
                cout << " - " << (theta + dth * (size - 1)) * (180.0 / Pi) << ", azimuth = " << 90.0 - phi * (180.0 / Pi);
            }
            
            k = prop_rk4(solution, break_check, length);
            if(!break_check){
                for(int n_bnc = 1; n_bnc <= bncs; n_bnc++){
                    set_refl(solution, k);
                    k = prop_rk4(solution, break_check, length);
                    if(break_check){
                        break;
                    }
                }
            }
            
            if(break_check){
                arrival_rng = rng_max;
                arrival_az_dev = 100.0;
            } else {
                arrival_rng = globe::gc_dist(src[1], src[2], solution[k][1], solution[k][2]);
                        
                arrival_az_dev = rcvr_az - globe::bearing(src[1], src[2], solution[k][1], solution[k][2]) * (180.0 / Pi);
                while(arrival_az_dev > 180.0){
                    arrival_az_dev -= 360.0;
                }
                while(arrival_az_dev < -180.0){
                    arrival_az_dev += 360.0;
                }
            }
            MPI_Barrier(comm);       
            MPI_Gather(&arrival_rng, 1, MPI_DOUBLE, rngs, 1, MPI_DOUBLE, 0, comm);
            MPI_Gather(&arrival_az_dev, 1, MPI_DOUBLE, az_devs, 1, MPI_DOUBLE, 0, comm);
            MPI_Barrier(comm);
            
            if (rank == 0){
                if(verbose && verbose_opt == rcvr_id){
                    if (fabs(rngs[0] - rng_max) > 1.0e-6){
                        cout << '\t' << "Arrival ranges: " << rngs[0];
                    } else {
                        cout << '\t' << "Arrival ranges: -";
                    }

                    for(int n = 1; n < size; n++){
                        if (fabs(rngs[n] - rng_max) > 1.0e-6){
                            cout << ", " << rngs[n];
                        } else {
                            cout << ", - ";
                        }
                    }
                    cout << '\n';
                }
                
                if(((rngs[0] - rcvr_rng) * (prev_rng - rcvr_rng) <= 0.0) && (rngs[0] < rng_max) && (prev_rng < rng_max) && (th - dth > th_min)){
                    if(iterations == 0){
                        results_buffer[2] = theta;
                    }
                    success[0] = 1;
                        
                    if(fabs(prev_az_dev) < az_err_lim){
                        if(verbose && verbose_opt == rcvr_id){
                            cout << '\t' << '\t' << '\t' << "Azimuth deviation = " << prev_az_dev << ", less than " << az_err_lim;
                            cout << " degrees: estimates acceptable." << '\n' << '\n';
                        }

                        results_buffer[0] = theta - dth;
                        results_buffer[1] = phi;
                        success[1] = 1;

                    } else {
                        if(verbose && verbose_opt == rcvr_id){
                            cout << '\n' << '\t' << '\t' <<'\t' << "Azimuth deviation = " << prev_az_dev << ", greater than " << az_err_lim << " degrees: compensating and searching inclinations again." << '\n';
                            cout << '\t' << '\t' << '\t' << "Launch azimuth correction: " << 90.0 - phi * (180.0 / Pi) << " --> " << 90.0 - ((phi * (180.0 / Pi) - prev_az_dev * 0.9)) << '\n' << '\n';
                        }

                        results_buffer[0] = max(th_min - 5.0 * (Pi / 180.0), th_min + dth * size);
                        results_buffer[1] = phi - prev_az_dev * (Pi / 180.0) * 0.9;
                        success[1] = 0;

                        iterations++;                        
                    }
                } else {
                    for(int k = 0; k < size - 1; k++){
                        if(((rngs[k] - rcvr_rng) * (rngs[k + 1] - rcvr_rng) <= 0.0) && (rngs[k] < rng_max) && (rngs[k + 1] < rng_max) && (th + (dth * k) > th_min)){
                            if(iterations == 0){
                                results_buffer[2] = theta + (k + 1) * dth;
                            }
                            success[0] = 1;
                        
                            if(fabs(az_devs[k]) < az_err_lim){
                                if(verbose && verbose_opt == rcvr_id){
                                    cout << '\t' << '\t' << '\t' << "Azimuth deviation = " << az_devs[k] << ", less than " << az_err_lim;
                                    cout << " degrees: estimates acceptable." << '\n' << '\n';
                                }

                                results_buffer[0] = theta + k * dth;
                                results_buffer[1] = phi;
                                success[1] = 1;
                            } else {
                                if(verbose && verbose_opt == rcvr_id){
                                    cout << '\n' << '\t' << '\t' <<'\t' << "Azimuth deviation = " << az_devs[k] << ", greater than " << az_err_lim << " degrees: compensating and searching inclinations again." << '\n';
                                    cout << '\t' << '\t' << '\t' << "Launch azimuth correction: " << 90.0 - phi * (180.0 / Pi) << " --> " << 90.0 - ((phi * (180.0 / Pi) - az_devs[k] * 0.9)) << '\n' << '\n';
                                }

                                results_buffer[0] = max(th_min - 5.0 * (Pi / 180.0), th_min + dth * size);
                                results_buffer[1] = phi - az_devs[k] * (Pi / 180.0) * 0.9;
                                success[1] = 0;

                                iterations++;
                            }
                            break;
                        }
                    }
                }
                prev_rng = rngs[size - 1];
                prev_az_dev = az_devs[size - 1];
            }
            MPI_Barrier(comm);
            MPI_Bcast(success, 2, MPI_INT, 0, comm);
            MPI_Bcast(results_buffer, 3, MPI_DOUBLE, 0, comm);
            MPI_Barrier(comm);

            if (success[1] == 1){
                th_est = results_buffer[0];
                ph_est = results_buffer[1];
                th_next = results_buffer[2];
                break;
            } else if (success[0] == 1) {
                th_min = results_buffer[0];
                phi = results_buffer[1];
                break;
            }

            if(iterations >= 3){
                dth = dth_big / (iterations + 2);
            }
        }

        MPI_Barrier(comm);
        if(th_max_reached){
            th_next = th_max;
            break;
        }
        
        if (success[1] == 1){
            break;
        }
    }
    
    delete_solution(solution, length);
    if(verbose && verbose_opt == rcvr_id && rank == 0 && !success[1]){
        cout << '\t' << '\t' << "Reached maximum inclination angle or iteration limit." << '\n';
    }

    return success[1];
}


bool geoac::find_eigenray(double src[3], double rcvr[2], double & th_est, double & ph_est, double freq, int bnc_cnt, int iterate_limit, char title[], int rcvr_id){
	bool break_check, success=false;
    char output_buffer [512];
    ofstream raypath;
    double D, attenuation, travel_time_sum, r_max, inclination, back_az, back_az_dev, dr, dr_prev = 10000.0;
    
	double dth_lim = dth_big * 0.9, dph_lim = 0.5 * (Pi / 180.0), step_sc = 1.0;
    long double lat, lon, r_grnd, c_src, c_grnd, dzg_dlat, dzg_dlon, ds_norm, ds_dth, ds_dph, dlat, dlon, dlat_dth, dlon_dth, dlat_dph, dlon_dph;
    long double det, dth, dph;
    
    int	k, length = int(s_max / ds_min);

    calc_amp = true;
    configure();

    double** solution;
    build_solution(solution, length);

    theta =	th_est;
    phi = ph_est;
    
    if(verbose && verbose_opt == rcvr_id){
        cout << '\t' << '\t' << "Searching for exact eigenray using auxiliary parameters." << '\n';
    }
	for(int n = 0; n <= iterate_limit; n++){
        if(n == iterate_limit){
            if(verbose && verbose_opt == rcvr_id){
                cout << '\t' <<'\t' << '\t' << "Search for exact eigenray maxed out iterations.  No eigneray idenfied." << '\n';
            }
            break;
        }
        
        if(verbose && verbose_opt == rcvr_id){
            cout << '\t' << '\t' << "Calculating ray path: " << theta * (180.0 / Pi) << " degrees inclination, ";
            cout << 90.0 - phi * (180.0 / Pi) << " degrees azimuth";
        }
        
		set_initial(solution, src[0], src[1], src[2]);
        k = prop_rk4(solution, break_check, length);
        if(break_check){
            if(verbose&& verbose_opt == rcvr_id){
                cout << '\t' << "Ray path left propagation region, reversing step and adjusting step scaling." << '\n';
            }

            theta -= dth * step_sc;
            phi -= dph * step_sc;
            step_sc /= 2.0;
        } else {
            for(int n_bnc = 1; n_bnc <= bnc_cnt; n_bnc++){
                set_refl(solution, k);
                k = prop_rk4(solution, break_check, length);
                if(break_check){
                    if(verbose && verbose_opt == rcvr_id){
                        cout << '\t' << "Ray path left propagation region, reversing step and adjusting step scaling." << '\n';
                    }

                    theta -= dth * step_sc;
                    phi -= dph * step_sc;
                    step_sc /= 2.0;
                    break;
                }
            } 
        }

        if(!break_check){
    		// Determine arrival location and check if it's within the defined tolerance
            dr = globe::gc_dist(solution[k][1], solution[k][2], rcvr[0], rcvr[1]);
            if(verbose && verbose_opt == rcvr_id){
                cout << '\t' << '\t' << "Arrival at (" << setprecision(8) << solution[k][1] * (180.0 / Pi);
                cout << ", " << solution[k][2] * (180.0 / Pi) << "), distance to receiver = " << dr << " km." << '\n';
            }
        
            if(dr < tolerance) {
                sprintf(output_buffer, "%s%i.dat", title, eigenray_cnt);
                raypath.open(output_buffer);

                raypath << "# lat [deg]";
                raypath << '\t' << "lon [deg]";
                raypath << '\t' << "z [km]";
                raypath << '\t' << "trans. coeff. [dB]";
                raypath << '\t' << "absorption [dB]";
                raypath << '\t' << "time [s]";
                raypath << '\n';
                
                attenuation = 0.0;
                travel_time_sum = 0.0;
                r_max = 0.0;
            
                set_initial(solution, src[0], src[1], src[2]);
                k = prop_rk4(solution, break_check, length);

                for(int m = 1; m < k; m++){
                    travel_time(travel_time_sum, solution, m - 1, m);
                    atten(attenuation, solution, m - 1, m, freq);
                    r_max = max (r_max, solution[m][0] - globe::r0);
                
                    if(m == 1 || m % 15 == 0){
                        raypath << setprecision(8) << solution[m][1] * (180.0 / Pi);
                        raypath << '\t' << setprecision(8) << solution[m][2] * (180.0 / Pi);
                        raypath << '\t' << solution[m][0] - globe::r0;
                        raypath << '\t' << 20.0 * log10(amp(solution, m));
                        raypath << '\t' << -2.0 * attenuation;
                        raypath << '\t' << travel_time_sum << '\n';
                    }
                }
                for(int n_bnc = 1; n_bnc <= bnc_cnt; n_bnc++){
                    set_refl(solution, k);
                
                    k = prop_rk4(solution, break_check, length);
                    for(int m = 1; m < k; m++){
                        travel_time(travel_time_sum, solution, m - 1, m);
                        atten(attenuation, solution, m - 1, m, freq);
                        r_max = max(r_max, solution[m][0] - globe::r0);
                    
                        if(m == 1 || m % 15 == 0){
                            raypath << setprecision(8) << solution[m][1] * (180.0 / Pi);
                            raypath << '\t' << setprecision(8) << solution[m][2] * (180.0 / Pi);
                            raypath << '\t' << solution[m][0] - globe::r0;
                            raypath << '\t' << 20.0 * log10(amp(solution, m));
                            raypath << '\t' << -2.0 * attenuation;
                            raypath << '\t' << travel_time_sum << '\n';
                        }
                    }
                }
                raypath.close();
                        
                inclination = - asin(atmo::c(solution[k][0], solution[k][1], solution[k][2]) / atmo::c(src[0], src[1], src[2]) * solution[k][3]) * 180.0 / Pi;
                back_az = 90.0 - atan2(-solution[k][4], -solution[k][5]) * 180.0 / Pi;
                while(back_az < -180.0) back_az += 360.0;
                while(back_az >  180.0) back_az -= 360.0;
            
                back_az_dev = ((Pi / 2.0 - atan2(-solution[k][4], -solution[k][5]) ) - globe::bearing(rcvr[0], rcvr[1], src[1], src[2]) ) * (180.0 / Pi);
                if(back_az_dev >  180.0) back_az_dev -= 360.0;
                if(back_az_dev < -180.0) back_az_dev += 360.0;
            
                if(verbose && verbose_opt == rcvr_id){
                    cout << '\t' << '\t' << "Eigenray-" << eigenray_cnt << ":" << '\n';
                    cout << '\t' << '\t' << '\t' << "inclination [deg] = " << theta * (180.0 / Pi) << '\n';
                    cout << '\t' << '\t' << '\t' << "azimuth [deg] = " << 90.0 - phi * (180.0 / Pi) << '\n';
                    cout << '\t' << '\t' << '\t' << "bounces [-] = " << bnc_cnt << '\n';
                    cout << '\t' << '\t' << '\t' << "latitude [deg] = " << setprecision(8) << solution[k][1] * 180.0 / Pi << '\n';
                    cout << '\t' << '\t' << '\t' << "longitude [deg] = " << setprecision(8) << solution[k][2] * 180.0 / Pi << '\n';
                    cout << '\t' << '\t' << '\t' << "time [s] = " << travel_time_sum << '\n';
                    cout << '\t' << '\t' << '\t' << "celerity [km/s] = " << globe::gc_dist(solution[k][1], solution[k][2], src[1], src[2]) / travel_time_sum << '\n';
                    cout << '\t' << '\t' << '\t' << "turning height [km] = " << r_max << '\n';
                    cout << '\t' << '\t' << '\t' << "arrival inclination [deg] = " << inclination << '\n';
                    cout << '\t' << '\t' << '\t' << "back azimuth [deg] = " << back_az << '\n';
                    cout << '\t' << '\t' << '\t' << "attenuation (geometric) [dB] = " << 20.0 * log10(geoac::amp(solution,k)) << '\n';
                    cout << '\t' << '\t' << '\t' << "absorption [dB] = " << -2.0 * attenuation << '\n' << '\n';
                }

                eig_results << setprecision(8) << theta * (180.0 / Pi);
                eig_results << '\t' << setprecision(8) << 90.0 - phi * (180.0 / Pi);
                eig_results << '\t' << bnc_cnt;
                eig_results << '\t' << setprecision(8) << solution[k][1] * 180.0 / Pi;
                eig_results << '\t' << setprecision(8) << solution[k][2] * 180.0 / Pi;
                eig_results << '\t' << travel_time_sum;
                eig_results << '\t' << globe::gc_dist(solution[k][1], solution[k][2], src[1], src[2]) / travel_time_sum;
                eig_results << '\t' << r_max;
                eig_results << '\t' << inclination;
                eig_results << '\t' << back_az;
                eig_results << '\t' << 20.0 * log10(geoac::amp(solution,k));
                eig_results << '\t' << -2.0 * attenuation;
                eig_results << '\n';

                th_est = theta;
                ph_est = phi;

                eigenray_cnt++;
                success = true;
                break;
            } else if (n > 0 && dr > dr_prev){
                if(verbose && verbose_opt == rcvr_id){
                    cout << '\t' << '\t' <<  '\t' << "Distance to receiver increased, reversing step and adjusting step scaling." << '\n';
                }

                theta -= dth * step_sc;
                phi -= dph * step_sc;
                step_sc /= 2.0;
            } else {
                step_sc = min(1.0, step_sc * 1.25);
            
                dlat = rcvr[0] - solution[k][1];
                dlon = rcvr[1] - solution[k][2];
            
                r_grnd = topo::z(solution[k][1], solution[k][2]);
                c_grnd = atmo::c(r_grnd, solution[k][1], solution[k][2]);
                c_src = atmo::c(src[0] + globe::r0, src[1], src[2]);
            
                dzg_dlat = topo::dz(solution[k][1], solution[k][2], 0) / r_grnd;
                dzg_dlon = topo::dz(solution[k][1], solution[k][2], 1) / (r_grnd * cos(solution[k][1]));
                ds_norm = solution[k][3] - dzg_dlat * solution[k][4] -  dzg_dlon  * solution[k][5];
            
                ds_dth = - c_src / c_grnd * (solution[k][6] -  dzg_dlat * r_grnd * solution[k][7] -  dzg_dlon * (r_grnd * cos(solution[k][1])) * solution[k][8])  / ds_norm;
                ds_dph = - c_src / c_grnd * (solution[k][12] - dzg_dlat * r_grnd * solution[k][13] - dzg_dlon * (r_grnd * cos(solution[k][1])) * solution[k][14]) / ds_norm;

                dlat_dth = solution[k][7]  + solution[k][4] / r_grnd * ds_dth;  dlon_dth = solution[k][8]  + solution[k][5] / (r_grnd * cos(solution[k][1])) * ds_dth;
                dlat_dph = solution[k][13] + solution[k][4] / r_grnd * ds_dph;  dlon_dph = solution[k][14] + solution[k][5] / (r_grnd * cos(solution[k][1])) * ds_dph;

                det = pow(1.0 + damping, 2) * dlat_dth * dlon_dph - dlat_dph * dlon_dth;
            
                dth = ((1.0 + damping) * dlon_dph * dlat - dlat_dph * dlon) / det;
                dph = ((1.0 + damping) * dlat_dth * dlon - dlon_dth * dlat) / det;
            
                theta += dth * step_sc;
                phi += dph * step_sc;

                dr_prev = dr;
            }
            clear_solution(solution, k);
        }
        if(sqrt(dth * dth + dph * dph) * step_sc < 1.0e-12){
            if (verbose && verbose_opt == rcvr_id){
                cout << '\t' << '\t' <<  '\t' << "Step size too small, near-critical ray path likely." << '\n' << '\n';
            }
            break;
        }
	}
    delete_solution(solution, length);
    return success;
}

#endif /* _GEOAC_EIGENRAY_SPH_MPI_CPP_ */
