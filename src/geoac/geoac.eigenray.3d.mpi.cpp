#ifndef _GEOAC_EIGENRAY_3D_MPI_CPP_
#define _GEOAC_EIGENRAY_3D_MPI_CPP_

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

// NOTE: All angles are assumed to be radians (convert deg --> rad in main.cpp file)

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
    double rcvr_rng = sqrt(pow(rcvr[0] - src[0], 2) + pow(rcvr[1] - src[1], 2));
    double rcvr_az = atan2(rcvr[1] - src[1], rcvr[0] - src[0]);
    
    if(verbose && verbose_opt == rcvr_id && rank == 0){
        cout << '\t' << "Estimating eigenray angles for source-receiver separated by " << rcvr_rng << " km, and azimuth ";
        cout << 90.0 - rcvr_az * (180.0 / Pi) << " degrees from N.  Inclination limits: " << th_min * 180.0 / Pi << ", " << th_max * 180.0 / Pi << "." << '\n';
    }
    
    int	iterations = 0, k, length = s_max * int(1.0 / (ds_min * 10)), success[2];
    double dth = dth_big, dph = 100.0, arrival_rng, arrival_az_dev, results_buffer[3], rngs[size], az_devs[size], prev_rng, prev_az_dev;
    bool break_check, th_max_reached;
    MPI_Status mpi_status;
    
    calc_amp = false;
    configure();
    
    double** solution;
    build_solution(solution, length);
    
    phi = rcvr_az;
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
                cout << " - " << (theta + dth * (size - 1)) * (180.0 / Pi) << ", azimuth = " << 90.0 - phi * (180.0 / Pi);}
            
            k = prop_rk4(solution, break_check);
            if(!break_check){
                for(int n_bnc = 1; n_bnc <= bncs; n_bnc++){
                    set_refl(solution, k); k = prop_rk4(solution, break_check);
                    if(break_check) break;
                }
            }
            
            if(break_check){
                arrival_rng = rcvr_rng;
            } else {
                arrival_rng = sqrt(pow(solution[k][0] - src[0], 2) + pow(solution[k][1] - src[1], 2));
                arrival_az_dev = rcvr_az - atan2(solution[k][1] - src[1], solution[k][0] - src[0]);
                while(arrival_az_dev > Pi){
                    arrival_az_dev -= 2.0 * Pi;
                }
                while(arrival_az_dev < -Pi){
                    arrival_az_dev += 2.0 * Pi;
                }
            }
            MPI_Barrier(comm);
            
            MPI_Gather(&arrival_rng, 1, MPI_DOUBLE, rngs, 1, MPI_DOUBLE, 0, comm);
            MPI_Gather(&arrival_az_dev, 1, MPI_DOUBLE, az_devs, 1, MPI_DOUBLE, 0, comm);
            
            if (rank == 0){
                if(verbose && verbose_opt == rcvr_id){
                    if (fabs(rngs[0] - rcvr_rng) > 1.0e-6){
                        cout << '\t' << "Arrival ranges: " << rngs[0];
                    } else {
                        cout << '\t' << "Arrival ranges: -";
                    }
                    
                    for(int n = 1; n < size; n++){
                        if (fabs(rngs[n] - rcvr_rng) > 1.0e-6){
                            cout << ", " << rngs[n];
                        } else {
                            cout << ", - ";
                        }
                    }
                    cout << "." << '\n';
                }
                
                if((theta > th_min) && (prev_rng - rcvr_rng) * (rngs[0] - rcvr_rng) < 0.0) {
                    if(iterations == 0){
                        results_buffer[2] = theta;
                    }
                    success[0] = 1;
                        
                    if(fabs(prev_az_dev) < az_err_lim){
                        if(verbose && verbose_opt == rcvr_id){
                            cout << '\t' << '\t' << '\t' << "Identified eigenray estimate.  ";
                            cout << "Azimuth deviation = " << prev_az_dev << ".  Estimate acceptable." << '\n' << '\n';
                        }
                        results_buffer[0] = theta;
                        results_buffer[1] = phi;
                        success[1] = 1;
                            
                    } else {
                        if(verbose && verbose_opt == rcvr_id){
                            cout << '\t' << '\t' << '\t' << "Identified eigenray estimate.  ";
                            cout << "Azimuth deviation = " << prev_az_dev << ".  Shifting azimuth and searching inclinations again." << '\n' << '\n';
                        }
                        results_buffer[0] = max(theta - 10.0, th_min);
                        results_buffer[1] = phi + 0.75 * prev_az_dev;
                        success[1] = 0;
                    }
                }
                if(success[0] != 1){
                    for(int k = 0; k < size - 1; k++){
                        if((rngs[k] - rcvr_rng) * (rngs[k + 1] - rcvr_rng) < 0.0){
                            if(iterations == 0){
                                results_buffer[2] = theta + (k + 1) * dth;
                            }
                            success[0] = 1;
                            
                            if(fabs(az_devs[k]) < az_err_lim){
                                if(verbose && verbose_opt == rcvr_id){
                                    cout << '\t' << '\t' << '\t' << "Identified eigenray estimate.  ";
                                    cout << "Azimuth deviation = " << az_devs[k] << ".  Estimate acceptable." << '\n' << '\n';
                                }
                                results_buffer[0] = theta + k * dth;
                                results_buffer[1] = phi;
                                success[1] = 1;
                                break;
                            } else {
                                if(verbose && verbose_opt == rcvr_id){
                                    cout << '\t' << '\t' << '\t' << "Identified eigenray estimate.  ";
                                    cout << "Azimuth deviation = " << az_devs[k] << ".  Shifting azimuth and searching inclinations again." << '\n' << '\n';
                                }
                                results_buffer[0] = max(theta - 10.0, th_min);
                                results_buffer[1] = phi + 0.75 * az_devs[k];
                                success[1] = 0;
                                break;
                            }
                        }
                    }}
                prev_rng = rngs[size - 1];
                prev_az_dev = az_devs[size - 1];
            }
            MPI_Barrier(comm);
            
            // CPU-0 sends out success state and results summaries
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
        }
        
        MPI_Barrier(comm);
        if (success[1] == 1){
            break;
        }
        if(th_max_reached){
            th_next = th_max;
            if(rank == 0){
                cout << '\n';
            }
            break;
        }
        
        iterations++;
        dth = dth_big / iterations;
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
    double D, attenuation, travel_time_sum, z_max, inclination, back_az, back_az_dev, dr, dr_prev = 10000.0;
    
    double dth_lim = dth_big * 0.9, dph_lim = 0.5 * (Pi / 180.0), step_sc = 1.0;
    long double x, y, c_src, c_grnd, dzg_dx, dzg_dy, ds_dth, ds_dph, dx, dy, dx_dth, dy_dth, dx_dph, dy_dph;
	long double det, dth, dph;
    
    int	k, length = s_max * int(1.0 / (ds_min * 10));

    calc_amp = true;
    configure();

    double** solution;
    build_solution(solution, length);
    
    theta =	th_est;
    phi = 	ph_est;

    if(verbose && verbose_opt == rcvr_id){
        cout << '\t' << '\t' << "Searching for exact eigenray using auxiliary parameters." << '\n';
    }
	for(int n = 0; n <= iterate_limit; n++){
        if(n == iterate_limit){
            if(verbose && verbose_opt == rcvr_id){
                cout << '\t' <<'\t' << '\t' << "Search for exact eigenray maxed out iterations.  ";
                cout << "No eigneray idenfied." << '\n' << '\n';
            }
            break;
        }
        
        set_initial(solution, src[0], src[1], src[2]);
        if(verbose && verbose_opt == rcvr_id){
            cout << '\t' << '\t' << "Calculating ray path: " << theta * (180.0 / Pi);
            cout << " degrees inclination, " << 90.0 - phi * (180.0 / Pi) << " degrees azimuth";
        }

        k = prop_rk4(solution, break_check);
        if(break_check){
            if(verbose && verbose_opt == rcvr_id){
                cout << '\t' << "Ray path left propagation region." << '\n';
            }
            break;
        }
        
        for(int n_bnc = 1; n_bnc <= bnc_cnt; n_bnc++){
            set_refl(solution,k); k = prop_rk4(solution, break_check);
            if(break_check){
                if(verbose && verbose_opt == rcvr_id){
                    cout << '\t' << "Ray path left propagation region." << '\n';
                }
                break;
            }
        }
        if(break_check){
            break;
        }
		// Compute arrival location and check if it's within the defined tolerance
        x = solution[k][0];     dx = rcvr[0] - x;
        y = solution[k][1];     dy = rcvr[1] - y;
        dr = sqrt(pow(dx, 2) + pow(dy, 2));
        if(verbose && verbose_opt == rcvr_id){
            cout << '\t' << '\t' << "Arrival after " << bnc_cnt << " reflections at (" << x << ", " << y << "), ";
            cout << "distance to receiver = " << dr << " km." << '\n';
        }
            
        if(dr < tolerance){
            sprintf(output_buffer, "%s%i.dat", title, eigenray_cnt);
            raypath.open(output_buffer);
            
            raypath << "# x [km]";
            raypath << '\t' << "y [km]";
            raypath << '\t' << "z [km]";
            raypath << '\t' << "geo. atten. [dB]";
            raypath << '\t' << "absorption [dB]";
            raypath << '\t' << "time [s]" << '\n';
            
            attenuation = 0.0;
            travel_time_sum = 0.0;
            z_max = 0.0;
            
            set_initial(solution, src[0], src[1], src[2]);
            k = prop_rk4(solution, break_check);
            
            for(int m = 1; m < k ; m++){
                travel_time(travel_time_sum, solution, m - 1,m);
                atten(attenuation, solution, m - 1, m, freq);
                z_max = max(z_max, solution[m][2]);
                
                if(m == 1 || m % 15 == 0){
                    raypath << solution[m][0];
                    raypath << '\t' << solution[m][1];
                    raypath << '\t' << max(solution[m][2],0.0);
                    raypath << '\t' << 10.0 * log10(amp(solution, m));
                    raypath << '\t' << -attenuation;
                    raypath << '\t' << travel_time_sum << '\n';

                }
            }
            for(int n_bnc = 1; n_bnc <= bnc_cnt; n_bnc++){
                set_refl(solution, k);
                k = prop_rk4(solution, break_check);

                for(int m = 1; m < k; m++){
                    travel_time(travel_time_sum, solution, m - 1,m);
                    atten(attenuation, solution, m - 1, m, freq);
                    z_max = max(z_max, solution[m][2]);
                    
                    if(m == 1 || m % 15 == 0){
                        raypath << solution[m][0];
                        raypath << '\t' << solution[m][1];
                        raypath << '\t' << max(solution[m][2], topo::z(solution[m][0], solution[m][1]));
                        raypath << '\t' << 10.0 * log10(amp(solution,m));
                        raypath << '\t' << attenuation;
                        raypath << '\t' << travel_time_sum << '\n';
                    }
                }
                
            }
            raypath.close();
            
            inclination = - asin(atmo::c(solution[k][0], solution[k][1], solution[k][2]) / atmo::c(src[0], src[1], src[2]) * solution[k][5]) * (180.0 / Pi);
            back_az = 90.0 - atan2(-solution[k][4], -solution[k][3]) * (180.0 / Pi);
            while(back_az > 180.0)      back_az -= 360.0;
            while(back_az < -180.0)     back_az += 360.0;
            
            back_az_dev = back_az - (90.0 - atan2(src[1] - rcvr[1], src[0] - rcvr[0]) * 180.0/Pi);
            while(back_az_dev > 180.0)  back_az_dev -= 360.0;
            while(back_az_dev < -180.0) back_az_dev += 360.0;
            
            if(verbose && verbose_opt == rcvr_id){
                cout << '\t' << '\t' << "Eigenray Identified:" << '\n';
                cout << '\t' << '\t' << '\t' << "inclination [deg] = " << setprecision(8) << theta * (180.0 / Pi)  << '\n';
                cout << '\t' << '\t' << '\t' << "azimuth [deg] = " << setprecision(8) << 90.0 - phi * (180.0 / Pi) << '\n';
                cout << '\t' << '\t' << '\t' << "bounces [-] = " << bnc_cnt << '\n';
                cout << '\t' << '\t' << '\t' << "x [km E-W] = " << solution[k][0] << '\n';
                cout << '\t' << '\t' << '\t' << "y [km N-S] = " << solution[k][1] << '\n';
                cout << '\t' << '\t' << '\t' << "time [s] = " << travel_time_sum << '\n';
                cout << '\t' << '\t' << '\t' << "celerity [km/s] = " << sqrt(pow(solution[k][0] - src[0], 2) + pow(solution[k][1] - src[1], 2)) / travel_time_sum << '\n';
                cout << '\t' << '\t' << '\t' << "turning height [km] = " << z_max << '\n';
                cout << '\t' << '\t' << '\t' << "arrival inclination [deg] = " << inclination << '\n';
                cout << '\t' << '\t' << '\t' << "back azimuth = " << back_az  << '\n';
                cout << '\t' << '\t' << '\t' << "attenuation (geometric) [dB] = " << 10.0 * log10(geoac::amp(solution, k)) << '\n';
                cout << '\t' << '\t' << '\t' << "absorption [dB] = " << -attenuation << '\n' << '\n';
            } else {
                cout << '\t' << "Eigenray identified:" << '\t' << "theta, phi = " << setprecision(8) << theta * (180.0 / Pi) << ", " << 90.0 - phi * (180.0 / Pi) << " degrees." << '\n';
            }
            
            eig_results << setprecision(8) << theta * (180.0 / Pi);
            eig_results << '\t' << setprecision(8) << 90.0 - phi * (180.0 / Pi);
            eig_results << '\t' << bnc_cnt;
            eig_results << '\t' << solution[k][0];
            eig_results << '\t' << solution[k][1];
            eig_results << '\t' << travel_time_sum;
            eig_results << '\t' << sqrt(pow(solution[k][0]- src[0], 2) + pow(solution[k][1] - src[1], 2)) / travel_time_sum;
            eig_results << '\t' << z_max;
            eig_results << '\t' << inclination;
            eig_results << '\t' << back_az;
            eig_results << '\t' << 10.0 * log10(geoac::amp(solution, k));
            eig_results << '\t' << -attenuation;
            eig_results << '\n';

            th_est = theta;
            ph_est = phi;
            
            eigenray_cnt++;
            success = true;
            break;
        } else if(n > 0 && dr > dr_prev){
            // If the range to the receiver has increased, undo the previous changes to theta and phi,
            // half the step scalar and repeat the step using the new scaled increments
            theta -= dth * step_sc;
            phi -= dph * step_sc;
            step_sc /= 2.0;
            if(sqrt(pow(dth, 2) + pow(dph, 2)) * step_sc < 1.0e-12){
                if (verbose && verbose_opt == rcvr_id){
                    cout << '\t' << '\t' <<  '\t' << "Step size too small, psuedo-critical ray path likely." << '\n' << '\n';
                }
                break;
            }
        } else {
            step_sc = min(1.0, step_sc * 1.5);
            
            c_src = atmo::c(src[0], src[1], src[2]);
            c_grnd = atmo::c(solution[k][0], solution[k][1], topo::z(solution[k][0], solution[k][1]));

            if (is_topo){
                dzg_dx = topo::dz(solution[k][0], solution[k][1], 0);
                dzg_dy = topo::dz(solution[k][0], solution[k][1], 1);
                
                ds_dth = - c_src / c_grnd * (solution[k][8] -  dzg_dx * solution[k][6] -  dzg_dy * solution[k][7]) /  (solution[k][5] - dzg_dx * solution[k][3] - dzg_dy * solution[k][4]);
                ds_dph = - c_src / c_grnd * (solution[k][14] - dzg_dx * solution[k][12] - dzg_dy * solution[k][13]) / (solution[k][5] - dzg_dx * solution[k][3] - dzg_dy * solution[k][4]);
            } else {
                ds_dth = - c_src / c_grnd * solution[k][8] / solution[k][5];
                ds_dph = - c_src / c_grnd * solution[k][14] / solution[k][5];
            }
            dx_dth = solution[k][6]  + solution[k][3] * ds_dth; dy_dth = solution[k][7]  + solution[k][4] * ds_dth;
            dx_dph = solution[k][12] + solution[k][3] * ds_dph; dy_dph = solution[k][13] + solution[k][4] * ds_dph;

            det = pow(1.0 + damping, 2) * dx_dth * dy_dph - dx_dph * dy_dth;
            
            dth = 1.0 / det * ((1.0 + damping) * dy_dph * dx - dx_dph * dy);
            dph = 1.0 / det * ((1.0 + damping) * dx_dth * dy - dy_dth * dx);
           
            theta += dth * step_sc;
            phi += dph * step_sc;
            
            dr_prev = dr;
        }
        clear_solution(solution, k);
	}
    delete_solution(solution, length);
    return success;
}

#endif /* _GEOAC_EIGENRAY_3D_MPI_CPP_ */
