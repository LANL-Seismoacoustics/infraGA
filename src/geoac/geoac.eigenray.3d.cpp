#ifndef _GEOAC_EIGENRAY_3D_CPP_
#define _GEOAC_EIGENRAY_3D_CPP_

#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "geoac.params.h"
#include "geoac.eqset.h"
#include "geoac.interface.h"
#include "geoac.eigenray.h"

#include "../atmo/atmo_state.h"
#include "../util/rk4solver.h"

// NOTE: All angles are assumed to be radians (convert deg->rad in main.cpp file)

using namespace std;

bool geoac::verbose = false;
int geoac::eigenray_cnt = 0;

double geoac::damping = 1.0e-3;
double geoac::tolerance = 0.1;

double geoac::dth_big = 0.1 * (Pi / 180.0);
double geoac::dth_sml = 0.001 * (Pi / 180.0);

ofstream geoac::eig_results;

double geoac::mod_dth(double dr, double dr_dtheta){ return dth_big - (dth_big - dth_sml) * exp(- 1.0 / 2.0 * pow(dr / dr_dtheta, 2));}

bool geoac::est_eigenray(double src[3], double rcvr[2], double th_min, double th_max, double & th_est, double & ph_est, double & th_next, int bncs, double az_err_lim){
    double rcvr_rng = sqrt(pow(rcvr[0] - src[0], 2) + pow(rcvr[1] - src[1], 2));
    double rcvr_az = atan2(rcvr[1] - src[1], rcvr[0] - src[0]);
    
    if(verbose){
        cout << '\t' << "Estimating eigenray angles for source-receiver separated by " << rcvr_rng << " km, and azimuth " << 90.0 - rcvr_az * (180.0 / Pi);
        cout << " degrees from N.  Inclination limits: " << th_min * 180.0 / Pi << ", " << th_max * 180.0 / Pi << "." << '\n';
    }
    
    int	iterations = 0, k, length = int(s_max / ds_min);
    double r, r_prev, dth = dth_big, dph = 100.0;
    bool break_check, th_max_reached;

    calc_amp = false;
    configure();

    double** solution;
    build_solution(solution, length);
    
    phi = rcvr_az;

    th_max_reached = false;
    while(fabs(dph * (180.0 / Pi)) > az_err_lim && iterations < 5){
        r = rcvr_rng;
        r_prev = rcvr_rng;

        for(double th = th_min; th <= th_max; th += dth){
            if(th + dth >= th_max){
                th_max_reached = true;
            }

            theta =	th;
            set_initial(solution, src[0], src[1], src[2]);

            k = prop_rk4(solution, break_check, length);
            if(!break_check){
                for(int n_bnc = 1; n_bnc <= bncs; n_bnc++){
                    set_refl(solution,k);
                    k = prop_rk4(solution, break_check, length);
                    if(break_check){
                        break;
                    }
                }
            }
        
            r = sqrt(pow(solution[k][0] - src[0], 2) + pow(solution[k][1] - src[1], 2));
            if(break_check){
                r = rng_max;
                if(verbose){
                    cout << '\t' << '\t' << "Ray launched with inclination " << theta * (180.0 / Pi) << ", azimuth " << 90.0 - phi * (180.0 / Pi) << " doesn't return to the ground" << '\n';
                }
            } else {
                if(verbose){
                    cout << '\t' << '\t' << "Ray launched with inclination " << th * (180.0 / Pi) << ", azimuth " << 90.0 - phi * (180.0 / Pi) << " arrives at range " << r;
                    cout << " km after " << bncs << " bounce(s)." << '\t' << "Exact arrival at " << solution[k][0] << " km East, " << solution[k][1] << " km North" << '\n';
                }
            }
       
            if(((r - rcvr_rng) * (r_prev - rcvr_rng) <= 0.0) && (r < rng_max) && (r_prev < rng_max) && (th > th_min)){
                if(iterations==0){
                    th_next = th;
                }
                
                dph  = atan2(solution[k][1] - src[1], solution[k][0] - src[0]) * (180.0 / Pi);
                dph -= atan2(rcvr[1] - src[1], rcvr[0] - src[0]) * (180.0 / Pi);
                while(dph >  180.0){
                    dph -= 360.0;
                }
                while(dph < -180.0){
                    dph += 360.0;
                }
                
                if(fabs(dph) < az_err_lim){
                    th_est = theta - dth;
                    ph_est = phi;

                    if(verbose){
                        cout << '\t' << '\t' << "Azimuth deviation = " << dph << ", less than " << az_err_lim << " degrees: estimates acceptable." << '\n' << '\n';
                    }

                    delete_solution(solution, length);
                    return true;
                } else {
                    if(verbose){
                        cout << '\n' << '\t' << '\t' << "Azimuth deviation = " << dph << ", greater than " << az_err_lim << " degrees: compensating and searching inclinations again." << '\n';
                        cout << '\t' << '\t' << "Launch azimuth correction: " << 90.0 - phi * (180.0 / Pi) << " --> " << 90.0 - (phi * (180.0 / Pi) - dph * 0.9) << '\n' << '\n';
                    }
                    phi -= dph * (Pi / 180.0) * 0.9;
                    th_min = max(th - 5.0 * (Pi / 180.0), th_min);
                }
                iterations++;
                break;
            }

            if(iterations >= 2){
                dth = mod_dth(r - rcvr_rng, (r - r_prev) / (2.0 * dth));
            }

            r_prev = r;
        }

        if(th_max_reached){
            th_next = th_max;
            break;
        }
    }
    delete_solution(solution, length);
    
    if(verbose) cout << '\t' << '\t' << "Reached maximum inclination angle or iteration limit." << '\n';
    return false;
}

void geoac::find_eigenray(double src[3], double rcvr[2], double th_est, double ph_est, double freq, int bnc_cnt, int iterate_limit, char title[]){
	bool break_check;
    char output_buffer [512];
    ofstream raypath;
    double D, attenuation, travel_time_sum, z_max, inclination, back_az, back_az_dev, dr, dr_prev = 10000.0, step_sc = 1.0;

    long double x, y, c_src, c_grnd, dzg_dx, dzg_dy, ds_dth, ds_dph, dx, dy, dx_dth, dy_dth, dx_dph, dy_dph;
	long double det, dth, dph;
    
    int	k, length = int(s_max / ds_min);

    calc_amp = true;
    configure();
    
    double** solution;
    build_solution(solution, length);
    
    theta =	th_est;
    phi = ph_est;
                    
    if(verbose){
        cout << '\t' << '\t' << "Searching for exact eigenray using auxiliary parameters." << '\n';
    }

	for(int n = 0; n <= iterate_limit; n++){
        if(n == iterate_limit){
            if(verbose){
                cout << '\t' <<'\t' << '\t' << "Search for exact eigenray maxed out iterations.  No eigneray idenfied." << '\n' << '\n';
            }
            break;
        }
        
        set_initial(solution, src[0], src[1], src[2]);
        if(verbose){
            cout << '\t' << '\t' << "Calculating ray path: " << theta * (180.0 / Pi) << " degrees inclination, " << 90.0 - phi * (180.0 / Pi) << " degrees azimuth";
        }
        k = prop_rk4(solution, break_check, length);
        if(break_check){
            if(verbose){
                cout << '\t' << "Ray path left propagation region, reversing step and adjusting step scaling." << '\n';
            }

            theta -= dth * step_sc;
            phi -= dph * step_sc;
            step_sc /= 2.0;
        }
        
        if(!break_check){
            for(int n_bnc = 1; n_bnc <= bnc_cnt; n_bnc++){
                set_refl(solution,k);
                k = prop_rk4(solution, break_check, length);
                if(break_check){
                    if(verbose){
                        cout << '\t' << "Ray path left propagation region, reversing step and adjusting step scaling." << '\n';
                    }

                    theta -= dth * step_sc;
                    phi -= dph * step_sc;
                    step_sc /= 2.0;            
                }
                break;
            }
        }
        
        if(!break_check){
    		// Compute arrival location and check if it's within the defined tolerance
            x = solution[k][0];     dx = rcvr[0] - x;
            y = solution[k][1];     dy = rcvr[1] - y;
            dr = sqrt(pow(dx, 2) + pow(dy, 2));
            if(verbose){
                cout << '\t' << '\t' << "Arrival after " << bnc_cnt << " reflections at (" << x << ", " << y << "), distance to receiver = " << dr << " km." << '\n';
            }

            if(dr < tolerance){
                sprintf(output_buffer, "%s.eigenray-%i.dat", title, eigenray_cnt);
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
                k = prop_rk4(solution, break_check, length);

                for(int m = 1; m < k ; m++){
                    travel_time(travel_time_sum, solution, m - 1,m);
                    atten(attenuation, solution, m - 1, m, freq);
                    z_max = max(z_max, solution[m][2]);
                
                    if(m == 1 || m % 15 == 0){
                        raypath << solution[m][0];
                        raypath << '\t' << solution[m][1];
                        raypath << '\t' << max(solution[m][2], topo::z(solution[m][0], solution[m][1]));
                        raypath << '\t' << 20.0 * log10(amp(solution, m));
                        raypath << '\t' << -attenuation;
                        raypath << '\t' << travel_time_sum << '\n';
                    }
                }
                for(int n_bnc = 1; n_bnc <= bnc_cnt; n_bnc++){
                    set_refl(solution, k);
                    k = prop_rk4(solution, break_check, length);

                    for(int m = 1; m < k; m++){
                        travel_time(travel_time_sum, solution, m - 1,m);
                        atten(attenuation, solution, m - 1, m, freq);
                        z_max = max(z_max, solution[m][2]);

                        if(m == 1 || m % 15 == 0){
                            raypath << solution[m][0];
                            raypath << '\t' << solution[m][1];
                            raypath << '\t' << solution[m][2];
                            raypath << '\t' << 20.0 * log10(amp(solution,m));
                            raypath << '\t' << attenuation;
                            raypath << '\t' << travel_time_sum << '\n';
                        }
                    }          
                }
                raypath.close();
            
                inclination = - asin(atmo::c(solution[k][0], solution[k][1], solution[k][2]) / atmo::c(src[0], src[1], src[2]) * solution[k][5]) * (180.0 / Pi);

                back_az = 90.0 - atan2(-solution[k][4], -solution[k][3]) * (180.0 / Pi);
                while(back_az > 180.0){
                    back_az -= 360.0;
                }
                while(back_az < -180.0){
                    back_az += 360.0;
                }

                back_az_dev = back_az - (90.0 - atan2(src[1] - rcvr[1], src[0] - rcvr[0]) * (180.0 / Pi));
                while(back_az_dev > 180.0){
                    back_az_dev -= 360.0;
                }
                while(back_az_dev < -180.0){
                    back_az_dev += 360.0;
                }

                if(verbose){
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
                    cout << '\t' << '\t' << '\t' << "attenuation (geometric) [dB] = " << 20.0 * log10(geoac::amp(solution, k)) << '\n';
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
                eig_results << '\t' << 20.0 * log10(geoac::amp(solution, k));
                eig_results << '\t' << -attenuation;
                eig_results << '\n';
            
                eigenray_cnt++;
            
                break;
            } else if(n > 0 && dr > dr_prev){
                // If the range to the receiver has increased, undo the previous changes to theta and phi,
                // half the step scalar and repeat the step using the new scaled increments
                theta -= dth * step_sc;
                phi -= dph * step_sc;
                step_sc /= 2.0;
            } else {
                step_sc = min(1.0, step_sc * 1.25);
            
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

                det = pow(1.0 + damping, 2) * (dx_dth * dy_dph - dx_dph * dy_dth);
            
                dth = ((1.0 + damping) * dy_dph * dx - dx_dph * dy) / det;
                dph = ((1.0 + damping) * dx_dth * dy - dy_dth * dx) / det;
            
                theta += dth * step_sc;
                phi += dph * step_sc;
            
                dr_prev = dr;
            }
            clear_solution(solution, k);
        }
        if(sqrt(dth * dth + dph * dph) * step_sc < 1.0e-12){
            if(verbose){
                cout << '\t' << '\t' <<  '\t' << "Step size too small, psuedo-critical ray path likely." << '\n' << '\n';
            }
            break;
        }
    }
    
    delete_solution(solution, length);
}

#endif /* _GEOAC_EIGENRAY_3D_CPP_ */
