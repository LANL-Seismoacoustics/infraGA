#ifndef _GEOAC_EQSET_3D_CPP_
#define _GEOAC_EQSET_3D_CPP_

#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>

#include <fftw3.h>

#include "geoac.params.h"
#include "geoac.eqset.h"
#include "../atmo/atmo_state.h"
#include "../atmo/atmo_io.3d.strat.h"
#include "../util/waveforms.h"

using namespace std;

//--------------------------------------------//
//-------Propagate in 3D Cartesian Space -----//
//--------------------------------------------//
void geoac::set_system(){
    dim = 3;
    is_strat = true;
}

//-----------------------------------------------------------//
//-------Structure containing source functions which are-----//
//-------called multiple times in the solver-----------------//
//-----------------------------------------------------------//
namespace geoac {
    struct src_refs{
        double src_loc[3]; double c0;               // Source location and thermodynamic sound speed at that location
    
        double c; double dc[5]; double ddc[3][2];   // Thermodynamic sound speed, first order derivatives (x,y,z,th,az), and second order (x, y, or z and inc or az)
        double u; double du[5]; double ddu[3][2];   // E-W wind component, first order derivatives (x,y,z,th,az), and second order (x, y, or z and inc or az)
        double v; double dv[5]; double ddv[3][2];   // N-S wind component, first order derivatives (x,y,z,th,az), and second order (x, y, or z and inc or az)
        double w; double dw[5]; double ddw[3][2];   // Vertical wind component, first order derivatives (x,y,z,th,az), and second order (x, y, or z and inc or az)
    
        double nu0; double nu_mag; double dnu_mag[2];                       // Eikonal vector magnitude at the source, along the ray path, and variation with respect to inc or az
        double cg[3]; double cg_mag; double dcg[3][2]; double dcg_mag[2];   // Group velocity components, magnitude, and variations of each
    };

    struct src_refs refs = {
        {0.0, 0.0, 0.0}, 0.0,
        0.0, {0.0, 0.0, 0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
        0.0, {0.0, 0.0, 0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
        0.0, {0.0, 0.0, 0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
        0.0, {0.0, 0.0, 0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
        0.0, 0.0, {0.0, 0.0},
        {0.0, 0.0, 0.0}, 0.0, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}, {0.0, 0.0}
    };
}

//----------------------------------------------------------------------//
//-------Fill in solution[0][n] with the appropriate initial values-----//
//----------------------------------------------------------------------//
void geoac::set_initial(double ** & solution, double x0, double y0, double z0){
    refs.src_loc[0] = x0;
    refs.src_loc[1] = y0;
    refs.src_loc[2] = z0;
	refs.c0 = atmo::c(x0, y0, z0);
    
    double mach_comps[3] =  { atmo::u(x0, y0, z0) / refs.c0,    atmo::v(x0, y0, z0) / refs.c0,      atmo::w(x0, y0, z0) / refs.c0};
    double nu0[3] =         { cos(theta) * cos(phi),            cos(theta) * sin(phi),              sin(theta)};
    double mu0_th[3] =      {-sin(theta) * cos(phi),           -sin(theta) * sin(phi),              cos(theta)};
    double mu0_ph[3] =      {-cos(theta) * sin(phi),            cos(theta) * cos(phi),              0.0};
    
    double mach_scalar = 1.0 + (nu0[0] * mach_comps[0] + nu0[1] * mach_comps[1] + nu0[2] * mach_comps[2]);
    refs.nu0 = 1.0 / mach_scalar;
    
	for(int n_eq = 0; n_eq < eq_cnt; n_eq++){
		switch(n_eq){
			case(0):    solution[0][n_eq] =  x0;   break;  // x(0) = x0
            case(1):    solution[0][n_eq] =  y0;   break;  // y(0) = y0
            case(2):    solution[0][n_eq] =  z0;   break;  // z(0) = z0
                
			case(3):                                                            // nu_x(0) = nu0_x / M
			case(4):                                                            // nu_y(0) = nu0_y / M
			case(5):    solution[0][n_eq] = nu0[n_eq - 3] / mach_scalar;  break;  // nu_z(0) = nu0_z / M
                
            case(6):                                        // Xt(0) = 0
			case(7):                                        // Yt(0) = 0
			case(8):                                        // Zt(0) = 0
			case(12):                                       // Xp(0) = 0
			case(13):                                       // Yp(0) = 0
			case(14):   solution[0][n_eq] = 0.0;    break;  // Zp(0) = 0
                
			case(9):			// mu_x_th(0) = d/d lt (nu_x)
			case(10):			// mu_y_th(0) = d/d lt (nu_y)
			case(11):			// mu_z_th(0) = d/d lt (nu_z)
                solution[0][n_eq] =  mu0_th[n_eq - 9] / mach_scalar - nu0[n_eq - 9] / pow(mach_scalar, 2) * (mu0_th[0] * mach_comps[0] + mu0_th[1] * mach_comps[1] + mu0_th[2] * mach_comps[2]);
				break;
                
			case(15):			// mu_x_ph(0) = d/d lp (nu_x)
            case(16):			// mu_y_ph(0) = d/d lp (nu_y)
            case(17):			// mu_z_ph(0) = d/d lp (nu_z)
                solution[0][n_eq] =  mu0_ph[n_eq - 15] / mach_scalar - nu0[n_eq - 15] / pow(mach_scalar, 2) * (mu0_ph[0] * mach_comps[0] + mu0_ph[1] * mach_comps[1] + mu0_ph[2] * mach_comps[2]);
                break;
                
			default:
                cout << "Unexpected index in Initial_Cond.  Model includes 18 variables." << '\n';
		}
	}    
}

//--------------------------------------------------------------------------//
//-------Taylor series fit to more accurately deterine intercept values-----//
//--------------------------------------------------------------------------//
void geoac::approx_intercept(double ** solution, int k, double* & prev){
    double dz_grnd = topo::z(solution[k - 1][0], solution[k - 1][1]) - solution[k - 1][2];  // Set dz from z_{k-1} to ground
    double dz_k = (solution[k - 1][2] - solution[k - 2][2])                                 // Set effective dz for step = z_{k - 2} - z_{k - 1} with ground slope
                    - (solution[k - 1][0] - solution[k - 2][0]) * topo::dz(solution[k - 1][0], solution[k-1][1], 0)
                    - (solution[k - 1][1] - solution[k - 2][1]) * topo::dz(solution[k - 1][0], solution[k-1][1], 1);
    
    for(int n_eq = 0; n_eq < eq_cnt; n_eq++){ prev[n_eq] = solution[k - 1][n_eq] + (solution[k - 1][n_eq] - solution[k - 2][n_eq]) / dz_k * dz_grnd;}
}

//-------------------------------------------------------------------------//
//-------Fill in solution[0][n] with the appropriate reflection values-----//
//-------------------------------------------------------------------------//
void geoac::set_refl(double** & solution, int k){
	double* prev = new double [eq_cnt];
	approx_intercept(solution, k, prev);

    double zg, c_grnd, a, b, norm1, norm2, C1x, C1y, C1z, C2x, C2y;

    double da_dth, db_dth, da_dph, db_dph;
    double dC1x_dth, dC1y_dth, dC1z_dth, dC2x_dth, dC2y_dth, ds0_dth;
    double dC1x_dph, dC1y_dph, dC1z_dph, dC2x_dph, dC2y_dph, ds0_dph;
    double dnux_ds, dnuy_ds, dnuz_ds;
    
    zg = topo::z(prev[0], prev[1]);
    c_grnd = atmo::c(prev[0], prev[1], zg);

    a = topo::dz(prev[0], prev[1], 0);
    b = topo::dz(prev[0], prev[1], 1);
    norm1 = 1.0 + pow(a, 2) + pow(b, 2);
    
    C1x = (1.0 - pow(a, 2) + pow(b, 2)) / norm1;    C2x = 2.0 * a / norm1;
    C1y = (1.0 + pow(a, 2) - pow(b, 2)) / norm1;    C2y = 2.0 * b / norm1;
    C1z = (1.0 - pow(a, 2) - pow(b, 2)) / norm1;
        
    if(calc_amp){
        norm2 = prev[5] - a * prev[3] - b * prev[4];

        da_dth =  (prev[6] * (prev[5] - b * prev[4]) - prev[3] * (prev[8] - b * prev[7])) / norm2 * topo::ddz(prev[0], prev[1], 0, 0)
                + (prev[7] * (prev[5] - a * prev[3]) - prev[4] * (prev[8] - a * prev[6])) / norm2 * topo::ddz(prev[0], prev[1], 0, 1);
        db_dth =  (prev[6] * (prev[5] - b * prev[4]) - prev[3] * (prev[8] - b * prev[7])) / norm2 * topo::ddz(prev[0], prev[1], 1, 0)
                + (prev[7] * (prev[5] - a * prev[3]) - prev[4] * (prev[8] - a * prev[6])) / norm2 * topo::ddz(prev[0], prev[1], 1, 1);
            
        da_dph =  (prev[12] * (prev[5] - b * prev[4]) - prev[3] * (prev[14] - b * prev[13])) / norm2 * topo::ddz(prev[0], prev[1], 0, 0)
                + (prev[13] * (prev[5] - a * prev[3]) - prev[4] * (prev[14] - a * prev[12])) / norm2 * topo::ddz(prev[0], prev[1], 0, 1);
        db_dph =  (prev[12] * (prev[5] - b * prev[4]) - prev[3] * (prev[14] - b * prev[13])) / norm2 * topo::ddz(prev[0], prev[1], 1, 0)
                + (prev[13] * (prev[5] - a * prev[3]) - prev[4] * (prev[14] - a * prev[12])) / norm2 * topo::ddz(prev[0], prev[1], 1, 1);
            
        dC1x_dth = -4.0 * a / pow(norm1, 2) * ((1.0 + pow(b, 2)) * da_dth - a * b * db_dth);                    dC1x_dph = -4.0 * a / pow(norm1, 2) * ((1.0 + pow(b, 2)) * da_dph - a * b * db_dph);
        dC1y_dth = -4.0 * b / pow(norm1, 2) * ((1.0 + pow(a, 2)) * db_dth - a * b * da_dth);                    dC1y_dph = -4.0 * b / pow(norm1, 2) * ((1.0 + pow(a, 2)) * db_dph - a * b * da_dph);
        dC1z_dth = -4.0 / pow(norm1, 2) * (a * da_dth + b * db_dth);                                            dC1z_dph = -4.0 / pow(norm1, 2) * (a * da_dph + b * db_dph);
        dC2x_dth =  2.0 / pow(norm1, 2) * ((1.0 - pow(a, 2) + pow(b, 2)) * da_dth - 2.0 * a * b * db_dth);      dC2x_dph =  2.0 / pow(norm1, 2) * ((1.0 - pow(a, 2) + pow(b, 2)) * da_dph - 2.0 * a * b * db_dph);
        dC2y_dth =  2.0 / pow(norm1, 2) * ((1.0 + pow(a, 2) - pow(b, 2)) * db_dth - 2.0 * a * b * da_dth);      dC2y_dph =  2.0 / pow(norm1, 2) * ((1.0 + pow(a, 2) - pow(b, 2)) * db_dph - 2.0 * a * b * da_dph);

        ds0_dth = - refs.c0 / c_grnd * (prev[8] - a * prev[6] - b * prev[7]) / norm2;
        ds0_dph = - refs.c0 / c_grnd * (prev[14] - a * prev[12] - b * prev[13]) / norm2;

        dnux_ds = - 1.0 / c_grnd * (refs.c0 / c_grnd * atmo::dc(prev[0], prev[1], zg, 0) + prev[3] * atmo::du(prev[0], prev[1], zg, 0) + prev[4] * atmo::dv(prev[0], prev[1], zg, 0) + prev[5] * atmo::dw(prev[0], prev[1], zg, 0));
        dnuy_ds = - 1.0 / c_grnd * (refs.c0 / c_grnd * atmo::dc(prev[0], prev[1], zg, 1) + prev[3] * atmo::du(prev[0], prev[1], zg, 1) + prev[4] * atmo::dv(prev[0], prev[1], zg, 1) + prev[5] * atmo::dw(prev[0], prev[1], zg, 1));
        dnuz_ds = - 1.0 / c_grnd * (refs.c0 / c_grnd * atmo::dc(prev[0], prev[1], zg, 2) + prev[3] * atmo::du(prev[0], prev[1], zg, 2) + prev[4] * atmo::dv(prev[0], prev[1], zg, 2) + prev[5] * atmo::dw(prev[0], prev[1], zg, 2));
    }
        
	for(int n_eq = 0; n_eq < eq_cnt; n_eq++){
		switch(n_eq){
			case(0):                                            // x(0) = x(s0)
			case(1):    solution[0][n_eq] = prev[n_eq]; break;  // y(0) = y(s0)

            case(2):    solution[0][n_eq] = zg; break;          // z(0) = z_g(x0,y0)

				
               
			case(3):                                                                                                // nux(0) = C1x * nux(s0) + C2x * (nuz(s0) - nuy(s0) * dzg/dy)
			case(6):                                                                                                // Xt(0) =  C1x * Xt(s0)  + C2x * (Zt(s0)  - Yt(s0)  * dzg/dy)
            case(12):   solution[0][n_eq] = C1x * prev[n_eq] + C2x * (prev[n_eq + 2] - prev[n_eq + 1] * b); break;  // Xp(0) =  C1x * Xp(s0)  + C2x * (Zp(s0)  - Yp(s0)  * dzg/dy)
				
                
            case(4):                                                                                                // nuy(0) = C1y * nuy(s0) + C2y * (nuz(s0) - nux(s0) * dzg/dx)
			case(7):                                                                                                // Yt(0) =  C1y * Yt(s0)  + C2y * (Zt(s0)  - Xt(s0)  * dzg/dx)
            case(13):   solution[0][n_eq] = C1y * prev[n_eq] + C2y * (prev[n_eq + 1] - prev[n_eq - 1] * a); break;  // Yp(0) =  C1y * Yp(s0)  + C2y * (Zp(s0)  - Xp(s0)  * dzg/dx)
                
			case(5):                                                                                                // nu_z =  -C1z * nuz(s0) + nux(s0) * C2x + nuy(s0) * C2y
			case(8):                                                                                                // Zt(0) = -C1z * Zt(s0)  + Xt(s0)  * C2x + Yt(s0)  * C2y
            case(14):   solution[0][n_eq] = -C1z * prev[n_eq] + prev[n_eq - 2] * C2x + prev[n_eq - 1] * C2y; break; // Zp(0) = -C1z * Zp(s0)  + Xp(s0)  * C2x + Xp(s0)  * C2y
                
            // mu_xt(0) =  C1x * mu_xt(s0) + ... (see notes)
            // mu_yt(0) =  C1y * mu_yt(s0) + ... (see notes)
            // mu_zt(0) = -C1z * mu_zt(s0) + ... (see notes)
            case(9):    solution[0][n_eq] =  C1x * prev[9] + C2x * (prev[11] - b * prev[10] - prev[4] * db_dth) + dC1x_dth * prev[3] + dC2x_dth * (prev[5] - b * prev[4]) + ((C1x - 1.0) * dnux_ds + C2x * (dnuz_ds - b * dnuy_ds)) * ds0_dth; break;
            case(10):   solution[0][n_eq] =  C1y * prev[10] + C2y * (prev[11] - a * prev[9] - prev[3] * da_dth) + dC1y_dth * prev[4] + dC2y_dth * (prev[5] - a * prev[3]) + ((C1y - 1.0) * dnuy_ds + C2y * (dnuz_ds - a * dnux_ds)) * ds0_dth; break;
            case(11):   solution[0][n_eq] = -C1z * prev[11] + C2x * prev[9] + C2y * prev[10] - dC1z_dth * prev[5] + dC2x_dth * prev[3] + dC2y_dth * prev[4] + (-(C1z + 1.0) * dnuz_ds + C2x * dnux_ds + C2y * dnuy_ds) * ds0_dth; break;

            // mu_xp(0) =  C1x * mu_xp(s0) + ... (see notes)
            // mu_yp(0) =  C1y * mu_yp(s0) + ... (see notes)
            // mu_zp(0) = -C1z * mu_zp(s0) + ... (see notes)
            case(15):   solution[0][n_eq] =  C1x * prev[15] + C2x * (prev[17] - b * prev[16] - prev[4] * db_dph) + dC1x_dph * prev[3] + dC2x_dph * (prev[5] - b * prev[4]) + ((C1x - 1.0) * dnux_ds + C2x * (dnuz_ds - b * dnuy_ds)) * ds0_dph; break;
            case(16):   solution[0][n_eq] =  C1y * prev[16] + C2y * (prev[17] - a * prev[15] - prev[3] * da_dph) + dC1y_dph * prev[4] + dC2y_dph * (prev[5] - a * prev[3]) + ((C1y - 1.0) * dnuy_ds + C2y * (dnuz_ds - a * dnux_ds)) * ds0_dph; break;
            case(17):   solution[0][n_eq] = -C1z * prev[17] + C2x * prev[15] + C2y * prev[16] - dC1z_dph * prev[5] + dC2x_dph * prev[3] + dC2y_dph * prev[4] + (-(C1z + 1.0) * dnuz_ds + C2x * dnux_ds + C2y * dnuy_ds) * ds0_dph; break;
		}
	}
	delete [] prev;
}

//-----------------------------------------------------------------------------------//
//-------Vary the solver step size, currently uses smaller steps near the ground-----//
//-----------------------------------------------------------------------------------//
double geoac::set_ds(double* current_values){
    double result = ds_max;
    if (current_values[2] < topo::z_bndlyr){
        result -= (ds_max - ds_min) * exp(-(current_values[2] - topo::z(current_values[0], current_values[1])) / 0.2);
    }
    return result;
}

//---------------------------------------//
//-------Update the source functions-----//
//---------------------------------------//
void geoac::update_refs(double ray_length, double* current_values){
    double x = current_values[0],		y = current_values[1], 	z = current_values[2];
	double nu[3] = {current_values[3], 	current_values[4], 		current_values[5]};

    double c_temp, dc_temp, ddc_temp;
    double dd_u_temp[6], dd_v_temp[6], dd_w_temp[6];
    double dd_c[3][3], dd_u[3][3], dd_v[3][3], dd_w[3][3];
    double X_th[3], X_ph[3], mu_th[3], mu_ph[3];

    refs.nu_mag = sqrt(pow(nu[0], 2) + pow(nu[1], 2) + pow(nu[2], 2));
    if(!calc_amp){
        double z_eval = interp::in_interval(current_values[2], atmo::c_spline.x_vals[0], atmo::c_spline.x_vals[atmo::c_spline.length - 1]);
        interp::eval_all(z_eval, atmo::c_spline, refs.c, dc_temp);
        for(int n = 0; n < 2; n++){
            refs.dc[n] = 0.0;
        }
        refs.dc[2] = dc_temp;

        atmo::u_spline.accel = atmo::c_spline.accel;
        atmo::v_spline.accel = atmo::c_spline.accel;
        atmo::calc_uvw(x, y, z, refs.u, refs.v, refs.w, refs.du, refs.dv, refs.dw);
        
        refs.cg[0] =  refs.c * nu[0] / refs.nu_mag + refs.u;
        refs.cg[1] =  refs.c * nu[1] / refs.nu_mag + refs.v;
        refs.cg[2] =  refs.c * nu[2] / refs.nu_mag + refs.w;
        refs.cg_mag = sqrt(pow(refs.cg[0], 2) + pow(refs.cg[1], 2) + pow(refs.cg[2], 2));
    } else {
        double z_eval = interp::in_interval(z, atmo::c_spline.x_vals[0], atmo::c_spline.x_vals[atmo::c_spline.length - 1]);
        interp::eval_all(z_eval, atmo::c_spline, refs.c, dc_temp, ddc_temp);
        for(int n = 0; n < 3; n++){
            refs.dc[n] = 0.0;
            for(int m = 0; m < 3; m++){
                dd_c[n][m] = 0.0;
            }
        }
        refs.dc[2] = dc_temp;
        dd_c[2][2] = ddc_temp;

        atmo::u_spline.accel = atmo::c_spline.accel;
        atmo::v_spline.accel = atmo::c_spline.accel;
        atmo::calc_uvw(x, y, z, refs.u, refs.v, refs.w, refs.du, refs.dv, refs.dw, dd_u_temp, dd_v_temp, dd_w_temp);
        
        refs.cg[0] =  refs.c * nu[0] / refs.nu_mag + refs.u;
        refs.cg[1] =  refs.c * nu[1] / refs.nu_mag + refs.v;
        refs.cg[2] =  refs.c * nu[2] / refs.nu_mag + refs.w;
        refs.cg_mag = sqrt(pow(refs.cg[0], 2) + pow(refs.cg[1], 2) + pow(refs.cg[2], 2));

        for(int n = 0; n < 3; n++){
            dd_u[n][n] = dd_u_temp[n];
            dd_v[n][n] = dd_v_temp[n];
            dd_w[n][n] = dd_w_temp[n];
        }
        dd_u[0][1] = dd_u_temp[3];  dd_u[1][0] = dd_u[0][1];
        dd_u[0][2] = dd_u_temp[4];  dd_u[2][0] = dd_u[0][2];
        dd_u[1][2] = dd_u_temp[5];  dd_u[2][1] = dd_u[1][2];
        
        dd_v[0][1] = dd_v_temp[3];  dd_v[1][0] = dd_v[0][1];
        dd_v[0][2] = dd_v_temp[4];  dd_v[2][0] = dd_v[0][2];
        dd_v[1][2] = dd_v_temp[5];  dd_v[2][1] = dd_v[1][2];
        
        dd_w[0][1] = dd_w_temp[3];  dd_w[1][0] = dd_w[0][1];
        dd_w[0][2] = dd_w_temp[4];  dd_w[2][0] = dd_w[0][2];
        dd_w[1][2] = dd_w_temp[5];  dd_w[2][1] = dd_w[1][2];

        X_th[0]  = current_values[6];       X_th[1]  = current_values[7];       X_th[2]  = current_values[8];
        mu_th[0] = current_values[9];		mu_th[1] = current_values[10];		mu_th[2] = current_values[11];
        X_ph[0]  = current_values[12];      X_ph[1]  = current_values[13];      X_ph[2]  = current_values[14];
        mu_ph[0] = current_values[15];		mu_ph[1] = current_values[16];		mu_ph[2] = current_values[17];
        
        refs.dc[3] = 0.0;  refs.du[3] = 0.0;  refs.dv[3] = 0.0;  refs.dw[3] = 0.0;
        refs.dc[4] = 0.0;  refs.du[4] = 0.0;  refs.dv[4] = 0.0;  refs.dw[4] = 0.0;
        for(int n = 0; n < 3; n++){
            refs.dc[3] += X_th[n] * refs.dc[n]; refs.du[3] += X_th[n] * refs.du[n]; refs.dv[3] += X_th[n] * refs.dv[n]; refs.dw[3] += X_th[n] * refs.dw[n];
            refs.dc[4] += X_ph[n] * refs.dc[n]; refs.du[4] += X_ph[n] * refs.du[n]; refs.dv[3] += X_ph[n] * refs.dv[n]; refs.dw[4] += X_ph[n] * refs.dw[n];
                
            refs.ddc[n][0] = 0.0;   refs.ddc[n][1] = 0.0;
            refs.ddu[n][0] = 0.0;   refs.ddu[n][1] = 0.0;
            refs.ddv[n][0] = 0.0;   refs.ddv[n][1] = 0.0;
            refs.ddw[n][0] = 0.0;   refs.ddw[n][1] = 0.0;
                
            for(int m = n; m < 3; m++){
                refs.ddc[n][0] += X_th[m] * dd_c[n][m];     refs.ddc[n][1] += X_ph[m] * dd_c[n][m];
                refs.ddu[n][0] += X_th[m] * dd_u[n][m];     refs.ddu[n][1] += X_ph[m] * dd_u[n][m];
                refs.ddv[n][0] += X_th[m] * dd_v[n][m];     refs.ddv[n][1] += X_ph[m] * dd_v[n][m];
                refs.ddw[n][0] += X_th[m] * dd_w[n][m];     refs.ddw[n][1] += X_ph[m] * dd_w[n][m];
            }
        }
        
        refs.dnu_mag[0] = (nu[0] * mu_th[0] + nu[1] * mu_th[1] + nu[2] * mu_th[2]) / refs.nu_mag;
        refs.dnu_mag[1] = (nu[0] * mu_ph[0] + nu[1] * mu_ph[1] + nu[2] * mu_ph[2]) / refs.nu_mag;
            
        for (int n = 0; n < 3; n++){
            refs.dcg[n][0] = nu[n] / refs.nu_mag * refs.dc[3] + refs.c * mu_th[n] / refs.nu_mag - refs.c * nu[n] / pow(refs.nu_mag,2) * refs.dnu_mag[0];
            refs.dcg[n][1] = nu[n] / refs.nu_mag * refs.dc[4] + refs.c * mu_ph[n] / refs.nu_mag - refs.c * nu[n] / pow(refs.nu_mag,2) * refs.dnu_mag[1];
        }
            
        refs.dcg[0][0] += refs.du[3]; refs.dcg[0][1] += refs.du[4];
        refs.dcg[1][0] += refs.dv[3]; refs.dcg[1][1] += refs.dv[4];
        refs.dcg[2][0] += refs.dw[3]; refs.dcg[2][1] += refs.dw[4];
            
        refs.dcg_mag[0] = (refs.cg[0] * refs.dcg[0][0] + refs.cg[1] * refs.dcg[1][0] + refs.cg[2] * refs.dcg[2][0]) / refs.cg_mag;
        refs.dcg_mag[1] = (refs.cg[0] * refs.dcg[0][1] + refs.cg[1] * refs.dcg[1][1] + refs.cg[2] * refs.dcg[2][1]) / refs.cg_mag;
    }
}

//-----------------------------------------------------------//
//-------Evaluate the Source Equation For Specific Index-----//
//-----------------------------------------------------------//
double geoac::eval_src_eq(double ray_length, double* current_values, int n_eq){
    double result = 0.0;
    
    double nu[3], mu_th[3], mu_ph[3];
    nu[0] = current_values[3]; nu[1] = current_values[4]; nu[2] = current_values[5];
    if (calc_amp){
        mu_th[0] = current_values[9];   mu_th[1] = current_values[10];  mu_th[2] = current_values[11];
        mu_ph[0] = current_values[15];  mu_ph[1] = current_values[16];  mu_ph[2] = current_values[17];
    }
    
	switch(n_eq){
		case(0):	// d x / d s
        case(1):    // d y / d s
        case(2):    // d z / d s
			result = refs.cg[n_eq] / refs.cg_mag;
            break;
            
		case(3): 	// d nu_x /ds
		case(4): 	// d nu_y / ds
		case(5): 	// d nu_z / ds
            if (is_topo || n_eq == 5){      result = -1.0 / refs.cg_mag * (refs.nu_mag * refs.dc[n_eq - 3]
                                                                             + nu[0] * refs.du[n_eq - 3] + nu[1] * refs.dv[n_eq - 3] + nu[2] * refs.dw[n_eq - 3]);}
            break;
            
		case(6):	// d X_th/ds
        case(7):	// d Y_th/ds
        case(8): 	// d Z_th/ds
			result = refs.dcg[n_eq - 6][0] / refs.cg_mag
                    - refs.cg[n_eq - 6] / pow(refs.cg_mag, 2) * refs.dcg_mag[0];
            break;
            
		case(9): 	// d mu_x_th /ds
		case(10): 	// d mu_y_th /ds
		case(11): 	// d mu_z_th /ds
            if (is_topo || n_eq == 11){ result = refs.dcg_mag[0] / pow(refs.cg_mag, 2) * (refs.nu_mag * refs.dc[n_eq - 9] + nu[0] * refs.du[n_eq - 9] + nu[1] * refs.dv[n_eq - 9] + nu[2] * refs.dw[n_eq - 9])
                                                    - 1.0 / refs.cg_mag * (refs.dnu_mag[0] * refs.dc[n_eq - 9] + refs.nu_mag * refs.ddc[n_eq - 9][0] + mu_th[0] * refs.du[n_eq - 9] + mu_th[1] * refs.dv[n_eq - 9] + mu_th[2] * refs.dw[n_eq - 9]
                                                                                + nu[0] * refs.ddu[n_eq - 9][0] + nu[1] * refs.ddv[n_eq - 9][0] + nu[2] * refs.ddw[n_eq - 9][0]);}
            break;
            
        case(12):  	// d X_ph/ds
        case(13):  	// d Y_ph/ds
        case(14):	// d Z_ph/ds
            result = refs.dcg[n_eq - 12][1] / refs.cg_mag
                    - refs.cg[n_eq - 12] / pow(refs.cg_mag, 2) * refs.dcg_mag[1];
            break;

		case(15): 	// d mu_x_ph/ds
		case(16): 	// d mu_y_ph/ds
        case(17): 	// d mu_z_ph/ds
            if (is_topo || n_eq == 17){     result = refs.dcg_mag[1] / pow(refs.cg_mag, 2) * (refs.nu_mag * refs.dc[n_eq - 15] + nu[0] * refs.du[n_eq - 15] + nu[1] * refs.dv[n_eq - 15] + nu[2] * refs.dw[n_eq - 15])
                                                        - 1.0 / refs.cg_mag * (refs.dnu_mag[1] * refs.dc[n_eq - 15] + refs.nu_mag * refs.ddc[n_eq - 15][1] + mu_ph[0] * refs.du[n_eq - 15] + mu_ph[1] * refs.dv[n_eq - 15] + mu_ph[2] * refs.dw[n_eq - 15]
                                                                                    + nu[0] * refs.ddu[n_eq - 15][1] + nu[1] * refs.ddv[n_eq - 15][1] + nu[2] * refs.ddw[n_eq - 15][1]);}
            break;
            
	}
    return result;
}
                                      

//-------------------------------------------------------------------//
//-------Calculate the Hamiltonian (Eikonal) To Check For Errors-----//
//-------------------------------------------------------------------//
double geoac::eval_eikonal(double ray_length, double* current_values){
	double x = current_values[0],  y = current_values[1], z = current_values[2];
	double nu[3] = {current_values[3], current_values[4], current_values[5]};
	
	return sqrt(nu[0] * nu[0] + nu[1] * nu[1] + nu[2] * nu[2]) - refs.c0 / atmo::c(x, y, z)
                    + (atmo::u(x, y, z) * nu[0] + atmo::v(x ,y ,z) * nu[1] + atmo::w(x, y, z) * nu[2]) / atmo::c(x, y, z);
}

double geoac::eval_eikonal(double** solution, int index){
	double x = solution[index][0],          y = solution[index][1],         z = solution[index][2];
	double nu[3] = {solution[index][3],     solution[index][4],             solution[index][5]};
    
	return sqrt(nu[0] * nu[0] + nu[1] * nu[1] + nu[2] * nu[2]) - refs.c0 / atmo::c(x, y, z)
                    + (atmo::u(x, y, z) * nu[0] + atmo::v(x ,y ,z) * nu[1] + atmo::w(x, y, z) * nu[2]) / atmo::c(x, y, z);
}

double geoac::eval_eikonal_deriv(double** solution, int index){
	double  x =           solution[index][0],y = solution[index][1], z = solution[index][2],
            X_th[3] =    {solution[index][6],    solution[index][7],     solution[index][8]},
            X_ph[3] =    {solution[index][12],   solution[index][13],    solution[index][14]},
            nu[3] =      {solution[index][3],    solution[index][4],     solution[index][5]},
            mu_th[3] =   {solution[index][9],    solution[index][10],    solution[index][11]},
            mu_ph[3] =   {solution[index][15],   solution[index][16],    solution[index][17]};
    
    double mag_nu = sqrt(nu[0] * nu[0] + nu[1] * nu[1] + nu[2] * nu[2]);
    
    double dc_dlt = X_th[0] * atmo::dc(x,y,z,0) + X_th[1] * atmo::dc(x,y,z,1) + X_th[2]*atmo::dc(x,y,z,2);    double dc_dlp = X_ph[0] * atmo::dc(x,y,z,0) + X_ph[1] * atmo::dc(x,y,z,1) + X_ph[2] * atmo::dc(x,y,z,2);
    double du_dlt = X_th[0] * atmo::du(x,y,z,0) + X_th[1] * atmo::du(x,y,z,1) + X_th[2]*atmo::du(x,y,z,2);    double du_dlp = X_ph[0] * atmo::du(x,y,z,0) + X_ph[1] * atmo::du(x,y,z,1) + X_ph[2] * atmo::du(x,y,z,2);
    double dv_dlt = X_th[0] * atmo::dv(x,y,z,0) + X_th[1] * atmo::dv(x,y,z,1) + X_th[2]*atmo::dv(x,y,z,2);    double dv_dlp = X_ph[0] * atmo::dv(x,y,z,0) + X_ph[1] * atmo::dv(x,y,z,1) + X_ph[2] * atmo::dv(x,y,z,2);
    double dw_dlt = X_th[0] * atmo::dw(x,y,z,0) + X_th[1] * atmo::dw(x,y,z,1) + X_th[2]*atmo::dw(x,y,z,2);    double dw_dlp = X_ph[0] * atmo::dw(x,y,z,0) + X_ph[1] * atmo::dw(x,y,z,1) + X_ph[2] * atmo::dw(x,y,z,2);
    
    double Resid_th = (nu[0] * mu_th[0] + nu[1] * mu_th[1] + nu[2] * mu_th[2]) / mag_nu + mag_nu / atmo::c(x,y,z) * dc_dlt
                            + 1.0 / atmo::c(x,y,z) * (mu_th[0] * atmo::u(x,y,z) + mu_th[1] * atmo::v(x,y,z) + mu_th[2] * atmo::w(x,y,z) + nu[0] * du_dlt + nu[1] * dv_dlt + nu[2] * dw_dlt);
    double Resid_ph = (nu[0] * mu_ph[0] + nu[1] * mu_ph[1] + nu[2] * mu_ph[2]) / mag_nu + mag_nu / atmo::c(x,y,z) * dc_dlp
                            + 1.0 / atmo::c(x,y,z) * (mu_ph[0] * atmo::u(x,y,z) + mu_ph[1] * atmo::v(x,y,z) + mu_ph[2] * atmo::w(x,y,z) + nu[0] * du_dlp + nu[1] * dv_dlp + nu[2] * dw_dlp);
    
    return sqrt(Resid_th * Resid_th + Resid_ph * Resid_ph);
}

//---------------------------------------------------------//
//--------Check if the ray has left the propagation--------//
//-------------region or returned to the ground------------//
//---------------------------------------------------------//
bool geoac::break_check(double ** & solution, int k){
    if(sqrt(pow(solution[k][0] - refs.src_loc[0], 2)
            + pow(solution[k][1] - refs.src_loc[1], 2)) >  rng_max){    return true;}
	if(solution[k][0] > x_max || solution[k][0] < x_min){               return true;}
    if(solution[k][1] > y_max || solution[k][1] < y_min){               return true;}
    if(solution[k][2] > alt_max){                                       return true;}
    
    return false;
}

bool geoac::ground_check(double ** solution, int k){
    if (solution[k][2] <= topo::z_max){
        if(solution[k][2] < topo::z(solution[k][0], solution[k][1])){ return true;}
    }
    return false;
}

bool geoac::reflect_check(double ** solution, int k){
    return false;
}


//----------------------------------------------------------------------------------//
//-------Calculate the travel time from source to location or between locations-----//
//----------------------------------------------------------------------------------//
double geoac::travel_time(double ** solution, int k){
    double dx, dy, dz, ds, x, y, z, nu[3], nu_mag, cg[3], cg_mag;
    double sndspd;
	double traveltime = 0;
	
	for (int n = 0; n < k; n++){
		dx = solution[n+1][0] - solution[n][0];	x = solution[n][0] + dx/2.0;
        dy = solution[n+1][1] - solution[n][1]; y = solution[n][1] + dy/2.0;
		dz = solution[n+1][2] - solution[n][2]; z = solution[n][2] + dz/2.0;
		ds = sqrt(dx * dx + dy * dy + dz * dz);
		
		nu[0] = solution[n][3] + (solution[n+1][3] - solution[n][3])/2.0;
		nu[1] = solution[n][4] + (solution[n+1][4] - solution[n][4])/2.0;
		nu[2] = solution[n][5] + (solution[n+1][5] - solution[n][5])/2.0;
		nu_mag = sqrt(nu[0] * nu[0] + nu[1] * nu[1] + nu[2] * nu[2]);
        
        sndspd = atmo::c(x, y, z);
		cg[0] = sndspd * nu[0] / nu_mag + atmo::u(x, y, z);
		cg[1] = sndspd * nu[1] / nu_mag + atmo::v(x, y, z);
		cg[2] = sndspd * nu[2] / nu_mag + atmo::w(x, y, z);
		cg_mag = sqrt(pow(cg[0], 2) + pow(cg[1], 2) + pow(cg[2], 2));
        
		traveltime += ds/cg_mag;
    }
	return traveltime;
}

void geoac::travel_time(double & t, double ** solution, int k1, int k2){
	double dx, dy, dz, ds, x, y, z, nu[3], nu_mag, cg[3], cg_mag;
    double sndspd;
    
	for (int n = k1; n < k2; n++){
		dx = solution[n+1][0] - solution[n][0];	x = solution[n][0] + dx/2.0;
        dy = solution[n+1][1] - solution[n][1]; y = solution[n][1] + dy/2.0;
		dz = solution[n+1][2] - solution[n][2]; z = solution[n][2] + dz/2.0;
		ds = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
        
		nu[0] = solution[n][3] + (solution[n+1][3] - solution[n][3])/2.0;
		nu[1] = solution[n][4] + (solution[n+1][4] - solution[n][4])/2.0;
		nu[2] = solution[n][5] + (solution[n+1][5] - solution[n][5])/2.0;
		nu_mag = sqrt(pow(nu[0], 2) + pow(nu[1], 2) + pow(nu[2], 2));
        
        sndspd = atmo::c(x, y, z);
		cg[0] = sndspd * nu[0] / nu_mag + atmo::u(x, y, z);
		cg[1] = sndspd * nu[1] / nu_mag + atmo::v(x, y, z);
		cg[2] = sndspd * nu[2] / nu_mag + atmo::w(x, y, z);
		cg_mag = sqrt(pow(cg[0], 2) + pow(cg[1], 2) + pow(cg[2], 2));
        
		t += ds/cg_mag;
	}
}

//------------------------------------------------//
//-------Calculate the travel time variance-------//
//------------------------------------------------//
void geoac::travel_time_var(double ** solution, int k, double & tt, double & tt_var_incl, double & tt_var_az){
    double dx, dy, dz, ds, x, y, z, nu[3], nu_mag, sndspd, cg[3], cg_mag;
    double X[3], mu[3], sndspd_scalar, dcg[3];

	tt = 0.0;
    tt_var_incl = 0.0;
    tt_var_az = 0.0;
	for (int n = 0; n < k; n++){
        // Calculate  ds along with x, y, z at mid point
		dx = solution[n + 1][0] - solution[n][0]; x = solution[n][0] + dx / 2.0;
        dy = solution[n + 1][1] - solution[n][1]; y = solution[n][1] + dy / 2.0;
		dz = solution[n + 1][2] - solution[n][2]; z = solution[n][2] + dz / 2.0;
		ds = sqrt(dx * dx + dy * dy + dz * dz);
		
        // Calculate c_g to define travel time contribution
		nu[0] = solution[n][3] + (solution[n + 1][3] - solution[n][3]) / 2.0;
		nu[1] = solution[n][4] + (solution[n + 1][4] - solution[n][4]) / 2.0;
		nu[2] = solution[n][5] + (solution[n + 1][5] - solution[n][5]) / 2.0;
		nu_mag = sqrt(nu[0] * nu[0] + nu[1] * nu[1] + nu[2] * nu[2]);
        
        sndspd = atmo::c(x, y, z);
		cg[0] = sndspd * nu[0] / nu_mag + atmo::u(x, y, z);
		cg[1] = sndspd * nu[1] / nu_mag + atmo::v(x, y, z);
		cg[2] = sndspd * nu[2] / nu_mag + atmo::w(x, y, z);
		cg_mag = sqrt(pow(cg[0], 2) + pow(cg[1], 2) + pow(cg[2], 2));

		tt += ds / cg_mag;

        // Calculate the inclination variation
		X[0] = solution[n][6] + (solution[n + 1][6] - solution[n][6]) / 2.0;
		X[1] = solution[n][7] + (solution[n + 1][7] - solution[n][7]) / 2.0;
		X[2] = solution[n][8] + (solution[n + 1][8] - solution[n][8]) / 2.0;

		mu[0] = solution[n][9]  + (solution[n + 1][9]  - solution[n][9]) / 2.0;
		mu[1] = solution[n][10] + (solution[n + 1][10] - solution[n][10]) / 2.0;
		mu[2] = solution[n][11] + (solution[n + 1][11] - solution[n][11]) / 2.0;

        sndspd_scalar  = X[0] * atmo::dc(x, y, z, 0);
        sndspd_scalar += X[1] * atmo::dc(x, y, z, 1);
        sndspd_scalar += X[2] * atmo::dc(x, y, z, 2);
        sndspd_scalar -= sndspd / pow(nu_mag, 2) * (mu[0] * nu[0] + mu[1] * nu[1] + mu[2] * nu[2]);

        dcg[0] = nu[0] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[0];
        dcg[1] = nu[1] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[1];
        dcg[2] = nu[2] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[2];

        dcg[0] += X[0] * atmo::du(x, y, z, 0) + X[1] * atmo::du(x, y, z, 1) + X[2] * atmo::du(x, y, z, 2);
        dcg[1] += X[0] * atmo::dv(x, y, z, 0) + X[1] * atmo::dv(x, y, z, 1) + X[2] * atmo::dv(x, y, z, 2);
        dcg[2] += X[0] * atmo::dw(x, y, z, 0) + X[1] * atmo::dw(x, y, z, 1) + X[2] * atmo::dw(x, y, z, 2);

        tt_var_incl += (dcg[0] * cg[0] + dcg[1] * cg[1] * dcg[2] * cg[2]) / pow(cg_mag, 3) * ds;

        // Repeat for azimuth variation
		X[0] = solution[n][12] + (solution[n + 1][12] - solution[n][12]) / 2.0;
		X[1] = solution[n][13] + (solution[n + 1][13] - solution[n][13]) / 2.0;
		X[2] = solution[n][14] + (solution[n + 1][14] - solution[n][14]) / 2.0;

		mu[0] = solution[n][15] + (solution[n + 1][15] - solution[n][15]) / 2.0;
		mu[1] = solution[n][16] + (solution[n + 1][16] - solution[n][16]) / 2.0;
		mu[2] = solution[n][17] + (solution[n + 1][17] - solution[n][17]) / 2.0;

        sndspd_scalar  = X[0] * atmo::dc(x, y, z, 0);
        sndspd_scalar += X[1] * atmo::dc(x, y, z, 1);
        sndspd_scalar += X[2] * atmo::dc(x, y, z, 2);
        sndspd_scalar -= sndspd / pow(nu_mag, 2) * (mu[0] * nu[0] + mu[1] * nu[1] + mu[2] * nu[2]);

        dcg[0] = nu[0] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[0];
        dcg[1] = nu[1] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[1];
        dcg[2] = nu[2] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[2];

        dcg[0] += X[0] * atmo::du(x, y, z, 0) + X[1] * atmo::du(x, y, z, 1) + X[2] * atmo::du(x, y, z, 2);
        dcg[1] += X[0] * atmo::dv(x, y, z, 0) + X[1] * atmo::dv(x, y, z, 1) + X[2] * atmo::dv(x, y, z, 2);
        dcg[2] += X[0] * atmo::dw(x, y, z, 0) + X[1] * atmo::dw(x, y, z, 1) + X[2] * atmo::dw(x, y, z, 2);

        tt_var_az += (dcg[0] * cg[0] + dcg[1] * cg[1] * dcg[2] * cg[2]) / pow(cg_mag, 3) * ds;
    }
}


void geoac::travel_time_var(double & tt, double & tt_var_incl, double & tt_var_az, double** solution, int k1, int k2){
   double dx, dy, dz, ds, x, y, z, nu[3], nu_mag, sndspd, cg[3], cg_mag;
    double X[3], mu[3], sndspd_scalar, dcg[3];

	for (int n = k1; n < k2; n++){
        // Calculate  ds along with x, y, z at mid point
		dx = solution[n + 1][0] - solution[n][0]; x = solution[n][0] + dx / 2.0;
        dy = solution[n + 1][1] - solution[n][1]; y = solution[n][1] + dy / 2.0;
		dz = solution[n + 1][2] - solution[n][2]; z = solution[n][2] + dz / 2.0;
		ds = sqrt(dx * dx + dy * dy + dz * dz);
		
        // Calculate c_g to define travel time contribution
		nu[0] = solution[n][3] + (solution[n + 1][3] - solution[n][3]) / 2.0;
		nu[1] = solution[n][4] + (solution[n + 1][4] - solution[n][4]) / 2.0;
		nu[2] = solution[n][5] + (solution[n + 1][5] - solution[n][5]) / 2.0;
		nu_mag = sqrt(nu[0] * nu[0] + nu[1] * nu[1] + nu[2] * nu[2]);
        
        sndspd = atmo::c(x, y, z);
		cg[0] = sndspd * nu[0] / nu_mag + atmo::u(x, y, z);
		cg[1] = sndspd * nu[1] / nu_mag + atmo::v(x, y, z);
		cg[2] = sndspd * nu[2] / nu_mag + atmo::w(x, y, z);
		cg_mag = sqrt(pow(cg[0], 2) + pow(cg[1], 2) + pow(cg[2], 2));

		tt += ds / cg_mag;

        // Calculate the inclination variation
		X[0] = solution[n][6] + (solution[n + 1][6] - solution[n][6]) / 2.0;
		X[1] = solution[n][7] + (solution[n + 1][7] - solution[n][7]) / 2.0;
		X[2] = solution[n][8] + (solution[n + 1][8] - solution[n][8]) / 2.0;

		mu[0] = solution[n][9]  + (solution[n + 1][9]  - solution[n][9]) / 2.0;
		mu[1] = solution[n][10] + (solution[n + 1][10] - solution[n][10]) / 2.0;
		mu[2] = solution[n][11] + (solution[n + 1][11] - solution[n][11]) / 2.0;

        sndspd_scalar  = X[0] * atmo::dc(x, y, z, 0);
        sndspd_scalar += X[1] * atmo::dc(x, y, z, 1);
        sndspd_scalar += X[2] * atmo::dc(x, y, z, 2);
        sndspd_scalar -= sndspd / pow(nu_mag, 2) * (mu[0] * nu[0] + mu[1] * nu[1] + mu[2] * nu[2]);

        dcg[0] = nu[0] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[0];
        dcg[1] = nu[1] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[1];
        dcg[2] = nu[2] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[2];

        dcg[0] += X[0] * atmo::du(x, y, z, 0) + X[1] * atmo::du(x, y, z, 1) + X[2] * atmo::du(x, y, z, 2);
        dcg[1] += X[0] * atmo::dv(x, y, z, 0) + X[1] * atmo::dv(x, y, z, 1) + X[2] * atmo::dv(x, y, z, 2);
        dcg[2] += X[0] * atmo::dw(x, y, z, 0) + X[1] * atmo::dw(x, y, z, 1) + X[2] * atmo::dw(x, y, z, 2);

        tt_var_incl -= (dcg[0] * cg[0] + dcg[1] * cg[1] * dcg[2] * cg[2]) / pow(cg_mag, 3) * ds;

        // Repeat for azimuth variation
		X[0] = solution[n][12] + (solution[n + 1][12] - solution[n][12]) / 2.0;
		X[1] = solution[n][13] + (solution[n + 1][13] - solution[n][13]) / 2.0;
		X[2] = solution[n][14] + (solution[n + 1][14] - solution[n][14]) / 2.0;

		mu[0] = solution[n][15] + (solution[n + 1][15] - solution[n][15]) / 2.0;
		mu[1] = solution[n][16] + (solution[n + 1][16] - solution[n][16]) / 2.0;
		mu[2] = solution[n][17] + (solution[n + 1][17] - solution[n][17]) / 2.0;

        sndspd_scalar  = X[0] * atmo::dc(x, y, z, 0);
        sndspd_scalar += X[1] * atmo::dc(x, y, z, 1);
        sndspd_scalar += X[2] * atmo::dc(x, y, z, 2);
        sndspd_scalar -= sndspd / pow(nu_mag, 2) * (mu[0] * nu[0] + mu[1] * nu[1] + mu[2] * nu[2]);

        dcg[0] = nu[0] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[0];
        dcg[1] = nu[1] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[1];
        dcg[2] = nu[2] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[2];

        dcg[0] += X[0] * atmo::du(x, y, z, 0) + X[1] * atmo::du(x, y, z, 1) + X[2] * atmo::du(x, y, z, 2);
        dcg[1] += X[0] * atmo::dv(x, y, z, 0) + X[1] * atmo::dv(x, y, z, 1) + X[2] * atmo::dv(x, y, z, 2);
        dcg[2] += X[0] * atmo::dw(x, y, z, 0) + X[1] * atmo::dw(x, y, z, 1) + X[2] * atmo::dw(x, y, z, 2);

        tt_var_az -= (dcg[0] * cg[0] + dcg[1] * cg[1] * dcg[2] * cg[2]) / pow(cg_mag, 3) * ds;
    }
}


//--------------------------------------------------------------------------//
//------Integrate the Sutherland Bass Attenuation Through the Ray Path------//
//--------------------------------------------------------------------------//
double geoac::atten(double ** solution, int k, double freq){
	double dx, dy, dz, ds, x, y, z, atten = 0.0;
	for (int n = 0; n < k; n++){
		dx = solution[n + 1][0] - solution[n][0]; x = solution[n][0] + dx / 2.0;
        dy = solution[n + 1][1] - solution[n][1]; y = solution[n][1] + dy / 2.0;
        dz = solution[n + 1][2] - solution[n][2]; z = solution[n][2] + dz / 2.0;
		ds = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
        
		atten += atmo::SB_alpha(x, y, z, freq) * ds;
	}
    return atten;
}

void geoac::atten(double & atten, double ** solution, int start, int end, double freq){
	double dx, dy, dz, ds, x, y, z;
	for (int n = start; n < end; n++){
		dx = solution[n + 1][0] - solution[n][0]; x = solution[n][0] + dx / 2.0;
        dy = solution[n + 1][1] - solution[n][1]; y = solution[n][1] + dy / 2.0;
        dz = solution[n + 1][2] - solution[n][2]; z = solution[n][2] + dz / 2.0;
		ds = sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
        
		atten += atmo::SB_alpha(x, y, z, freq) * ds;
	}
}

//-----------------------------------------------------------------------------------//
//-------Calculate the Jacobian determinant and from it the ampltude coefficient-----//
//-----------------------------------------------------------------------------------//
double geoac::jacobian(double ** solution, int k){
    double sndspd, nu_mag, cg[3], cg_mag;
    double dxds, dyds, dzds, dxdth, dydth, dzdth, dxdph, dydph, dzdph;
    
    sndspd = atmo::c(solution[k][0], solution[k][1], solution[k][2]);
    nu_mag = sqrt(pow(solution[k][3], 2) + pow(solution[k][4], 2) + pow(solution[k][5], 2));

	cg[0] = sndspd * solution[k][3] / nu_mag + atmo::u(solution[k][0], solution[k][1], solution[k][2]);
    cg[1] = sndspd * solution[k][4] / nu_mag + atmo::v(solution[k][0], solution[k][1], solution[k][2]);
    cg[2] = sndspd * solution[k][5] / nu_mag + atmo::w(solution[k][0], solution[k][1], solution[k][2]);
	cg_mag = sqrt(pow(cg[0], 2) + pow(cg[1], 2) + pow(cg[2], 2));
    
    dxds = cg[0] / cg_mag;  dyds = cg[1] / cg_mag;  dzds = cg[2] / cg_mag;
	dxdth = solution[k][6];         dydth = solution[k][7];         dzdth = solution[k][8];
	dxdph = solution[k][12];        dydph = solution[k][13];        dzdph = solution[k][14];
    
	return	dxds * (dydth * dzdph - dydph * dzdth)
            - dxdth * (dyds * dzdph - dzds * dydph)
                + dxdph * (dyds * dzdth - dzds * dydth);
}

double geoac::amp(double ** solution, int k){
    double sndspd, winds[3], winds0[3], nu_mag, nu_mag0, cg[3], cg0[3], cg_mag, cg_mag0, num, den;
    
    sndspd = atmo::c(solution[k][0], solution[k][1], solution[k][2]);
    winds[0] = atmo::u(solution[k][0], solution[k][1], solution[k][2]); winds0[0] = atmo::u(refs.src_loc[0], refs.src_loc[1], refs.src_loc[2]),
    winds[1] = atmo::v(solution[k][0], solution[k][1], solution[k][2]); winds0[1] = atmo::v(refs.src_loc[0], refs.src_loc[1], refs.src_loc[2]),
    winds[2] = atmo::w(solution[k][0], solution[k][1], solution[k][2]); winds0[2] = atmo::w(refs.src_loc[0], refs.src_loc[1], refs.src_loc[2]),

    nu_mag = sqrt(pow(solution[k][3], 2) + pow(solution[k][4], 2) + pow(solution[k][5], 2));
    nu_mag0 = 1.0 - (solution[k][3] * winds0[0] + solution[k][4] * winds0[1] + solution[k][5] * winds0[2]) / refs.c0;
    
    for (int n = 0; n < 3; n++){ cg[n] =  sndspd * solution[k][3 + n] / nu_mag + winds[n];}
    cg0[0] = refs.c0 * cos(theta) * cos(phi) + winds0[0];
    cg0[1] = refs.c0 * cos(theta) * sin(phi) + winds0[1];
    cg0[2] = refs.c0 * sin(theta) + winds0[2];
    
	cg_mag = sqrt(pow(cg[0], 2) + pow(cg[1], 2) + pow(cg[2], 2));
    cg_mag0 = sqrt(pow(cg0[0], 2) + pow(cg0[1], 2) + pow(cg0[2], 2));
	
	num = atmo::rho(solution[k][0], solution[k][1], solution[k][2]) * nu_mag * pow(sndspd, 3) * cg_mag0 * cos(theta);
	den = atmo::rho(refs.src_loc[0], refs.src_loc[1], refs.src_loc[2]) * nu_mag0* pow(refs.c0,3) * cg_mag * jacobian(solution, k);
    
	return sqrt( fabs(num / den));
}


//------------------------------------------------//
//-------Calculate Weakly Non-Linear Waveform-----//
//------------------------------------------------//
double geoac::wnl_wvfrm(double** solution, double** & u, int k1, int k2, double s_prev, double c0, double cg0, double nu0, double rho0, double D0, double p0, double output_step){
    // NOTE: function assumes specified u(t, s1) is the scaled pressure waveform.

    int fft_len = wvfrm::len / 2 + 1;
    double s, ds, s_tot, tukey_param;
    double beta0, dt, t0, Uf_max, beta_curr, beta_next;

    char output_buffer [512];
    ofstream file_out;

    // Set up the fft time series, spectrum, and plans
    double* fft_time;           fft_time = new double [wvfrm::len];
	fftw_complex* fft_spec;     fft_spec = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * fft_len);
	fftw_plan fwd_plan;         fwd_plan = fftw_plan_dft_r2c_1d(wvfrm::len, fft_time, fft_spec, FFTW_MEASURE);
    fftw_plan bwd_plan;         bwd_plan = fftw_plan_dft_c2r_1d(wvfrm::len, fft_spec, fft_time, FFTW_MEASURE);

    // Define previous u(t), and current, squared, and two-steps ago FFT, U(f).
    double* u_sqr;      u_sqr = new double [wvfrm::len];
    double* du_sqr;     du_sqr = new double [wvfrm::len];
    double* tukey_win;  tukey_win = new double [wvfrm::len];

    double* filter1;     filter1 = new double [fft_len];
    double* filter2;     filter2 = new double [fft_len];
    
    double** U;         U = new double* [fft_len];
    double** A;         A = new double* [fft_len];
    double** B;         B = new double* [fft_len];
    double** C;         C = new double* [fft_len];

    for(int n = 0; n < fft_len; n++){
        U[n] = new double [2];
        A[n] = new double [2];
        B[n] = new double [2];
        C[n] = new double [2];
    }

    // Convert density [g / cm^3] to [kg / m^3] by factor of 1.0e3 and sound speed [km/s] to [m/s]
    // by the same factor in order to balance p0 being specified in [Pa] = [kg / (m s^2)]; beta in air is 1.2
    beta0 = 1.2 * p0 / ((rho0 * 1.0e3) * pow(c0 * 1.0e3, 2)) * (nu0 * c0) / (cg0 * atmo::c(refs.src_loc[0], refs.src_loc[1], refs.src_loc[2]));

    struct interp::linear_spline_1D ray_x;     interp::prep(ray_x, k2 - k1);
    struct interp::linear_spline_1D ray_y;     interp::prep(ray_y, k2 - k1);
    struct interp::linear_spline_1D ray_z;     interp::prep(ray_z, k2 - k1);

    struct interp::linear_spline_1D D;               interp::prep(D, k2 - k1);
    struct interp::linear_spline_1D beta;            interp::prep(beta, k2 - k1);
    struct interp::linear_spline_1D beta_caustic;    interp::prep(beta_caustic, k2 - k1);

    s = 0.0;
    for (int k = 0; k < (k2 - k1); k++) {
        double c = atmo::c(solution[k1 + k][0], solution[k1 + k][1], solution[k1 + k][2]);
        double rho = atmo::rho(solution[k1 + k][0], solution[k1 + k][1], solution[k1 + k][2]);

        double nu = sqrt(pow(solution[k1 + k][3], 2) + pow(solution[k1 + k][4], 2) + pow(solution[k1 + k][5], 2));

        double cg_x = c * solution[k1 + k][3] / nu * atmo::u(solution[k1 + k][0], solution[k1 + k][1], solution[k1 + k][2]);
        double cg_y = c * solution[k1 + k][4] / nu * atmo::v(solution[k1 + k][0], solution[k1 + k][1], solution[k1 + k][2]);
        double cg_z = c * solution[k1 + k][5] / nu * atmo::w(solution[k1 + k][0], solution[k1 + k][1], solution[k1 + k][2]);
        double cg = sqrt(pow(cg_x, 2) + pow(cg_y, 2) + pow(cg_z, 2));

        ray_x.x_vals[k] = s;    ray_x.f_vals[k] = solution[k1 + k][0];
        ray_y.x_vals[k] = s;    ray_y.f_vals[k] = solution[k1 + k][1];
        ray_z.x_vals[k] = s;    ray_z.f_vals[k] = solution[k1 + k][2];

        D.x_vals[k] = s;        D.f_vals[k] = jacobian(solution, k1 + k);
        beta.x_vals[k] = s;     beta.f_vals[k] = beta0 * sqrt(fabs((D0 * rho0 * pow(cg0, 3)) / (D.f_vals[k] * rho * pow(c, 3)) * (c * pow(nu, 3)) / (c0 * pow(nu0, 3))));

        beta_caustic.x_vals[k] = s;
        if(k < 1){  beta_caustic.f_vals[k] = 0.0;}
        else {      beta_caustic.f_vals[k] = beta0 * sqrt(fabs((D0 * rho0 * pow(cg0, 3)) / ((D.f_vals[k] - D.f_vals[k - 1]) / (D.x_vals[k] - D.x_vals[k - 1]) * rho * pow(c, 3)) * (c * pow(nu, 3)) / (c0 * pow(nu0, 3))));}


        if ((k + 1) < (k2 - k1)){
            s += sqrt(pow(solution[k1 + k + 1][0] - solution[k1 + k][0], 2)
                    + pow(solution[k1 + k + 1][1] - solution[k1 + k][1], 2)
                    + pow(solution[k1 + k + 1][2] - solution[k1 + k][2], 2));
        }
    }
    beta_caustic.f_vals[0] = beta_caustic.f_vals[1];
    s_tot = s;

    interp::set(ray_x);
    interp::set(ray_y);
    interp::set(ray_z);
    interp::set(D);
    interp::set(beta);
    interp::set(beta_caustic);

    // Identify waveform starting time and sample rate
    dt = u[1][0] - u[0][0];
    t0 = u[0][0];

    // Define the Tukey time series window to keep energy at the beginning
    // and end of the record from causing problems
    tukey_param = 0.8;
    for(int n = 0; n < wvfrm::len; n++){
        if(n < tukey_param * (wvfrm::len - 1) / 2.0){
            tukey_win[n] = 1.0 / 2.0 * (1.0 + cos(Pi * (2.0 * n / (tukey_param * (wvfrm::len - 1)) - 1.0)));
        } else if (n > (wvfrm::len - 1) * (1.0 - tukey_param / 2.0)){
            tukey_win[n] = 1.0 / 2.0 * (1.0 + cos(Pi * (2.0 * n / (tukey_param * (wvfrm::len - 1)) - 2.0 / tukey_param + 1.0)));
        } else {
            tukey_win[n] = 1.0;
        }
    }


    // Define low-pass filters to prevent Gibbs (ringing) phenomena
    // Note: filter1 is applied to the non-linear contribution and damps out high
    // frequency there, while filter2 is applied to the entire signal to
    // prevent ringing if energy gets through filter1 but avoids filtering out
    // the original signal
    for(int n = 0; n < fft_len; n++){
        filter1[n] = pow(1.0 + pow(n / (wvfrm::filt1_g * fft_len), wvfrm::filt1_n1), -wvfrm::filt1_n2);
        filter2[n] = pow(1.0 + pow(n / (wvfrm::filt2_g * fft_len), wvfrm::filt2_n1), -wvfrm::filt2_n2);
    }

    // Compute u(s0, t) = p(s0, t) / p0 and U(s0, f) = FFT{u(s0, t)} to initialize
    for(int n = 0; n < wvfrm::len; n++){
        fft_time[n] = u[n][1];
    }
    fftw_execute(fwd_plan);
    for(int n = 0; n < fft_len; n++){
        U[n][0] = fft_spec[n][0];
        U[n][1] = fft_spec[n][1];
    }

    Uf_max = 0.0;
    for(int n = 0; n < fft_len; n++){
        Uf_max = max(Uf_max, (1.0 / dt * n / wvfrm::len) * sqrt(pow(U[n][0], 2) + pow(U[n][1], 2)));
    }

    // Use Heun's Method (similar to RK2) to compute the waveform along the ray path
    s = 0.0;
    ds = ds_wvfrm;
    while(s < s_tot){
        // Note: it is assumed that u(s, t) and U(s, f) are known to start step
        if ((output_step < 1.0e4) && (fabs(output_step * ceil((s + s_prev + ds) / output_step) - output_step * ceil((s + s_prev) / output_step)) > 0.0)){
            sprintf(output_buffer, "sc_spectrum-%04g.dat", output_step * ceil((s + s_prev) / output_step));
            file_out.open(output_buffer);
            for (int n = 0; n < fft_len; n++){
                file_out << 1.0 / dt * n / wvfrm::len << '\t' << U[n][0] << '\t' << U[n][1] << '\n';
            }
            file_out.close();

            sprintf(output_buffer, "sc_waveform-%04g.dat", output_step * ceil((s + s_prev) / output_step));
            file_out.open(output_buffer);
            for (int n = 0; n < wvfrm::len; n++){
                file_out << u[n][0] << '\t' << u[n][1] << '\n';
            }
            file_out.close();
        }

        // Compute beta, the step size, and beta_next
        // Note: near caustics, use transformed Burger's equation with beta_caustic to avoid the singularity
        if(fabs(interp::eval_f(s, D)) >= 0.5){
            beta_curr = interp::eval_f(s, beta);
        } else {
            beta_curr = interp::eval_f(s, beta_caustic);
        }

        // Use current value of beta to define step size for stability
        ds = ds_wvfrm / (Pi * Uf_max * beta_curr);
        ds = min(ds, 1.0);
        ds = max(ds, 1.0e-5);

        if(fabs(interp::eval_f(min(s + ds, ray_z.x_vals[(k2 - k1) - 1]), D)) >= 0.5){
            beta_next = interp::eval_f(min(s + ds, ray_z.x_vals[(k2 - k1) - 1]), beta);
        } else {
            beta_next = interp::eval_f(min(s + ds, ray_z.x_vals[(k2 - k1) - 1]), beta_caustic);
        }

        if (fabs(10.0 * ceil((s + s_prev + ds) / 10.0) - 10.0 * ceil((s + s_prev) / 10.0)) > 0.0){
            cout << '\t' << "Computing waveform at ray length: " << 10.0 * ceil((s + s_prev) / 10.0) << '\t' << "altitude: " << interp::eval_f(s, ray_z) << '\t' << "scaled non-linearity factor: " << beta_curr << '\n';
        }


        // Compute A = FFT{u^2}
        for(int n = 0; n < wvfrm::len; n++){
            u_sqr[n] = pow(u[n][1], 2);
            fft_time[n] = u_sqr[n];
        }
        fftw_execute(fwd_plan);
        for(int n = 0; n < fft_len; n++){
            A[n][0] = fft_spec[n][0];
            A[n][1] = fft_spec[n][1];
        }

        // Compute B = FFT{u d/dt(u^2)} and C = FFT{(d/dt(u^2))^2}
        du_sqr[0] = (u_sqr[1] - u_sqr[0]) / (u[1][0] - u[0][0]);
        for(int n = 1; n < wvfrm::len - 1; n++){
            du_sqr[n] = (u_sqr[n + 1] - u_sqr[n - 1]) / (u[n + 1][0] - u[n - 1][0]);
        }
        du_sqr[wvfrm::len - 1] = (u_sqr[wvfrm::len - 1] - u_sqr[wvfrm::len - 2]) / (u[wvfrm::len - 1][0] - u[wvfrm::len - 2][0]);

        for(int n = 0; n < wvfrm::len; n++){
            fft_time[n] = u[n][1] * du_sqr[n];
        }
        fftw_execute(fwd_plan);
        for(int n = 0; n < fft_len; n++){
            B[n][0] = fft_spec[n][0];
            B[n][1] = fft_spec[n][1];
        }

        for(int n = 0; n < wvfrm::len; n++){
            fft_time[n] = pow(du_sqr[n], 2);
        }
        fftw_execute(fwd_plan);
        for(int n = 0; n < fft_len; n++){
            C[n][0] = fft_spec[n][0];
            C[n][1] = fft_spec[n][1];
        }

        // Step through frequencies and update the spectra through the step
        // U(s + ds, f) = U(s, f) + i pi f ds (...) (see notes)

        // Check for caustic between steps and apply phase shift if necessary (Hilbert transform)
        if(s + ds <= s_tot) {
            if(interp::eval_f(s, D) * interp::eval_f(s + ds, D) <= 0.0){
                cout << '\t' << "Caustic encountered." << '\n';
                for(int n = 0; n < fft_len; n++){
                    double temp = U[n][0];
                    U[n][0] = - U[n][1];
                    U[n][1] = temp;
                }
            }
        }

        // Apply absorption and apply Heun's Method step using symmetry of rft
        for(int n = 0; n < fft_len; n++){
            // Absorpotion computed in dB as: 10 * log10(p(s + ds)) - 10 * log10(p(s)) = alpha * ds [dB] -> p(s + ds) = p(s) * 10^(alpha(s) * ds / 10) with c/c0 factor from Lonzaga paper
            double freq = 1.0 / dt * n / wvfrm::len;
            double absorp = pow(10.0, - atmo::c(interp::eval_f(s, ray_x), interp::eval_f(s, ray_y), interp::eval_f(s, ray_z)) / c0 * atmo::SB_alpha(interp::eval_f(s, ray_x), interp::eval_f(s, ray_y), interp::eval_f(s, ray_z), freq) * ds / 10.0);

            U[n][0] -= Pi * freq * (ds * (beta_curr + beta_next) / 2.0 * A[n][1] + pow(ds, 2) / 2.0 * beta_curr * beta_next * B[n][1] + pow(ds, 3) / 8.0 * pow(beta_curr, 2) * beta_next * C[n][1]) * filter1[n];
            U[n][1] += Pi * freq * (ds * (beta_curr + beta_next) / 2.0 * A[n][0] + pow(ds, 2) / 2.0 * beta_curr * beta_next * B[n][0] + pow(ds, 3) / 8.0 * pow(beta_curr, 2) * beta_next * C[n][0]) * filter1[n];

            U[n][0] *= absorp * filter2[n];
            U[n][1] *= absorp * filter2[n];
        }

        // Reverse the FFT to obtain the waveform at s + ds, check the new maximum scaled waveform, and increment s -> d + ds
        // Use Tukey window to to maintain
        for(int n = 0; n < fft_len; n++){
            fft_spec[n][0] = U[n][0];
            fft_spec[n][1] = U[n][1];
        }
        fftw_execute(bwd_plan);

        for(int n = 0; n < wvfrm::len; n++){
            u[n][1] = fft_time[n]  * tukey_win[n] / wvfrm::len;
        }

        Uf_max = 0.0;
        for(int n = 0; n < fft_len; n++){
            Uf_max = max(Uf_max, (1.0 / dt * n / wvfrm::len) * sqrt(pow(U[n][0], 2) + pow(U[n][1], 2)));
        }

        s += ds;
    }

    // Clear all memory used in analysis
    interp::clear(ray_x);
    interp::clear(ray_y);
    interp::clear(ray_z);

    interp::clear(D);
    interp::clear(beta);
    interp::clear(beta_caustic);

    delete [] u_sqr;
    delete [] du_sqr;
    delete [] tukey_win;
    
    delete [] filter1;
    delete [] filter2;

    for (int n = 0; n < fft_len; n++){
        delete U[n];    delete A[n];
        delete B[n];    delete C[n];
    }

    delete [] U;        delete [] A;
    delete [] B;        delete [] C;
    
    fftw_destroy_plan(bwd_plan);
    fftw_destroy_plan(fwd_plan);

    fftw_free(fft_spec);
    delete [] fft_time;

    return s_tot;
}



#endif /* _GEOAC_EQSET_3D_CPP_ */
