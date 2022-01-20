#ifndef _GEOAC_EQSET_SPH_RNGDEP_CPP_
#define _GEOAC_EQSET_SPH_RNGDEP_CPP_

#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>

#include <fftw3.h>

#include "geoac.params.h"
#include "geoac.eqset.h"

#include "../atmo/atmo_state.h"
#include "../atmo/atmo_io.sph.rngdep.h"

#include "../util/interpolation.h"
#include "../util/globe.h"
#include "../util/waveforms.h"


using namespace std;

//--------------------------------------------------//
//-------Propagate in 3D, assume non-stratified-----//
//--------------------------------------------------//
void geoac::set_system(){
    dim = 3;
    is_strat = false;
}

//-----------------------------------------------------------//
//-------Structure containing source functions which are-----//
//-------called multiple times in the solver-----------------//
//-----------------------------------------------------------//
namespace geoac {
    struct src_refs{
        double src_loc[3]; double c0;                   // Location of the source (r, theta, phi) and sound speed at the location
        
        double c;   double dc[5];   double ddc[3][2];   // Thermodynamic sound speed and its derivatives
        double w;	double dw[5];   double ddw[3][2];   // Vertical wind speed and its derivatives
        double v;   double dv[5];   double ddv[3][2];   // N-S wind speed (theta component, positive values toward south pole)
        double u;   double du[5];   double ddu[3][2];   // E-W wind speed (phi component, positive values to the west)
        
        double nu0; double nu_mag;  double dnu_mag[2];  // Eikonal vector magnitude at the source, along the ray, and derivatives
        double c_gr[3]; double c_gr_mag;                // Group velocity along ray paths and its magnitude
        double dc_gr[3][2]; double dc_gr_mag[2];        // Group velocity derivatives
        
        double G[3]; double dG[3][2];                   // Scalar coefficients for coordiante transformation (spherical -> 1.0, 1.0/r, 1.0/(r*sin(theta)))
        double T[3]; double dT[3][2];                   // Extra terms which keep the eikonal vector direction constant as the unit vectors vary with position
    };
    
    struct src_refs refs = {
        {0.0, 0.0, 0.0}, 0.0,
        0.0, {0.0, 0.0, 0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
        0.0, {0.0, 0.0, 0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
        0.0, {0.0, 0.0, 0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
        0.0, {0.0, 0.0, 0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
        0.0, 0.0, {0.0, 0.0},
        {0.0, 0.0, 0.0}, 0.0, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}, {0.0, 0.0},
        {0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}},
        {0.0, 0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}}
    };
}
//----------------------------------------------------------------------//
//-------Fill in solution[0][n] with the appropriate initial values-----//
//----------------------------------------------------------------------//
void geoac::set_initial(double ** & solution, double z_src, double lat_src, double lon_src){
    refs.src_loc[0] = globe::r0 + z_src;
    refs.src_loc[1] = lat_src;
    refs.src_loc[2] = lon_src;
	refs.c0 = atmo::c(refs.src_loc[0], refs.src_loc[1], refs.src_loc[2]);
    
    double mach_comps[3] = {atmo::w(refs.src_loc[0], refs.src_loc[1], refs.src_loc[2]) / refs.c0,
                            atmo::v(refs.src_loc[0], refs.src_loc[1], refs.src_loc[2]) / refs.c0,
                            atmo::u(refs.src_loc[0], refs.src_loc[1], refs.src_loc[2]) / refs.c0};
    
    double nu0[3] =    {sin(theta),  cos(theta) * sin(phi),  cos(theta) * cos(phi)};
    double mu0_lt[3] = {cos(theta), -sin(theta) * sin(phi), -sin(theta) * cos(phi)};
    double mu0_lp[3] = {0.0,         cos(theta) * cos(phi), -cos(theta) * sin(phi)};
    double mach_sc = 1.0 + (nu0[0] * mach_comps[0] + nu0[1] * mach_comps[1] + nu0[2] * mach_comps[2]);
    refs.nu0 = 1.0 / mach_sc;
    
	for(int n_eq = 0; n_eq < eq_cnt; n_eq++){
		switch(n_eq){
			case(0):    solution[0][n_eq] =  refs.src_loc[0];   break;  // r(0) = r0 (relative to earth radius)
            case(1):    solution[0][n_eq] =  refs.src_loc[1];   break;  // t(0) = theta0
            case(2):    solution[0][n_eq] =  refs.src_loc[2];   break;  // p(0) = phi0
                
			case(3):                                                        // nu_r(0)
			case(4):                                                        // nu_t(0)
			case(5):    solution[0][n_eq] = nu0[n_eq - 3] / mach_sc;  break;  // nu_p(0)
				
            case(6):                                        // Rt(0) = 0
			case(7):                                        // Tht(0) = 0
			case(8):                                        // Pht(0) = 0
			case(12):                                       // Rp(0) = 0
			case(13):                                       // Thp(0) = 0
			case(14):   solution[0][n_eq] =  0.0;   break;  // Php(0) = 0
                
            case(9):                                                                                                                                                                                        // mu_r_lt(0) = d/d lt (nu_r)
            case(10):                                                                                                                                                                                       // mu_t_lt(0) = d/d lt (nu_t)
            case(11):   solution[0][n_eq] =  mu0_lt[n_eq - 9] / mach_sc - nu0[n_eq - 9] / pow(mach_sc, 2) * (mu0_lt[0] * mach_comps[0] + mu0_lt[1] * mach_comps[1] + mu0_lt[2] * mach_comps[2]);    break;  // mu_p_lt(0) = d/d lt (nu_p)
                
            case(15):                                                                                                                                                                                       // mu_r_lp(0) = d/d lp (nu_r)
            case(16):                                                                                                                                                                                       // mu_t_lp(0) = d/d lp (nu_t)
            case(17):   solution[0][n_eq] =  mu0_lp[n_eq - 15] / mach_sc - nu0[n_eq - 15] / pow(mach_sc, 2) * (mu0_lp[0] * mach_comps[0] + mu0_lp[1] * mach_comps[1] + mu0_lp[2] * mach_comps[2]);  break;  // mu_p_lp(0) = d/d lp (nu_p)
                
            default:    cout << "Unexpected index in Initial_Cond.  Model includes 18 variables." << '\n'; break;
		}
	}
}

//--------------------------------------------------------------------------//
//-------Taylor series fit to more accurately deterine intercept values-----//
//--------------------------------------------------------------------------//
void geoac::approx_intercept(double ** solution, int k, double* & prev){
	double dr_grnd = topo::z(solution[k - 1][1], solution[k - 1][2]) - solution[k - 1][0];                  // Set dr from r_{k-1} to ground
	double dr_k = (solution[k - 1][0] - solution[k - 2][0])                                                 // Set effective dr with ground slope
                    - (solution[k - 1][1] - solution[k - 2][1]) * 1.0 / solution[k - 1][0] * topo::dz(solution[k - 1][1], solution[k-1][2], 0)
                    - (solution[k - 1][2] - solution[k - 2][2]) * 1.0 / (solution[k - 1][0] * cos(solution[k - 1][1])) * topo::dz(solution[k - 1][1], solution[k-1][2], 1);
    
	for(int n_eq = 0; n_eq < eq_cnt; n_eq++){ prev[n_eq] = solution[k - 1][n_eq] + (solution[k - 1][n_eq] - solution[k - 2][n_eq]) / dr_k * dr_grnd;}
}

//-------------------------------------------------------------------------//
//-------Fill in solution[0][n] with the appropriate reflection values-----//
//-------------------------------------------------------------------------//
void geoac::set_refl(double** & solution, int k_end){
    double* prev = new double [eq_cnt];
    approx_intercept(solution, k_end, prev);
    
    double rg, c_grnd, a, b, norm1, norm2, C1r, C1th, C1ph, C2th, C2ph;
    
    double da_dlth, db_dlth, da_dlph, db_dlph;
    double dC1r_dlth, dC1th_dlth, dC1ph_dlth, dC2th_dlth, dC2ph_dlth, ds0_dlth;
    double dC1r_dlph, dC1th_dlph, dC1ph_dlph, dC2th_dlph, dC2ph_dlph, ds0_dlph;
    double dnu_r_ds, dnu_th_ds, dnu_ph_ds;
    
    rg = topo::z(prev[1], prev[2]);
    c_grnd = atmo::c(rg, prev[1], prev[2]);
    
    a = 1.0 / rg * topo::dz(prev[1], prev[2], 0);
    b = 1.0 / (rg * cos(prev[1])) * topo::dz(prev[1], prev[2], 1);
    
    norm1 = 1.0 + pow(a, 2) + pow(b, 2);
    
    C1r =  (1.0 - pow(a, 2) - pow(b, 2)) / norm1;
    C1th = (1.0 - pow(a, 2) + pow(b, 2)) / norm1;   C2th = 2.0 * a / norm1;
    C1ph = (1.0 + pow(a, 2) - pow(b, 2)) / norm1;   C2ph = 2.0 * b / norm1;
    
    if(calc_amp){
        norm2 = prev[3] - a * prev[4] - b * prev[5];
        
        da_dlth = - prev[6] / rg * a + 1.0 / rg * ((prev[7] * (prev[3] - b * prev[5]) - prev[4] * (prev[6] / rg  - cos(prev[1]) * b * prev[8])) / norm2 * topo::ddz(prev[1], prev[2], 0, 0)
                                                   + (prev[8] * (prev[3] - a * prev[4]) - prev[5] * (prev[6] / (rg * cos(prev[1])) - a / cos(prev[1]) * prev[7])) / norm2 * topo::ddz(prev[1], prev[2], 0, 1));
        
        da_dlph = - prev[12] / rg * a + 1.0 / rg * ((prev[13] * (prev[3] - b * prev[5]) - prev[4] * (prev[12] / rg - cos(prev[1]) * b * prev[14])) / norm2 * topo::ddz(prev[1], prev[2], 0, 0)
                                                    + (prev[14] * (prev[3] - a * prev[4]) - prev[5] * (prev[12] / (rg * cos(prev[1])) - a / cos(prev[1]) * prev[13])) / norm2 * topo::ddz(prev[1], prev[2], 0, 1));
        
        db_dlth = - (prev[6] / rg - tan(prev[4]) * prev[7]) * b + 1.0 / (rg * cos(prev[1])) * ((prev[7] * (prev[3] - b * prev[5]) - prev[4] * (prev[6] / rg - cos(prev[1]) * b * prev[8])) / norm2 * topo::ddz(prev[1], prev[2], 1, 0)
                                                                                               + (prev[8] * (prev[3] - a * prev[4]) - prev[5] * (prev[6] / (rg * cos(prev[1])) - a / cos(prev[1]) * prev[7])) / norm2 * topo::ddz(prev[1], prev[2], 1, 1));
        
        db_dlph = - (prev[12] / rg - tan(prev[4]) * prev[13]) * b + 1.0 / (rg * cos(prev[1])) * ((prev[13] * (prev[3] - b * prev[5]) - prev[4] * (prev[12] / rg - cos(prev[1]) * b * prev[14])) / norm2 * topo::ddz(prev[1], prev[2], 1, 0)
                                                                                                 + (prev[14] * (prev[3] - a * prev[4]) - prev[5] * (prev[12] / (rg * cos(prev[1])) - a / cos(prev[1]) * prev[13])) / norm2 * topo::ddz(prev[1], prev[2], 1, 1));
        
        dC1r_dlth =  -4.0 / pow(norm1, 2) * (a * da_dlth + b * db_dlth);                                        dC1r_dlph =  -4.0 / pow(norm1, 2) * (a * da_dlph + b * db_dlph);
        dC1th_dlth = -4.0 * a / pow(norm1, 2) * ((1.0 + pow(b, 2)) * da_dlth - a * b * db_dlth);                dC1th_dlph = -4.0 * a / pow(norm1, 2) * ((1.0 + pow(b, 2)) * da_dlph - a * b * db_dlph);
        dC1ph_dlth = -4.0 * b / pow(norm1, 2) * ((1.0 + pow(a, 2)) * db_dlth - a * b * da_dlth);                dC1ph_dlph = -4.0 * b / pow(norm1, 2) * ((1.0 + pow(a, 2)) * db_dlph - a * b * da_dlph);
        dC2th_dlth =  2.0 / pow(norm1, 2) * ((1.0 - pow(a, 2) + pow(b, 2)) * da_dlth - 2.0 * a * b * db_dlth);  dC2th_dlph =  2.0 / pow(norm1, 2) * ((1.0 - pow(a, 2) + pow(b, 2)) * da_dlph - 2.0 * a * b * db_dlph);
        dC2ph_dlth =  2.0 / pow(norm1, 2) * ((1.0 + pow(a, 2) - pow(b, 2)) * db_dlth - 2.0 * a * b * db_dlth);  dC2ph_dlph =  2.0 / pow(norm1, 2) * ((1.0 + pow(a, 2) - pow(b, 2)) * db_dlph - 2.0 * a * b * da_dlph);
        
        ds0_dlth = -refs.c0 / c_grnd * (prev[6] - rg * a * prev[7] - rg * cos(prev[1]) * b * prev[8]) / norm2;
        ds0_dlph = -refs.c0 / c_grnd * (prev[12] - rg * a * prev[13] - rg * cos(prev[1]) * b * prev[14]) / norm2;
        
        dnu_r_ds  = refs.c0 / c_grnd * atmo::dc(rg, prev[1], prev[2], 0) + prev[3] * atmo::dw(rg, prev[1], prev[2], 0) + prev[4] * atmo::dv(rg, prev[1], prev[2], 0) + prev[5] * atmo::du(rg, prev[1], prev[2], 0);
        dnu_th_ds = refs.c0 / c_grnd * atmo::dc(rg, prev[1], prev[2], 1) + prev[3] * atmo::dw(rg, prev[1], prev[2], 1) + prev[4] * atmo::dv(rg, prev[1], prev[2], 1) + prev[5] * atmo::du(rg, prev[1], prev[2], 1);
        dnu_ph_ds = refs.c0 / c_grnd * atmo::dc(rg, prev[1], prev[2], 2) + prev[3] * atmo::dw(rg, prev[1], prev[2], 2) + prev[4] * atmo::dv(rg, prev[1], prev[2], 2) + prev[5] * atmo::du(rg, prev[1], prev[2], 2);
        
        dnu_r_ds  += pow(c_grnd, 2) / refs.c0 * 1.0 / rg * (pow(prev[4], 2) + pow(prev[5], 2));
        dnu_th_ds -= pow(c_grnd, 2) / refs.c0 * (prev[3] * prev[4] + tan(prev[1]) * pow(prev[5], 2));
        dnu_ph_ds -= pow(c_grnd, 2) / refs.c0 * prev[5] * (prev[3] * cos(prev[1]) + prev[4] * sin(prev[1]));
        
        dnu_r_ds *= - 1.0 / c_grnd;
        dnu_th_ds *= - 1.0 / (c_grnd * rg);
        dnu_ph_ds *= - 1.0 / (c_grnd * rg * cos(prev[1]));
    }
    
    for(int n_eq = 0; n_eq < eq_cnt; n_eq++){
        switch(n_eq){
            case(0):    solution[0][n_eq] = rg;             break;  // r = r_grnd
                
            case(1):
            case(2):    solution[0][n_eq] =  prev[n_eq];    break;  // theta, phi are continuous
                
            case(3):    solution[0][n_eq] =  -C1r * prev[3] + C2th * prev[4] + C2ph * prev[5];   break;  // nu_r
            case(4):    solution[0][n_eq] =  C1th * prev[4] + C2th * (prev[3] - b * prev[5]);    break;  // nu_th
            case(5):	solution[0][n_eq] =  C1ph * prev[5] + C2ph * (prev[3] - a * prev[4]);    break;  // nu_ph
                
            case(6):                                                                                                                                // R_lth
            case(12):   solution[0][n_eq] =  -C1r * prev[n_eq] + rg * C2th * prev[n_eq + 1] + rg * cos(prev[1]) * C2ph * prev[n_eq + 2];    break;  // R_lph
                
            case(7):                                                                                                                        // Theta_lth
            case(13):   solution[0][n_eq] =  C1th * prev[n_eq] + C2th * (prev[n_eq - 1] / rg  - b * cos(prev[1]) * prev[n_eq + 1]); break;  // Theta_lph
                
            case(8):                                                                                                                                        // Phi_lth
            case(14):   solution[0][n_eq] =  C1ph * prev[n_eq] + C2ph * (prev[n_eq - 2] / (rg * cos(prev[1])) - a / cos(prev[1]) * prev[n_eq - 1]); break;  // Phi_lph
                
                
            case(9):    solution[0][n_eq] = -C1r * prev[9] + C2th * prev[10] + C2ph * prev[11] - dC1r_dlth * prev[3] + dC2th_dlth * prev[4] + dC2ph_dlth * prev[5] + (-(C1r + 1.0) * dnu_r_ds + C2th * dnu_th_ds + C2ph * dnu_ph_ds) * ds0_dlth;                break;  // mu_r_lth
            case(10):	solution[0][n_eq] =  C1th * prev[10] + C2th * (prev[9] - b * prev[11] - db_dlth * prev[5]) + dC1th_dlth * prev[4] + dC2th_dlth * (prev[3] - b * prev[5]) + ((C1th - 1.0) * dnu_th_ds + C2th * (dnu_r_ds - b * dnu_ph_ds)) * ds0_dlth;   break;  // mu_th_lth
            case(11):   solution[0][n_eq] =  C1ph * prev[11] + C2ph * (prev[9] - a * prev[10] - da_dlth * prev[4]) + dC1ph_dlth * prev[5] + dC2ph_dlth * (prev[3] - a * prev[4]) + ((C1ph - 1.0) * dnu_ph_ds + C2ph * (dnu_r_ds - a * dnu_th_ds)) * ds0_dlth;   break;  // mu_ph_lth
                
            case(15):   solution[0][n_eq] = -C1r * prev[15] + C2th * prev[16] + C2ph * prev[17] - dC1r_dlph * prev[3] + dC2th_dlph * prev[4] + dC2ph_dlph * prev[5] + (-(C1r + 1.0) * dnu_r_ds + C2th * dnu_th_ds + C2ph * dnu_ph_ds) * ds0_dlph;               break;  // mu_r_lph
            case(16):	solution[0][n_eq] =  C1th * prev[16] + C2th * (prev[15] - b * prev[17] - db_dlph * prev[5]) + dC1th_dlph * prev[4] + dC2th_dlph * (prev[3] - b * prev[5]) + ((C1th - 1.0) * dnu_th_ds + C2th * (dnu_r_ds - b * dnu_ph_ds)) * ds0_dlph;  break;  // mu_th_lplh
            case(17):   solution[0][n_eq] =  C1ph * prev[17] + C2ph * (prev[15] - a * prev[16] - da_dlph * prev[4]) + dC1ph_dlph * prev[5] + dC2ph_dlph * (prev[3] - a * prev[4]) + ((C1ph - 1.0) * dnu_ph_ds + C2ph * (dnu_r_ds - a * dnu_th_ds)) * ds0_dlph;  break;  // mu_ph_lph
        }
    }
    
    delete [] prev;
}


//-----------------------------------------------------------------------------------//
//-------Vary the solver step size, currently uses smaller steps near the ground-----//
//-----------------------------------------------------------------------------------//
double geoac::set_ds(double* current_values){
    double result = ds_max;
    if (current_values[0] - globe::r0 < topo::z_bndlyr){
        result -= (ds_max - ds_min) * exp(-(current_values[0] - topo::z(current_values[1], current_values[2])) / 0.2);
    }
	return result;
}

//---------------------------------------//
//-------Update the source functions-----//
//---------------------------------------//
void geoac::update_refs(double ray_length, double* current_values){
    double r = current_values[0],       t = current_values[1],  p = current_values[2];
	double nu[3] = {current_values[3],  current_values[4],      current_values[5]};
    refs.nu_mag = sqrt(pow(nu[0], 2) + pow(nu[1], 2) + pow(nu[2], 2));

    double R_lt[3], R_lp[3], mu_lt[3], mu_lp[3];
    
    double dc_temp[3], dd_c_temp[6], dd_u_temp[6], dd_v_temp[6], dd_w_temp[6];
    double dd_c[3][3], dd_u[3][3], dd_v[3][3], dd_w[3][3];

    if(!calc_amp){
        interp::eval_all(t, p, r - globe::r0, atmo::c_spline, refs.c, dc_temp);
        refs.dc[1] = dc_temp[0];
        refs.dc[2] = dc_temp[1];
        refs.dc[0] = dc_temp[2];

        for(int n = 0; n < 3; n++){
            atmo::u_spline.accel[n] = atmo::c_spline.accel[n];
            atmo::v_spline.accel[n] = atmo::c_spline.accel[n];
        }
        atmo::calc_uvw(r, t, p, refs.u, refs.v, refs.w, refs.du, refs.dv, refs.dw);
    
        refs.c_gr[0] = refs.c * nu[0] / refs.nu_mag + refs.w;
        refs.c_gr[1] = refs.c * nu[1] / refs.nu_mag + refs.v;
        refs.c_gr[2] = refs.c * nu[2] / refs.nu_mag + refs.u;
        refs.c_gr_mag = sqrt(pow(refs.c_gr[0], 2) + pow(refs.c_gr[1], 2) + pow(refs.c_gr[2], 2));
        
        refs.G[0] = 1.0;                    refs.T[0] = 1.0 / r * (nu[1] * refs.c_gr[1] + nu[2] * refs.c_gr[2]);
        refs.G[1] = 1.0 / r;                refs.T[1] = (nu[0] * refs.v - nu[1] * refs.w) - nu[0] * refs.c_gr[1] + nu[2] * refs.c_gr[2] * tan(t);
        refs.G[2] = 1.0 / (r * cos(t));     refs.T[2] = (nu[0] * refs.u - nu[2] * refs.w) * cos(t) + (nu[1] * refs.u - nu[2] * refs.v) * sin(t) - refs.c_gr[2] * (nu[0] * cos(t) + nu[1] * sin(t));
        
    } else {
        interp::eval_all(t, p, r - globe::r0, atmo::c_spline, refs.c, dc_temp, dd_c_temp);
        refs.dc[1] = dc_temp[0];        dd_c[1][1] = dd_c_temp[0];
        refs.dc[2] = dc_temp[1];        dd_c[2][2] = dd_c_temp[1];
        refs.dc[0] = dc_temp[2];        dd_c[0][0] = dd_c_temp[2];

        dd_c[1][2] = dd_c_temp[3];      dd_c[2][1] = dd_c[1][2];
        dd_c[0][1] = dd_c_temp[4];      dd_c[1][0] = dd_c[0][1];
        dd_c[0][2] = dd_c_temp[5];      dd_c[2][0] = dd_c[0][2];

        for(int n = 0; n < 3; n++){
            atmo::u_spline.accel[n] = atmo::c_spline.accel[n];
            atmo::v_spline.accel[n] = atmo::c_spline.accel[n];
        }
        atmo::calc_uvw(r, t, p, refs.u, refs.v, refs.w, refs.du, refs.dv, refs.dw, dd_u_temp, dd_v_temp, dd_w_temp);

        refs.c_gr[0] =  refs.c * nu[0] / refs.nu_mag + refs.w;
        refs.c_gr[1] =  refs.c * nu[1] / refs.nu_mag + refs.v;
        refs.c_gr[2] =  refs.c * nu[2] / refs.nu_mag + refs.u;
        refs.c_gr_mag = sqrt(pow(refs.c_gr[0], 2) + pow(refs.c_gr[1], 2) + pow(refs.c_gr[2], 2));
        
        refs.G[0] = 1.0;                    refs.T[0] = 1.0 / r * (nu[1] * refs.c_gr[1] + nu[2] * refs.c_gr[2]);
        refs.G[1] = 1.0 / r;                refs.T[1] = (nu[0] * refs.v - nu[1] * refs.w) - nu[0] * refs.c_gr[1] + nu[2] * refs.c_gr[2] * tan(t);
        refs.G[2] = 1.0 / (r * cos(t));     refs.T[2] = (nu[0] * refs.u - nu[2] * refs.w) * cos(t) + (nu[1] * refs.u - nu[2] * refs.v) * sin(t) - refs.c_gr[2] * (nu[0] * cos(t) + nu[1] * sin(t));

        for(int n = 0; n < 3; n++){
            dd_u[n][n] = dd_u_temp[n];
            dd_v[n][n] = dd_v_temp[n];
            dd_w[n][n] = dd_w_temp[n];
        }
        dd_w[0][1] = dd_w_temp[3];  dd_w[1][0] = dd_w[0][1];
        dd_w[0][2] = dd_w_temp[4];  dd_w[2][0] = dd_w[0][2];
        dd_w[1][2] = dd_w_temp[5];  dd_w[2][1] = dd_w[1][2];
        
        dd_v[0][1] = dd_v_temp[3];  dd_v[1][0] = dd_v[0][1];
        dd_v[0][2] = dd_v_temp[4];  dd_v[2][0] = dd_v[0][2];
        dd_v[1][2] = dd_v_temp[5];  dd_v[2][1] = dd_v[1][2];
        
        dd_u[0][1] = dd_u_temp[3];  dd_u[1][0] = dd_u[0][1];
        dd_u[0][2] = dd_u_temp[4];  dd_u[2][0] = dd_u[0][2];
        dd_u[1][2] = dd_u_temp[5];  dd_u[2][1] = dd_u[1][2];
        
        R_lt[0]  = current_values[6];   R_lt[1]  = current_values[7];   R_lt[2]  = current_values[8];
		mu_lt[0] = current_values[9];   mu_lt[1] = current_values[10];  mu_lt[2] = current_values[11];
        R_lp[0]  = current_values[12];  R_lp[1]  = current_values[13];  R_lp[2]  = current_values[14];
		mu_lp[0] = current_values[15];  mu_lp[1] = current_values[16];  mu_lp[2] = current_values[17];
    
        refs.dc[3] = 0.0;   refs.dw[3] = 0.0;   refs.dv[3] = 0.0;  refs.du[3] = 0.0;
        refs.dc[4] = 0.0;   refs.dw[4] = 0.0;   refs.dv[4] = 0.0;  refs.du[4] = 0.0;
        for(int n = 0; n < 3; n++){
            refs.ddc[n][0] = 0.0;   refs.ddw[n][0] = 0.0;   refs.ddv[n][0] = 0.0;   refs.ddu[n][0] = 0.0;
            refs.ddc[n][1] = 0.0;   refs.ddw[n][1] = 0.0;   refs.ddv[n][1] = 0.0;   refs.ddu[n][1] = 0.0;
        }
        
        for(int n = 0; n < 3; n++){
            refs.dc[3] += R_lt[n] * refs.dc[n];   refs.dc[4] += R_lp[n] * refs.dc[n];
            refs.dw[3] += R_lt[n] * refs.dw[n];   refs.dw[4] += R_lp[n] * refs.dw[n];
            refs.dv[3] += R_lt[n] * refs.dv[n];   refs.dv[4] += R_lp[n] * refs.dv[n];
            refs.du[3] += R_lt[n] * refs.du[n];   refs.du[4] += R_lp[n] * refs.du[n];
            
            for(int m = 0; m < 3; m++){
                refs.ddc[m][0] += R_lt[n] * dd_c[m][n];     refs.ddc[m][1] += R_lp[n] * dd_c[m][n];
                refs.ddw[m][0] += R_lt[n] * dd_w[m][n];     refs.ddw[m][1] += R_lp[n] * dd_w[m][n];
                refs.ddv[m][0] += R_lt[n] * dd_v[m][n];     refs.ddv[m][1] += R_lp[n] * dd_v[m][n];
                refs.ddu[m][0] += R_lt[n] * dd_u[m][n];     refs.ddu[m][1] += R_lp[n] * dd_u[m][n];
            }
        }
        
		refs.dnu_mag[0] = (nu[0] * mu_lt[0] + nu[1] * mu_lt[1] + nu[2] * mu_lt[2]) / refs.nu_mag;
		refs.dnu_mag[1] = (nu[0] * mu_lp[0] + nu[1] * mu_lp[1] + nu[2] * mu_lp[2]) / refs.nu_mag;
        
        refs.dc_gr[0][0] = nu[0] / refs.nu_mag * refs.dc[3] + refs.c * mu_lt[0] / refs.nu_mag - refs.c * nu[0] / pow(refs.nu_mag, 2) * refs.dnu_mag[0] + refs.dw[3];
        refs.dc_gr[1][0] = nu[1] / refs.nu_mag * refs.dc[3] + refs.c * mu_lt[1] / refs.nu_mag - refs.c * nu[1] / pow(refs.nu_mag, 2) * refs.dnu_mag[0] + refs.dv[3];
        refs.dc_gr[2][0] = nu[2] / refs.nu_mag * refs.dc[3] + refs.c * mu_lt[2] / refs.nu_mag - refs.c * nu[2] / pow(refs.nu_mag, 2) * refs.dnu_mag[0] + refs.du[3];
        
        refs.dc_gr[0][1] = nu[0] / refs.nu_mag * refs.dc[4] + refs.c * mu_lp[0] / refs.nu_mag - refs.c * nu[0] / pow(refs.nu_mag, 2) * refs.dnu_mag[1] + refs.dw[4];
        refs.dc_gr[1][1] = nu[1] / refs.nu_mag * refs.dc[4] + refs.c * mu_lp[1] / refs.nu_mag - refs.c * nu[1] / pow(refs.nu_mag, 2) * refs.dnu_mag[1] + refs.dv[4];
        refs.dc_gr[2][1] = nu[2] / refs.nu_mag * refs.dc[4] + refs.c * mu_lp[2] / refs.nu_mag - refs.c * nu[2] / pow(refs.nu_mag, 2) * refs.dnu_mag[1] + refs.du[4];
        
        refs.dc_gr_mag[0] = (refs.c_gr[0] * refs.dc_gr[0][0] + refs.c_gr[1] * refs.dc_gr[1][0] + refs.c_gr[2] * refs.dc_gr[2][0]) / refs.c_gr_mag;
        refs.dc_gr_mag[1] = (refs.c_gr[0] * refs.dc_gr[0][1] + refs.c_gr[1] * refs.dc_gr[1][1] + refs.c_gr[2] * refs.dc_gr[2][1]) / refs.c_gr_mag;
        
        refs.dG[0][0] = 0.0;                                                                        refs.dG[0][1] = 0.0;
        refs.dG[1][0] = -R_lt[0] / (pow(r, 2));                                                     refs.dG[1][1] = -R_lp[0] / (pow(r,2));
        refs.dG[2][0] = -R_lt[0] / (pow(r, 2) * cos(t)) + sin(t) / (r * pow(cos(t),2)) * R_lt[1];   refs.dG[2][1] = -R_lp[0] / (pow(r,2) * cos(t)) + sin(t) / (r * pow(cos(t),2)) * R_lp[1];
        
        refs.dT[0][0] = 0.0;
        refs.dT[1][0] = (mu_lt[0] * refs.v + nu[0] * refs.dv[3] - mu_lt[1] * refs.w - nu[1] * refs.dw[3]);
        refs.dT[2][0] = (mu_lt[0] * refs.u + nu[0] * refs.du[3] - mu_lt[2] * refs.w - nu[2] * refs.dw[3]) * cos(t) - (nu[0] * refs.u - nu[2] * refs.w) * R_lt[1] * sin(t)
                            + (mu_lt[1] * refs.u + nu[1] * refs.du[3] - mu_lt[2] * refs.v - nu[2] * refs.dv[3]) * sin(t) + (nu[1] * refs.u - nu[2] * refs.v) * R_lt[1] * cos(t);
        
        refs.dT[0][0] += -R_lt[0] / pow(r, 2) * (nu[1] * refs.c_gr[1] + nu[2] * refs.c_gr[2]) + 1.0 / r * (mu_lt[1] * refs.c_gr[1] + nu[1] * refs.dc_gr[1][0] + mu_lt[2] * refs.c_gr[2] + nu[2] * refs.dc_gr[2][0]);
        refs.dT[1][0] += -mu_lt[0] * refs.c_gr[1] - nu[0] * refs.dc_gr[1][0] + mu_lt[2] * refs.c_gr[2] * tan(t) + nu[2] * refs.dc_gr[2][0] * tan(t) + nu[2] * refs.c_gr[2] * R_lt[1]/pow(cos(t), 2);
        refs.dT[2][0] += -refs.dc_gr[2][0] * (nu[0] * cos(t) + nu[1] * sin(t)) - refs.c_gr[2] * (mu_lt[0] * cos(t) - nu[0] * R_lt[1] * sin(t) + mu_lt[1] * sin(t) + nu[1] * R_lt[1] * cos(t));
        
        refs.dT[0][1] = 0.0;
        refs.dT[1][1] = (mu_lp[0] * refs.v + nu[0] * refs.dv[4] - mu_lp[1] * refs.w - nu[1] * refs.dw[4]);
        refs.dT[2][1] = (mu_lp[0] * refs.u + nu[0] * refs.du[4] - mu_lp[2] * refs.w - nu[2] * refs.dw[4]) * cos(t) - (nu[0] * refs.u - nu[2] * refs.w) * R_lp[1] * sin(t)
                            + (mu_lp[1] * refs.u + nu[1] * refs.du[4] - mu_lp[2] * refs.v - nu[2] * refs.dv[4]) * sin(t) + (nu[1] * refs.u - nu[2] * refs.v) * R_lp[1] * cos(t);
        
        refs.dT[0][1] += -R_lp[0] / pow(r, 2) * (nu[1] * refs.c_gr[1] + nu[2] * refs.c_gr[2]) + 1.0 / r * (mu_lp[1] * refs.c_gr[1] + nu[1] * refs.dc_gr[1][1] + mu_lp[2] * refs.c_gr[2] + nu[2] * refs.dc_gr[2][1]);
        refs.dT[1][1] += -mu_lp[0] * refs.c_gr[1] - nu[0] * refs.dc_gr[1][1] + mu_lp[2] * refs.c_gr[2] * tan(t) + nu[2] * refs.dc_gr[2][1] * tan(t) + nu[2] * refs.c_gr[2] * R_lp[1]/pow(cos(t), 2);
        refs.dT[2][1] += -refs.dc_gr[2][1] * (nu[0] * cos(t) + nu[1] * sin(t)) - refs.c_gr[2] * (mu_lp[0] * cos(t) - nu[0] * R_lp[1] * sin(t) + mu_lp[1] * sin(t) + nu[1] * R_lp[1] * cos(t));
	}
}

//-----------------------------------------------------------//
//-------Evaluate the Source Equation For Specific Index-----//
//-----------------------------------------------------------//
double geoac::eval_src_eq(double ray_length, double* current_values, int n_eq){
	double result;
    
	// Set variables used in all equations
    double nu[3] =      {current_values[3], 	current_values[4], 		current_values[5]};
    double mu_lth[3] =  {current_values[9], 	current_values[10],		current_values[11]};
    double mu_lph[3] =  {current_values[15], 	current_values[16],		current_values[17]};
    
	switch(n_eq){
		case(0):	// d r / d s
        case(1):    // d theta / d s
        case(2):    // d phi / d s
			result = refs.G[n_eq] * refs.c_gr[n_eq] / refs.c_gr_mag;
			break;
            
		case(3): 	// d nu_r /ds
		case(4): 	// d nu_theta / ds
		case(5): 	// d nu_phi / ds
			result = -refs.G[n_eq - 3] / refs.c_gr_mag * (refs.nu_mag * refs.dc[n_eq - 3]
                                                          + nu[0] * refs.dw[n_eq - 3] + nu[1] * refs.dv[n_eq - 3] + nu[2] * refs.du[n_eq - 3] + refs.T[n_eq - 3]);
			break;
            
            
		case(6):	// d R_lt/ds
        case(7):	// d Theta_lt/ds
        case(8): 	// d Phi_lt/ds
			result = refs.dG[n_eq - 6][0] * refs.c_gr[n_eq - 6] / refs.c_gr_mag + refs.G[n_eq - 6] * refs.dc_gr[n_eq - 6][0] / refs.c_gr_mag
                        - refs.G[n_eq - 6] * refs.c_gr[n_eq - 6] / pow(refs.c_gr_mag, 2) * refs.dc_gr_mag[0];
			break;
            
            
		case(9): 	// d mu_r_lt /ds
		case(10): 	// d mu_theta_lt /ds
		case(11): 	// d mu_phi_lt /ds
			result = -refs.dG[n_eq - 9][0] / refs.c_gr_mag * (refs.nu_mag * refs.dc[n_eq - 9] + nu[0] * refs.dw[n_eq - 9] + nu[1] * refs.dv[n_eq - 9] + nu[2] * refs.du[n_eq - 9] + refs.T[n_eq - 9])
                        + refs.G[n_eq - 9] / pow(refs.c_gr_mag, 2) * refs.dc_gr_mag[0] * (refs.nu_mag * refs.dc[n_eq - 9] + nu[0] * refs.dw[n_eq - 9] + nu[1] * refs.dv[n_eq - 9] + nu[2] * refs.du[n_eq - 9])
                            - refs.G[n_eq - 9] / refs.c_gr_mag * (refs.dnu_mag[0] * refs.dc[n_eq - 9] + refs.nu_mag * refs.ddc[n_eq - 9][0] + mu_lth[0] * refs.dw[n_eq - 9] + mu_lth[1] * refs.dv[n_eq - 9] + mu_lth[2] * refs.du[n_eq - 9]
                                                                                                                        + nu[0] * refs.ddw[n_eq - 9][0] + nu[1] * refs.ddv[n_eq - 9][0] + nu[2] * refs.ddu[n_eq - 9][0] + refs.dT[n_eq - 9][0]);
			break;
            
            
		case(12):  	// dR_lp/ds
		case(13):  	// dTheta_lp/ds
        case(14):	// dPhi_lp/ds
			result = refs.dG[n_eq - 12][1] * refs.c_gr[n_eq - 12] / refs.c_gr_mag + refs.G[n_eq - 12] * refs.dc_gr[n_eq - 12][1] / refs.c_gr_mag
                        - refs.G[n_eq - 12] * refs.c_gr[n_eq - 12] / pow(refs.c_gr_mag, 2) * refs.dc_gr_mag[1];
			break;
            
		case(15): 	// d mu_r_lp/ds
		case(16): 	// d mu_t_lp/ds
		case(17): 	// d mu_p_lp/ds
			result = -refs.dG[n_eq - 15][1] / refs.c_gr_mag * (refs.nu_mag * refs.dc[n_eq - 15] + nu[0] * refs.dw[n_eq - 15] + nu[1] * refs.dv[n_eq - 15] + nu[2] * refs.du[n_eq - 15] + refs.T[n_eq - 15])
                        + refs.G[n_eq - 15] / pow(refs.c_gr_mag, 2) * refs.dc_gr_mag[1] * (refs.nu_mag * refs.dc[n_eq - 15] + nu[0] * refs.dw[n_eq - 15] + nu[1] * refs.dv[n_eq - 15] + nu[2] * refs.du[n_eq - 15])
                            - refs.G[n_eq - 15] / refs.c_gr_mag * (refs.dnu_mag[1] * refs.dc[n_eq - 15] + refs.nu_mag * refs.ddc[n_eq - 15][1] + mu_lph[0] * refs.dw[n_eq - 15] + mu_lph[1] * refs.dv[n_eq - 15] + mu_lph[2] * refs.du[n_eq - 15]
                                                                                                                        + nu[0] * refs.ddw[n_eq - 15][1] + nu[1] * refs.ddv[n_eq - 15][1] + nu[2] * refs.ddu[n_eq - 15][1] + refs.dT[n_eq - 15][1]);
            break;
	}
	return result;
}

//-------------------------------------------------------------------//
//-------Calculate the Hamiltonian (Eikonal) To Check For Errors-----//
//-------------------------------------------------------------------//
double geoac::eval_eikonal(double ray_length, double* current_values){
	double r = current_values[0],  t = current_values[1], p = current_values[2];
	double nu[3] = {current_values[3], current_values[4], current_values[5]};
	
	return sqrt(nu[0] * nu[0] + nu[1] * nu[1] + nu[2] * nu[2]) - refs.c0 / atmo::c(r, t, p)
                + (atmo::w(r, t, p) * nu[0] + atmo::v(r, t, p) * nu[1] + atmo::u(r, t, p) * nu[2]) / atmo::c(r, t, p);
}

double geoac::eval_eikonal(double** solution, int k){
	double r = solution[k][0],      t = solution[k][1],     p = solution[k][2];
	double nu[3] = {solution[k][3], solution[k][4],         solution[k][5]};
    
	return sqrt(nu[0] * nu[0] + nu[1] * nu[1] + nu[2] * nu[2]) - refs.c0 / atmo::c(r, t, p)
                + (atmo::w(r, t, p) * nu[0] + atmo::v(r, t, p) * nu[1] + atmo::u(r, t, p) * nu[2]) / atmo::c(r, t, p);
}

double geoac::eval_eikonal_deriv(double** solution, int k){
	double  r = solution[k][0],             t = solution[k][1], p = solution[k][2],
    nu[3] =      {solution[k][3],   solution[k][4],     solution[k][5]},
    R_lth[3] =   {solution[k][6],   solution[k][7],     solution[k][8]},
    mu_lth[3] =  {solution[k][9],   solution[k][10],    solution[k][11]},
    R_lph[3] =   {solution[k][12],  solution[k][13],    solution[k][14]},
    mu_lph[3] =  {solution[k][15],  solution[k][16],    solution[k][17]};
    
    double mag_nu = sqrt(nu[0] * nu[0] + nu[1] * nu[1] + nu[2] * nu[2]);
    
    double dc_dlth = R_lth[0] * atmo::dc(r, t, p, 0) + R_lth[1] * atmo::dc(r, t, p, 1) + R_lth[2] * atmo::dc(r, t, p, 2),   dc_dlph = R_lph[0] * atmo::dc(r, t, p, 0) + R_lph[1] * atmo::dc(r, t, p, 1) + R_lph[2] * atmo::dc(r, t, p, 2);
    double du_dlth = R_lth[0] * atmo::du(r, t, p, 0) + R_lth[1] * atmo::du(r, t, p, 1) + R_lth[2] * atmo::du(r, t, p, 2),   du_dlph = R_lph[0] * atmo::du(r, t, p, 0) + R_lph[1] * atmo::du(r, t, p, 1) + R_lph[2] * atmo::du(r, t, p, 2);
    double dv_dlth = R_lth[0] * atmo::dv(r, t, p, 0) + R_lth[1] * atmo::dv(r, t, p, 1) + R_lth[2] * atmo::dv(r, t, p, 2),   dv_dlph = R_lph[0] * atmo::dv(r, t, p, 0) + R_lph[1] * atmo::dv(r, t, p, 1) + R_lph[2] * atmo::dv(r, t, p, 2);
    double dw_dlth = R_lth[0] * atmo::dw(r, t, p, 0) + R_lth[1] * atmo::dw(r, t, p, 1) + R_lth[2] * atmo::dw(r, t, p, 2),   dw_dlph = R_lph[0] * atmo::dw(r, t, p, 0) + R_lph[1] * atmo::dw(r, t, p, 1) + R_lph[2] * atmo::dw(r, t, p, 2);
    
    double resid_lth = (nu[0] * mu_lth[0] + nu[1] * mu_lth[1] + nu[2] * mu_lth[2]) / mag_nu + mag_nu / atmo::c(r, t, p) * dc_dlth
                            + 1.0 / atmo::c(r, t, p) * (mu_lth[0] * atmo::w(r, t, p) + mu_lth[1] * atmo::v(r, t, p) + mu_lth[2] * atmo::u(r, t, p) + nu[0] * dw_dlth + nu[1] * dv_dlth + nu[2] * du_dlth);
    double resid_lph = (nu[0] * mu_lph[0] + nu[1] * mu_lph[1] + nu[2] * mu_lph[2]) / mag_nu + mag_nu / atmo::c(r, t, p) * dc_dlph
                            + 1.0 / atmo::c(r, t, p) * (mu_lph[0] * atmo::w(r, t, p) + mu_lph[1] * atmo::v(r, t, p) + mu_lph[2] * atmo::u(r, t, p) + nu[0] * dw_dlph + nu[1] * dv_dlph + nu[2] * du_dlph);
    
    return sqrt(pow(resid_lth, 2) + pow(resid_lph, 2));
}

//--------------------------------------------------------------------------//
//-------Check if ray has left propagation region or returned to ground-----//
//--------------------------------------------------------------------------//
bool geoac::break_check(double ** & solution, int k){
    if (solution[k][2] > Pi) {       solution[k][2] -= 2.0 * Pi;}
    else if (solution[k][2] < - Pi){ solution[k][2] += 2.0 * Pi;}
    
    if(solution[k][1] < -Pi / 2.0 || solution[k][1] > Pi / 2.0){ return true;}
    
    if(solution[k][0] - globe::r0 > alt_max){                                          return true;}
    if(globe::gc_dist(solution[k][1], solution[k][2], refs.src_loc[1], refs.src_loc[2]) > rng_max){   return true;}
    if(solution[k][1] < lat_min || solution[k][1] > lat_max){                                       return true;}
    if(solution[k][2] < lon_min || solution[k][2] > lon_max){                                       return true;}
    return false;
}

bool geoac::ground_check(double ** solution, int k){
    if (solution[k][0] - globe::r0 <= topo::z_max){
        if(solution[k][0] < topo::z(solution[k][1], solution[k][2])){
            return true;
        }
    }
    return false;
}

//----------------------------------------------------------------------------------//
//-------Calculate the travel time from source to location or between locations-----//
//----------------------------------------------------------------------------------//
double geoac::travel_time(double ** solution, int k){
    double dr, dt, dp, ds, r, t, p, nu[3], nu_mag, c_prop[3], c_prop_mag;
    double sndspd, traveltime = 0.0;
	
	for (int n = 0; n < k; n++){
		dr = solution[n + 1][0] - solution[n][0];   r = solution[n][0] + dr / 2.0;
		dt = solution[n + 1][1] - solution[n][1];   t = solution[n][1] + dt / 2.0;
		dp = solution[n + 1][2] - solution[n][2];   p = solution[n][2] + dp / 2.0;

        if(dp > 2.0 * Pi){
            dp -= 2.0 * Pi;
        } else if (dp < -2.0 * Pi){
            dp += 2.0 * Pi;
        }

        ds = sqrt(pow(dr, 2) + pow(r * dt, 2) + pow(r * cos(t) * dp, 2));
        
		nu[0] = solution[n][3] + (solution[n + 1][3] - solution[n][3]) / 2.0;
		nu[1] = solution[n][4] + (solution[n + 1][4] - solution[n][4]) / 2.0;
		nu[2] = solution[n][5] + (solution[n + 1][5] - solution[n][5]) / 2.0;
		nu_mag = sqrt(pow(nu[0], 2) + pow(nu[1], 2) + pow(nu[2], 2));
        
        sndspd = atmo::c(r,t,p);
		c_prop[0] = sndspd * nu[0] / nu_mag + atmo::w(r, t, p);
		c_prop[1] = sndspd * nu[1] / nu_mag + atmo::v(r, t, p);
		c_prop[2] = sndspd * nu[2] / nu_mag + atmo::u(r, t, p);
		c_prop_mag = sqrt(pow(c_prop[0], 2) + pow(c_prop[1], 2) + pow(c_prop[2], 2));
        
		traveltime += ds / c_prop_mag;	// Add contribution to the travel time
	}
    
	return traveltime;
}

void geoac::travel_time(double & time, double ** solution, int k1, int k2){
    double dr, dt, dp, ds, r, t, p, nu[3], nu_mag, c_prop[3], c_prop_mag, sndspd;
    
	for (int n = k1; n < k2; n++){
		dr = solution[n + 1][0] - solution[n][0];   r = solution[n][0] + dr / 2.0;
		dt = solution[n + 1][1] - solution[n][1];   t = solution[n][1] + dt / 2.0;
		dp = solution[n + 1][2] - solution[n][2];   p = solution[n][2] + dp / 2.0;

        if(dp > 2.0 * Pi){
            dp -= 2.0 * Pi;
        } else if (dp < -2.0 * Pi){
            dp += 2.0 * Pi;
        }

        ds = sqrt(pow(dr, 2) + pow(r * dt, 2) + pow(r * cos(t) * dp, 2));
        
		nu[0] = solution[n][3] + (solution[n + 1][3] - solution[n][3]) / 2.0;
		nu[1] = solution[n][4] + (solution[n + 1][4] - solution[n][4]) / 2.0;
		nu[2] = solution[n][5] + (solution[n + 1][5] - solution[n][5]) / 2.0;
		nu_mag = sqrt(pow(nu[0], 2) + pow(nu[1], 2) + pow(nu[2], 2));
        
        sndspd = atmo::c(r, t, p);
		c_prop[0] = sndspd * nu[0] / nu_mag + atmo::w(r, t, p);
		c_prop[1] = sndspd * nu[1] / nu_mag + atmo::v(r, t, p);
		c_prop[2] = sndspd * nu[2] / nu_mag + atmo::u(r, t, p);
		c_prop_mag = sqrt(pow(c_prop[0], 2) + pow(c_prop[1], 2) + pow(c_prop[2], 2));
		time += ds / c_prop_mag;	// Add contribution to the travel time
	}
}


//------------------------------------------------//
//-------Calculate the travel time variance-------//
//------------------------------------------------//
void geoac::travel_time_var(double ** solution, int k, double & tt, double & tt_var_incl, double & tt_var_az){
    double dr, dt, dp, ds, r, t, p, nu[3], nu_mag, sndspd, cg[3], cg_mag;
    double R[3], mu[3], sndspd_scalar, dcg[3];

	tt = 0.0;
    tt_var_incl = 0.0;
    tt_var_az = 0.0;
	for (int n = 0; n < k; n++){
        // Calculate  ds along with r, t, p at mid point
		dr = solution[n + 1][0] - solution[n][0];   r = solution[n][0] + dr / 2.0;
		dt = solution[n + 1][1] - solution[n][1];   t = solution[n][1] + dt / 2.0;
		dp = solution[n + 1][2] - solution[n][2];   p = solution[n][2] + dp / 2.0;

        if(dp > 2.0 * Pi){
            dp -= 2.0 * Pi;
        } else if (dp < -2.0 * Pi){
            dp += 2.0 * Pi;
        }

        ds = sqrt(pow(dr, 2) + pow(r * dt, 2) + pow(r * cos(t) * dp, 2));
		
        // Calculate c_g to define travel time contribution
		nu[0] = solution[n][3] + (solution[n + 1][3] - solution[n][3]) / 2.0;
		nu[1] = solution[n][4] + (solution[n + 1][4] - solution[n][4]) / 2.0;
		nu[2] = solution[n][5] + (solution[n + 1][5] - solution[n][5]) / 2.0;
		nu_mag = sqrt(nu[0] * nu[0] + nu[1] * nu[1] + nu[2] * nu[2]);
        
        sndspd = atmo::c(r, t, p);
		cg[0] = sndspd * nu[0] / nu_mag + atmo::w(r, t, p);
		cg[1] = sndspd * nu[1] / nu_mag + atmo::v(r, t, p);
		cg[2] = sndspd * nu[2] / nu_mag + atmo::u(r, t, p);
		cg_mag = sqrt(pow(cg[0], 2) + pow(cg[1], 2) + pow(cg[2], 2));

		tt += ds / cg_mag;

        // Calculate the inclination variation
		R[0] = solution[n][6] + (solution[n + 1][6] - solution[n][6]) / 2.0;
		R[1] = solution[n][7] + (solution[n + 1][7] - solution[n][7]) / 2.0;
		R[2] = solution[n][8] + (solution[n + 1][8] - solution[n][8]) / 2.0;

		mu[0] = solution[n][9]  + (solution[n + 1][9]  - solution[n][9]) / 2.0;
		mu[1] = solution[n][10] + (solution[n + 1][10] - solution[n][10]) / 2.0;
		mu[2] = solution[n][11] + (solution[n + 1][11] - solution[n][11]) / 2.0;

        sndspd_scalar  = R[0] * atmo::dc(r, t, p, 0);
        sndspd_scalar += R[1] * atmo::dc(r, t, p, 1);
        sndspd_scalar += R[2] * atmo::dc(r, t, p, 2);
        sndspd_scalar -= sndspd / pow(nu_mag, 2) * (mu[0] * nu[0] + mu[1] * nu[1] + mu[2] * nu[2]);

        dcg[0] = nu[0] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[0];
        dcg[1] = nu[1] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[1];
        dcg[2] = nu[2] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[2];

        dcg[0] += R[0] * atmo::dw(r, t, p, 0) + R[1] * atmo::dw(r, t, p, 1) + R[2] * atmo::dw(r, t, p, 2);
        dcg[1] += R[0] * atmo::dv(r, t, p, 0) + R[1] * atmo::dv(r, t, p, 1) + R[2] * atmo::dv(r, t, p, 2);
        dcg[2] += R[0] * atmo::du(r, t, p, 0) + R[1] * atmo::du(r, t, p, 1) + R[2] * atmo::du(r, t, p, 2);

        tt_var_incl += (dcg[0] * cg[0] + dcg[1] * cg[1] * dcg[2] * cg[2]) / pow(cg_mag, 3) * ds;

        // Repeat for azimuth variation
		R[0] = solution[n][12] + (solution[n + 1][12] - solution[n][12]) / 2.0;
		R[1] = solution[n][13] + (solution[n + 1][13] - solution[n][13]) / 2.0;
		R[2] = solution[n][14] + (solution[n + 1][14] - solution[n][14]) / 2.0;

		mu[0] = solution[n][15] + (solution[n + 1][15] - solution[n][15]) / 2.0;
		mu[1] = solution[n][16] + (solution[n + 1][16] - solution[n][16]) / 2.0;
		mu[2] = solution[n][17] + (solution[n + 1][17] - solution[n][17]) / 2.0;

        sndspd_scalar  = R[0] * atmo::dc(r, t, p, 0);
        sndspd_scalar += R[1] * atmo::dc(r, t, p, 1);
        sndspd_scalar += R[2] * atmo::dc(r, t, p, 2);
        sndspd_scalar -= sndspd / pow(nu_mag, 2) * (mu[0] * nu[0] + mu[1] * nu[1] + mu[2] * nu[2]);

        dcg[0] = nu[0] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[0];
        dcg[1] = nu[1] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[1];
        dcg[2] = nu[2] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[2];

        dcg[0] += R[0] * atmo::dw(r, t, p, 0) + R[1] * atmo::dw(r, t, p, 1) + R[2] * atmo::dw(r, t, p, 2);
        dcg[1] += R[0] * atmo::dv(r, t, p, 0) + R[1] * atmo::dv(r, t, p, 1) + R[2] * atmo::dv(r, t, p, 2);
        dcg[2] += R[0] * atmo::du(r, t, p, 0) + R[1] * atmo::du(r, t, p, 1) + R[2] * atmo::du(r, t, p, 2);

        tt_var_az += (dcg[0] * cg[0] + dcg[1] * cg[1] * dcg[2] * cg[2]) / pow(cg_mag, 3) * ds;
    }
}


void geoac::travel_time_var(double & tt, double & tt_var_incl, double & tt_var_az, double** solution, int k1, int k2){
   double dr, dt, dp, ds, r, t, p, nu[3], nu_mag, sndspd, cg[3], cg_mag;
    double R[3], mu[3], sndspd_scalar, dcg[3];

	for (int n = k1; n < k2; n++){
        // Calculate  ds along with r, t, p at mid point
		dr = solution[n + 1][0] - solution[n][0];   r = solution[n][0] + dr / 2.0;
		dt = solution[n + 1][1] - solution[n][1];   t = solution[n][1] + dt / 2.0;
		dp = solution[n + 1][2] - solution[n][2];   p = solution[n][2] + dp / 2.0;
        
        if(dp > 2.0 * Pi){
            dp -= 2.0 * Pi;
        } else if (dp < -2.0 * Pi){
            dp += 2.0 * Pi;
        }

        ds = sqrt(pow(dr, 2) + pow(r * dt, 2) + pow(r * cos(t) * dp, 2));
		
        // Calculate c_g to define travel time contribution
		nu[0] = solution[n][3] + (solution[n + 1][3] - solution[n][3]) / 2.0;
		nu[1] = solution[n][4] + (solution[n + 1][4] - solution[n][4]) / 2.0;
		nu[2] = solution[n][5] + (solution[n + 1][5] - solution[n][5]) / 2.0;
		nu_mag = sqrt(nu[0] * nu[0] + nu[1] * nu[1] + nu[2] * nu[2]);
        
        sndspd = atmo::c(r, t, p);
		cg[0] = sndspd * nu[0] / nu_mag + atmo::w(r, t, p);
		cg[1] = sndspd * nu[1] / nu_mag + atmo::v(r, t, p);
		cg[2] = sndspd * nu[2] / nu_mag + atmo::u(r, t, p);
		cg_mag = sqrt(pow(cg[0], 2) + pow(cg[1], 2) + pow(cg[2], 2));

		tt += ds / cg_mag;

        // Calculate the inclination variation
		R[0] = solution[n][6] + (solution[n + 1][6] - solution[n][6]) / 2.0;
		R[1] = solution[n][7] + (solution[n + 1][7] - solution[n][7]) / 2.0;
		R[2] = solution[n][8] + (solution[n + 1][8] - solution[n][8]) / 2.0;

		mu[0] = solution[n][9]  + (solution[n + 1][9]  - solution[n][9]) / 2.0;
		mu[1] = solution[n][10] + (solution[n + 1][10] - solution[n][10]) / 2.0;
		mu[2] = solution[n][11] + (solution[n + 1][11] - solution[n][11]) / 2.0;

        sndspd_scalar  = R[0] * atmo::dc(r, t, p, 0);
        sndspd_scalar += R[1] * atmo::dc(r, t, p, 1);
        sndspd_scalar += R[2] * atmo::dc(r, t, p, 2);
        sndspd_scalar -= sndspd / pow(nu_mag, 2) * (mu[0] * nu[0] + mu[1] * nu[1] + mu[2] * nu[2]);

        dcg[0] = nu[0] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[0];
        dcg[1] = nu[1] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[1];
        dcg[2] = nu[2] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[2];

        dcg[0] += R[0] * atmo::dw(r, t, p, 0) + R[1] * atmo::dw(r, t, p, 1) + R[2] * atmo::dw(r, t, p, 2);
        dcg[1] += R[0] * atmo::dv(r, t, p, 0) + R[1] * atmo::dv(r, t, p, 1) + R[2] * atmo::dv(r, t, p, 2);
        dcg[2] += R[0] * atmo::du(r, t, p, 0) + R[1] * atmo::du(r, t, p, 1) + R[2] * atmo::du(r, t, p, 2);

        tt_var_incl -= (dcg[0] * cg[0] + dcg[1] * cg[1] * dcg[2] * cg[2]) / pow(cg_mag, 3) * ds;

        // Repeat for azimuth variation
		R[0] = solution[n][12] + (solution[n + 1][12] - solution[n][12]) / 2.0;
		R[1] = solution[n][13] + (solution[n + 1][13] - solution[n][13]) / 2.0;
		R[2] = solution[n][14] + (solution[n + 1][14] - solution[n][14]) / 2.0;

		mu[0] = solution[n][15] + (solution[n + 1][15] - solution[n][15]) / 2.0;
		mu[1] = solution[n][16] + (solution[n + 1][16] - solution[n][16]) / 2.0;
		mu[2] = solution[n][17] + (solution[n + 1][17] - solution[n][17]) / 2.0;

        sndspd_scalar  = R[0] * atmo::dc(r, t, p, 0);
        sndspd_scalar += R[1] * atmo::dc(r, t, p, 1);
        sndspd_scalar += R[2] * atmo::dc(r, t, p, 2);
        sndspd_scalar -= sndspd / pow(nu_mag, 2) * (mu[0] * nu[0] + mu[1] * nu[1] + mu[2] * nu[2]);

        dcg[0] = nu[0] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[0];
        dcg[1] = nu[1] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[1];
        dcg[2] = nu[2] / nu_mag * sndspd_scalar + sndspd / nu_mag * mu[2];

        dcg[0] += R[0] * atmo::dw(r, t, p, 0) + R[1] * atmo::dw(r, t, p, 1) + R[2] * atmo::dw(r, t, p, 2);
        dcg[1] += R[0] * atmo::dv(r, t, p, 0) + R[1] * atmo::dv(r, t, p, 1) + R[2] * atmo::dv(r, t, p, 2);
        dcg[2] += R[0] * atmo::du(r, t, p, 0) + R[1] * atmo::du(r, t, p, 1) + R[2] * atmo::du(r, t, p, 2);

        tt_var_az -= (dcg[0] * cg[0] + dcg[1] * cg[1] * dcg[2] * cg[2]) / pow(cg_mag, 3) * ds;
    }
}



//--------------------------------------------------------------------------//
//------Integrate the Sutherland Bass Attenuation Through the Ray Path------//
//--------------------------------------------------------------------------//
double geoac::atten(double ** solution, int k, double freq){
	double dr, dt, dp, ds, r, t, p;
    double atten = 0.0;
	for (int n = 0; n < k; n++){
		dr = solution[n + 1][0] - solution[n][0];   r = solution[n][0] + dr / 2.0;
		dt = solution[n + 1][1] - solution[n][1];   t = solution[n][1] + dt / 2.0;
		dp = solution[n + 1][2] - solution[n][2];   p = solution[n][2] + dp / 2.0;

        if(dp > 2.0 * Pi){
            dp -= 2.0 * Pi;
        } else if (dp < -2.0 * Pi){
            dp += 2.0 * Pi;
        }

        ds = sqrt(pow(dr, 2) + pow(r * dt, 2) + pow(r * cos(t) * dp, 2));
        
		atten += atmo::SB_alpha(r, t, p, freq) * ds;	// Add contribution to the travel time
	}
    return atten;
}

void geoac::atten(double & atten, double ** solution, int k1, int k2, double freq){
	double dr, dt, dp, ds, r, t, p;
	for (int n = k1; n < k2; n++){
		dr = solution[n + 1][0] - solution[n][0];   r = solution[n][0] + dr / 2.0;
		dt = solution[n + 1][1] - solution[n][1];   t = solution[n][1] + dt / 2.0;
		dp = solution[n + 1][2] - solution[n][2];   p = solution[n][2] + dp / 2.0;

        if(dp > 2.0 * Pi){
            dp -= 2.0 * Pi;
        } else if (dp < -2.0 * Pi){
            dp += 2.0 * Pi;
        }

        ds = sqrt(pow(dr, 2) + pow(r * dt, 2) + pow(r * cos(t) * dp, 2));
        
		atten += atmo::SB_alpha(r, t, p, freq)*ds;	// Add contribution to the travel time
	}
}


//-----------------------------------------------------------------------------------//
//-------Calculate the Jacobian determinant and from it the ampltude coefficient-----//
//-----------------------------------------------------------------------------------//
double geoac::jacobian(double ** solution, int k){
    double r, t, p;
    double nu[3], nu_mag;
    double c_prop[3], c_prop_mag;
    double dr_ds, dt_ds, dp_ds, dr_dlth, dt_dlth, dp_dlth, dr_dlph, dt_dlph, dp_dlph;
    
	r = solution[k][0]; nu[0] = solution[k][3];
    t = solution[k][1]; nu[1] = solution[k][4];
    p = solution[k][2];	nu[2] = solution[k][5];
    nu_mag = sqrt(pow(nu[0], 2) + pow(nu[1], 2) + pow(nu[2], 2));
    
	c_prop[0] = atmo::c(r, t, p) * nu[0] / nu_mag + atmo::w(r, t, p);
    c_prop[1] = atmo::c(r, t, p) * nu[1] / nu_mag + atmo::v(r, t, p);
    c_prop[2] = atmo::c(r, t, p) * nu[2] / nu_mag + atmo::u(r, t, p);
	c_prop_mag = sqrt(pow(c_prop[0], 2) + pow(c_prop[1], 2) + pow(c_prop[2], 2));
    
	dr_ds = c_prop[0] / c_prop_mag; dt_ds = 1.0 / r * c_prop[1] / c_prop_mag;   dp_ds = 1.0 / (r * cos(t)) * c_prop[2] / c_prop_mag;
	dr_dlth = solution[k][6];       dt_dlth = solution[k][7];                   dp_dlth = solution[k][8];
	dr_dlph = solution[k][12];      dt_dlph = solution[k][13];                  dp_dlph = solution[k][14];
    
	return	pow(r, 2) * cos(t) * (dr_ds * (dt_dlth * dp_dlph - dt_dlph * dp_dlth)
                                    - dr_dlth * (dt_ds * dp_dlph - dp_ds * dt_dlph)
                                        + dr_dlph * (dt_ds * dp_dlth - dp_ds * dt_dlth));
}

double geoac::amp(double ** solution, int k){
    double  r, t, p, r0, t0, p0;
    double nu[3], nu0[3], nu_mag, nu_mag0;
    double c_prop0[3], c_prop[3], c_prop_mag0, c_prop_mag;
    double num, den;
    
    r = solution[k][0]; r0 = refs.src_loc[0];   nu[0] = solution[k][3]; nu0[0] = refs.nu0 * sin(theta);
    t = solution[k][1]; t0 = refs.src_loc[1];   nu[1] = solution[k][4]; nu0[1] = refs.nu0 * cos(theta) * sin(phi);
    p = solution[k][2]; p0 = refs.src_loc[2];   nu[2] = solution[k][5]; nu0[2] = refs.nu0 * cos(theta) * cos(phi);
    
    nu_mag = (refs.c0 - nu[0] * atmo::w(r,t,p) - nu[1] * atmo::v(r,t,p) - nu[2] * atmo::u(r,t,p)) / atmo::c(r,t,p);
    nu_mag0 = refs.nu0;
    
    c_prop0[0] = refs.c0 * nu0[0] / nu_mag0 + atmo::w(r0, t0, p0);  c_prop[0] = atmo::c(r, t, p) * nu[0] / nu_mag + atmo::w(r, t, p);
    c_prop0[1] = refs.c0 * nu0[1] / nu_mag0 + atmo::v(r0, t0, p0);  c_prop[1] = atmo::c(r, t, p) * nu[1] / nu_mag + atmo::v(r, t, p);
    c_prop0[2] = refs.c0 * nu0[2] / nu_mag0 + atmo::u(r0, t0, p0);  c_prop[2] = atmo::c(r, t, p) * nu[2] / nu_mag + atmo::u(r, t, p);
    
    c_prop_mag =  sqrt(pow(c_prop[0],2) +  pow(c_prop[1],2) +  pow(c_prop[2],2));
    c_prop_mag0 = sqrt(pow(c_prop0[0],2) + pow(c_prop0[1],2) + pow(c_prop0[2],2));
    
	num = atmo::rho(r, t, p) *  nu_mag * pow(atmo::c(r, t, p), 3)  *  c_prop_mag0 * cos(theta);
	den = atmo::rho(r0, t0, p0) * nu_mag0 * pow(atmo::c(r0, t0, p0), 3) * c_prop_mag  * jacobian(solution, k);
    
	return sqrt(fabs(num / den));
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

    struct interp::linear_spline_1D ray_r;          interp::prep(ray_r, k2 - k1);
    struct interp::linear_spline_1D ray_th;         interp::prep(ray_th, k2 - k1);
    struct interp::linear_spline_1D ray_ph;         interp::prep(ray_ph, k2 - k1);

    struct interp::linear_spline_1D D;              interp::prep(D, k2 - k1);
    struct interp::linear_spline_1D beta;           interp::prep(beta, k2 - k1);
    struct interp::linear_spline_1D beta_caustic;   interp::prep(beta_caustic, k2 - k1);

    s = 0.0;
    for (int k = 0; k < (k2 - k1); k++) {
        double c = atmo::c(solution[k1 + k][0], solution[k1 + k][1], solution[k1 + k][2]);
        double rho = atmo::rho(solution[k1 + k][0], solution[k1 + k][1], solution[k1 + k][2]);

        double nu = sqrt(pow(solution[k1 + k][3], 2) + pow(solution[k1 + k][4], 2) + pow(solution[k1 + k][5], 2));

        double cg_r  = c * solution[k1 + k][3] / nu * atmo::w(solution[k1 + k][0], solution[k1 + k][1], solution[k1 + k][2]);
        double cg_th = c * solution[k1 + k][4] / nu * atmo::v(solution[k1 + k][0], solution[k1 + k][1], solution[k1 + k][2]);
        double cg_ph = c * solution[k1 + k][5] / nu * atmo::u(solution[k1 + k][0], solution[k1 + k][1], solution[k1 + k][2]);
        double cg = sqrt(pow(cg_r, 2) + pow(cg_th, 2) + pow(cg_ph, 2));

        ray_r.x_vals[k] = s;    ray_r.f_vals[k] = solution[k1 + k][0];
        ray_th.x_vals[k] = s;   ray_th.f_vals[k] = solution[k1 + k][1];
        ray_ph.x_vals[k] = s;   ray_ph.f_vals[k] = solution[k1 + k][2];

        D.x_vals[k] = s;        D.f_vals[k] = jacobian(solution, k1 + k);
        beta.x_vals[k] = s;     beta.f_vals[k] = beta0 * sqrt(fabs((D0 * rho0 * pow(cg0, 3)) / (D.f_vals[k] * rho * pow(c, 3)) * (c * pow(nu, 3)) / (c0 * pow(nu0, 3))));

        beta_caustic.x_vals[k] = s;
        if(k < 1){  beta_caustic.f_vals[k] = 0.0;}
        else {      beta_caustic.f_vals[k] = beta0 * sqrt(fabs((D0 * rho0 * pow(cg0, 3)) / ((D.f_vals[k] - D.f_vals[k - 1]) / (D.x_vals[k] - D.x_vals[k - 1]) * rho * pow(c, 3)) * (c * pow(nu, 3)) / (c0 * pow(nu0, 3))));}


        if ((k + 1) < (k2 - k1)){
            s += sqrt(pow(solution[k1 + k + 1][0] - solution[k1 + k][0], 2)
                    + pow(solution[k1 + k][0] * (solution[k1 + k + 1][1] - solution[k1 + k][1]), 2)
                    + pow(solution[k1 + k][0] * cos(solution[k1 + k][1]) * (solution[k1 + k + 1][2] - solution[k1 + k][2]), 2));
        }
    }
    beta_caustic.f_vals[0] = beta_caustic.f_vals[1];
    s_tot = s;

    interp::set(ray_r);
    interp::set(ray_th);
    interp::set(ray_ph);
    interp::set(D);
    interp::set(beta);

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
        if(fabs(interp::eval_f(s, D)) >= 0.5){
            beta_curr = interp::eval_f(s, beta);
        } else {
            beta_curr = interp::eval_f(s, beta_caustic);
        }

        // Use current value of beta to define step size for stability
        ds = ds_wvfrm / (Pi * Uf_max * beta_curr);
        ds = min(ds, 1.0);
        ds = max(ds, 1.0e-6);

        if(fabs(interp::eval_f(min(s + ds, ray_r.x_vals[(k2 - k1) - 1]), D)) >= 0.5){
            beta_next = interp::eval_f(min(s + ds, ray_r.x_vals[(k2 - k1) - 1]), beta);
        } else {
            beta_next = interp::eval_f(min(s + ds, ray_r.x_vals[(k2 - k1) - 1]), beta_caustic);
        }

        if (fabs(10.0 * ceil((s + s_prev + ds) / 10.0) - 10.0 * ceil((s + s_prev) / 10.0)) > 0.0){
            cout << '\t' << "Generating waveform at ray length: " << 10.0 * ceil((s + s_prev) / 10.0) << '\t' << "altitude: " << interp::eval_f(s, ray_r) - globe::r0<< '\t' << "scaled non-linearity factor: " << beta_curr << '\n';
        }

        // Copy the squared waveform into u_sqr and compute derivative
        // values for du_sqr = d/dt (u^2(t)) (use central difference away from edges)
        for(int n = 0; n < wvfrm::len; n++){
            u_sqr[n] = pow(u[n][1], 2);
        }

        du_sqr[0] = (u_sqr[1] - u_sqr[0]) / (u[1][0] - u[0][0]);
        for(int n = 1; n < wvfrm::len - 1; n++){
            du_sqr[n] = (u_sqr[n + 1] - u_sqr[n - 1]) / (u[n + 1][0] - u[n - 1][0]);
        }
        du_sqr[wvfrm::len - 1] = (u_sqr[wvfrm::len - 1] - u_sqr[wvfrm::len - 2]) / (u[wvfrm::len - 1][0] - u[wvfrm::len - 2][0]);

        // Compute A = FFT{u^2}
        for(int n = 0; n < wvfrm::len; n++){
            fft_time[n] = u_sqr[n];
        }
        fftw_execute(fwd_plan);
        for(int n = 0; n < fft_len; n++){
            A[n][0] = fft_spec[n][0];
            A[n][1] = fft_spec[n][1];
        }

        // Compute B = FFT{u d/dt(u^2)}
        for(int n = 0; n < wvfrm::len; n++){
            fft_time[n] = u[n][1] * du_sqr[n];
        }
        fftw_execute(fwd_plan);
        for(int n = 0; n < fft_len; n++){
            B[n][0] = fft_spec[n][0];
            B[n][1] = fft_spec[n][1];
        }

        // Compute C = FFT{(d/dt(u^2))^2}
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
        if(s + ds < s_tot) {
            if(interp::eval_f(s, D) * interp::eval_f(s + ds, D) < 0.0){
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
            // Absorpotion computed in dB as: 10 * log10(p(s + ds)) - 10 * log10(p(s)) = alpha * ds [dB] -> p(s + ds) = p(s) * 10^(alpha * ds / 20) with c/c0 factor from Lonzaga paper
            double freq = 1.0 / dt * n / wvfrm::len;
            double absorp = pow(10.0, - atmo::c(interp::eval_f(s, ray_r), interp::eval_f(s, ray_th), interp::eval_f(s, ray_ph)) / c0 * atmo::SB_alpha(interp::eval_f(s, ray_r), interp::eval_f(s, ray_th), interp::eval_f(s, ray_ph), freq) * ds / 10.0);

            U[n][0] -= Pi * freq * (ds * (beta_curr + beta_next) / 2.0 * A[n][1] + pow(ds, 2) / 2.0 * beta_curr * beta_next * B[n][1] + pow(ds, 3) / 8.0 * pow(beta_curr, 2) * beta_next * C[n][1]) * filter1[n];
            U[n][1] += Pi * freq * (ds * (beta_curr + beta_next) / 2.0 * A[n][0] + pow(ds, 2) / 2.0 * beta_curr * beta_next * B[n][0] + pow(ds, 3) / 8.0 * pow(beta_curr, 2) * beta_next * C[n][0]) * filter1[n];

            U[n][0] *= absorp * filter2[n];
            U[n][1] *= absorp * filter2[n];
        }

        // Reverse the FFT to obtain the waveform at s + ds, check the new maximum scaled waveform, and increment s -> d + ds
        for(int n = 0; n < fft_len; n++){
            fft_spec[n][0] = U[n][0];
            fft_spec[n][1] = U[n][1];
        }
        fftw_execute(bwd_plan);

        for(int n = 0; n < wvfrm::len; n++){
            u[n][1] = fft_time[n] * tukey_win[n] / wvfrm::len;
        }

        Uf_max = 0.0;
        for(int n = 0; n < fft_len; n++){
            Uf_max = max(Uf_max, (1.0 / dt * n / wvfrm::len) * sqrt(pow(U[n][0], 2) + pow(U[n][1], 2)));
        }

        s += ds;
    }


    // Clear all memory used in analysis
    interp::clear(ray_r);
    interp::clear(ray_th);
    interp::clear(ray_ph);

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


#endif /* _GEOAC_EQSET_SPH_RNGDEP_CPP_ */
