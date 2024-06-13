#ifndef _GEOAC_EQSET_2D_CPP_
#define _GEOAC_EQSET_2D_CPP_

#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>

#include <fftw3.h>

#include "geoac.params.h"
#include "geoac.eqset.h"
#include "../atmo/atmo_state.h"
#include "../atmo/atmo_io.2d-rngdep.h"
#include "../util/interpolation.h"
#include "../util/waveforms.h"

using namespace std;

//----------------------------------------------//
//-------Propagate in 2D, assume stratified-----//
//----------------------------------------------//
void geoac::set_system(){
    dim = 2;
    is_strat = false;
}

//-----------------------------------------------------------//
//-----------Quantities used multiple times during-----------//
//-------ray path computations are computed separately-------//
//-----------------------------------------------------------//
namespace geoac {
    struct src_refs{
        double z0;              // Altitude of the source
        double c_eff0;          // Effective sound speed at the source
        double c_eff;           // Effective sound speed
        double dc_eff[2];       // First order derivatives of c_eff
        double ddc_eff[2][2];   // Second order derivatives of c_eff
    };

    struct src_refs refs = {0.0, 0.0, 0.0, {0.0, 0.0}, {{0.0, 0.0}, {0.0, 0.0}}};
}
//----------------------------------------------------------------------//
//-------Fill in solution[0][n] with the appropriate initial values-----//
//----------------------------------------------------------------------//
void geoac::set_initial(double ** & solution, double r0, double z0){
    refs.z0 = z0;
    refs.c_eff0 = atmo::c(0.0, 0.0, z0) + atmo::u(0.0, 0.0, z0) * cos(phi) + atmo::v(0.0, 0.0, z0) * sin(phi);

	for(int index = 0; index < eq_cnt; index++){
		switch(index){
			case(0):    solution[0][index] = r0; break;                     // r(0) = r0
			case(1):    solution[0][index] = z0; break;                     // z(0) = z0
                
            case(2):
            case(7):    solution[0][index] =  cos(theta); break;        // nu_r(0), mu_z(0) = cos(theta)
                
			case(3):    solution[0][index] =  sin(theta); break;        // nu_z(0) = sin(theta)
                
			case(4):
			case(5):    solution[0][index] =  0.0; break;               // R(0), Z(0) = 0.0
                
            case(6):    solution[0][index] =  -sin(theta); break;       // mu_r(0) = -sin(theta)

            default:    cout << "Unexpected index in Initial_Cond.  Model includes 8 variables." << '\n'; break;
		}
	}
}

//----------------------------------------------------------------------//
//-------Taylor series to more accurately deterine intercept values-----//
//----------------------------------------------------------------------//
void geoac::approx_intercept(double ** solution, int k, double* & prev){
	double dz_grnd = topo::z(solution[k-1][0]) - solution[k-1][1];                               // set dz from z_{k-1} to ground
	double dz_k = (solution[k-1][1] - solution[k-2][1])
                    - topo::dz(solution[k - 1][0]) * (solution[k-1][0] - solution[k-2][0]);  // set effective dz for step = z_k - z_{k-1} with ground slope
    
	for(int n_eq = 0; n_eq < eq_cnt; n_eq++){ prev[n_eq] = solution[k-1][n_eq] + (solution[k-1][n_eq] - solution[k-2][n_eq]) / dz_k * dz_grnd;}
}

//-------------------------------------------------------------------------//
//-------Fill in solution[0][n] with the appropriate reflection values-----//
//-------------------------------------------------------------------------//
void geoac::set_refl(double** & solution, int k_end){

	double* prev = new double [eq_cnt];
	approx_intercept(solution, k_end, prev);

    double zg, c_eff, c_eff_diff, C1, C2, dC1, dC2, ds0_dtheta;
    double a, norm, da_dtheta, dnur_ds, dnuz_ds;

    zg = topo::z(prev[0]);
    c_eff = atmo::c(0.0, 0.0, zg) + atmo::u(0.0, 0.0, zg) * cos(phi) + atmo::v(0.0, 0.0, zg) * sin(phi);
    c_eff_diff = atmo::dc(0.0, 0.0, zg, 2) + atmo::du(0.0, 0.0, zg, 2) * cos(phi) + atmo::dv(0.0, 0.0, zg, 2) * sin(phi);
    if (is_topo){
        a = topo::dz(prev[0]);
        norm = 1.0 + pow(a, 2);
        da_dtheta = (prev[4] * prev[3] - prev[5] * prev[2]) / (prev[3] - a * prev[2]) * topo::ddz(prev[0]);

        C1 = (1.0 - pow(a,2)) / norm;   dC1 = -2.0 * da_dtheta / norm * C2;
        C2 = 2.0 * a / norm;            dC2 = 2.0 * da_dtheta / norm * C1;

        ds0_dtheta = - refs.c_eff0 / c_eff * (prev[5] - a * prev[4]) / (prev[3] - a * prev[2]);
    } else {
        C1 = 1.0;  dC1 = 0.0;
        C2 = 0.0;  dC2 = 0.0;
        ds0_dtheta = - refs.c_eff0 / c_eff * prev[5] / prev[3];
    }

    dnur_ds = 0.0;
    dnuz_ds = - refs.c_eff0 / pow(c_eff, 2) * c_eff_diff;

	for(int n_eq = 0; n_eq < eq_cnt; n_eq++){
		switch(n_eq){
			case(0):    solution[0][n_eq] = prev[n_eq];     break;  // r(s0 + 0^+) = r0
			case(1):    solution[0][n_eq] = zg;             break;  // z(s0 + 0^+) = z_ground(r0)

            case(2):                                                                        // nu_r(s0 + 0^+) = C1 nu_r0 + C2 nu_z0
            case(4):    solution[0][n_eq] = C1 * prev[n_eq] + C2 * prev[n_eq + 1];   break; // R(s0 + 0^+) = C1 R0 + C2 Z0

            case(3):                                                                        // nu_z(s0 + 0^+)=	-C1 nu_z0 + C2 nu_r0
			case(5):    solution[0][n_eq] = -C1 * prev[n_eq] + C2 * prev[n_eq - 1];  break; // Z(s0 + 0^+) = 	-C1 Z0 + C2 R0

            case(6):    solution[0][n_eq] = C1 * prev[n_eq] + C2 * prev[n_eq + 1] + dC1 * prev[2] + dC2 * prev[3] + ((C1 - 1.0) * dnur_ds + C2 * dnuz_ds) * ds0_dtheta;  break;
            case(7):    solution[0][n_eq] = - C1 * prev[n_eq] + C2 * prev[n_eq - 1] -  dC1 * prev[3] + dC2 * prev[2] + (-(C1 + 1.0) * dnuz_ds + C2 * dnur_ds) * ds0_dtheta;  break;
                // mu_r(s0 + 0^+) = C1 mu_r0  + C2 mu_z0 + dC1 nu_r0 + dC2 nu_z0 + ((C1 - 1) d nu_r/ds + C2 d nu_z/ds) ds0/dtheta
                // mu_z(s0 + 0^+) = -C1 mu_z0  + C2 mu_r0 - dC1 nu_z0 + dC2 nu_r0 + (- (C1 + 1) d nu_z/ds + C2 d nu_r/ds) ds0/dtheta
                
            default:
				cout << "Unexpected index in Initial_Cond.  Model includes 6 variables." << '\n';
		}        
	}
	delete [] prev;
}

//-----------------------------------------------------------------------------------//
//-------Vary the solver step size, currently uses smaller steps near the ground-----//
//-----------------------------------------------------------------------------------//
double geoac::set_ds(double* current_values){
    double result = ds_max;
    if (current_values[1] < topo::z_bndlyr){
        result -= (ds_max - ds_min) * exp(-(current_values[1] - topo::z(current_values[0])) / 0.2);
    }
	return result;
}

//---------------------------------------//
//-------Update the source functions-----//
//---------------------------------------//
void geoac::update_refs(double ray_length, double* current_values){
    double z_eval = interp::in_interval(current_values[1], atmo::c_spline.x_vals[0], atmo::c_spline.x_vals[atmo::c_spline.length_x - 1]);
    double c, dc, ddc, u, du, ddu, v, dv, ddv;

    if(calc_amp){
        interp::eval_all(z_eval, atmo::c_spline, c, dc, ddc);
        interp::eval_all(z_eval, atmo::u_spline, u, du, ddu);
        interp::eval_all(z_eval, atmo::v_spline, v, dv, ddv);

        refs.c_eff = c + u * cos(phi) + v * sin(phi);
        refs.dc_eff = dc + du * cos(phi) + dv * sin(phi);
        refs.ddc_eff = ddc + ddu * cos(phi) + ddv * sin(phi);
    } else {
        interp::eval_all(z_eval, atmo::c_spline, c, dc);
        interp::eval_all(z_eval, atmo::u_spline, u, du);
        interp::eval_all(z_eval, atmo::v_spline, v, dv);

        refs.c_eff = c + u * cos(phi) + v * sin(phi);                                      
        refs.dc_eff = dc + du * cos(phi) + dv * sin(phi);
    }
}

//-----------------------------------------------------------//
//-------Evaluate the Source Equation For Specific Index-----//
//-----------------------------------------------------------//
double geoac::eval_src_eq(double ray_length, double* current_values, int eq_n){
	double result;

	double nu_r = current_values[2],    nu_z = current_values[3];
    double c0 = refs.c_eff0,            c = refs.c_eff;
    double dc_dr = refs.dc_eff[0],      dc_dz = refs.dc_eff[1];

    double R, Z, mu_r, mu_z;
    double ddc_ddr, ddc_ddz, ddc_drdz;
    double dc_dtheta, ddc_drdtheta, ddc_dzdtheta;

    if(calc_amp){
        R = current_values[4];  mu_r =  current_values[6];
        Z = current_values[5];  mu_z =  current_values[7];

        ddc_ddr = 0.0;
        ddc_ddz = refs.ddc_eff;
        ddc_drdz = 0.0;
        
        dc_dtheta = dc_dr * R + dc_dz * Z;
        ddc_drdtheta = ddc_ddr * R + ddc_drdz * Z;
        ddc_dzdtheta = ddc_drdz * R + ddc_ddz * Z;
    }

	switch(eq_n){
		case(0):    result = c / c0 * nu_r;     break;
		case(1):    result = c / c0 * nu_z;     break;
            
		case(2): 	result = -c0 / pow(c, 2) * dc_dr;    break;
		case(3):	result = -c0 / pow(c, 2) * dc_dz;    break;
            
		case(4):	result = dc_dtheta / c0 * nu_r + c / c0 * mu_r; break;
		case(5): 	result = dc_dtheta / c0 * nu_z + c / c0 * mu_z; break;
            
		case(6): 	result = 2.0 * c0 / pow(c,3) * dc_dr * dc_dtheta - c0/pow(c,2) * ddc_drdtheta;  break;
		case(7): 	result = 2.0 * c0 / pow(c,3) * dc_dz * dc_dtheta - c0/pow(c,2) * ddc_dzdtheta;  break;
	}
	return result;
}

//-------------------------------------------------------------------//
//-------Calculate the Hamiltonian (Eikonal) To Check For Errors-----//
//-------------------------------------------------------------------//
double geoac::eval_eikonal(double ray_length, double* current_values){
    double c = atmo::c(0.0, 0.0, current_values[1]) + atmo::u(0.0, 0.0, current_values[1]) * cos(phi) + atmo::v(0.0, 0.0, current_values[1]) * sin(phi);
    return sqrt(pow(current_values[2], 2) + pow(current_values[3], 2)) - refs.c_eff0/c;

}
double geoac::eval_eikonal(double ** solution, int k){
    double c = atmo::c(0.0, 0.0, solution[k][1]) + atmo::u(0.0, 0.0, solution[k][1]) * cos(phi) + atmo::v(0.0, 0.0, solution[k][1]) * sin(phi);
    return sqrt(pow(solution[k][2], 2) + pow(solution[k][3], 2)) - refs.c_eff0 / c;
}

//--------------------------------------------------------------------------//
//-------Check if ray has left propagation region or returned to ground-----//
//--------------------------------------------------------------------------//
bool geoac::break_check(double ** & solution, int k){
    if((solution[k][0] > rng_max) || (solution[k][1] > alt_max)){
        return true;
    }
	return false;
}

bool geoac::ground_check(double ** solution, int k){
    if (solution[k][1] <= topo::z_max){
        if(solution[k][1] < topo::z(solution[k][0])){
            return true;
        }
    }
    return false;
}

//----------------------------------------------------------------------------------//
//-------Calculate the travel time from source to location or between locations-----//
//----------------------------------------------------------------------------------//
double geoac::travel_time(double ** solution, int k){
	double dr, dz, ds, z_avg, c_eff;
	double traveltime = 0.0;

	for (int n = 0; n < k; n++){
		dr = solution[n+1][0] - solution[n][0];
		dz = solution[n+1][1] - solution[n][1];
		z_avg = solution[n][1] + dz/2.0;

		c_eff = atmo::c(0,0,z_avg) + atmo::u(0,0,z_avg) * cos(phi) + atmo::v(0,0,z_avg) * sin(phi);
		ds = sqrt(pow(dr,2) + pow(dz,2));
		traveltime += ds/c_eff;
	}
	return traveltime;	
}

void geoac::travel_time(double & traveltime, double ** solution, int k1, int k2){
	double dr, dz, ds, z_avg, c_eff;
    
	for (int n = k1; n < k2; n++){
		dr = solution[n+1][0] - solution[n][0];
		dz = solution[n+1][1] - solution[n][1];
		z_avg = solution[n][1] + dz/2.0;
        
		c_eff = atmo::c(0,0,z_avg) + atmo::u(0,0,z_avg) * cos(phi) + atmo::v(0,0,z_avg) * sin(phi);
		ds = sqrt(pow(dr,2) + pow(dz,2));
		traveltime += ds/c_eff;
	}
}

//--------------------------------------------------------------------------//
//------Integrate the Sutherland Bass Attenuation Through the Ray Path------//
//--------------------------------------------------------------------------//
double geoac::atten(double ** solution, int k, double freq){
	double dr, dz, ds, z_avg;
    double atten = 0.0;
	for (int n = 0; n < k; n++){
		dr = solution[n+1][0] - solution[n][0];
		dz = solution[n+1][1] - solution[n][1];
		z_avg = solution[n][1] + dz/2.0;
        
        ds = sqrt(pow(dr, 2) + pow(dz, 2));
		atten += atmo::SB_alpha(0.0, 0.0, z_avg, freq) * ds;
	}
    return atten;
}

void geoac::atten(double & atten, double ** solution, int k1, int k2, double freq){
	double dr, dz, ds, z_avg;
	for (int n = k1; n < k2; n++){
		dr = solution[n + 1][0] - solution[n][0];
		dz = solution[n + 1][1] - solution[n][1];
		z_avg = solution[n][1] + dz / 2.0;

        ds = sqrt(pow(dr, 2) + pow(dz, 2));
		atten += atmo::SB_alpha(0.0, 0.0, z_avg, freq) * ds;
	}
}

//-----------------------------------------------------------------------------------//
//-------Calculate the Jacobian determinant and from it the ampltude coefficient-----//
//-----------------------------------------------------------------------------------//
double geoac::jacobian(double ** solution, int k){
	double r = solution[k][0];
    double drds = atmo::c(0.0, 0.0, solution[k][1]) / refs.c_eff0 * solution[k][2], drdtheta = solution[k][4];
    double dzds = atmo::c(0.0, 0.0, solution[k][1]) / refs.c_eff0 * solution[k][3], dzdtheta = solution[k][5];
		
	return r * (drds * dzdtheta - dzds * drdtheta);
}


double geoac::amp(double ** solution, int k){
	double z = solution[k][1];
	
	double num = atmo::rho(0.0, 0.0, z) * atmo::c(0.0, 0.0, z) * cos(theta);
	double den = atmo::rho(0.0, 0.0, refs.z0) * refs.c_eff0 * jacobian(solution, k);
	
	return sqrt(fabs(num/den));
}

//------------------------------------------------//
//-------Calculate the Out of Plane Deviation-----//
//------------------------------------------------//
double geoac::est_dev(double ** solution, int k){
	double ds, z_avg, c_eff, v_perp;
	double dev = 0.0;
    
	for (int n = 0; n < k; n++){
		ds = sqrt(pow(solution[n + 1][0] - solution[n][0],2) + pow(solution[n + 1][1] - solution[n][1],2));
		z_avg = (solution[n + 1][1] + solution[n][1])/2.0;

		c_eff = atmo::c(0, 0, z_avg) + atmo::u(0,0,z_avg) * cos(phi) + atmo::v(0,0,z_avg) * sin(phi);
        v_perp = - atmo::u(0, 0, z_avg) * sin(phi) + atmo::v(0, 0, z_avg) * cos(phi);
        
		dev += v_perp / c_eff * 1.0 / sqrt(1.0  + 2.0 * v_perp / c_eff + (pow(atmo::u(0,0,z_avg),2) + pow(atmo::v(0,0,z_avg),2)) / pow(c_eff,2)) * ds;
	}
	return dev;
}

//------------------------------------------------//
//-------Calculate Weakly Non-Linear Waveform-----//
//------------------------------------------------//
double geoac::wnl_wvfrm(double** solution, double** & u, int k1, int k2, double s_prev, double c0, double rho0, double D0, double p0, double output_step){
    // NOTE: function assumes specified u(t, s1) is the scaled pressure waveform.

    int fft_len = wvfrm::len / 2 + 1, step_cnt;
    double s, ds, s_tot, tukey_param;
    double c, beta0, dt, t0, Uf_max, beta_curr, beta_next;

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
    beta0 = 1.2 * p0 / ((rho0 * 1.0e3) * pow(c0 * 1.0e3, 2));

    struct interp::linear_spline_1D alt;            interp::prep(alt, k2 - k1);
    struct interp::linear_spline_1D D;              interp::prep(D, k2 - k1);
    struct interp::linear_spline_1D beta;           interp::prep(beta, k2 - k1);
    struct interp::linear_spline_1D beta_caustic;   interp::prep(beta_caustic, k2 - k1);

    s = 0.0;
    for (int k = 0; k < (k2 - k1); k++) {
        alt.x_vals[k] = s;
        alt.f_vals[k] = solution[k1 + k][1];

        D.x_vals[k] = s;
        D.f_vals[k] = jacobian(solution, k1 + k);

        c = atmo::c(0.0, 0.0, solution[k1 + k][1]) + atmo::u(0.0, 0.0, solution[k1 + k][1]) * cos(phi) + atmo::v(0.0, 0.0, solution[k1 + k][1]) * sin(phi);

        beta.x_vals[k] = s;
        beta.f_vals[k] = beta0 / c0 * sqrt(fabs((D0 * rho0 * c0) / (D.f_vals[k] * atmo::rho(0.0, 0.0, solution[k1 + k][1]) * c)));

        beta_caustic.x_vals[k] = s;
        if(k < 1){  beta_caustic.f_vals[k] = 0.0;}
        else {      beta_caustic.f_vals[k] = beta0 / c0 * sqrt(fabs((D0 * rho0 * c0) / ((D.f_vals[k] - D.f_vals[k - 1]) / (D.x_vals[k] - D.x_vals[k - 1]) * atmo::rho(0.0, 0.0, solution[k1 + k][1]) * c)));}

        if ((k + 1) < (k2 - k1)){
            s += sqrt(pow(solution[k1 + k + 1][0] - solution[k1 + k][0], 2) + pow(solution[k1 + k + 1][1] - solution[k1 + k][1], 2));
        }
    }
    beta_caustic.f_vals[0] = beta_caustic.f_vals[1];
    s_tot = s;

    interp::set(alt);
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
    step_cnt = 0;
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
        ds = max(ds, 1.0e-7);

        if(fabs(interp::eval_f(min(s + ds, alt.x_vals[(k2 - k1) - 1]), D)) >= 0.5){
            beta_next = interp::eval_f(min(s + ds, alt.x_vals[(k2 - k1) - 1]), beta);
        } else {
            beta_next = interp::eval_f(min(s + ds, alt.x_vals[(k2 - k1) - 1]), beta_caustic);
        }

        if (fabs(10.0 * ceil((s + s_prev + ds) / 10.0) - 10.0 * ceil((s + s_prev) / 10.0)) > 0.0){
            cout << '\t' << "Generating waveform at ray length: " << 10.0 * ceil((s + s_prev) / 10.0) << '\t' << "altitude: " << interp::eval_f(s, alt) << '\t' << "scaled non-linearity factor: " << beta_curr << '\n';
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
            // Absorpotion computed in dB as: 10 * log10(p(s + ds)) - 10 * log10(p(s)) = alpha * ds [dB] -> p(s + ds) = p(s) * 10^(alpha * ds / 10) with c/c0 factor from Lonzaga paper
            double freq = 1.0 / dt * n / wvfrm::len;
            double absorp = pow(10.0, - c / c0 * atmo::SB_alpha(0.0, 0.0, interp::eval_f(s, alt), freq) * ds / 10.0);

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

        step_cnt++;
        s += ds;
    }

    // Clear all memory used in analysis
    interp::clear(alt);
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



#endif /* _GEOAC_EQSET_2D_CPP_ */
