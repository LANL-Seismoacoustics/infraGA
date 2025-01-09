# ifndef _ATMO_STATE_ABS_CART_CPP_
# define _ATMO_STATE_ABS_CART_CPP_

#include <iostream>
#include <math.h>

#include "../geoac/geoac.params.h"
#include "atmo_state.h"

double atmo::tweak_abs = 0.3; // tweak absorption alpha by this factor

double atmo::SB_alpha(double x, double y, double z, double freq){
    // Expressions based on Bass and Sutherland, JASA 2004
    // Computes alpha(freq) for given G2S output
    // Subroutine can easily be modified to include dispersion effects
    // In that case, make alpha a complex array and
    // use the real part as absorption and the imaginary part for dispersion

    int m;
    double T_o, P_o, S;
    double X[7], X_ON, Z_rot[2], Z_rot_;
    double sigma, nn, chi, cchi, mu, nu, mu_o;
    double beta_0, beta_1, beta_2, alpha_1, alpha_2;
    double a_cl, a_rot, a_vib;
    double T_z, P_z, c_snd_z;
    double A1, A2, B, C, D, E, F, G, H, I, J, K, L, ZZ, hu;
    double f_vib[4], a_vib_c[4], Cp_R[4], Cv_R[4], theta[4], C_R, A_max, Tr;
       
    // Atmospheric composition constants
    mu_o  = 18.192e-6;      // Reference viscosity [kg / (m * s)]
    T_o   = 293.15;         // Reference temperature [K]
    P_o   = 101325 ;        // Reference pressure [Pa]
    S     = 117.0;          // Sutherland constant [K]
    
    Cv_R[0] = 5.0 / 2.0;        // Heat capacity|volume (O2)
    Cv_R[1] = 5.0 / 2.0;        // Heat capacity|volume (N2)
    Cv_R[2] = 3.0;              // Heat capacity|volume (CO2)
    Cv_R[3] = 3.0;              // Heat capacity|volume (O3)

    Cp_R[0] = 7.0 / 2.0;        // Heat capacity|pressure (O2)
    Cp_R[1] = 7.0 / 2.0;        // Heat capacity|pressure (N2)
    Cp_R[2] = 4.0;              // Heat capacity|pressure (CO2)
    Cp_R[3] = 4.0;              // Heat capacity|pressure (O3)

    theta[0]= 2239.1;           // Charact. temperature (O2)
    theta[1]= 3352.0;           // Charact. temperature (N2)
    theta[2]= 915.0;            // Charact. temperature (CO2)
    theta[3]= 1037.0;           // Charact. temperature (O3)
    
    T_z     = pow(c(x,y,z) * 1.0e3, 2) / (R * gam);                     // Temperature at location [K]
    P_z     = rho(x,y,z) * 1.0e3 * pow(c(x,y,z) * 1.0e3, 2) / gam ;     // Pressure at location [Pa]
    c_snd_z = c(x,y,z) * 1.0e3;                                         // Sound speed at location [m/s]
    
    mu      = mu_o * sqrt(T_z / T_o) * ((1.0 + S / T_o) / (1.0 + S / T_z)); // Viscosity [kg/(m*s)]
    nu      = (8.0 * Pi * freq * mu) / (3.0 * P_z);                         // Nondimensional frequency
        
    //-------- Gas fraction polynomial fits -----------------------------------
    // O2 profile
    if (z > 90.0)   X[0] = pow(10.0, 49.296 - (1.5524 * z) + (1.8714e-2 * pow(z, 2)) - (1.1069e-4 * pow(z, 3)) + (3.199e-7 * pow(z, 4)) - (3.6211e-10 * pow(z, 5)));
    else            X[0] = pow(10.0, -0.67887);
        
    // N2 profile
    if (z > 76.0)   X[1] = pow(10.0, (1.3972e-1) - (5.6269e-3 * z) + (3.9407e-5 * pow(z, 2)) - (1.0737e-7 * pow(z, 3)));
    else            X[1] = pow(10.0, -0.10744);
        
    // C02 profile
    X[2]  = pow(10.0, -3.3979);
        
    // O3 profile
    if (z > 80.0)   X[3] = pow(10.0, -4.234 - (3.0975e-2 * z));
    else            X[3] = pow(10.0, -19.027 + (1.3093 * z) - (4.6496e-2 * pow(z, 2)) + (7.8543e-4 * pow(z, 3)) - (6.5169e-6 * pow(z, 4)) + (2.1343e-8 * pow(z, 5)));
        
    // O profile
    if (z > 95.0)   X[4] = pow(10.0, -3.2456 + (4.6642e-2 * z) - (2.6894e-4 * pow(z, 2)) + (5.264e-7 * pow(z, 3)));
    else            X[4] = pow(10.0, -11.195 + (1.5408e-1 * z) - (1.4348e-3 * pow(z, 2)) + (1.0166e-5 * pow(z, 3)));
        
    // N profile
    X[5]  = pow(10.0, -53.746 + (1.5439 * z) - (1.8824e-2 * pow(z, 2)) + (1.1587e-4 * pow(z, 3)) - (3.5399e-7 * pow(z, 2)) + (4.2609e-10 * pow(z, 5)));
        
    // H20 profile
    if (z > 30.0)           X[6] = pow(10.0, -4.2563 + (7.6245e-2 * z) - (2.1824e-3 * pow(z, 2)) - (2.3010e-6 * pow(z, 3)) + (2.4265e-7 * pow(z, 4)) - (1.2500e-9 * pow(z, 5)));
    else{  if (z > 100.0)   X[6] = pow(10.0, -0.62534 - (8.3665e-2 * z));
           else             X[6] = pow(10.0, -1.7491 + (4.4986e-2 * z) - (6.8549e-2 * pow(z, 2)) + (5.4639e-3 * pow(z, 3)) - (1.5539e-4 * pow(z, 4)) + (1.5063e-6 * pow(z, 5)));
    }
        
    X_ON = (X[0] + X[1]) / 0.9903;
        
    //-------- Rotational collision number-------------------------------------
    Z_rot[0] = 54.1 * exp(-17.3 * (pow(T_z, -1.0 / 3.0)));   // O2
    Z_rot[1] = 63.3 * exp(-16.7 * (pow(T_z, -1.0 / 3.0)));   // N2
    Z_rot_   = 1.0 / ((X[1] / Z_rot[1]) + (X[0] / Z_rot[0]));
        
    //-------- Nondimensional atmospheric quantities---------------------------
    sigma = 5.0 / sqrt(21.0);
    nn = (4.0 / 5.0) * sqrt(3.0 / 7.0) * Z_rot_;
    chi = 3.0 * nn * nu / 4.0;
    cchi = 2.36 * chi;
        
    //---------Classical + rotational loss/dispersion--------------------------
    beta_0  = 2.0 * Pi * freq / c_snd_z;
    beta_1  = sqrt(0.5 * (sqrt(1.0 + pow(nu, 2)) + 1.0) / (1.0 + pow(nu, 2)));
    beta_2  = sqrt((1.0 + pow(cchi, 2)) / (1.0 + pow((sigma * cchi), 2)));
    alpha_1 = beta_0 * sqrt(0.5 * (sqrt(1.0 + pow(nu, 2)) - 1.0) / (1.0 + pow(nu, 2)));
    alpha_2 = beta_0 * (((sigma / 2.0 - 1.0 / (2.0 * sigma)) * chi) / (sqrt((1.0 + pow(cchi, 2)) * (1.0 + pow(sigma * cchi, 2)))));

    a_cl    = alpha_1 * beta_2;
    a_rot = alpha_2 * beta_1 * X_ON;

    //---------Vibrational relaxation-------------------------------------------
    Tr = pow(T_z / T_o, -1.0 / 3.0) - 1.0;
    A1 = (X[0] + X[1]) * 24.0 * exp(-9.16 * Tr);
    A2 = (X[4] + X[5]) * 2400.0;
    B  = 40400.0 * exp(10.0 * Tr);
    C  = 0.02 * exp(-11.2 * Tr);
    D  = 0.391 * exp(8.41 * Tr);
    E  = 9.0 * exp(-19.9 * Tr);
    F  = 60000.0;
    G  = 28000.0 * exp(-4.17 * Tr);
    H  = 22000.0 * exp(-7.68 * Tr);
    I  = 15100.0 * exp(-10.4 * Tr);
    J  = 11500.0 * exp(-9.17 * Tr);
    K  = (8.48e08) * exp(9.17 * Tr);
    L  = exp(-7.72 * Tr);
    ZZ = H * X[2] + I * (X[0] + 0.5 * X[4]) + J * (X[1] + 0.5 * X[5]) + K * (X[6] + X[3]);
    hu = 100.0 * (X[3] + X[6]);

    f_vib[0] = (P_z / P_o) * (mu_o / mu) * (A1 + A2 + B * hu * (C + hu) * (D + hu));
    f_vib[1] = (P_z / P_o) * (mu_o / mu) * (E + F * X[3] + G * X[6]);
    f_vib[2] = (P_z / P_o) * (mu_o / mu) * ZZ;
    f_vib[3] = (P_z / P_o) * (mu_o / mu) * 1.2e5 * L;
        
    a_vib = 0.0;
    for (m = 0; m < 4; m++){
        C_R = ((pow(theta[m] / T_z, 2)) * exp(-theta[m] / T_z)) / (pow(1.0 - exp(-theta[m] / T_z), 2));
        A_max = (X[m] * (Pi / 2) * C_R) / (Cp_R[m] * (Cv_R[m] + C_R));
        a_vib += (A_max / c_snd_z) * ((2.0 * (pow(freq, 2)) / f_vib[m]) / (1.0 + pow(freq / f_vib[m], 2)));
    }

    // Value computed is in Np/m (nepers per meter) and is expect to be in dB/km for use here.
    // Need to scale by 8.685889e3 to convert Np/m to dB/km; also scales by a "tweak" factor to adjust absorption
    return (a_cl * 1.003 + a_rot + a_vib) * tweak_abs * 8.685889e3;
}

#endif  /*_ATMO_STATE_ABS_CART_CPP_*/
