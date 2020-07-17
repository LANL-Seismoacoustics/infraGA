# ifndef _ATMO_STATE_ABS_SPH_CPP_
# define _ATMO_STATE_ABS_SPH_CPP_

#include <math.h>

#include "atmo_state.h"

#include "../util/globe.h"
#include "../geoac/geoac.params.h"

using namespace std;
double atmo::tweak_abs = 0.3; // tweak absorption alpha by this factor

double atmo::SB_alpha(double r, double t, double p, double freq){
    // Expressions based on Bass and Sutherland, JASA 2004
    // Computes alpha(freq) for given G2S output
    // Subroutine can easily be modified to include dispersion effects
    // In that case, make alpha a complex array and
    // use the real part as absorption and the imaginary part for dispersion

    int m;
    double r_earth = globe::r0;
    double T_o, P_o, S;
    double X[7], X_ON, Z_rot[2], Z_rot_;
    double sigma, nn, chi, cchi, mu, nu, mu_o;
    double beta_0, beta_1, beta_2, alpha_1, alpha_2;
    double a_cl, a_rot, a_diff, a_vib;
    double T_z, P_z, c_snd_z;
    double A1, A2, B, C, D, E, F, G, H, I, J, K, L, ZZ, hu;
    double f_vib[4], a_vib_c[4], Cp_R[4], Cv_R[4], Theta[4], C_R, A_max, Tr;
        
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
    Theta[0]= 2239.1;           // Charact. temperature (O2)
    Theta[1]= 3352.0;           // Charact. temperature (N2)
    Theta[2]= 915.0;            // Charact. temperature (CO2)
    Theta[3]= 1037.0;           // Charact. temperature (O3)
    
    T_z     = pow(c(r, t, p) * 1.0e3, 2) / (R * gam);                   // Temperature at location [K]
    P_z     = rho(r, t, p) * 1.0e3 * pow(c(r, t, p) * 1.0e3, 2) / gam;  // Pressure at location [Pa]
    c_snd_z = c(r, t, p) * 1.0e3;                                       // Sound speed at location [km/s]
    
    mu      = mu_o * sqrt(T_z / T_o) * ((1.0 + S / T_o) / (1.0 + S / T_z)); // Viscosity [kg/(m*s)]
    nu      = (8.0 * Pi * freq * mu) / (3.0 * P_z);                         // Nondimensional frequency
        
    //-------- Gas fraction polynomial fits -----------------------------------
    // O2 profile
    if ((r - r_earth) > 90.)    X[0] = pow(10.0,49.296-(1.5524*(r - r_earth)) + (1.8714E-2 * pow(r - r_earth, 2)) - (1.1069E-4 * pow(r - r_earth, 3))
                                            + (3.199E-7 * pow(r - r_earth, 4)) - (3.6211E-10 * pow(r - r_earth, 5)));
    else                        X[0] = pow(10.0, -0.67887);
        
    // N2 profile
    if ((r - r_earth) > 76.)    X[1] = pow(10.0, (1.3972E-1) - (5.6269E-3 * (r - r_earth)) + (3.9407E-5 * pow(r - r_earth, 2)) - (1.0737E-7 * pow(r - r_earth, 3)));
    else                        X[1] = pow(10.0, -0.10744);
        
    // C02 profile
    X[2]  = pow(10, -3.3979);
        
    // O3 profile
    if ((r - r_earth) > 80. )   X[3] = pow(10.0, -4.234 - (3.0975E-2 * (r - r_earth)));
    else                        X[3] = pow(10.0, -19.027 + (1.3093 * (r - r_earth)) - (4.6496E-2 * pow(r - r_earth, 2)) + (7.8543E-4 * pow(r - r_earth, 3))
                                                                - (6.5169E-6 * pow(r - r_earth, 4)) + (2.1343E-8 * pow(r - r_earth,5)));
        
    // O profile
    if ((r - r_earth) > 95. )   X[4] = pow(10.0, -3.2456 + (4.6642E-2 * (r - r_earth)) - (2.6894E-4 * pow(r - r_earth, 2)) + (5.264E-7 * pow(r - r_earth, 3)));
    else                        X[4] = pow(10.0, -11.195 + (1.5408E-1 * (r - r_earth)) - (1.4348E-3 * pow(r - r_earth, 2)) + (1.0166E-5 * pow(r - r_earth, 3)));
        
    // N profile
    X[5]  = pow(10.0, -53.746 + (1.5439 * (r - r_earth)) - (1.8824E-2 * pow(r - r_earth, 2)) + (1.1587E-4 * pow(r - r_earth, 3))
                                    - (3.5399E-7 * pow(r - r_earth, 4)) + (4.2609E-10 * pow(r - r_earth, 5)));
        
    // H20 profile
    if ((r - r_earth) > 30. )           X[6] = pow(10.0, -4.2563 + (7.6245E-2 * (r - r_earth)) - (2.1824E-3 * pow(r - r_earth, 2)) - (2.3010E-6 * pow(r - r_earth, 3))
                                                   + (2.4265E-7 * pow(r - r_earth, 4)) - (1.2500E-09 * pow(r - r_earth, 5)));
    else{  if ((r - r_earth) > 100.)    X[6] = pow(10.0, -0.62534 - (8.3665E-2 * (r - r_earth)));
           else                         X[6] = pow(10.0, -1.7491 + (4.4986E-2 * (r - r_earth)) - (6.8549E-2 * pow(r - r_earth, 2)) + (5.4639E-3 * pow(r - r_earth, 3))
                                                   - (1.5539E-4 * pow(r - r_earth, 4)) + (1.5063E-06 * pow(r - r_earth, 5)));
    }
        
    X_ON = (X[0] + X[1]) / 0.9903;
        
    //-------- Rotational collision number-------------------------------------
    Z_rot[0] = 54.1 * exp(-17.3 * (pow(T_z, -1.0/3.0)));   // O2
    Z_rot[1] = 63.3 * exp(-16.7 * (pow(T_z, -1.0/3.0)));   // N2
    Z_rot_   = 1.0 / ((X[1] / Z_rot[1]) + (X[0] / Z_rot[0]));
        
    //-------- Nondimensional atmospheric quantities---------------------------
    sigma = 5.0 / sqrt(21.0);
    nn = (4.0 / 5.0) * sqrt(3.0 / 7.0) * Z_rot_;
    chi=3.0 * nn * nu/4.0;
    cchi=2.36 * chi;
        
    //---------Classical + rotational loss/dispersion--------------------------
    beta_0  = 2.0 * Pi * freq / c_snd_z;
    beta_1  = beta_0 * sqrt(0.5 * (sqrt(1.0 + pow(nu, 2)) + 1.0) / (1.0 + pow(nu, 2)));
    beta_2  = beta_0 * sqrt((1.0 + pow(chi, 2))/(1 + pow((sigma * chi), 2)));
    alpha_1 = beta_0 * sqrt(0.5 * (sqrt(1.0 + pow(nu, 2)) - 1.0) / (1.0 + pow(nu, 2)));
    alpha_2 = beta_0 * (((sigma / 2.0 - 1.0 / (2.0 * sigma)) * chi) / (sqrt((1.0 + pow(chi, 2)) * (1.0 + pow(sigma * chi, 2)))));
    //a_cl    = alpha_1*(beta_2/beta_0);
    //a_rot = alpha_2*(beta_1/beta_0)*X_ON;
        
    a_cl    = (2.0 * Pi * freq / c_snd_z) * sqrt(0.5 * (sqrt(1.0 + pow(nu, 2)) - 1.0) * (1.0 + pow(cchi, 2)) / ((1.0 + pow(nu, 2)) * (1.0 + pow(sigma * cchi, 2))));
    a_rot   = (2.0 * Pi * freq / c_snd_z) * X_ON * ((pow(sigma, 2) - 1.0) * chi / (2 * sigma)) * sqrt(0.5 * (sqrt(1.0 + pow(nu, 2)) + 1.0) / ((1.0 + pow(nu, 2)) * (1.0 + pow(cchi, 2))));
    a_diff  = 0.003 * a_cl;
        
    //---------Vibrational relaxation-------------------------------------------
    Tr = pow(T_z / T_o, -1.0/3.0) - 1.0;
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
    K  = (8.48E08)*exp(9.17 * Tr);
    L  = exp(-7.72 * Tr);
    ZZ = H * X[2] + I * (X[0] + 0.5 * X[4]) + J * (X[1] + 0.5 * X[5]) + K * (X[6] + X[3]);
    hu = 100.0 * (X[3] + X[6]);
    f_vib[0] = (P_z / P_o) * (mu_o / mu) * (A1 + A2 + B * hu * (C + hu) * (D + hu));
    f_vib[1] = (P_z / P_o) * (mu_o / mu) * (E + F * X[3] + G * X[6]);
    f_vib[2] = (P_z / P_o) * (mu_o / mu) * ZZ;
    f_vib[3] = (P_z / P_o) * (mu_o / mu) * (1.2E5) * L;
        
    a_vib = 0.0;
    for (m=0; m<4; m++){
        C_R          = ((pow(Theta[m] / T_z, 2)) * exp(-Theta[m] / T_z)) / (pow(1.0 - exp(-Theta[m] / T_z), 2));
        A_max        = (X[m] * (Pi / 2) * C_R) / (Cp_R[m] * (Cv_R[m] + C_R));
        a_vib_c[m]   = (A_max / c_snd_z) * ((2 * (pow(freq, 2)) / f_vib[m]) / (1.0 + pow(freq / f_vib[m], 2)));
        a_vib        += a_vib_c[m];
    }
        
    // Value computed is in Np/m (nepers per meter) and is expect to be in dB/km for use here.
    // Need to scale by 8.685889e3 to convert Np/m to dB/km; also scales by a "tweak" factor to adjust absorption
    return (a_cl + a_rot + a_diff + a_vib) * tweak_abs * 8.685889e3;
}


#endif  /*_ATMO_STATE_ABS_SPH_CPP_*/
