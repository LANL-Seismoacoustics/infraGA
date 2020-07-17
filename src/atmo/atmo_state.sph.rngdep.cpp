# ifndef _ATMO_STATE_SPH_RNGDEP_CPP_
# define _ATMO_STATE_SPH_RNGDEP_CPP_

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <sstream>
#include <math.h>

#include "atmo_state.h"
#include "atmo_io.sph.rngdep.h"

#include "../util/interpolation.h"
#include "../util/fileIO.h"
#include "../util/globe.h"

#include "../geoac/geoac.params.h"

using namespace std;

//----------------------------//
//-----Physical Constants-----//
//----------------------------//
double atmo::gam = 1.4;
double atmo::R = 287.058;


//---------------------------------------------------//
//---------Functions Defining the Topography---------//
//--------------For 3D Global Propgation-------------//
//---------------------------------------------------//
//-------- NOTE: z0, z_max, and z_bndlyr are --------//
//---defined by topogrpahy (relative to ellipsoid)---//
//---------------------------------------------------//
double topo::z0 = 0.0;
double topo::z_max;
double topo::z_bndlyr;
bool topo::use_BLw = false;

void topo::set_bndlyr(){
    if (geoac::is_topo){
        z_max = 0.0;
        for (int n1 = 0; n1 < spline.length_x; n1++){
        for (int n2 = 0; n2 < spline.length_y; n2++){
            z_max = max(z_max, spline.f_vals[n1][n2]);
        }}
        z_bndlyr = z_max + 2.0;
        z_max += 0.01;
    } else {
        z_max = z0;
        z_bndlyr = z_max + 2.0;
    }
}

double topo::z(double lat, double lon){
    double result = globe::r0;
    if(geoac::is_topo){
        result += interp::eval_f(interp::in_interval(lat, spline.x_vals[0], spline.x_vals[spline.length_x - 1]),
                                 interp::in_interval(lon, spline.y_vals[0], spline.y_vals[spline.length_y - 1]), spline);
    } else {
        result += z0;
    }
    return result;
}

double topo::dz(double lat, double lon, int n){
    double result = 0.0;
    if (geoac::is_topo){
        result += interp::eval_df(interp::in_interval(lat, spline.x_vals[0], spline.x_vals[spline.length_x - 1]),
                                  interp::in_interval(lon, spline.y_vals[0], spline.y_vals[spline.length_y - 1]), n, spline);
    }
    return result;
}

double topo::ddz(double lat, double lon, int n1, int n2){
    double result = 0.0;
    if (geoac::is_topo){
        result += interp::eval_ddf(interp::in_interval(lat, spline.x_vals[0], spline.x_vals[spline.length_x - 1]),
                                   interp::in_interval(lon, spline.y_vals[0], spline.y_vals[spline.length_y - 1]), n1, n2, spline);
    }
    return result;
}
double topo::dddz(double lat, double lon, int n1, int n2, int n3){
    double result = 0.0;
    if (geoac::is_topo){
        return interp::eval_dddf(interp::in_interval(lat, spline.x_vals[0], spline.x_vals[spline.length_x - 1]),
                                 interp::in_interval(lon, spline.y_vals[0], spline.y_vals[spline.length_y - 1]), n1, n2, n3, spline);
    }
    return result;
}


double topo::bndlyr_param = 0.1;
double topo::vert_wind_grad = 0.5;
double topo::bndlyr_sc(double r, double lat, double lon){
    double X, result = 1.0;
    if (r - globe::r0 < topo::z_bndlyr) {
        X = exp(-(r - z(lat, lon)) / bndlyr_param);
        result = (1.0 - X) / (1.0 + X);
    }
    return result;
}
double topo::bndlyr_dsc(double r, double lat, double lon, int n){
    double X, dX, result = 0.0;
    if (r - globe::r0 < topo::z_bndlyr) {
        X = exp(-(r - z(lat, lon)) / bndlyr_param);
        if(n > 0){  dX = X / bndlyr_param * dz(lat, lon, n - 1);}
        else{       dX = - X / bndlyr_param;}

        result = - 2.0 * dX / pow(1.0 + X, 2);
    }
    return result;
}
double topo::bndlyr_ddsc(double r, double lat, double lon, int n1, int n2){
    double X, dX1, dX2, ddX, result = 0.0;
    
    if (r - globe::r0 < topo::z_bndlyr) {
        X = exp(-(r - z(lat, lon)) / bndlyr_param);
        
        if(n1 > 0){ dX1 = X / bndlyr_param * dz(lat, lon, n1 - 1);} else{ dX1 = - X / bndlyr_param;}
        if(n2 > 0){ dX2 = X / bndlyr_param * dz(lat, lon, n2 - 1);} else{ dX2 = - X / bndlyr_param;}
        
        if(n1 > 0 && n2 > 0){       ddX = 1.0 / bndlyr_param * (ddz(lat, lon, n1 - 1, n2 - 1) + 1.0 / bndlyr_param * dz(lat, lon, n1 - 1) * dz(lat, lon, n2 - 1)) * X;}
        else if(n1 > 0 && n2 == 0){ ddX = - X / pow(bndlyr_param, 2) * dz(lat, lon, n1 - 1);}
        else if(n1 == 0 && n2 > 0){ ddX = - X / pow(bndlyr_param, 2) * dz(lat, lon, n2 - 1);}
        else {                      ddX = X / pow(bndlyr_param, 2);}
        
        result = - 2.0 * (ddX / pow(1.0 + X, 2) - 2.0 * dX1 * dX2 / pow(1.0 + X, 3));
    }
    return result;
}


//--------------------------------------------------//
//---------Functions defining the atmosphere--------//
//------for global range dependent propagation------//
//--------------------------------------------------//
//---------Note that while the interpolated---------//
//-------specifications are in (lat, lon, z),-------//
//----atmospheric functions are in (r, lat, lon)----//
//--------------------------------------------------//
//-------Derivative indices are wrapped by 1:-------//
//--------------x, y, z <--> r, th, ph--------------//
//---------------0, 1, 2 <--> 1, 2, 0---------------//
//--------------------------------------------------//

double atmo::rho(double r, double lat, double lon){
    double lat_eval = interp::in_interval(lat, rho_spline.x_vals[0], rho_spline.x_vals[rho_spline.length_x - 1]);
    double lon_eval = interp::in_interval(lon, rho_spline.y_vals[0], rho_spline.y_vals[rho_spline.length_y - 1]);
    double z_eval = interp::in_interval(r - globe::r0, rho_spline.z_vals[0], rho_spline.z_vals[rho_spline.length_z - 1]);
    return interp::eval_f(lat_eval, lon_eval, z_eval, rho_spline);
}

double atmo::c(double r, double lat, double lon){
    double lat_eval = interp::in_interval(lat, c_spline.x_vals[0], c_spline.x_vals[c_spline.length_x - 1]);
    double lon_eval = interp::in_interval(lon, c_spline.y_vals[0], c_spline.y_vals[c_spline.length_y - 1]);
    double z_eval = interp::in_interval(r - globe::r0, c_spline.z_vals[0], c_spline.z_vals[c_spline.length_z - 1]);

    return interp::eval_f(lat_eval, lon_eval, z_eval, c_spline);
}

double atmo::dc(double r, double lat, double lon, int n){
    double lat_eval = interp::in_interval(lat, c_spline.x_vals[0], c_spline.x_vals[c_spline.length_x - 1]);
    double lon_eval = interp::in_interval(lon, c_spline.y_vals[0], c_spline.y_vals[c_spline.length_y - 1]);
    double z_eval = interp::in_interval(r - globe::r0, c_spline.z_vals[0], c_spline.z_vals[c_spline.length_z - 1]);

    int n_eval;
    if(n == 0){
        n_eval = 2;
    } else {
        n_eval = n - 1;
    }
    
    return interp::eval_df(lat_eval, lon_eval, z_eval, n_eval, c_spline);
}

double atmo::ddc(double r, double lat, double lon, int n1, int n2){
    double lat_eval = interp::in_interval(lat, c_spline.x_vals[0], c_spline.x_vals[c_spline.length_x - 1]);
    double lon_eval = interp::in_interval(lon, c_spline.y_vals[0], c_spline.y_vals[c_spline.length_y - 1]);
    double z_eval = interp::in_interval(r - globe::r0, c_spline.z_vals[0], c_spline.z_vals[c_spline.length_z - 1]);

    int n1_eval, n2_eval;
    if(n1 == 0){
        n1_eval = 2;
    } else {
        n1_eval = n1 - 1;
    }
    
    if(n2 == 0){
        n2_eval = 2;
    } else {
        n2_eval = n2 - 1;
    }

    return interp::eval_ddf(lat_eval, lon_eval, z_eval, n1_eval, n2_eval, c_spline);
}

double atmo::u(double r, double lat, double lon){
    double lat_eval = interp::in_interval(lat, u_spline.x_vals[0], u_spline.x_vals[u_spline.length_x - 1]);
    double lon_eval = interp::in_interval(lon, u_spline.y_vals[0], u_spline.y_vals[u_spline.length_y - 1]);
    double z_eval = interp::in_interval(r - globe::r0, u_spline.z_vals[0], u_spline.z_vals[u_spline.length_z - 1]);
    
    return interp::eval_f(lat_eval, lon_eval, z_eval, u_spline) * topo::bndlyr_sc(r, lat, lon);
}
double atmo::du(double r, double lat, double lon, int n){
    double result;
    double lat_eval = interp::in_interval(lat, u_spline.x_vals[0], u_spline.x_vals[u_spline.length_x - 1]);
    double lon_eval = interp::in_interval(lon, u_spline.y_vals[0], u_spline.y_vals[u_spline.length_y - 1]);
    double z_eval  =  interp::in_interval(r - globe::r0, u_spline.z_vals[0], u_spline.z_vals[u_spline.length_z - 1]);
    
    int n_eval;
    if(n == 0){
        n_eval = 2;
    } else {
        n_eval = n - 1;
    }
    
    result  = interp::eval_df(lat_eval, lon_eval, z_eval, n_eval, u_spline) * topo::bndlyr_sc(r, lat, lon);
    result += interp::eval_f(lat_eval, lon_eval, z_eval, u_spline) * topo::bndlyr_dsc(r, lat, lon, n);
    return result;
}
double atmo::ddu(double r, double lat, double lon, int n1, int n2){
    double result;
    double lat_eval = interp::in_interval(lat, u_spline.x_vals[0], u_spline.x_vals[u_spline.length_x - 1]);
    double lon_eval = interp::in_interval(lon, u_spline.y_vals[0], u_spline.y_vals[u_spline.length_y - 1]);
    double z_eval = interp::in_interval(r - globe::r0, u_spline.z_vals[0], u_spline.z_vals[u_spline.length_z - 1]);

    int n1_eval, n2_eval;
    if(n1 == 0){
        n1_eval = 2;
    } else {
        n1_eval = n1 - 1;
    }

    if(n2 == 0){
        n2_eval = 2;
    } else {
        n2_eval = n2 - 1;
    }
    
    result  = interp::eval_ddf(lat_eval, lon_eval, z_eval, n1_eval, n2_eval, u_spline) * topo::bndlyr_sc(r, lat, lon);
    result += interp::eval_f(lat_eval, lon_eval, z_eval, u_spline) * topo::bndlyr_ddsc(r, lat, lon, n1, n2);
    if(n1_eval == n2_eval){
        result += 2.0 * interp::eval_df(lat_eval, lon_eval, z_eval, n1_eval, u_spline) * topo::bndlyr_dsc(r, lat, lon, n1);
    } else {
        result += interp::eval_df(lat_eval, lon_eval, z_eval, n1_eval, u_spline) * topo::bndlyr_dsc(r, lat, lon, n2);
        result += interp::eval_df(lat_eval, lon_eval, z_eval, n2_eval, u_spline) * topo::bndlyr_dsc(r, lat, lon, n1);
    }
    return result;
}

double atmo::v(double r, double lat, double lon){
    double lat_eval = interp::in_interval(lat, v_spline.x_vals[0], v_spline.x_vals[v_spline.length_x - 1]);
    double lon_eval = interp::in_interval(lon, v_spline.y_vals[0], v_spline.y_vals[v_spline.length_y - 1]);
    double z_eval = interp::in_interval(r - globe::r0, v_spline.z_vals[0], v_spline.z_vals[v_spline.length_z - 1]);
    
    return interp::eval_f(lat_eval, lon_eval, z_eval, v_spline) * topo::bndlyr_sc(r, lat, lon);
}
double atmo::dv(double r, double lat, double lon, int n){
    double result;
    double lat_eval = interp::in_interval(lat, v_spline.x_vals[0], v_spline.x_vals[v_spline.length_x - 1]);
    double lon_eval = interp::in_interval(lon, v_spline.y_vals[0], v_spline.y_vals[v_spline.length_y - 1]);
    double z_eval = interp::in_interval(r - globe::r0, v_spline.z_vals[0], v_spline.z_vals[v_spline.length_z - 1]);
 
    int n_eval;
    if(n == 0){
        n_eval = 2;
    } else {
        n_eval = n - 1;
    }
    
    result  = interp::eval_df(lat_eval, lon_eval, z_eval, n_eval, v_spline) * topo::bndlyr_sc(r, lat, lon);
    result += interp::eval_f(lat_eval, lon_eval, z_eval, v_spline) * topo::bndlyr_dsc(r, lat, lon, n);
    return result;
}
double atmo::ddv(double r, double lat, double lon, int n1, int n2){
    double result;
    double lat_eval = interp::in_interval(lat, v_spline.x_vals[0], v_spline.x_vals[v_spline.length_x - 1]);
    double lon_eval = interp::in_interval(lon, v_spline.y_vals[0], v_spline.y_vals[v_spline.length_y - 1]);
    double z_eval = interp::in_interval(r - globe::r0, v_spline.z_vals[0], v_spline.z_vals[v_spline.length_z - 1]);

    int n1_eval, n2_eval;
    if(n1 == 0){
        n1_eval = 2;
    } else {
        n1_eval = n1 - 1;
    }
    
    if(n2 == 0){
        n2_eval = 2;
    } else {
        n2_eval = n2 - 1;
    }
    
    result  = interp::eval_ddf(lat_eval, lon_eval, z_eval, n1_eval, n2_eval, v_spline) * topo::bndlyr_sc(r, lat, lon);
    result += interp::eval_f(lat_eval, lon_eval, z_eval, v_spline) * topo::bndlyr_ddsc(r, lat, lon, n1, n2);
    if(n1_eval == n2_eval){
        result += 2.0 * interp::eval_df(lat_eval, lon_eval, z_eval, n1_eval, v_spline) * topo::bndlyr_dsc(r, lat, lon, n1);
    } else {
        result += interp::eval_df(lat_eval, lon_eval, z_eval, n1_eval, v_spline) * topo::bndlyr_dsc(r, lat, lon, n2);
        result += interp::eval_df(lat_eval, lon_eval, z_eval, n2_eval, v_spline) * topo::bndlyr_dsc(r, lat, lon, n1);
    }
    return result;
}

double atmo::w(double r, double lat, double lon){
    double z_eval, lat_eval, lon_eval, zg, dzg[2];
    double X, S, dS[2], w0, sg_ratio, result = 0.0;

    if(geoac::is_topo && (r - globe::r0 < topo::z_bndlyr) && topo::use_BLw){
        lat_eval = interp::in_interval(lat, v_spline.x_vals[0], v_spline.x_vals[v_spline.length_x - 1]);
        lon_eval = interp::in_interval(lon, v_spline.y_vals[0], v_spline.y_vals[v_spline.length_y - 1]);
        z_eval = interp::in_interval(r - globe::r0, v_spline.z_vals[0], v_spline.z_vals[v_spline.length_z - 1]);

        zg = topo::z(lat, lon);
        dzg[0] = topo::dz(lat, lon, 0);
        dzg[1] = topo::dz(lat, lon, 1);

        sg_ratio = topo::bndlyr_param / topo::vert_wind_grad;
        
        X = exp(-(r - zg) / topo::bndlyr_param);
        S = (1.0 - X) / (1.0 + X);

        dS[0] = -2.0 / topo::bndlyr_param * (X / pow(1.0 + X, 2)) * dzg[0];
        dS[1] = -2.0 / topo::bndlyr_param * (X / pow(1.0 + X, 2)) * dzg[1];

        w0  = interp::eval_f(lat_eval, lon_eval, z_eval, v_spline) * dS[0];
        w0 += interp::eval_f(lat_eval, lon_eval, z_eval, u_spline) * dS[1] / cos(lat);
        w0 *= -2.0 * topo::bndlyr_param / (globe::r0 + zg);

        result  = S * X * w0;
        result /= (X + (sg_ratio / 2.0) * (1.0 - pow(X, 2)));
    }
    return result;

}

double atmo::dw(double r, double lat, double lon, int n){
    double z_eval, lat_eval, lon_eval, X, dX[3], ddX[2], S, dS[3], ddS[2], N, D, sg_ratio, result = 0.0;
    double zg, dzg[2], ddzg[2], u0, v0, w0, dw0, dN, dD;

    if (geoac::is_topo && (r - globe::r0 < topo::z_bndlyr) && topo::use_BLw){
        lat_eval = interp::in_interval(lat, v_spline.x_vals[0], v_spline.x_vals[v_spline.length_x - 1]);
        lon_eval = interp::in_interval(lon, v_spline.y_vals[0], v_spline.y_vals[v_spline.length_y - 1]);
        z_eval = interp::in_interval(r - globe::r0, v_spline.z_vals[0], v_spline.z_vals[v_spline.length_z - 1]);

        sg_ratio = topo::bndlyr_param / topo::vert_wind_grad;

        // Compute the topogrpahy and derivatives       
        zg = topo::z(lat, lon);
        dzg[0] = topo::dz(lat, lon, 0);
        dzg[1] = topo::dz(lat, lon, 1);
        
        // Evaluate the u and v splines
        u0 = interp::eval_f(lat_eval, lon_eval, z_eval, u_spline);
        v0 = interp::eval_f(lat_eval, lon_eval, z_eval, v_spline);

        // Compute the boundary layer (BL) scaling        
        X = exp(-(r - zg) / topo::bndlyr_param);
        for(int j = 0; j < 3; j++){
            if(j > 0){  dX[j] = X / topo::bndlyr_param * dzg[j - 1];}
            else {      dX[j] = -X / topo::bndlyr_param;}
        }
        S = (1.0 - X) / (1.0 + X);
        for(int j = 0; j < 3; j++){
            dS[j] = -2.0 * dX[j] / pow(1.0 + X, 2);
        }

        // ddS is dlatdn and dlondn
        if(n > 0){ 
            ddX[0] = X / topo::bndlyr_param * (topo::ddz(lat, lon, 0, n - 1) + 1.0 / topo::bndlyr_param * dzg[0] * dzg[n - 1]);
            ddX[1] = X / topo::bndlyr_param * (topo::ddz(lat, lon, 1, n - 1) + 1.0 / topo::bndlyr_param * dzg[1] * dzg[n - 1]);
        } else {
            ddX[0] = - X / pow(topo::bndlyr_param, 2) * dzg[0];
            ddX[1] = - X / pow(topo::bndlyr_param, 2) * dzg[1];
        }

        ddS[0] = -2.0 * (ddX[0] / pow(1.0 + X, 2) - 2.0 * dX[1] * dX[n] / pow(1.0 + X, 3));
        ddS[1] = -2.0 * (ddX[1] / pow(1.0 + X, 2) - 2.0 * dX[2] * dX[n] / pow(1.0 + X, 3));
     
        w0  = v0 * dS[1] + u0 * dS[2] / cos(lat);
        w0 *= -2.0 * topo::bndlyr_param / (globe::r0 + zg);

        N = S * X * w0;
        D = X + (sg_ratio / 2.0) * (1.0 - pow(X, 2));

        dw0 = v0 * ddS[0] + u0 * ddS[1] / cos(lat);
        if (n == 0) {
            dw0 += interp::eval_df(lat_eval, lon_eval, z_eval, 2, v_spline) * dS[1];
            dw0 += interp::eval_df(lat_eval, lon_eval, z_eval, 2, u_spline) * dS[2] / cos(lat);
        } else if (n == 1){
            dw0 += u0 * dS[2] * sin(lat) / pow(cos(lat), 2);
        }
        dw0 *= -2.0 * topo::bndlyr_param / (globe::r0 + zg);
            
        dN = (dS[n] * X + S * dX[n]) * w0 + S * X * dw0;
        dD = (1.0 - sg_ratio * X) * dX[n];
        result = (dN * D - N * dD) / pow(D, 2);
    }
    return result;
}


double atmo::ddw(double r, double lat, double lon, int n1, int n2){
    double z_eval, lat_eval, lon_eval, X, S, N, D, sg_ratio, result = 0.0;
    double zg, dzg[2], u0, v0, w0, du0[3], dv0[3],  dw0[2], dX[3], ddX[3][3], dddX, dS[3],  dN[2], dD[2];
    double ddzg[2][2], ddu0[3][3], ddv0[3][3], ddS[3][3], dddS[2], ddN, ddD, ddw0;

    if (geoac::is_topo && (r - globe::r0 < topo::z_bndlyr) && topo::use_BLw){
        lat_eval = interp::in_interval(lat, u_spline.x_vals[0], u_spline.x_vals[u_spline.length_x - 1]);
        lon_eval = interp::in_interval(lon, u_spline.y_vals[0], u_spline.y_vals[u_spline.length_y - 1]);
        z_eval = interp::in_interval(r - globe::r0, u_spline.z_vals[0], u_spline.z_vals[u_spline.length_z - 1]);
        sg_ratio = topo::bndlyr_param / topo::vert_wind_grad;

        // Compute the topogrpahy and derivatives       
        zg = topo::z(lat, lon);
        for(int j = 0; j < 2; j++){
            dzg[j] = topo::dz(lat, lon, j);
            for(int k = j; k < 2; k++){
                ddzg[j][k] = topo::ddz(lat, lon, j, k);
                ddzg[k][j] = ddzg[j][k];
            }
        }

        // Evaluate the u and v splines
        u0 = interp::eval_f(lat_eval, lon_eval, z_eval, u_spline);
        v0 = interp::eval_f(lat_eval, lon_eval, z_eval, v_spline);
        for(int n = 0; n < 3; n++){
            if(n == 0){
                du0[n] = interp::eval_df(lat_eval, lon_eval, z_eval, 2, u_spline);
                dv0[n] = interp::eval_df(lat_eval, lon_eval, z_eval, 2, v_spline);
            } else {
                du0[n] = 0.0;
                dv0[n] = 0.0;
            }
            for(int m = 0; m < 3; m++){
                if((n == 0) && (m == 0)){
                    ddu0[n][m] = interp::eval_ddf(lat_eval, lon_eval, z_eval, 2, 2, u_spline);
                    ddv0[n][m] = interp::eval_ddf(lat_eval, lon_eval, z_eval, 2, 2, v_spline);
                } else {
                    ddu0[n][m] = 0.0;
                    ddu0[n][m] = 0.0;
                }
            }
        }

       // Compute the boundary layer (BL) scaling and derivatives
        X = exp(-(r - zg) / topo::bndlyr_param);
        for(int n = 0; n < 3; n++){
            if(n > 0){  dX[n] = X / topo::bndlyr_param * dzg[n - 1];}
            else {      dX[n] = -X / topo::bndlyr_param;}
            for(int m = 0; m < 3; m++){
                if((n > 0) && (m > 0)){         ddX[n][m] = X / topo::bndlyr_param * (ddzg[n - 1][m - 1] + 1.0 / topo::bndlyr_param * dzg[n - 1] * dzg[m - 1]);}
                else if((n == 0) && (m > 0)){   ddX[n][m] = - X / pow(topo::bndlyr_param, 2) * dzg[m - 1];}
                else if((n > 0) && (m == 0)){   ddX[n][m] = - X / pow(topo::bndlyr_param, 2) * dzg[n - 1];}
                else {                          ddX[n][m] = X / pow(topo::bndlyr_param, 2);}
            }
        }
        S = (1.0 - X) / (1.0 + X);
        for(int n = 0; n < 3; n++){
            dS[n] = -2.0 * dX[n] / pow(1.0 + X, 2);
            for(int m = 0; m < 3; m++){
                ddS[n][m] = -2.0 * (ddX[n][m] / pow(1.0 + X, 2) - 2.0 * dX[n] * dX[m] / pow(1.0 + X, 3));
            }
        }

        // Compute w0, N, and D along with their n1 and n2 derivatives
        w0  = v0 * dS[1] + u0 * dS[2] / cos(lat);
        w0 *= -2.0 * topo::bndlyr_param / (globe::r0 + zg);

        sg_ratio = topo::bndlyr_param / topo::vert_wind_grad;
        
        N = S * X * w0;
        D = X + (sg_ratio / 2.0) * (1.0 - pow(X, 2));

        dw0[0] = v0 * ddS[1][n1] + u0 * ddS[2][n1] / cos(lat);
        if (n1 == 0){
            dw0[0] += dv0[0] * dS[1] + du0[0] * dS[2] / cos(lat);
        } else if (n1 == 1){
            dw0[0] += sin(lat) / pow(cos(lat), 2) * u0 * dS[2];
        }
        dw0[0] *= -2.0 * topo::bndlyr_param / (globe::r0 + zg);
        dN[0] = (dS[n1] * X + S * dX[n1]) * w0 + S * X * dw0[0];
        dD[0] = (1.0 - sg_ratio * X) * dX[n1];

        if(n2 == n1){
            dw0[1] = dw0[0];
            dN[1] = dN[0];
            dD[1] = dD[0];
        } else {
            dw0[1] = v0 * ddS[1][n2] + u0 * ddS[2][n2] / cos(lat);
            if (n2 == 0){
                dw0[1] += dv0[0] * dS[1] + du0[0] * dS[2] / cos(lat);
            } else if (n2 == 1){
                dw0[1] += sin(lat) / pow(cos(lat), 2) * u0 * dS[2];
            }
            dw0[1] *= -2.0 * topo::bndlyr_param / (globe::r0 + zg);
            dN[1] = (dS[n2] * X + S * dX[n2]) * w0 + S * X * dw0[1];
            dD[1] = (1.0 - sg_ratio * X) * dX[n2];
        }
        
        if (n1 > 0 && n2 > 0){
            for(int k = 0; k < 2; k++){
                dddX  = topo::dddz(lat, lon, k, n1 - 1, n2 - 1);
                dddX += (ddzg[n1 - 1][n2 - 1] * dzg[k] + ddzg[k][n2 - 1] * dzg[n1 - 1] + ddzg[k][n1 - 1] * dzg[n2 - 1]) / topo::bndlyr_param;
                dddX += dzg[k] * dzg[n1 - 1] * dzg[n2 - 1] / pow(topo::bndlyr_param, 2);
                dddX *= X / topo::bndlyr_param;

                dddS[k]  = dddX / pow(1.0 + X, 2);
                dddS[k] -= 2.0 * (ddX[n1][n2] * dX[k + 1] + ddX[n1][k + 1] * dX[n2] + ddX[k + 1][n2] * dX[n1]) / pow(1.0 + X, 3);
                dddS[k] += 6.0 * dX[n1] * dX[n2] * dX[k + 1] / pow(1.0 + X, 4);
                dddS[k] *= -2.0;
            }

            ddw0 = dddS[0] * v0 + dddS[1] * u0 / cos(lat);
            if ((n1 == 1) && (n2 == 1)){
                ddw0 += 2.0 * sin(lat) / pow(cos(lat), 2) * ddS[1][2] * u0;
                ddw0 += (pow(sin(lat), 2) + 1.0) / pow(cos(lat), 3) * dS[2] * u0;
            } else if (((n1 == 1) && (n2 == 2)) || ((n1 == 2) && (n2 == 1))){            
                ddw0 += sin(lat) / pow(cos(lat), 2) * ddS[2][2] * u0;
            }
        } else if (n1 == 0 && n2 > 0){
            for(int k = 0; k < 2; k++){
                dddX  = ddzg[k][n2 - 1] + dzg[k] * dzg[n2 - 1] / topo::bndlyr_param;
                dddX *= -X / pow(topo::bndlyr_param, 2);

                dddS[k]  = dddX / pow(1.0 + X, 2);
                dddS[k] -= 2.0 * (ddX[n1][n2] * dX[k + 1] + ddX[n1][k + 1] * dX[n2] + ddX[k + 1][n2] * dX[n1]) / pow(1.0 + X, 3);
                dddS[k] += 6.0 * dX[n1] * dX[n2] * dX[k + 1] / pow(1.0 + X, 4);
                dddS[k] *= -2.0;
            }
            ddw0  = dddS[0] * v0 + ddS[n2][1] * dv0[0];
            ddw0 += (dddS[1] * u0 + ddS[n2][2] * du0[0]) / cos(lat);;
            if (n2 == 1){
                ddw0 += sin(lat) / pow(cos(lat), 2) * (ddS[0][2] * u0 + dS[2] * du0[0]);
            }
        } else if (n1 > 0 && n2 == 0){
            for(int k = 0; k < 2; k++){
                dddX  = ddzg[k][n1 - 1] + dzg[k] * dzg[n1 - 1] / topo::bndlyr_param;
                dddX *= -X / pow(topo::bndlyr_param, 2);

                dddS[k]  = dddX / pow(1.0 + X, 2);
                dddS[k] -= 2.0 * (ddX[n1][n2] * dX[k + 1] + ddX[n1][k + 1] * dX[n2] + ddX[k + 1][n2] * dX[n1]) / pow(1.0 + X, 3);
                dddS[k] += 6.0 * dX[n1] * dX[n2] * dX[k + 1] / pow(1.0 + X, 4);
                dddS[k] *= -2.0;
            }
            ddw0  = dddS[0] * v0 + ddS[n1][1] * dv0[0];
            ddw0 += (dddS[1] * u0 + ddS[n1][2] * du0[0]) / cos(lat);;
            if (n1 == 1){
                ddw0 += sin(lat) / pow(cos(lat), 2) * (ddS[0][2] * u0 + dS[2] * du0[0]);
            }

        } else if (n1 == 0 && n2 == 0){
            for(int k = 0; k < 2; k++){
                dddX = X / pow(topo::bndlyr_param, 3) * dzg[k];

                dddS[k]  = dddX / pow(1.0 + X, 2);
                dddS[k] -= 2.0 * (ddX[n1][n2] * dX[k + 1] + ddX[n1][k + 1] * dX[n2] + ddX[k + 1][n2] * dX[n1]) / pow(1.0 + X, 3);
                dddS[k] += 6.0 * dX[n1] * dX[n2] * dX[k + 1] / pow(1.0 + X, 4);
                dddS[k] *= -2.0;
            }

            ddw0  = dddS[0] * v0 + 2.0 * ddS[0][1] * dv0[0] + dS[1] * interp::eval_ddf(lat_eval, lon_eval, z_eval, 2, 2, v_spline);
            ddw0 += (dddS[1] * u0 + 2.0 * ddS[0][2] * du0[0] + dS[2] * interp::eval_ddf(lat_eval, lon_eval, z_eval, 2, 2, u_spline)) / cos(lat);
        }
        ddw0 *= -2.0 * topo::bndlyr_param / (globe::r0 + zg);

        ddN  = ddS[n1][n2] * X * w0 + S * ddX[n1][n2] * w0 + S * X * ddw0;
        ddN += dS[n1] * dX[n2] * w0 + dS[n1] * X * dw0[1] + S * dX[n1] * dw0[1];
        ddN += dS[n2] * dX[n1] * w0 + dS[n2] * X * dw0[0] + S * dX[n2] * dw0[0];
        ddD = -sg_ratio * dX[n1] * dX[n2] + (1.0 - sg_ratio * X) * ddX[n1][n2];

        result = ((ddN * D - dN[0] * dD[1] - dN[1] * dD[0] - N * ddD) * pow(D, 2) + 2.0 * N * D * dD[0] * dD[1]) / pow(D, 4);
    }
    return result;
}

//--------------------------------------------//
//-------------Functions To Define------------//
//------------Winds Simultaneously------------//
//--------------------------------------------//
void atmo::calc_uvw(double r, double lat, double lon, double & u, double & v, double & w, double du[3], double dv[3], double dw[3]){
    double lat_eval, lon_eval, z_eval, zg, dzg[2], ddzg[2][2], ddzg_temp[3];
    double u0, v0, du0[3], dv0[3], du0_temp[3], dv0_temp[3], ddu0_temp[6], ddv0_temp[6];
    double X, dX[3], ddX[3][3], sc0, dsc0[3], ddsc0[3][3];
    double w0, sg_ratio, N, D, dw0, dN, dD;
    
    // Compute z_eval, zg(x, y), the zg(x, y) derivatives
    lat_eval = interp::in_interval(lat, v_spline.x_vals[0], v_spline.x_vals[v_spline.length_x - 1]);
    lon_eval = interp::in_interval(lon, v_spline.y_vals[0], v_spline.y_vals[v_spline.length_y - 1]);
    z_eval = interp::in_interval(r - globe::r0, v_spline.z_vals[0], v_spline.z_vals[v_spline.length_z - 1]);

    if (geoac::is_topo && (r - globe::r0 < topo::z_bndlyr)){
        interp::eval_all(lat, lon, topo::spline, zg, dzg, ddzg_temp);
        zg += globe::r0;
        for(int n = 0; n < 2; n++){
            for(int m = 0; m < 2; m++){
                ddzg[n][m] = ddzg_temp[n + m];
            }
        }

        for(int n = 0; n < 2; n++){
            dzg[n] = topo::dz(lat, lon, n);
            for(int m = n; m < 2; m++){
                ddzg[n][m] = topo::ddz(lat, lon, n, m);
                ddzg[m][n] = ddzg[n][m];
            }
        }
    } else {
        zg = topo::z(lat, lon);
        for(int n = 0; n < 2; n++){
            dzg[n] = 0.0;
            for(int m = 0; m < 2; m++){
                ddzg[n][m] = 0.0;
            }
        }
    }
    
    // Evaluate the u and v splines and derivatives
    interp::eval_all(lat_eval, lon_eval, z_eval, atmo::u_spline, u0, du0_temp);
    interp::eval_all(lat_eval, lon_eval, z_eval, atmo::v_spline, v0, dv0_temp);
    
    du0[1] = du0_temp[0];    dv0[1] = dv0_temp[0];
    du0[2] = du0_temp[1];    dv0[2] = dv0_temp[1];
    du0[0] = du0_temp[2];    dv0[0] = dv0_temp[2];

    // Compute the boundary layer (BL) scaling
    X = exp(-(r - zg) / topo::bndlyr_param);
    for(int n = 0; n < 3; n++){
        if(n > 0){  dX[n] = X / topo::bndlyr_param * dzg[n - 1];}
        else {      dX[n] = -X / topo::bndlyr_param;}
        for(int m = 0; m < 3; m++){
            if((n > 0) && (m > 0)){         ddX[n][m] = X / topo::bndlyr_param * (ddzg[n - 1][m - 1] + 1.0 / topo::bndlyr_param * dzg[n - 1] * dzg[m - 1]);}
            else if((n == 0) && (m > 0)){   ddX[n][m] = - X / pow(topo::bndlyr_param, 2) * dzg[m - 1];}
            else if((n > 0) && (m == 0)){   ddX[n][m] = - X / pow(topo::bndlyr_param, 2) * dzg[n - 1];}
            else {                          ddX[n][m] = X / pow(topo::bndlyr_param, 2);}
        }
    }
    
    sc0 = (1.0 - X) / (1.0 + X);
    for(int n = 0; n < 3; n++){
        dsc0[n] = -2.0 * dX[n] / pow(1.0 + X, 2);
        for(int m = 0; m < 3; m++){
            ddsc0[n][m] = -2.0 * (ddX[n][m] / pow(1.0 + X, 2) - 2.0 * dX[n] * dX[m] / pow(1.0 + X, 3));
        }
    }
    
    // Combine u0, v0 with the BL scaling
    u = u0 * sc0;
    v = v0 * sc0;
    for(int n = 0; n < 3; n++){
        du[n] = du0[n] * sc0 + u0 * dsc0[n];
        dv[n] = dv0[n] * sc0 + v0 * dsc0[n];
    }

    // If using topography and within the BL, compute vertical winds
    if (geoac::is_topo && (r - globe::r0 < topo::z_bndlyr) && topo::use_BLw){
        w0  = v0 * dsc0[1] + u0 * dsc0[2] / cos(lat);
        w0 *= -2.0 * topo::bndlyr_param / (globe::r0 + zg);

        sg_ratio = topo::bndlyr_param / topo::vert_wind_grad;
        
        N = sc0 * X * w0;
        D = X + (sg_ratio / 2.0) * (1.0 - pow(X, 2));
        w = N / D;
        
        for(int n = 0; n < 3; n++){
            if (n == 0) {
                dw0  = dv0[0] * dsc0[1] + v0 * ddsc0[1][0];
                dw0 += (du0[0] * dsc0[2] + u0 * ddsc0[2][0]) / cos(lat);
            } else if (n == 1){
                dw0  = v0 * ddsc0[1][1];
                dw0 += u0 / cos(lat) * (ddsc0[1][2]  + dsc0[2] * tan(lat)) ;
            } else {
                dw0  = v0 * ddsc0[1][2];
                dw0 += u0 * ddsc0[2][2] / cos(lat);
            }
            dw0 *= -2.0 * topo::bndlyr_param / (globe::r0 + zg);

            dN = (dsc0[n] * X + sc0 * dX[n]) * w0 + sc0 * X * dw0;
            dD = (1.0 - sg_ratio * X) * dX[n];
            
            dw[n] = (dN * D - N * dD) / pow(D, 2);
        }
    } else {
        w = 0.0;
        for(int n = 0; n < 3; n++){ dw[n] = 0.0;}
    }
}

void atmo::calc_uvw(double r, double lat, double lon, double & u, double & v, double & w, double du[3], double dv[3], double dw[3], double ddu[6], double ddv[6], double ddw[6]){
    double lat_eval, lon_eval, z_eval, zg, dzg[2], ddzg[2][2], dddzg[2][2][2], ddzg_temp [3], dddzg_temp[4];
    double u0, v0, du0[3], dv0[3], ddu0[6], ddv0[6];
    double du0_temp[3], dv0_temp[3], ddu0_temp[6], ddv0_temp[6];
    double X, dX[3], ddX[3][3], dddX[2][3][3], sc0, dsc0[3], ddsc0[3][3], dddsc0[2][3][3];
    double w0, sg_ratio, N, D, dw0[3], dN[3], dD[3], ddw0, ddN, ddD;
    
    // Compute z_eval, zg(x, y), the zg(x, y) derivatives
    lat_eval = interp::in_interval(lat, v_spline.x_vals[0], v_spline.x_vals[v_spline.length_x - 1]);
    lon_eval = interp::in_interval(lon, v_spline.y_vals[0], v_spline.y_vals[v_spline.length_y - 1]);
    z_eval = interp::in_interval(r - globe::r0, v_spline.z_vals[0], v_spline.z_vals[v_spline.length_z - 1]);

    if (geoac::is_topo && (r - globe::r0 < topo::z_bndlyr)){
        interp::eval_all(lat, lon, topo::spline, zg, dzg, ddzg_temp, dddzg_temp);
        zg += globe::r0;
        for(int n = 0; n < 2; n++){
            for(int m = 0; m < 2; m++){
                ddzg[n][m] = ddzg_temp[n + m];
                for(int k = 0; k < 2; k++){
                    dddzg[n][m][k] = dddzg_temp[n + m + k];
                }
            }
        }
    } else {
        zg = topo::z(lat, lon);
        for(int n = 0; n < 2; n++){
            dzg[n] = 0.0;
            for(int m = 0; m < 2; m++){
                ddzg[n][m] = 0.0;
            }
        }
    }
    
    // Evaluate the u and v splines and derivatives
    interp::eval_all(lat_eval, lon_eval, z_eval, atmo::u_spline, u0, du0_temp, ddu0_temp);
    interp::eval_all(lat_eval, lon_eval, z_eval, atmo::v_spline, v0, dv0_temp, ddv0_temp);
    
    du0[1] = du0_temp[0];    dv0[1] = dv0_temp[0];
    du0[2] = du0_temp[1];    dv0[2] = dv0_temp[1];
    du0[0] = du0_temp[2];    dv0[0] = dv0_temp[2];

    ddu0[1] = ddu0_temp[0];    ddv0[1] = ddv0_temp[0];
    ddu0[2] = ddu0_temp[1];    ddv0[2] = ddv0_temp[1];
    ddu0[0] = ddu0_temp[2];    ddv0[0] = ddv0_temp[2];

    ddu0[5] = ddu0_temp[3];    ddv0[5] = ddv0_temp[3];
    ddu0[3] = ddu0_temp[4];    ddv0[3] = ddv0_temp[4];
    ddu0[4] = ddu0_temp[5];    ddv0[4] = ddv0_temp[5];

    // Compute the boundary layer (BL) scaling
    X = exp(-(r - zg) / topo::bndlyr_param);
    for(int n = 1; n < 3; n++){
        dX[n] = X / topo::bndlyr_param * dzg[n - 1];
        for(int m = 0; m < 3; m++){
            if(m == 0){
                ddX[n][m] = - X / pow(topo::bndlyr_param, 2) * dzg[n - 1];
                ddX[m][n] = ddX[n][m];

                for(int k = 0; k < 2; k++){
                    dddX[k][n][m]  = ddzg[n - 1][k] + dzg[n - 1] * dzg[k] / topo::bndlyr_param;
                    dddX[k][n][m] *= -X / pow(topo::bndlyr_param, 2);
                    dddX[k][m][n] = dddX[k][n][m];
                }
            } else {
                ddX[n][m] = ddzg[n - 1][m - 1] + dzg[n - 1] * dzg[m - 1] / topo::bndlyr_param;
                ddX[n][m] *= X / topo::bndlyr_param;
                ddX[m][n] = ddX[n][m];

                for(int k = 0; k < 2; k++){
                    dddX[k][n][m]  = dddzg[k][n - 1][m - 1];
                    dddX[k][n][m] += (ddzg[n - 1][m - 1] * dzg[k] + ddzg[k][m - 1] * dzg[n - 1] + ddzg[k][n - 1] * dzg[m - 1]) / topo::bndlyr_param;
                    dddX[k][n][m] += dzg[k] * dzg[n - 1] * dzg[m - 1] / pow(topo::bndlyr_param, 2);
                    dddX[k][n][m] *= X / topo::bndlyr_param;
                    dddX[k][m][n] = dddX[k][n][m];
                }
            }
        }
    }
    dX[0] = -X / topo::bndlyr_param;
    ddX[0][0] = X / pow(topo::bndlyr_param, 2);
    dddX[0][0][0] = X / pow(topo::bndlyr_param, 3) * dzg[0];
    dddX[1][0][0] = X / pow(topo::bndlyr_param, 3) * dzg[1];
    
    sc0 = (1.0 - X) / (1.0 + X);
    for(int n = 0; n < 3; n++){
        dsc0[n] = -2.0 * dX[n] / pow(1.0 + X, 2);
        for(int m = 0; m < 3; m++){
            ddsc0[n][m] = -2.0 * (ddX[n][m] / pow(1.0 + X, 2) - 2.0 * dX[n] * dX[m] / pow(1.0 + X, 3));
            for(int k = 0; k < 2; k++){
                dddsc0[k][n][m]  = dddX[k][n][m] / pow(1.0 + X, 2);
                dddsc0[k][n][m] -= 2.0 * (ddX[n][m] * dX[k + 1] + ddX[n][k + 1] * dX[m] + ddX[k + 1][m] * dX[n]) / pow(1.0 + X, 3);
                dddsc0[k][n][m] += 6.0 * dX[n] * dX[m] * dX[k + 1] / pow(1.0 + X, 4);
                dddsc0[k][n][m] *= -2.0;
            }
        }
    }
    
    // Combine u0, v0 with the BL scaling
    u = u0 * sc0;
    v = v0 * sc0;
    for(int n = 0; n < 3; n++){
        du[n] = du0[n] * sc0 + u0 * dsc0[n];
        dv[n] = dv0[n] * sc0 + v0 * dsc0[n];
        
        ddu[n] = ddu0[n] * sc0 + 2.0 * du0[n] * dsc0[n] + u0 * ddsc0[n][n];
        ddv[n] = ddv0[n] * sc0 + 2.0 * dv0[n] * dsc0[n] + v0 * ddsc0[n][n];
    }
    ddu[3] = ddu0[3] * sc0 + du0[0] * dsc0[1] + du0[1] * dsc0[0] + u0 * ddsc0[0][1];
    ddu[4] = ddu0[4] * sc0 + du0[0] * dsc0[2] + du0[2] * dsc0[0] + u0 * ddsc0[0][2];
    ddu[5] = ddu0[5] * sc0 + du0[1] * dsc0[2] + du0[2] * dsc0[1] + u0 * ddsc0[1][2];
    
    ddv[3] = ddv0[3] * sc0 + dv0[0] * dsc0[1] + dv0[1] * dsc0[0] + v0 * ddsc0[0][1];
    ddv[4] = ddv0[4] * sc0 + dv0[0] * dsc0[2] + dv0[2] * dsc0[0] + v0 * ddsc0[0][2];
    ddv[5] = ddv0[5] * sc0 + dv0[1] * dsc0[2] + dv0[2] * dsc0[1] + v0 * ddsc0[1][2];
    
    // If using topography and within the BL, compute vertical winds
    if (geoac::is_topo && (r - globe::r0 < topo::z_bndlyr) && topo::use_BLw){
        w0  = v0 * dsc0[1] + u0 * dsc0[2] / cos(lat);
        w0 *= -2.0 * topo::bndlyr_param / (globe::r0 + zg);

        sg_ratio = topo::bndlyr_param / topo::vert_wind_grad;
        
        N = sc0 * X * w0;
        D = X + (sg_ratio / 2.0) * (1.0 - pow(X, 2));
        w = N / D;
        
        for(int n = 0; n < 3; n++){
            dw0[n] = v0 * ddsc0[1][n] + u0 * ddsc0[2][n] / cos(lat);
            if(n == 0){
                dw0[n] += dv0[0] * dsc0[1] + du0[0] * dsc0[2] / cos(lat);
            } else if (n == 1){
                dw0[n] += u0 * dsc0[2] * sin(lat) / pow(cos(lat), 2);
            }
            dw0[n] *= -2.0 * topo::bndlyr_param / (globe::r0 + zg);
            
            dN[n] = (dsc0[n] * X + sc0 * dX[n]) * w0 + sc0 * X * dw0[n];
            dD[n] = (1.0 - sg_ratio * X) * dX[n];
            dw[n] = (dN[n] * D - N * dD[n]) / pow(D, 2);
            
            ddw0  = dddsc0[0][n][n] * v0 + dddsc0[1][n][n] * u0 / cos(lat);
            if(n == 0){
                ddw0 += 2.0 * ddsc0[1][0] * dv0[0] + dsc0[1] * ddv0[0];
                ddw0 += (2.0 * ddsc0[2][0] * du0[0] + dsc0[2] * ddu0[0]) / cos(lat);
            } else if(n == 1){
                ddw0 += 2.0 * sin(lat) / pow(cos(lat), 2) * ddsc0[1][2] * u0;
                ddw0 += (pow(sin(lat), 2) + 1.0) / pow(cos(lat), 3) * dsc0[2] * u0;
            }
            ddw0 *= -2.0 * topo::bndlyr_param / (globe::r0 + zg);

            ddN  = ddsc0[n][n] * X * w0 + sc0 * ddX[n][n] * w0 + sc0 * X * ddw0;
            ddN += 2.0 * (dsc0[n] * dX[n] * w0 + dsc0[n] * X * dw0[n] + sc0 * dX[n] * dw0[n]);
            ddD = -sg_ratio * pow(dX[n], 2) + (1.0 - sg_ratio * X) * ddX[n][n];
            ddw[n] = ((ddN * D - 2.0 * dN[n] * dD[n] - N * ddD) * pow(D, 2) + 2.0 * N * D * pow(dD[n], 2)) / pow(D, 4);
        }
        
        // ddw[3-5] are the drdt, drdp, and dtdp derivatives
        ddw0  = dddsc0[0][0][1] * v0 + ddsc0[1][1] * dv0[0] + (dddsc0[1][0][1] * u0 + ddsc0[1][2] * du0[0]) / cos(lat);
        ddw0 += sin(lat) / pow(cos(lat), 2) * (ddsc0[0][2] * u0 + dsc0[2] * du0[0]);    
        ddw0 *= -2.0 * topo::bndlyr_param / (globe::r0 + zg);

        ddN  = ddsc0[0][1] * X * w0 + sc0 * ddX[0][1] * w0 + sc0 * X * ddw0;
        ddN += dsc0[0] * dX[1] * w0 + dsc0[0] * X * dw0[1] + sc0 * dX[0] * dw0[1];
        ddN += dsc0[1] * dX[0] * w0 + dsc0[1] * X * dw0[0] + sc0 * dX[1] * dw0[0];       
        ddD = - sg_ratio * dX[0] * dX[1] + (1.0 - sg_ratio * X) * ddX[0][1];
        ddw[3] = ((ddN * D - dN[0] * dD[1] - dN[1] * dD[0] - N * ddD) * pow(D, 2) + 2.0 * N * D * dD[0] * dD[1]) / pow(D, 4);
        
        ddw0  = dddsc0[0][0][2] * v0 + ddsc0[1][2] * dv0[0] + (dddsc0[1][0][2] * u0 + ddsc0[2][2] * du0[0]) / cos(lat);
        ddw0 *= -2.0 * topo::bndlyr_param / (globe::r0 + zg);

        ddN  = ddsc0[0][2] * X * w0 + sc0 * ddX[0][2] * w0 + sc0 * X * ddw0;
        ddN += dsc0[0] * dX[2] * w0 + dsc0[0] * X * dw0[2] + sc0 * dX[0] * dw0[2];
        ddN += dsc0[2] * dX[0] * w0 + dsc0[2] * X * dw0[0] + sc0 * dX[2] * dw0[0];
        ddD = - sg_ratio * dX[0] * dX[2] + (1.0 - sg_ratio * X) * ddX[0][2];
        ddw[4] = ((ddN * D - dN[0] * dD[2] - dN[2] * dD[0] - N * ddD) * pow(D, 2) + 2.0 * N * D * dD[0] * dD[2]) / pow(D, 4);

        ddw0 = dddsc0[0][1][2] * v0 + dddsc0[1][1][2] * u0 / cos(lat);
        ddw0 += sin(lat) / pow(cos(lat), 2) * ddsc0[2][2] * u0;
        ddw0 *= -2.0 * topo::bndlyr_param / (globe::r0 + zg);

        ddN  = ddsc0[1][2] * X * w0 + sc0 * ddX[1][2] * w0 + sc0 * X * ddw0;
        ddN += dsc0[1] * dX[2] * w0 + dsc0[1] * X * dw0[2] + sc0 * dX[1] * dw0[2];
        ddN += dsc0[2] * dX[1] * w0 + dsc0[2] * X * dw0[1] + sc0 * dX[2] * dw0[1];
        ddD = - sg_ratio * dX[1] * dX[2] + (1.0 - sg_ratio * X) * ddX[1][2];
        ddw[5] = ((ddN * D - dN[1] * dD[2] - dN[2] * dD[1] - N * ddD) * pow(D, 2) + 2.0 * N * D * dD[1] * dD[2]) / pow(D, 4);
    } else {
        w = 0.0;
        for(int n = 0; n < 3; n++){ dw[n] = 0.0;}
        for(int n = 0; n < 6; n++){ ddw[n] = 0.0;}
    }
}





#endif /* _ATMO_STATE_SPH_RNGDEP_CPP_ */
