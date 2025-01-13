# ifndef _ATMO_STATE_3D_RNGDEP_CPP_
# define _ATMO_STATE_3D_RNGDEP_CPP_

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <sstream>
#include <math.h>

#include "atmo_state.h"
#include "atmo_io.3d.rngdep.h"
#include "../util/interpolation.h"
#include "../util/fileIO.h"
#include "../geoac/geoac.params.h"

using namespace std;

//----------------------------//
//-----Physical Constants-----//
//----------------------------//
double atmo::gam = 1.4;
double atmo::R = 287.058;

//---------------------------------------------------//
//---------Functions Defining the Topography---------//
//------------For 3D Cartesian Propgation------------//
//---------------------------------------------------//
double topo::z0 = 0.0;
double topo::z_max;
double topo::z_bndlyr;
bool topo::use_BLw = false;

void topo::set_bndlyr(){
    if (geoac::is_topo){
        z_max = 0.0;
        for (int nx = 0; nx < topo::spline.length_x; nx++){ for (int ny = 0; ny < topo::spline.length_y; ny++){ z_max = max(z_max, topo::spline.f_vals[nx][ny]);}}
        z_bndlyr = z_max + 2.0; z_max += 0.01;
    } else {
        z_max = z0; z_bndlyr = z_max + 2.0;
    }
}

double topo::z(double x, double y){
    double result = 0.0;
    if (geoac::is_topo){
        return interp::eval_f(interp::in_interval(x, spline.x_vals[0], spline.x_vals[spline.length_x - 1]),
                              interp::in_interval(y, spline.y_vals[0], spline.y_vals[spline.length_y - 1]), spline);
    }
    return result;
}

double topo::dz(double x, double y, int n){
    double result = 0.0;
    if (geoac::is_topo){
        result = interp::eval_df(interp::in_interval(x, spline.x_vals[0], spline.x_vals[spline.length_x - 1]),
                                 interp::in_interval(y, spline.y_vals[0], spline.y_vals[spline.length_y - 1]), n, spline);
    }
    return result;
}

double topo::ddz(double x, double y, int n1, int n2){
    double result = 0.0;
    if (geoac::is_topo){
        result = interp::eval_ddf(interp::in_interval(x, spline.x_vals[0], spline.x_vals[spline.length_x - 1]),
                                  interp::in_interval(y, spline.y_vals[0], spline.y_vals[spline.length_y - 1]), n1, n2, spline);
    }
    return result;
}

double topo::dddz(double x, double y, int n1, int n2, int n3){
    double result = 0.0;
    if (geoac::is_topo){
        return interp::eval_dddf(interp::in_interval(x, spline.x_vals[0], spline.x_vals[spline.length_x - 1]),
                                 interp::in_interval(y, spline.y_vals[0], spline.y_vals[spline.length_y - 1]), n1, n2, n3, spline);
    }
    return result;
}

double topo::bndlyr_param = 0.1;
double topo::vert_wind_grad = 1.0;
double topo::bndlyr_sc(double x_val, double y_val, double z_val){
    double X, result = 1.0;
    if (z_val < z_bndlyr){
        X = exp(-(z_val - z(x_val, y_val)) / bndlyr_param);
        result = (1.0 - X) / (1.0 + X);
    }
    return result;
}
double topo::bndlyr_dsc(double x_val, double y_val, double z_val, int n){
    double X, dX, result = 0.0;
    if (z_val < z_bndlyr){
        X = exp(-(z_val - z(x_val, y_val)) / bndlyr_param);
        if(n < 2){ 
            dX = X / bndlyr_param * dz(x_val, y_val, n);
        } else{
            dX = - X / bndlyr_param;
        }
        result = -2.0 * dX / pow(1.0 + X, 2);
    }
    return result;
}
double topo::bndlyr_ddsc(double x_val, double y_val, double z_val, int n1, int n2){
    double X, dX1, dX2, ddX, result = 0.0;
    if (z_val < z_bndlyr){
        X = exp(-(z_val - z(x_val, y_val)) / bndlyr_param);
        if(n1 < 2){
            dX1 = X / bndlyr_param * dz(x_val, y_val, n1);
        } else {
            dX1 = - X / bndlyr_param;
        }

        if(n2 < 2){
            dX2 = X / bndlyr_param * dz(x_val, y_val, n2);
        } else {
            dX2 = - X / bndlyr_param;
        }
    
        if(n1 < 2 && n2 < 2){
            ddX =  X / bndlyr_param * (ddz(x_val, y_val, n1, n2) + dz(x_val, y_val, n1) * dz(x_val, y_val, n2) / bndlyr_param);
        } else if(n1 < 2 && n2 == 2){ 
            ddX = -X / pow(bndlyr_param, 2) * dz(x_val, y_val, n1);
        } else if(n1 == 2 && n2 < 2){ 
            ddX = -X / pow(bndlyr_param, 2) * dz(x_val, y_val, n2);
        } else {
            ddX =  X / pow(bndlyr_param, 2);
        }
    
        result = -2.0 * (ddX / pow(1.0 + X, 2) - 2.0 * dX1 * dX2 / pow(1.0 + X, 3.0));
    }
    return result;
}


//---------------------------------------------------//
//---------Functions Defining the Atmosphere---------//
//------------For 3D Cartesian Propgation------------//
//---------------------------------------------------//
double atmo::rho(double x, double y, double z){
    double x_eval = interp::in_interval(x, rho_spline.x_vals[0], rho_spline.x_vals[rho_spline.length_x - 1]);
    double y_eval = interp::in_interval(y, rho_spline.y_vals[0], rho_spline.y_vals[rho_spline.length_y - 1]);
    double z_eval = interp::in_interval(z, rho_spline.z_vals[0], rho_spline.z_vals[rho_spline.length_z - 1]);
    
    return interp::eval_f(x_eval, y_eval, z_eval, rho_spline);
}

double atmo::c(double x, double y, double z){
    double x_eval = interp::in_interval(x, c_spline.x_vals[0], c_spline.x_vals[c_spline.length_x - 1]);
    double y_eval = interp::in_interval(y, c_spline.y_vals[0], c_spline.y_vals[c_spline.length_y - 1]);
    double z_eval = interp::in_interval(z, c_spline.z_vals[0], c_spline.z_vals[c_spline.length_z - 1]);
    
    return interp::eval_f(x_eval, y_eval, z_eval, c_spline);
}

double atmo::dc(double x, double y, double z, int n){
    double x_eval = interp::in_interval(x, c_spline.x_vals[0], c_spline.x_vals[c_spline.length_x - 1]);
    double y_eval = interp::in_interval(y, c_spline.y_vals[0], c_spline.y_vals[c_spline.length_y - 1]);
    double z_eval = interp::in_interval(z, c_spline.z_vals[0], c_spline.z_vals[c_spline.length_z - 1]);
    
    return interp::eval_df(x_eval, y_eval, z_eval, n, c_spline);
}

double atmo::ddc(double x, double y, double z, int n1, int n2){
    double x_eval = interp::in_interval(x, c_spline.x_vals[0], c_spline.x_vals[c_spline.length_x - 1]);
    double y_eval = interp::in_interval(y, c_spline.y_vals[0], c_spline.y_vals[c_spline.length_y - 1]);
    double z_eval = interp::in_interval(z, c_spline.z_vals[0], c_spline.z_vals[c_spline.length_z - 1]);
    
    return interp::eval_ddf(x_eval, y_eval, z_eval, n1, n2, c_spline);
}

double atmo::u(double x, double y, double z){
    double x_eval = interp::in_interval(x, u_spline.x_vals[0], u_spline.x_vals[u_spline.length_x - 1]);
    double y_eval = interp::in_interval(y, u_spline.y_vals[0], u_spline.y_vals[u_spline.length_y - 1]);
    double z_eval = interp::in_interval(z, u_spline.z_vals[0], u_spline.z_vals[u_spline.length_z - 1]);
    
    return interp::eval_f(x_eval, y_eval, z_eval, u_spline) * topo::bndlyr_sc(x, y, z);
}
double atmo::du(double x, double y, double z, int n){
    double x_eval = interp::in_interval(x, u_spline.x_vals[0], u_spline.x_vals[u_spline.length_x - 1]);
    double y_eval = interp::in_interval(y, u_spline.y_vals[0], u_spline.y_vals[u_spline.length_y - 1]);
    double z_eval = interp::in_interval(z, u_spline.z_vals[0], u_spline.z_vals[u_spline.length_z - 1]);
    
    double result;
    result  = interp::eval_df(x_eval, y_eval, z_eval, n, u_spline) * topo::bndlyr_sc(x, y, z);
    result += interp::eval_f(x_eval, y_eval, z_eval, u_spline) * topo::bndlyr_dsc(x, y, z, n);
    return result;
}
double atmo::ddu(double x, double y, double z, int n1, int n2){
    double x_eval = interp::in_interval(x, u_spline.x_vals[0], u_spline.x_vals[u_spline.length_x - 1]);
    double y_eval = interp::in_interval(y, u_spline.y_vals[0], u_spline.y_vals[u_spline.length_y - 1]);
    double z_eval = interp::in_interval(z, u_spline.z_vals[0], u_spline.z_vals[u_spline.length_z - 1]);

    double result;
    result  = interp::eval_ddf(x_eval, y_eval, z_eval, n1, n2, u_spline) * topo::bndlyr_sc(x, y, z);
    result += interp::eval_f(x_eval, y_eval, z_eval, u_spline) * topo::bndlyr_ddsc(x, y, z, n1, n2);
    if(n1 == n2){
        result += 2.0 * interp::eval_df(x_eval, y_eval, z_eval, n1, u_spline) * topo::bndlyr_dsc(x, y, z, n1);
    } else {
        result += interp::eval_df(x_eval, y_eval, z_eval, n1, u_spline) * topo::bndlyr_dsc(x, y, z, n2);
        result += interp::eval_df(x_eval, y_eval, z_eval, n2, u_spline) * topo::bndlyr_dsc(x, y, z, n1);
    }
    return result;
}

double atmo::v(double x, double y, double z){
    double x_eval = interp::in_interval(x, v_spline.x_vals[0], v_spline.x_vals[v_spline.length_x - 1]);
    double y_eval = interp::in_interval(y, v_spline.y_vals[0], v_spline.y_vals[v_spline.length_y - 1]);
    double z_eval = interp::in_interval(z, v_spline.z_vals[0], v_spline.z_vals[v_spline.length_z - 1]);
   
    return interp::eval_f(x_eval, y_eval, z_eval, v_spline) * topo::bndlyr_sc(x, y, z);
}
double atmo::dv(double x, double y, double z, int n){
    double x_eval = interp::in_interval(x, v_spline.x_vals[0], v_spline.x_vals[v_spline.length_x - 1]);
    double y_eval = interp::in_interval(y, v_spline.y_vals[0], v_spline.y_vals[v_spline.length_y - 1]);
    double z_eval = interp::in_interval(z, v_spline.z_vals[0], v_spline.z_vals[v_spline.length_z - 1]);
    
    double result;
    result  = interp::eval_df(x_eval, y_eval, z_eval, n, v_spline) * topo::bndlyr_sc(x, y, z);
    result += interp::eval_f(x_eval, y_eval, z_eval, v_spline) * topo::bndlyr_dsc(x, y, z, n);
    return result;
}
double atmo::ddv(double x, double y, double z, int n1, int n2){
    double x_eval = interp::in_interval(x, v_spline.x_vals[0], v_spline.x_vals[v_spline.length_x - 1]);
    double y_eval = interp::in_interval(y, v_spline.y_vals[0], v_spline.y_vals[v_spline.length_y - 1]);
    double z_eval = interp::in_interval(z, v_spline.z_vals[0], v_spline.z_vals[v_spline.length_z - 1]);
    
    double result;
    result  = interp::eval_ddf(x_eval, y_eval, z_eval, n1, n2, v_spline) * topo::bndlyr_sc(x, y, z);
    result += interp::eval_f(x_eval, y_eval, z_eval, v_spline) * topo::bndlyr_ddsc(x, y, z, n1, n2);
    if(n1 == n2){
        result += 2.0 * interp::eval_df(x_eval, y_eval, z_eval, n1, v_spline) * topo::bndlyr_dsc(x, y, z, n1);
    } else {
        result += interp::eval_df(x_eval, y_eval, z_eval, n1, v_spline) * topo::bndlyr_dsc(x, y, z, n2);
        result += interp::eval_df(x_eval, y_eval, z_eval, n2, v_spline) * topo::bndlyr_dsc(x, y, z, n1);
    }
    return result;
}

// NOTE: dw and ddw direct functions aren't returning the correct values for the BL winds
// when topography is included

double atmo::w(double x, double y, double z){
    double x_eval, y_eval, z_eval, X, dX[3], ddX[2], S, dS[3], ddS[2], N, D, sg_ratio, result = 0.0;
    double zg, dzg[2], ddzg[2], u0, v0, w0, dw0, dN, dD;
    
    if (geoac::is_topo && (z < topo::z_bndlyr) && topo::use_BLw){
        x_eval = interp::in_interval(x, u_spline.x_vals[0], u_spline.x_vals[u_spline.length_x - 1]);
        y_eval = interp::in_interval(y, u_spline.y_vals[0], u_spline.y_vals[u_spline.length_y - 1]);
        z_eval = interp::in_interval(z, u_spline.z_vals[0], u_spline.z_vals[u_spline.length_z - 1]);
        sg_ratio = topo::bndlyr_param / topo::vert_wind_grad;

        X = exp(-(z - topo::z(x, y)) / topo::bndlyr_param);
        S = (1.0 - X) / (1.0 + X);
        
        dS[0] = -2.0 / topo::bndlyr_param * (X / pow(1.0 + X, 2)) * topo::dz(x, y, 0);
        dS[1] = -2.0 / topo::bndlyr_param * (X / pow(1.0 + X, 2)) * topo::dz(x, y, 1);
        
        w0  = interp::eval_f(x_eval, y_eval, z_eval, u_spline) * dS[0];
        w0 += interp::eval_f(x_eval, y_eval, z_eval, v_spline) * dS[1];
        w0 *= -2.0 * topo::bndlyr_param;        
        
        result  = S * X * w0;
        result /= (X + (sg_ratio / 2.0) * (1.0 - pow(X, 2)));
    }
    return result;
}

double atmo::dw(double x, double y, double z, int n){
    double x_eval, y_eval, z_eval, X, dX[3], ddX[2], S, dS[3], ddS[2], N, D, sg_ratio, result = 0.0;
    double zg, dzg[2], ddzg[2], u0, v0, w0, dw0, dN, dD;
    
    if (geoac::is_topo && (z < topo::z_bndlyr) && topo::use_BLw){
        x_eval = interp::in_interval(x, u_spline.x_vals[0], u_spline.x_vals[u_spline.length_x - 1]);
        y_eval = interp::in_interval(y, u_spline.y_vals[0], u_spline.y_vals[u_spline.length_y - 1]);
        z_eval = interp::in_interval(z, u_spline.z_vals[0], u_spline.z_vals[u_spline.length_z - 1]);
        sg_ratio = topo::bndlyr_param / topo::vert_wind_grad;

        // Compute the topogrpahy and derivatives       
        zg = topo::z(x, y);
        dzg[0] = topo::dz(x, y, 0);
        dzg[1] = topo::dz(x, y, 1);

        // Evaluate the u and v splines
        u0 = interp::eval_f(x_eval, y_eval, z_eval, u_spline);
        v0 = interp::eval_f(x_eval, y_eval, z_eval, v_spline);
        
        // Compute the boundary layer (BL) scaling        
        X = exp(-(z - zg) / topo::bndlyr_param);
        for(int j = 0; j < 2; j++){
            dX[j] = X / topo::bndlyr_param * dzg[j];
        }
        dX[2] = -X / topo::bndlyr_param;

        S = (1.0 - X) / (1.0 + X);
        for(int j = 0; j < 3; j++){
            dS[j] = -2.0 * dX[j] / pow(1.0 + X, 2);
        }

        // ddS is dxdn and dydn
        if(n < 2){
            ddX[0] = X / topo::bndlyr_param * (topo::ddz(x, y, 0, n) + 1.0 / topo::bndlyr_param * dzg[0] * dzg[n]);
            ddX[1] = X / topo::bndlyr_param * (topo::ddz(x, y, 1, n) + 1.0 / topo::bndlyr_param * dzg[1] * dzg[n]);           
        } else {
            ddX[0] = - X / pow(topo::bndlyr_param, 2) * dzg[0];
            ddX[1] = - X / pow(topo::bndlyr_param, 2) * dzg[1];
        }
        
        ddS[0] = -2.0 * (ddX[0] / pow(1.0 + X, 2) - 2.0 * dX[0] * dX[n] / pow(1.0 + X, 3));
        ddS[1] = -2.0 * (ddX[1] / pow(1.0 + X, 2) - 2.0 * dX[1] * dX[n] / pow(1.0 + X, 3));

        // Calculate dw
        w0  = -2.0 * topo::bndlyr_param * (u0 * dS[0] + v0 * dS[1]);
        N = S * X * w0;
        D = X + (sg_ratio / 2.0) * (1.0 - pow(X, 2));
        
        dw0  = u0 * ddS[0] + v0 * ddS[1];
        if(n == 2){
            dw0 += interp::eval_df(x_eval, y_eval, z_eval, 2, u_spline) * dS[0];
            dw0 += interp::eval_df(x_eval, y_eval, z_eval, 2, v_spline) * dS[1];
        }
        dw0 *= -2.0 * topo::bndlyr_param;

        dN = (dS[n] * X + S * dX[n]) * w0 + S * X * dw0;
        dD = (1.0 - sg_ratio * X) * dX[n];

        result = (dN * D - N * dD) / pow(D, 2);
    }
    
    return result;
}

double atmo::ddw(double x, double y, double z, int n1, int n2){
    double x_eval, y_eval, z_eval, X, S, N, D, sg_ratio, result = 0.0;
    double zg, dzg[2], u0, v0, w0, du0 = 0.0, dv0 = 0.0, dw0[2], dX[3], ddX[3][3], dS[3], ddS[3][3], dN[2], dD[2];
    double ddzg[2][2], dddX, dddS[2], ddN, ddD, ddw0;
    
    if (geoac::is_topo && (z < topo::z_bndlyr) && topo::use_BLw){
        x_eval = interp::in_interval(x, u_spline.x_vals[0], u_spline.x_vals[u_spline.length_x - 1]);
        y_eval = interp::in_interval(y, u_spline.y_vals[0], u_spline.y_vals[u_spline.length_y - 1]);
        z_eval = interp::in_interval(z, u_spline.z_vals[0], u_spline.z_vals[u_spline.length_z - 1]);
        sg_ratio = topo::bndlyr_param / topo::vert_wind_grad;
        
        // Compute the topogrpahy and derivatives       
        zg = topo::z(x, y);
        for(int j = 0; j < 2; j++){
            dzg[j] = topo::dz(x, y, j);
            for(int k = j; k < 2; k++){
                ddzg[j][k] = topo::ddz(x, y, j, k);
                ddzg[k][j] = ddzg[j][k];
            }
        }

        // Evaluate the u and v splines and derivatives
        u0 = interp::eval_f(x_eval, y_eval, z_eval, u_spline);
        v0 = interp::eval_f(x_eval, y_eval, z_eval, v_spline);
        if((n1 == 2) || (n2 == 2)){
            du0 = interp::eval_df(x_eval, y_eval, z_eval, 2, u_spline);
            dv0 = interp::eval_df(x_eval, y_eval, z_eval, 2, v_spline);
        } else {
            du0 = 0.0;
            dv0 = 0.0;
        }
        
        // Compute the boundary layer (BL) scaling        
        X = exp(-(z - zg) / topo::bndlyr_param);
        for(int n = 0; n < 2; n++){
            dX[n] = X / topo::bndlyr_param * dzg[n];
            for(int m = n; m < 3; m++){
                if(m < 2){  ddX[n][m] = X / topo::bndlyr_param * (ddzg[n][m] + 1.0 / topo::bndlyr_param * dzg[n] * dzg[m]);}
                else {      ddX[n][m] = - X / pow(topo::bndlyr_param, 2) * dzg[n];}
                ddX[m][n] = ddX[n][m];
            }
        }
        dX[2] = -X / topo::bndlyr_param;
        ddX[2][2] = X / pow(topo::bndlyr_param, 2);

        S = (1.0 - X) / (1.0 + X);
        for(int n = 0; n < 3; n++){
            dS[n] = -2.0 * dX[n] / pow(1.0 + X, 2);
            for(int m = n; m < 3; m++){
                ddS[n][m] = -2.0 * (ddX[n][m] / pow(1.0 + X, 2) - 2.0 * dX[n] * dX[m] / pow(1.0 + X, 3));
                ddS[m][n] = ddS[n][m];
            }
        }
        
        // Calculate ddw
        w0  = u0 * dS[0] + v0 * dS[1];
        w0 *= -2.0 * topo::bndlyr_param;

        N = S * X * w0;
        D = X + (sg_ratio / 2.0) * (1.0 - pow(X, 2));
        
        dw0[0]  = u0 * ddS[0][n1] + v0 * ddS[1][n1];
        if (n1 == 2){ dw0[0] += du0 * dS[0] + dv0 * dS[1];}
        dw0[0] *= -2.0 * topo::bndlyr_param;

        dN[0] = (dS[n1] * X + S * dX[n1]) * w0 + S * X * dw0[0];
        dD[0] = (1.0 - sg_ratio * X) * dX[n1];
        
        if(n2 == n1){
            dw0[1] = dw0[0];
            dN[1] = dN[0];
            dD[1] = dD[0];
        } else {
            dw0[1] = u0 * ddS[0][n2] + v0 * ddS[1][n2];
            if (n2 == 2){ dw0[1] += du0 * dS[0] + dv0 * dS[1];}
            dw0[1] *= -2.0 * topo::bndlyr_param;
        
            dN[1] = (dS[n2] * X + S * dX[n2]) * w0 + S * X * dw0[1];
            dD[1] = (1.0 - sg_ratio * X) * dX[n2];
        }

        if (n1 < 2 && n2 < 2){
            for(int k = 0; k < 2; k++){
                dddX  = topo::dddz(x, y, n1, n2, k);
                dddX += (ddzg[n1][n2] * dzg[k] + ddzg[k][n2] * dzg[n1] + ddzg[k][n1] * dzg[n2]) / topo::bndlyr_param;
                dddX += dzg[k] * dzg[n1] * dzg[n2] / pow(topo::bndlyr_param, 2);
                dddX *= X / topo::bndlyr_param;

                dddS[k]  = dddX / pow(1.0 + X, 2);
                dddS[k] -= 2.0 * (ddX[n1][n2] * dX[k] + ddX[n1][k] * dX[n2] + ddX[k][n2] * dX[n1]) / pow(1.0 + X, 3);
                dddS[k] += 6.0 * dX[n1] * dX[n2] * dX[k] / pow(1.0 + X, 4);
                dddS[k] *= -2.0;
            }
            ddw0 = dddS[0] * u0 + dddS[1] * v0;
        } else if (n1 < 2 && n2 == 2){
            for(int k = 0; k < 2; k++){
                dddX  = ddzg[k][n1] + dzg[k] * dzg[n1] / topo::bndlyr_param;
                dddX *= -X / pow(topo::bndlyr_param, 2);

                dddS[k]  = dddX / pow(1.0 + X, 2);
                dddS[k] -= 2.0 * (ddX[n1][n2] * dX[k] + ddX[n1][k] * dX[n2] + ddX[k][n2] * dX[n1]) / pow(1.0 + X, 3);
                dddS[k] += 6.0 * dX[n1] * dX[n2] * dX[k] / pow(1.0 + X, 4);
                dddS[k] *= -2.0;
            }
            ddw0  = dddS[0] * u0 + ddS[0][n1] * du0;
            ddw0 += dddS[1] * v0 + ddS[1][n1] * dv0;
        } else if (n1 == 2 && n2 < 2){
            for(int k = 0; k < 2; k++){
                dddX  = ddzg[k][n2] + dzg[k] * dzg[n2] / topo::bndlyr_param;
                dddX *= -X / pow(topo::bndlyr_param, 2);

                dddS[k]  = dddX / pow(1.0 + X, 2);
                dddS[k] -= 2.0 * (ddX[n1][n2] * dX[k] + ddX[n1][k] * dX[n2] + ddX[k][n2] * dX[n1]) / pow(1.0 + X, 3);
                dddS[k] += 6.0 * dX[n1] * dX[n2] * dX[k] / pow(1.0 + X, 4);
                dddS[k] *= -2.0;
            }
            ddw0  = dddS[0] * u0 + ddS[0][n2] * du0;
            ddw0 += dddS[1] * v0 + ddS[1][n2] * dv0;
        } else if (n1 == 2 && n2 == 2){
            for(int k = 0; k < 2; k++){
                dddX = X / pow(topo::bndlyr_param, 3) * dzg[k];

                dddS[k]  = dddX / pow(1.0 + X, 2);
                dddS[k] -= 2.0 * (ddX[n1][n2] * dX[k] + ddX[n1][k] * dX[n2] + ddX[k][n2] * dX[n1]) / pow(1.0 + X, 3);
                dddS[k] += 6.0 * dX[n1] * dX[n2] * dX[k] / pow(1.0 + X, 4);
                dddS[k] *= -2.0;
            }

            ddw0  = dddS[0] * u0 + 2.0 * ddS[0][2] * du0 + dS[0] * interp::eval_ddf(x_eval, y_eval, z_eval, 2, 2, u_spline);
            ddw0 += dddS[1] * v0 + 2.0 * ddS[1][2] * dv0 + dS[1] * interp::eval_ddf(x_eval, y_eval, z_eval, 2, 2, v_spline);
        }
        ddw0 *= -2.0 * topo::bndlyr_param;    
        
        ddN  = ddS[n1][n2] * X * w0 + S * ddX[n1][n2] * w0 + S * X * ddw0;
        ddN += dS[n1] * dX[n2] * w0 + dS[n1] * X * dw0[1] + S * dX[n2] * dw0[0];
        ddN += dS[n2] * dX[n1] * w0 + dS[n2] * X * dw0[0] + S * dX[n1] * dw0[1];

        ddD = - sg_ratio * dX[n1] * dX[n2] + (1.0 - sg_ratio * X) * ddX[n1][n2];        

        result = ((ddN * D - dN[0] * dD[1] - dN[1] * dD[0] - N * ddD) * pow(D, 2) + 2.0 * N * D * dD[0] * dD[1]) / pow(D, 4);
    }
    return result;
}

//--------------------------------------------//
//-------------Functions To Define------------//
//------------Winds Simultaneously------------//
//--------------------------------------------//
void atmo::calc_uvw(double x, double y, double z, double & u, double & v, double & w, double du[3], double dv[3], double dw[3]){
    double x_eval, y_eval, z_eval, zg, dzg[2], ddzg[2][2], ddzg_temp [2];
    double u0, v0, du0[3], dv0[3];
    double X, dX[3], ddX[3][3], sc0, dsc0[3], ddsc0[3][3];
    double w0, sg_ratio, N, D, dw0, dN, dD;
    
    // Compute z_eval, zg(x, y), the zg(x, y) derivatives
    x_eval = interp::in_interval(x, u_spline.x_vals[0], u_spline.x_vals[u_spline.length_x - 1]);
    y_eval = interp::in_interval(y, u_spline.y_vals[0], u_spline.y_vals[u_spline.length_y - 1]);
    z_eval = interp::in_interval(z, u_spline.z_vals[0], u_spline.z_vals[u_spline.length_z - 1]);
    if (geoac::is_topo && (z < topo::z_bndlyr)){
        interp::eval_all(x, y, topo::spline, zg, dzg, ddzg_temp);
        for(int n = 0; n < 2; n++){
            for(int m = 0; m < 2; m++){
                ddzg[n][m] = ddzg_temp[n + m];
            }
        }
    } else {
        zg = topo::z(x, y);
        for(int n = 0; n < 2; n++){
            dzg[n] = 0.0;
            for(int m = 0; m < 2; m++){
                ddzg[n][m] = 0.0;
            }
        }
    }
    
    // Evaluate the u and v splines and first order derivatives
    interp::eval_all(x_eval, y_eval, z_eval, atmo::u_spline, u0, du0);
    interp::eval_all(x_eval, y_eval, z_eval, atmo::v_spline, v0, dv0);
    
    // Compute the boundary layer (BL) scaling
    X = exp(-(z - zg) / topo::bndlyr_param);
    for(int n = 0; n < 2; n++){
        dX[n] = X / topo::bndlyr_param * dzg[n];
        for(int m = n; m < 3; m++){
            if(m < 2){  ddX[n][m] = X / topo::bndlyr_param * (ddzg[n][m] + 1.0 / topo::bndlyr_param * dzg[n] * dzg[m]);}
            else {      ddX[n][m] = - X / pow(topo::bndlyr_param, 2) * dzg[n];}
            ddX[m][n] = ddX[n][m];
        }
    }
    dX[2] = -X / topo::bndlyr_param;
    ddX[2][2] = X / pow(topo::bndlyr_param, 2);

    sc0 = (1.0 - X) / (1.0 + X);
    for(int n = 0; n < 3; n++){
        dsc0[n] = -2.0 * dX[n] / pow(1.0 + X, 2);
        for(int m = n; m < 3; m++){
            ddsc0[n][m] = -2.0 * (ddX[n][m] / pow(1.0 + X, 2) - 2.0 * dX[n] * dX[m] / pow(1.0 + X, 3));
            ddsc0[m][n] = ddsc0[n][m];
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
    if (geoac::is_topo && (z < topo::z_bndlyr) && topo::use_BLw){
        w0 = -2.0 * topo::bndlyr_param * (u0 * dsc0[0] + v0 * dsc0[1]);
        sg_ratio = topo::bndlyr_param / topo::vert_wind_grad;
        
        N = sc0 * X * w0;
        D = X + (sg_ratio / 2.0) * (1.0 - pow(X, 2));
        
        w = N / D;
        
        for(int n = 0; n < 3; n++){
            if(n < 2){  dw0 = -2.0 * topo::bndlyr_param * (u0 * ddsc0[0][n] + v0 * ddsc0[1][n]);}
            else {      dw0 = -2.0 * topo::bndlyr_param * (du0[2] * dsc0[0] + u0 * ddsc0[0][2] + dv0[2] * dsc0[1] + v0 * ddsc0[1][2]);}
            
            dN = (dsc0[n] * X + sc0 * dX[n]) * w0 + sc0 * X * dw0;
            dD = (1.0 - sg_ratio * X) * dX[n];
            
            dw[n] = (dN * D - N * dD) / pow(D, 2);
        }
    } else {
        w = 0.0;
        for(int n = 0; n < 3; n++){ dw[n] = 0.0;}
    }
}

void atmo::calc_uvw(double x, double y, double z, double & u, double & v, double & w, double du[3], double dv[3], double dw[3], double ddu[6], double ddv[6], double ddw[6]){
    double x_eval, y_eval, z_eval, zg, dzg[2], ddzg[2][2], dddzg[2][2][2], ddzg_temp [3], dddzg_temp [4];
    double u0, v0, du0[3], dv0[3], ddu0[6], ddv0[6];
    double X, dX[3], ddX[3][3], dddX[2][3][3], sc0, dsc0[3], ddsc0[3][3], dddsc0[2][3][3];
    double w0, sg_ratio, N, D, dw0[3], dN[3], dD[3], ddw0, ddN, ddD;
    
    // Compute x, y, and z_eval, zg(x, y), the zg(x, y) derivatives
    x_eval = interp::in_interval(x, u_spline.x_vals[0], u_spline.x_vals[u_spline.length_x - 1]);
    y_eval = interp::in_interval(y, u_spline.y_vals[0], u_spline.y_vals[u_spline.length_y - 1]);
    z_eval = interp::in_interval(z, u_spline.z_vals[0], u_spline.z_vals[u_spline.length_z - 1]);
    if (geoac::is_topo && (z < topo::z_bndlyr)){
        interp::eval_all(x, y, topo::spline, zg, dzg, ddzg_temp, dddzg_temp);
        for(int n = 0; n < 2; n++){
            for(int m = 0; m < 2; m++){
                ddzg[n][m] = ddzg_temp[n + m];
                for(int k = 0; k < 2; k++){
                    dddzg[n][m][k] = dddzg_temp[n + m + k];
                }
            }
        }
    } else {
        zg = topo::z(x, y);
        for(int n = 0; n < 2; n++){
            dzg[n] = 0.0;
            for(int m = 0; m < 2; m++){
                ddzg[n][m] = 0.0;
            }
        }
        
    }
    // Evaluate the u and v splines and derivatives
    interp::eval_all(x_eval, y_eval, z_eval, atmo::u_spline, u0, du0, ddu0);
    interp::eval_all(x_eval, y_eval, z_eval, atmo::v_spline, v0, dv0, ddv0);
    
    // Compute the boundary layer (BL) scaling
    X = exp(-(z - zg) / topo::bndlyr_param);
    for(int n = 0; n < 2; n++){
        dX[n] = X / topo::bndlyr_param * dzg[n];
        for(int m = n; m < 3; m++){
            if(m < 2){
                ddX[n][m] = X / topo::bndlyr_param * (ddzg[n][m] + dzg[n] * dzg[m] / topo::bndlyr_param);
                for(int k = 0; k < 2; k++){
                    dddX[k][n][m]  = dddzg[n][m][k];
                    dddX[k][n][m] += (ddzg[n][m] * dzg[k] + ddzg[k][m] * dzg[n] + ddzg[k][n] * dzg[m]) / topo::bndlyr_param;
                    dddX[k][n][m] += dzg[k] * dzg[n] * dzg[m] / pow(topo::bndlyr_param, 2);
                    dddX[k][n][m] *= X / topo::bndlyr_param;
                } 
            } else {
                ddX[n][m] = - X / pow(topo::bndlyr_param, 2) * dzg[n];
                dddX[0][n][m] = -X / pow(topo::bndlyr_param, 2) * (ddzg[n][0] + dzg[n] * dzg[0] / topo::bndlyr_param);
                dddX[1][n][m] = -X / pow(topo::bndlyr_param, 2) * (ddzg[n][1] + dzg[n] * dzg[1] / topo::bndlyr_param);
            }
            ddX[m][n] = ddX[n][m];
            dddX[0][m][n] = dddX[0][n][m];
            dddX[1][m][n] = dddX[1][n][m];
        }
    }
    dX[2] = -X / topo::bndlyr_param;
    ddX[2][2] = X / pow(topo::bndlyr_param, 2);
    dddX[0][2][2] = X / pow(topo::bndlyr_param, 3) * dzg[0];
    dddX[1][2][2] = X / pow(topo::bndlyr_param, 3) * dzg[1];

    sc0 = (1.0 - X) / (1.0 + X);
    for(int n = 0; n < 3; n++){
        dsc0[n] = -2.0 * dX[n] / pow(1.0 + X, 2);
        for(int m = n; m < 3; m++){
            ddsc0[n][m] = -2.0 * (ddX[n][m] / pow(1.0 + X, 2) - 2.0 * dX[n] * dX[m] / pow(1.0 + X, 3));
            ddsc0[m][n] = ddsc0[n][m];

            for(int k = 0; k < 2; k++){
                dddsc0[k][n][m]  = dddX[k][n][m] / pow(1.0 + X, 2);
                dddsc0[k][n][m] -= 2.0 * (ddX[n][m] * dX[k] + ddX[n][k] * dX[m] + ddX[k][m] * dX[n]) / pow(1.0 + X, 3);
                dddsc0[k][n][m] += 6.0 * dX[n] * dX[m] * dX[k] / pow(1.0 + X, 4);
                dddsc0[k][n][m] *= -2.0;

                dddsc0[k][m][n] = dddsc0[k][n][m];
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
    if (geoac::is_topo && (z < topo::z_bndlyr) && topo::use_BLw){
        w0 = -2.0 * topo::bndlyr_param * (u0 * dsc0[0] + v0 * dsc0[1]);
        sg_ratio = topo::bndlyr_param / topo::vert_wind_grad;
        
        N = sc0 * X * w0;
        D = X + (sg_ratio / 2.0) * (1.0 - pow(X, 2));
        w = N / D;
        
        for(int n = 0; n < 3; n++){
            if(n < 2){  dw0[n] = -2.0 * topo::bndlyr_param * (u0 * ddsc0[0][n] + v0 * ddsc0[1][n]);}
            else {      dw0[n] = -2.0 * topo::bndlyr_param * (du0[2] * dsc0[0] + u0 * ddsc0[0][2] + dv0[2] * dsc0[1] + v0 * ddsc0[1][2]);}            
            dN[n] = (dsc0[n] * X + sc0 * dX[n]) * w0 + sc0 * X * dw0[n];
            dD[n] = (1.0 - sg_ratio * X) * dX[n];
            
            dw[n] = (dN[n] * D - N * dD[n]) / pow(D, 2);
            
            if(n < 2){
                ddw0 = -2.0 * topo::bndlyr_param * (dddsc0[0][n][n] * u0 + dddsc0[1][n][n] * v0);
            } else {
                ddw0  = dddsc0[0][2][2] * u0 + 2.0 * ddsc0[0][2] * du0[2] + dsc0[0] * ddu0[2];
                ddw0 += dddsc0[1][2][2] * v0 + 2.0 * ddsc0[1][2] * dv0[2] + dsc0[1] * ddv0[2];
                ddw0 *= -2.0 * topo::bndlyr_param;
            }

            ddN  = ddsc0[n][n] * X * w0 + sc0 * ddX[n][n] * w0 + sc0 * X * ddw0;
            ddN += 2.0 * (dsc0[n] * dX[n] * w0 + dsc0[n] * X * dw0[n] + sc0 * dX[n] * dw0[n]);
            ddD = - sg_ratio * pow(dX[n], 2) + (1.0 - sg_ratio * X) * ddX[n][n];
            
            ddw[n] = ((ddN * D - 2.0 * dN[n] * dD[n] - N * ddD) * pow(D, 2) + 2.0 * N * D * pow(dD[n], 2)) / pow(D, 4);
        }

        // ddw[3-5] are the dxdy, dxdz, and dydz derivatives ([0][1], [0][2],  and [1][2] indices)
        ddw0 = -2.0 * topo::bndlyr_param * (dddsc0[0][0][1] * u0 + dddsc0[1][0][1] * v0);
        ddN  = ddsc0[0][1] * X * w0 + sc0 * ddX[0][1] * w0 + sc0 * X * ddw0;
        ddN += dsc0[0] * dX[1] * w0 + dsc0[0] * X * dw0[1] + sc0 * dX[0] * dw0[1];
        ddN += dsc0[1] * dX[0] * w0 + dsc0[1] * X * dw0[0] + sc0 * dX[1] * dw0[0];
        ddD = - sg_ratio * dX[0] * dX[1] + (1.0 - sg_ratio * X) * ddX[0][1];
        ddw[3] = ((ddN * D - dN[0] * dD[1] - dN[1] * dD[0] - N * ddD) * pow(D, 2) + 2.0 * N * D * dD[0] * dD[1]) / pow(D, 4);
        
        ddw0  = dddsc0[0][0][2] * u0 + ddsc0[0][0] * du0[2];
        ddw0 += dddsc0[1][0][2] * v0 + ddsc0[1][0] * dv0[2];
        ddw0 *= -2.0 * topo::bndlyr_param;
        ddN  = ddsc0[0][2] * X * w0 + sc0 * ddX[0][2] * w0 + sc0 * X * ddw0;
        ddN += dsc0[0] * dX[2] * w0 + dsc0[0] * X * dw0[2] + sc0 * dX[0] * dw0[2];
        ddN += dsc0[2] * dX[0] * w0 + dsc0[2] * X * dw0[0] + sc0 * dX[2] * dw0[0];
        ddD = - sg_ratio * dX[0] * dX[2] + (1.0 - sg_ratio * X) * ddX[0][2];
        ddw[4] = ((ddN * D - dN[0] * dD[2] - dN[2] * dD[0] - N * ddD) * pow(D, 2) + 2.0 * N * D * dD[0] * dD[2]) / pow(D, 4);
        
        ddw0  = dddsc0[0][1][2] * u0 + ddsc0[0][1] * du0[2];
        ddw0 += dddsc0[1][1][2] * v0 + ddsc0[1][1] * dv0[2];
        ddw0 *= -2.0 * topo::bndlyr_param;
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

#endif /* _ATMO_STATE_3D_RNGDEP_CPP_ */
