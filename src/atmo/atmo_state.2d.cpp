# ifndef _ATMO_STATE_2D_CPP_
# define _ATMO_STATE_2D_CPP_

#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <math.h>

#include "atmo_state.h"
#include "atmo_io.2d.h"
#include "../util/interpolation.h"
#include "../util/fileIO.h"
#include "../geoac/geoac.params.h"

using namespace std;

//----------------------------//
//-----Physical Constants-----//
//----------------------------//
double atmo::gam = 1.4;
double atmo::R = 287.058;

double atmo::z_reflect = 100.0;

//---------------------------------------------------//
//---------Functions Defining the Topography---------//
//--------------For 2D Planar Propgation-------------//
//---------------------------------------------------//
double topo::z0 = 0.0;
double topo::z_max;
double topo::z_bndlyr;

void topo::set_bndlyr(){
    if (geoac::is_topo){
        z_max = 0.0;
        for (int nr = 0; nr < spline.length; nr++){
            z_max = max(z_max, spline.f_vals[nr]);
        }
        z_max += 0.1;
    } else {
        z_max = z0;
    }
    z_bndlyr = z_max + 2.0;    
}

double topo::z(double r){
    double result = z0;
    if (geoac::is_topo){
        result = interp::eval_f(interp::in_interval(r, spline.x_vals[0], spline.x_vals[spline.length - 1]), spline);
    }
    return result;
}

double topo::dz(double r){
    double result = 0.0;
    if (geoac::is_topo){
        result = interp::eval_df(interp::in_interval(r, spline.x_vals[0], spline.x_vals[spline.length - 1]), spline);
    }
    return result;
}

double topo::ddz(double r){ 
    double result = 0.0;
    if (geoac::is_topo){
        return interp::eval_ddf(interp::in_interval(r, spline.x_vals[0], spline.x_vals[spline.length - 1]), spline);
    }
    return result;
}

//---------------------------------------------------//
//---------Functions Defining the Atmosphere---------//
//--------------For 2D Planar Propgation-------------//
//---------------------------------------------------//
double atmo::rho(double x, double y, double z){
    double z_eval = interp::in_interval(z, rho_spline.x_vals[0], rho_spline.x_vals[rho_spline.length - 1]);
    return interp::eval_f(interp::in_interval(z_eval, rho_spline.x_vals[0], rho_spline.x_vals[rho_spline.length - 1]),rho_spline);
}


double atmo::c(double x, double y, double z){
    double z_eval = interp::in_interval(z, c_spline.x_vals[0], c_spline.x_vals[c_spline.length - 1]);
    return interp::eval_f(interp::in_interval(z_eval, c_spline.x_vals[0], c_spline.x_vals[c_spline.length - 1]), c_spline);
}

double atmo::dc(double x, double y, double z, int n){
    double z_eval, result = 0.0;
    if(n == 2){
        z_eval = interp::in_interval(z, c_spline.x_vals[0], c_spline.x_vals[c_spline.length - 1]);
        result = interp::eval_df(interp::in_interval(z_eval, c_spline.x_vals[0], c_spline.x_vals[c_spline.length - 1]), c_spline);
    }
    return result;
}

double atmo::ddc(double x, double y, double z, int n1, int n2){
    double z_eval, result = 0.0;
    if(n1 == 2 && n2 ==  2){
        z_eval = interp::in_interval(z, c_spline.x_vals[0], c_spline.x_vals[c_spline.length - 1]);
        result = interp::eval_ddf(interp::in_interval(z_eval, c_spline.x_vals[0], c_spline.x_vals[c_spline.length - 1]), c_spline);
    }
    return result;
}

double atmo::u(double x, double y, double z){
    double z_eval = interp::in_interval(z, u_spline.x_vals[0], u_spline.x_vals[u_spline.length - 1]);
    return interp::eval_f(interp::in_interval(z_eval, u_spline.x_vals[0], u_spline.x_vals[u_spline.length - 1]), u_spline);
}

double atmo::du(double x, double y, double z, int n){
    double z_eval, result = 0.0;
    if(n == 2){
        z_eval = interp::in_interval(z, u_spline.x_vals[0], u_spline.x_vals[u_spline.length - 1]);
        result = interp::eval_df(interp::in_interval(z_eval, u_spline.x_vals[0], u_spline.x_vals[u_spline.length - 1]), u_spline);
    }
    return result;
}

double atmo::ddu(double x, double y, double z, int n1, int n2){
    double z_eval, result = 0.0;
    if(n1 == 2 && n2 == 2){
        z_eval = interp::in_interval(z, u_spline.x_vals[0], u_spline.x_vals[u_spline.length - 1]);
        result = interp::eval_ddf(interp::in_interval(z_eval, u_spline.x_vals[0], u_spline.x_vals[u_spline.length - 1]), u_spline);
    }
    return result;
}

double atmo::v(double x, double y, double z){
    double z_eval = interp::in_interval(z, v_spline.x_vals[0], v_spline.x_vals[c_spline.length - 1]);
    return interp::eval_f(interp::in_interval(z_eval, v_spline.x_vals[0], v_spline.x_vals[v_spline.length - 1]), v_spline);
}

double atmo::dv(double x, double y, double z, int n){
    double z_eval, result = 0.0;
    if(n == 2){
        z_eval = interp::in_interval(z, v_spline.x_vals[0], v_spline.x_vals[v_spline.length - 1]);
        result = interp::eval_df(interp::in_interval(z_eval, v_spline.x_vals[0], v_spline.x_vals[v_spline.length - 1]), v_spline);
    }
    return result;
}

double atmo::ddv(double x, double y, double z, int n1, int n2){
    double z_eval, result = 0.0;
    if(n1 == 2 && n2  == 2){
        z_eval = interp::in_interval(z, v_spline.x_vals[0], v_spline.x_vals[v_spline.length - 1]);
        result = interp::eval_ddf(interp::in_interval(z_eval, v_spline.x_vals[0], v_spline.x_vals[v_spline.length - 1]), v_spline);
    }
    return result;
}

double atmo::w(double x, double y, double z){
    return 0.0;
}

double atmo::dw(double x, double y, double z, int n){
    return 0.0;
}

double atmo::ddw(double x, double y, double z, int n1, int n2){
    return 0.0;
}

#endif /* _ATMO_STATE_2D_CPP_ */
