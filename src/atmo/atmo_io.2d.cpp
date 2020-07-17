# ifndef _ATMO_IO_2D_CPP_
# define _ATMO_IO_2D_CPP_

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

//-----------------------------------------------//
//---------Define the Propagation Region---------//
//------------For 2D Planar Propgation-----------//
//-----------------------------------------------//
void geoac::set_limits(){
    if(geoac::is_topo){ rng_max = topo::spline.x_vals[topo::spline.length - 1];}
    else{               rng_max = 2000.0;}
    
    alt_max  =	atmo::c_spline.x_vals[atmo::c_spline.length - 1];
}

//-----------------------------------//
//-----Topography and Atmosphere-----//
//-----------Interpolations----------//
//-----------------------------------//
struct interp::natural_cubic_spline_1D topo::spline;

struct interp::natural_cubic_spline_1D atmo::rho_spline;
struct interp::natural_cubic_spline_1D atmo::c_spline;
struct interp::natural_cubic_spline_1D atmo::u_spline;
struct interp::natural_cubic_spline_1D atmo::v_spline;

//-------------------------------------------//
//--------------Set Up or Clear--------------//
//---------Topography and Atmosphere---------//
//--------------Interpolations---------------//
//-------------------------------------------// 
void set_region(char* atmo_file, char* atmo_format, bool invert_winds){
    cout << "Interpolating atmosphere data in '" << atmo_file << "' using format '" << atmo_format << "'..." << '\n';

    interp::prep(atmo::c_spline,    file_length(atmo_file));
    interp::prep(atmo::u_spline,    atmo::c_spline.length);
    interp::prep(atmo::v_spline,    atmo::c_spline.length);
    interp::prep(atmo::rho_spline,  atmo::c_spline.length);
    
    int n = 0;
    double temp;
    string line;
    ifstream file_in;
    
    file_in.open(atmo_file);
    while(!file_in.eof() && n < atmo::c_spline.length){
        getline (file_in, line);
        if(line.find("#") != 0){
            stringstream ss(line);
            if(strncmp(atmo_format, "zTuvdp", 6) == 0){
                ss >> atmo::c_spline.x_vals[n];    // Extract z_i value
                ss >> temp;                        // Extract T(z_i) but don't store it
                ss >> atmo::u_spline.f_vals[n];    // Extract u(z_i)
                ss >> atmo::v_spline.f_vals[n];    // Extract v(z_i)
                ss >> atmo::rho_spline.f_vals[n];  // Extract rho(z_i)
                ss >> atmo::c_spline.f_vals[n];    // Extract p(z_i) into c(z_i) and convert below
            } else if (strncmp(atmo_format, "zuvwTdp", 7) == 0){
                ss >> atmo::c_spline.x_vals[n];    // Extract z_i value
                ss >> atmo::u_spline.f_vals[n];    // Extract u(z_i)
                ss >> atmo::v_spline.f_vals[n];    // Extract v(z_i)
                ss >> temp;                        // Extract w(z_i) but don't store it
                ss >> temp;                        // Extract T(z_i) but don't store it
                ss >> atmo::rho_spline.f_vals[n];  // Extract rho(z_i)
                ss >> atmo::c_spline.f_vals[n];    // Extract p(z_i) into c(z_i) and convert below
            } else if (strncmp(atmo_format, "zcuvd", 5) == 0){
                ss >> atmo::c_spline.x_vals[n];    // Extract z_i value
                ss >> atmo::c_spline.f_vals[n];    // Extract c(z_i)
                ss >> atmo::u_spline.f_vals[n];    // Extract u(z_i)
                ss >> atmo::v_spline.f_vals[n];    // Extract v(z_i)
                ss >> atmo::rho_spline.f_vals[n];  // Extract rho(z_i)
            } else {
                cout << "Unrecognized profile option: " << atmo_format << ".  Valid options are: zTuvdp, zuvwTdp, or zcuvd" << '\n';
                break;
            }
       
            // Copy altitude values to other interpolations
            atmo::u_spline.x_vals[n] = atmo::c_spline.x_vals[n];
            atmo::v_spline.x_vals[n] = atmo::c_spline.x_vals[n];
            atmo::rho_spline.x_vals[n] = atmo::c_spline.x_vals[n];
        
            // Convert pressure and density to adiabatic sound speed unless c is specified and scale winds from m/s to km/s
            if (strncmp(atmo_format, "zTuvdp", 6) == 0 || strncmp(atmo_format, "zuvwTdp", 7) == 0){
                atmo::c_spline.f_vals[n] = sqrt(0.1 * atmo::gam * atmo::c_spline.f_vals[n] / atmo::rho_spline.f_vals[n]) / 1000.0;
            } else {
                atmo::c_spline.f_vals[n] /= 1000.0;
            }
        
            if(invert_winds){
                atmo::u_spline.f_vals[n] /= -1000.0;
                atmo::v_spline.f_vals[n] /= -1000.0;
            } else {
                atmo::u_spline.f_vals[n] /= 1000.0;
                atmo::v_spline.f_vals[n] /= 1000.0;
            }
            n++;
        }
    }
    file_in.close();

    interp::set(atmo::c_spline);    interp::set(atmo::u_spline);
    interp::set(atmo::rho_spline);  interp::set(atmo::v_spline);
    
    geoac::set_limits();
    cout << '\t' << "Propagation region limits:" << '\n';
    cout << '\t' << '\t' << "r = 0.0, " << geoac::rng_max << '\n';
    cout << '\t' << '\t' << "z = 0.0, " << geoac::alt_max << '\n';
    
    topo::set_bndlyr();
    // cout << '\t' << "Maximum topography height: " << topo::z_max << '\n';
    // cout << '\t' << "Boundary layer height: " << topo::z_bndlyr << '\n';
}


void set_region(char* atmo_file, char* topo_file, char* atmo_format, bool invert_winds){
    cout << "Interpolating atmosphere data in '" << atmo_file << "' and topography data in '" << topo_file << "'..." << '\n';

    ifstream file_in;
    
    interp::prep(topo::spline, file_length(topo_file));
    file_in.open(topo_file);
    for (int n = 0; n < topo::spline.length; n++){
        file_in >> topo::spline.x_vals[n];
        file_in >> topo::spline.f_vals[n];
    }
    file_in.close();

    interp::prep(atmo::c_spline, file_length(atmo_file));
    interp::prep(atmo::u_spline,    atmo::c_spline.length);
    interp::prep(atmo::v_spline,    atmo::c_spline.length);
    interp::prep(atmo::rho_spline,  atmo::c_spline.length);

    int n = 0;
    double temp;
    string line;
   
    file_in.open(atmo_file);
    while(!file_in.eof() && n < atmo::c_spline.length){
        getline (file_in, line);
        if(line.find("#") != 0){
            stringstream ss(line);
            if(strncmp(atmo_format, "zTuvdp", 6) == 0){
                ss >> atmo::c_spline.x_vals[n];    // Extract z_i value
                ss >> temp;                        // Extract T(z_i) but don't store it
                ss >> atmo::u_spline.f_vals[n];    // Extract u(z_i)
                ss >> atmo::v_spline.f_vals[n];    // Extract v(z_i)
                ss >> atmo::rho_spline.f_vals[n];  // Extract rho(z_i)
                ss >> atmo::c_spline.f_vals[n];    // Extract p(z_i) into c(z_i) and convert below
            } else if (strncmp(atmo_format, "zuvwTdp", 7) == 0){
                ss >> atmo::c_spline.x_vals[n];    // Extract z_i value
                ss >> atmo::u_spline.f_vals[n];    // Extract u(z_i)
                ss >> atmo::v_spline.f_vals[n];    // Extract v(z_i)
                ss >> temp;                        // Extract w(z_i) but don't store it
                ss >> temp;                        // Extract T(z_i) but don't store it
                ss >> atmo::rho_spline.f_vals[n];  // Extract rho(z_i)
                ss >> atmo::c_spline.f_vals[n];    // Extract p(z_i) into c(z_i) and convert below
            } else if (strncmp(atmo_format, "zcuvd", 5) == 0){
                ss >> atmo::c_spline.x_vals[n];    // Extract z_i value
                ss >> atmo::c_spline.f_vals[n];    // Extract c(z_i)
                ss >> atmo::u_spline.f_vals[n];    // Extract u(z_i)
                ss >> atmo::v_spline.f_vals[n];    // Extract v(z_i)
                ss >> atmo::rho_spline.f_vals[n];  // Extract rho(z_i)
            } else {
                cout << "Unrecognized profile option: " << atmo_format << ".  Valid options are: zTuvdp, zuvwTdp, or zcuvd" << '\n';
                break;
            }
        
            // Copy altitude values to other interpolations
            atmo::u_spline.x_vals[n] = atmo::c_spline.x_vals[n];
            atmo::v_spline.x_vals[n] = atmo::c_spline.x_vals[n];
            atmo::rho_spline.x_vals[n] = atmo::c_spline.x_vals[n];
        
            // Convert pressure and density to adiabatic sound speed unless c is specified and scale winds from m/s to km/s
            if (strncmp(atmo_format, "zTuvdp", 6) == 0 || strncmp(atmo_format, "zuvwTdp", 7) == 0){
                atmo::c_spline.f_vals[n] = sqrt(0.1 * atmo::gam * atmo::c_spline.f_vals[n] / atmo::rho_spline.f_vals[n]) / 1000.0;
            } else {
                atmo::c_spline.f_vals[n] /= 1000.0;
            }
        
            if(invert_winds){
                atmo::u_spline.f_vals[n] /= -1000.0;
                atmo::v_spline.f_vals[n] /= -1000.0;
            } else {
                atmo::u_spline.f_vals[n] /= 1000.0;
                atmo::v_spline.f_vals[n] /= 1000.0;
            }
            n++;
        }
    }
    file_in.close();

    interp::set(topo::spline);
    interp::set(atmo::c_spline);    interp::set(atmo::u_spline);
    interp::set(atmo::rho_spline);  interp::set(atmo::v_spline);

    geoac::set_limits();
    cout << '\t' << "Propagation region limits:" << '\n';
    cout << '\t' << '\t' << "r = " << topo::spline.x_vals[0] << ", " << topo::spline.x_vals[topo::spline.length - 1] << '\n';
    cout << '\t' << '\t' << "z = 0.0, " << geoac::alt_max << '\n' << '\n';

    topo::set_bndlyr();
    cout << '\t' << "Maximum topography height: " << topo::z_max << '\n';
    cout << '\t' << "Boundary layer height: " << topo::z_bndlyr << '\n';
}


void clear_region(){
    if(geoac::is_topo){                 interp::clear(topo::spline);}
    interp::clear(atmo::c_spline);      interp::clear(atmo::u_spline);
    interp::clear(atmo::rho_spline);    interp::clear(atmo::v_spline);
}


#endif /* _ATMO_IO_2D_CPP_ */
