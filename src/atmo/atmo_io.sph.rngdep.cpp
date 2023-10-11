# ifndef _ATMO_IO_SPH_RNGDEP_CPP_
# define _ATMO_IO_SPH_RNGDEP_CPP_

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

//-----------------------------------------------//
//---------Define the Propagation Region---------//
//------------For 3D Global Propgation-----------//
//-----------------------------------------------//
void geoac::set_limits(){
    double buffer = 1.0e-3 * (Pi / 180.0);

    lat_min = atmo::c_spline.x_vals[0] + buffer;  lat_max = atmo::c_spline.x_vals[atmo::c_spline.length_x - 1] - buffer;
    lon_min = atmo::c_spline.y_vals[0] + buffer;  lon_max = atmo::c_spline.y_vals[atmo::c_spline.length_y - 1] - buffer;
    
    if(is_topo){
        lat_min = max(lat_min, topo::spline.x_vals[0] + buffer); lat_max = min(lat_max, topo::spline.x_vals[topo::spline.length_x - 1] - buffer);
        lon_min = max(lon_min, topo::spline.y_vals[0] + buffer); lon_max = min(lon_max, topo::spline.y_vals[topo::spline.length_y - 1] - buffer);
        if(lat_max < lat_min || lon_max < lon_min){
            cout << '\n' << '\n' << "WARNING!!!  Specified atmosphere grid and topogrpahy grid have no overlap.  Check grid definitions." << '\n' << '\n';
        }
    }

    topo::z0 = atmo::c_spline.z_vals[0];
    alt_max = atmo::c_spline.z_vals[atmo::c_spline.length_z - 1];
}

//-----------------------------------//
//-----Topography and Atmosphere-----//
//-----------Interpolations----------//
//-----------------------------------//
struct interp::natural_cubic_spline_2D topo::spline;

struct interp::hybrid_spline_3D atmo::rho_spline;
struct interp::hybrid_spline_3D atmo::c_spline;
struct interp::hybrid_spline_3D atmo::u_spline;
struct interp::hybrid_spline_3D atmo::v_spline;

//-------------------------------------------//
//--------------Set Up or Clear--------------//
//---------Topography and Atmosphere---------//
//--------------Interpolations---------------//
//-------------------------------------------//
int set_region(char* atmo_prefix, char* atmo_locs_lat, char* atmo_locs_lon, char* atmo_format, bool invert_winds){
    cout << "Interpolating atmosphere data in '" << atmo_prefix << "'* using format '" << atmo_format << "'..." << '\n';
    if(invert_winds){
        cout << '\t' << "Inverting wind fields for back projection analysis..." << '\n';
    }

    int lat_cnt, lon_cnt, z_cnt;
    double temp;
    char output_buffer [512];
    string line;
    ifstream file_in;
    
    lat_cnt = file_length(atmo_locs_lat);
    lon_cnt = file_length(atmo_locs_lon);
    if(lat_cnt < 2 || lon_cnt < 2){
        cout << '\t' << "ERROR: Invalid grid specifications (check grid node files)" << '\n' << '\n';
        return 0;
    }

    sprintf(output_buffer, "%s%i.met", atmo_prefix, 0);
    file_in.open(output_buffer);
    if(!file_in.is_open()){
        cout << '\t' << "ERROR: Invalid atmospheric specification (" << output_buffer << ")" << '\n' << '\n';
        return 0;
    }
    file_in.close();

    z_cnt = file_length(output_buffer);
    interp::prep(atmo::c_spline, lat_cnt, lon_cnt, z_cnt);
    interp::prep(atmo::u_spline, lat_cnt, lon_cnt, z_cnt);
    interp::prep(atmo::v_spline, lat_cnt, lon_cnt, z_cnt);
    interp::prep(atmo::rho_spline, lat_cnt, lon_cnt, z_cnt);
    
    file_in.open(atmo_locs_lat);
    for (int n_lat = 0; n_lat < lat_cnt; n_lat++){
        file_in >> atmo::c_spline.x_vals[n_lat];
        
        atmo::c_spline.x_vals[n_lat] *= Pi / 180.0;
        atmo::u_spline.x_vals[n_lat] = atmo::c_spline.x_vals[n_lat];
        atmo::v_spline.x_vals[n_lat] = atmo::c_spline.x_vals[n_lat];
        atmo::rho_spline.x_vals[n_lat] = atmo::c_spline.x_vals[n_lat];
    }
    file_in.close();
    
    file_in.open(atmo_locs_lon);
    for (int n_lon = 0; n_lon < lon_cnt; n_lon++){
        file_in >> atmo::c_spline.y_vals[n_lon];
        
        atmo::c_spline.y_vals[n_lon] *= Pi / 180.0;
        atmo::u_spline.y_vals[n_lon] = atmo::c_spline.y_vals[n_lon];
        atmo::v_spline.y_vals[n_lon] = atmo::c_spline.y_vals[n_lon];
        atmo::rho_spline.y_vals[n_lon] = atmo::c_spline.y_vals[n_lon];
    }
    file_in.close();

    for(int n_lat = 0; n_lat < lat_cnt; n_lat++){
    for(int n_lon = 0; n_lon < lon_cnt; n_lon++){
        sprintf(output_buffer, "%s%i.met", atmo_prefix, n_lat * lon_cnt + n_lon);
        file_in.open(output_buffer);
        if(!file_in.is_open()){
            cout << '\t' << "ERROR: Invalid atmospheric specification (" << output_buffer << ")" << '\n' << '\n';
            return 0;
        }

        if((n_lat < 3) && (n_lon < 3)){
            cout << '\t' << "Setting grid node at (" << atmo::c_spline.x_vals[n_lat] * 180.0 / Pi << ", " << atmo::c_spline.y_vals[n_lon] * 180.0 / Pi;
            cout << ") with profile " << output_buffer << '\n';
        } else if((n_lat < 3) && (n_lon == 3)) {
            cout << '\t' << "..." << '\n';
        }
            
        int nz = 0;
        while(!file_in.eof() && nz < z_cnt){
            getline (file_in, line);
            if(line.find("#") != 0){
                stringstream ss(line);
                if(strncmp(atmo_format, "zTuvdp", 6) == 0){
                    ss >> atmo::c_spline.z_vals[nz];                // Extract z_i value
                    ss >> temp;                                     // Extract T(z_i) but don't store it
                    ss >> atmo::u_spline.f_vals[n_lat][n_lon][nz];        // Extract u(z_i)
                    ss >> atmo::v_spline.f_vals[n_lat][n_lon][nz];        // Extract v(z_i)
                    ss >> atmo::rho_spline.f_vals[n_lat][n_lon][nz];      // Extract rho(z_i)
                    ss >> atmo::c_spline.f_vals[n_lat][n_lon][nz];        // Extract p(z_i) into c(z_i) and convert below
                } else if (strncmp(atmo_format, "zuvwTdp", 7) == 0){
                    ss >> atmo::c_spline.z_vals[nz];                // Extract z_i value
                    ss >> atmo::u_spline.f_vals[n_lat][n_lon][nz];        // Extract u(z_i)
                    ss >> atmo::v_spline.f_vals[n_lat][n_lon][nz];        // Extract v(z_i)
                    ss >> temp;                                     // Extract w(z_i) but don't store it
                    ss >> temp;                                     // Extract T(z_i) but don't store it
                    ss >> atmo::rho_spline.f_vals[n_lat][n_lon][nz];      // Extract rho(z_i)
                    ss >> atmo::c_spline.f_vals[n_lat][n_lon][nz];        // Extract p(z_i) into c(z_i) and convert below
                } else if (strncmp(atmo_format, "zcuvd", 5) == 0){
                    ss >> atmo::c_spline.z_vals[nz];                // Extract z_i value
                    ss >> atmo::c_spline.f_vals[n_lat][n_lon][nz];        // Extract c(z_i)
                    ss >> atmo::u_spline.f_vals[n_lat][n_lon][nz];        // Extract u(z_i)
                    ss >> atmo::v_spline.f_vals[n_lat][n_lon][nz];        // Extract v(z_i)
                    ss >> atmo::rho_spline.f_vals[n_lat][n_lon][nz];      // Extract rho(z_i)
                } else {
                    cout << "Unrecognized profile option: " << atmo_format << ".  Valid options are: zTuvdp, zuvwTdp, or zcuvd" << '\n';
                    break;
                }
                
                // Copy altitude values to other interpolations
                atmo::u_spline.z_vals[nz] = atmo::c_spline.z_vals[nz];
                atmo::v_spline.z_vals[nz] = atmo::c_spline.z_vals[nz];
                atmo::rho_spline.z_vals[nz] = atmo::c_spline.z_vals[nz];
                
                // Convert pressure and density to adiabatic sound speed unless c is specified and scale winds from m/s to km/s
                if (strncmp(atmo_format, "zTuvdp", 6) == 0 || strncmp(atmo_format, "zuvwTdp", 7) == 0){
                    atmo::c_spline.f_vals[n_lat][n_lon][nz] = sqrt(0.1 * atmo::gam * atmo::c_spline.f_vals[n_lat][n_lon][nz] / atmo::rho_spline.f_vals[n_lat][n_lon][nz]) / 1000.0;
                } else {
                    atmo::c_spline.f_vals[n_lat][n_lon][nz] /= 1000.0;
                }
                
                if(invert_winds){
                    atmo::u_spline.f_vals[n_lat][n_lon][nz] /= -1000.0;
                    atmo::v_spline.f_vals[n_lat][n_lon][nz] /= -1000.0;
                } else {
                    atmo::u_spline.f_vals[n_lat][n_lon][nz] /= 1000.0;
                    atmo::v_spline.f_vals[n_lat][n_lon][nz] /= 1000.0;
                }
                nz++;
            }
        }
        file_in.close();
    }}
    cout << '\n';

    interp::set(atmo::c_spline);    interp::set(atmo::u_spline);
    interp::set(atmo::rho_spline);  interp::set(atmo::v_spline);
    
    geoac::set_limits();
    cout << '\t' << "Propagation region limits:" << '\n';
    cout << '\t' << '\t' << "latitude = " << geoac::lat_min * (180.0 / Pi) << ", " << geoac::lat_max * (180.0 / Pi) << '\n';
    cout << '\t' << '\t' << "longitude = " << geoac::lon_min * (180.0 / Pi) << ", " << geoac::lon_max * (180.0 / Pi) << '\n';
    cout << '\t' << '\t' << "altitude = " << topo::z0 << ", " << geoac::alt_max << '\n' << '\n';

    topo::set_bndlyr();
    return 1;
}


int set_region(char* atmo_prefix, char* atmo_locs_lat, char* atmo_locs_lon, char* topo_file, char* atmo_format, bool invert_winds){
    cout << "Interpolating atmosphere data in '" << atmo_prefix << "*' and topography data in '" << topo_file << "'..." << '\n';
    if(invert_winds){
        cout << '\t' << "Inverting wind fields for back projection analysis..." << '\n';
    }

    int lat_cnt, lon_cnt, z_cnt;
    double temp;
    char output_buffer [512];
    string line;
    ifstream file_in;
    
    file_2d_dims(topo_file, lat_cnt, lon_cnt);
    interp::prep(topo::spline, lat_cnt, lon_cnt);
    
    file_in.open(topo_file);
    if(!file_in.is_open()){
        cout << '\t' << "ERROR: Invalid terrain file (" << topo_file << ")" << '\n' << '\n';
        return 0;
    }
    for (int n_lat = 0; n_lat < topo::spline.length_x; n_lat++){
        for (int n_lon = 0; n_lon < topo::spline.length_y; n_lon++){
            file_in >> topo::spline.x_vals[n_lat];
            file_in >> topo::spline.y_vals[n_lon];
            file_in >> topo::spline.f_vals[n_lat][n_lon];
        }}
    file_in.close();
    
    for (int n_lat = 0; n_lat < topo::spline.length_x; n_lat++){ topo::spline.x_vals[n_lat] *= Pi / 180.0;}
    for (int n_lon = 0; n_lon < topo::spline.length_y; n_lon++){ topo::spline.y_vals[n_lon] *= Pi / 180.0;}
    
    lat_cnt = file_length(atmo_locs_lat);
    lon_cnt = file_length(atmo_locs_lon);
    if(lat_cnt < 2 || lon_cnt < 2){
        cout << '\t' << "ERROR: Invalid grid specifications (check grid node files)" << '\n' << '\n';
        return 0;
    }

    sprintf(output_buffer, "%s%i.met", atmo_prefix, 0);
    file_in.open(output_buffer);
    if(!file_in.is_open()){
        cout << '\t' << "ERROR: Invalid atmospheric specification (" << output_buffer << ")" << '\n' << '\n';
        return 0;
    }
    file_in.close();
    
    z_cnt = file_length(output_buffer);    
    interp::prep(atmo::c_spline, lat_cnt, lon_cnt, z_cnt);
    interp::prep(atmo::u_spline, lat_cnt, lon_cnt, z_cnt);
    interp::prep(atmo::v_spline, lat_cnt, lon_cnt, z_cnt);
    interp::prep(atmo::rho_spline, lat_cnt, lon_cnt, z_cnt);
    
    file_in.open(atmo_locs_lat);
    for (int n_lat = 0; n_lat < lat_cnt; n_lat++){
        file_in >> atmo::c_spline.x_vals[n_lat];
        
        atmo::c_spline.x_vals[n_lat] *= Pi / 180.0;
        atmo::u_spline.x_vals[n_lat] = atmo::c_spline.x_vals[n_lat];
        atmo::v_spline.x_vals[n_lat] = atmo::c_spline.x_vals[n_lat];
        atmo::rho_spline.x_vals[n_lat] = atmo::c_spline.x_vals[n_lat];
    }
    file_in.close();
    
    file_in.open(atmo_locs_lon);
    for (int n_lon = 0; n_lon < lon_cnt; n_lon++){
        file_in >> atmo::c_spline.y_vals[n_lon];
        
        atmo::c_spline.y_vals[n_lon] *= Pi / 180.0;
        atmo::u_spline.y_vals[n_lon] = atmo::c_spline.y_vals[n_lon];
        atmo::v_spline.y_vals[n_lon] = atmo::c_spline.y_vals[n_lon];
        atmo::rho_spline.y_vals[n_lon] = atmo::c_spline.y_vals[n_lon];
    }
    file_in.close();
    
    
    for(int n_lat = 0; n_lat < lat_cnt; n_lat++){
    for(int n_lon = 0; n_lon < lon_cnt; n_lon++){
        sprintf(output_buffer, "%s%i.met", atmo_prefix, n_lat * lon_cnt + n_lon);
        file_in.open(output_buffer);
        if(!file_in.is_open()){
            cout << '\t' << "ERROR: Invalid atmospheric specification (" << output_buffer << ")" << '\n' << '\n';
            return 0;
        }
    
        if((n_lat < 3) && (n_lon < 3)){
            cout << '\t' << "Setting grid node at (" << atmo::c_spline.x_vals[n_lat] * 180.0 / Pi << ", " << atmo::c_spline.y_vals[n_lon] * 180.0 / Pi;
            cout << ") with profile " << output_buffer << '\n';
        } else if((n_lat < 3) && (n_lon == 3)) {
            cout << '\t' << "..." << '\n';
        }
            
        int nz = 0;
        while(!file_in.eof() && nz < z_cnt){
            getline (file_in, line);
            if(line.find("#") != 0){
                stringstream ss(line);
                if(strncmp(atmo_format, "zTuvdp", 6) == 0){
                    ss >> atmo::c_spline.z_vals[nz];                // Extract z_i value
                    ss >> temp;                                     // Extract T(z_i) but don't store it
                    ss >> atmo::u_spline.f_vals[n_lat][n_lon][nz];        // Extract u(z_i)
                    ss >> atmo::v_spline.f_vals[n_lat][n_lon][nz];        // Extract v(z_i)
                    ss >> atmo::rho_spline.f_vals[n_lat][n_lon][nz];      // Extract rho(z_i)
                    ss >> atmo::c_spline.f_vals[n_lat][n_lon][nz];        // Extract p(z_i) into c(z_i) and convert below
                } else if (strncmp(atmo_format, "zuvwTdp", 7) == 0){
                    ss >> atmo::c_spline.z_vals[nz];                // Extract z_i value
                    ss >> atmo::u_spline.f_vals[n_lat][n_lon][nz];        // Extract u(z_i)
                    ss >> atmo::v_spline.f_vals[n_lat][n_lon][nz];        // Extract v(z_i)
                    ss >> temp;                                     // Extract w(z_i) but don't store it
                    ss >> temp;                                     // Extract T(z_i) but don't store it
                    ss >> atmo::rho_spline.f_vals[n_lat][n_lon][nz];      // Extract rho(z_i)
                    ss >> atmo::c_spline.f_vals[n_lat][n_lon][nz];        // Extract p(z_i) into c(z_i) and convert below
                } else if (strncmp(atmo_format, "zcuvd", 5) == 0){
                    ss >> atmo::c_spline.z_vals[nz];                // Extract z_i value
                    ss >> atmo::c_spline.f_vals[n_lat][n_lon][nz];        // Extract c(z_i)
                    ss >> atmo::u_spline.f_vals[n_lat][n_lon][nz];        // Extract u(z_i)
                    ss >> atmo::v_spline.f_vals[n_lat][n_lon][nz];        // Extract v(z_i)
                    ss >> atmo::rho_spline.f_vals[n_lat][n_lon][nz];      // Extract rho(z_i)
                } else {
                    cout << "Unrecognized profile option: " << atmo_format << ".  Valid options are: zTuvdp, zuvwTdp, or zcuvd" << '\n';
                    break;
                }
                
                // Copy altitude values to other interpolations
                atmo::u_spline.z_vals[nz] = atmo::c_spline.z_vals[nz];
                atmo::v_spline.z_vals[nz] = atmo::c_spline.z_vals[nz];
                atmo::rho_spline.z_vals[nz] = atmo::c_spline.z_vals[nz];
                
                // Convert pressure and density to adiabatic sound speed unless c is specified and scale winds from m/s to km/s
                if (strncmp(atmo_format, "zTuvdp", 6) == 0 || strncmp(atmo_format, "zuvwTdp", 7) == 0){
                    atmo::c_spline.f_vals[n_lat][n_lon][nz] = sqrt(0.1 * atmo::gam * atmo::c_spline.f_vals[n_lat][n_lon][nz] / atmo::rho_spline.f_vals[n_lat][n_lon][nz]) /1000.0;
                } else {
                    atmo::c_spline.f_vals[n_lat][n_lon][nz] /= 1000.0;
                }
                
                if(invert_winds){
                    atmo::u_spline.f_vals[n_lat][n_lon][nz] /= -1000.0;
                    atmo::v_spline.f_vals[n_lat][n_lon][nz] /= -1000.0;
                } else {
                    atmo::u_spline.f_vals[n_lat][n_lon][nz] /= 1000.0;
                    atmo::v_spline.f_vals[n_lat][n_lon][nz] /= 1000.0;
                }
                nz++;
            }
        }
        file_in.close();
    }}
    
    interp::set(topo::spline);
    interp::set(atmo::c_spline);    interp::set(atmo::u_spline);
    interp::set(atmo::rho_spline);  interp::set(atmo::v_spline);
    
    geoac::set_limits();
    cout << '\t' << "Propagation region limits:" << '\n';
    cout << '\t' << '\t' << "latitude = " << geoac::lat_min * (180.0 / Pi) << ", " << geoac::lat_max * (180.0 / Pi) << '\n';
    cout << '\t' << '\t' << "longitude = " << geoac::lon_min * (180.0 / Pi) << ", " << geoac::lon_max * (180.0 / Pi) << '\n';
    cout << '\t' << '\t' << "altitude = " << topo::z0 << ", " << geoac::alt_max << '\n' << '\n';

    
    topo::set_bndlyr();
    cout << '\t' << "Maximum topography height: " << topo::z_max << '\n';
    cout << '\t' << "Boundary layer height: " << topo::z_bndlyr << '\n';

    return 1;
}


void clear_region(){
    if(geoac::is_topo){
        interp::clear(topo::spline);
    }
    
    interp::clear(atmo::c_spline);      interp::clear(atmo::u_spline);
    interp::clear(atmo::rho_spline);    interp::clear(atmo::v_spline);
}

#endif /* _ATMO_IO_SPH_RNGDEP_CPP_ */
