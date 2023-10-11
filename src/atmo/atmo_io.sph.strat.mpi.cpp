# ifndef _ATMO_IO_SPH_STRAT_MPI_CPP_
# define _ATMO_IO_SPH_STRAT_MPI_CPP_

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <sstream>
#include <math.h>
#include <mpi.h>

#include "atmo_state.h"
#include "atmo_io.sph.strat.h"

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
    if(is_topo){
        lat_min = topo::spline.x_vals[0] + buffer;  lat_max = topo::spline.x_vals[topo::spline.length_x - 1] - buffer;
        lon_min = topo::spline.y_vals[0] + buffer;  lon_max = topo::spline.y_vals[topo::spline.length_y - 1] - buffer;
    } else{
        lat_min = -Pi / 2.0;    lat_max = Pi / 2.0;
        lon_min = -Pi;          lon_max = Pi;
    }
    
    topo::z0 = atmo::c_spline.x_vals[0];
    alt_max = atmo::c_spline.x_vals[atmo::c_spline.length - 1];
}

//-----------------------------------//
//-----Topography and Atmosphere-----//
//-----------Interpolations----------//
//-----------------------------------//
struct interp::natural_cubic_spline_2D  topo::spline;

struct interp::natural_cubic_spline_1D  atmo::rho_spline;
struct interp::natural_cubic_spline_1D  atmo::c_spline;
struct interp::natural_cubic_spline_1D  atmo::u_spline;
struct interp::natural_cubic_spline_1D  atmo::v_spline;

//-------------------------------------------//
//--------------Set Up or Clear--------------//
//---------Topography and Atmosphere---------//
//--------------Interpolations---------------//
//-------------------------------------------//
int set_region(char* atmo_file, char* atmo_format, bool invert_winds, int rank){
    int result = 1;
    if(rank == 0){
        cout << "Interpolating atmosphere data in '" << atmo_file << "' using format '" << atmo_format << "'..." << '\n';
        if(invert_winds){
            cout << '\t' << "Inverting wind fields for back projection analysis..." << '\n';
        }

        interp::prep(atmo::c_spline, file_length(atmo_file));
        interp::prep(atmo::u_spline,   atmo::c_spline.length);
        interp::prep(atmo::v_spline,   atmo::c_spline.length);
        interp::prep(atmo::rho_spline, atmo::c_spline.length);

        int n = 0;
        double temp;
        string line;
        ifstream file_in;

        file_in.open(atmo_file);
        if(!file_in.is_open()){
            cout << '\t' << "ERROR: Invalid atmospheric specification (" << atmo_file << ")" << '\n' << '\n';
            result = 0;
        } else {
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
        }
        file_in.close();
    }

    MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(result == 1){
        // Broadcast atmosphere splines
        int spline_len;
        double spline_vals [5];

        if(rank == 0) spline_len = atmo::c_spline.length;
        MPI_Bcast(&spline_len, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank != 0){
            interp::prep(atmo::c_spline, spline_len);
            interp::prep(atmo::u_spline, spline_len);
            interp::prep(atmo::v_spline, spline_len);
            interp::prep(atmo::rho_spline, spline_len);
        }
    
        for(int n = 0; n < spline_len; n++){
            if(rank == 0){
                spline_vals[0] = atmo::c_spline.x_vals[n];
                spline_vals[1] = atmo::c_spline.f_vals[n];
                spline_vals[2] = atmo::u_spline.f_vals[n];
                spline_vals[3] = atmo::v_spline.f_vals[n];
                spline_vals[4] = atmo::rho_spline.f_vals[n];
            }
            MPI_Bcast(&spline_vals, 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        
            if(rank != 0){
                atmo::c_spline.x_vals[n] = spline_vals[0];      atmo::c_spline.f_vals[n] = spline_vals[1];
                atmo::u_spline.x_vals[n] = spline_vals[0];      atmo::u_spline.f_vals[n] = spline_vals[2];
                atmo::v_spline.x_vals[n] = spline_vals[0];      atmo::v_spline.f_vals[n] = spline_vals[3];
                atmo::rho_spline.x_vals[n] = spline_vals[0];    atmo::rho_spline.f_vals[n] = spline_vals[4];
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    
        interp::set(atmo::c_spline);    interp::set(atmo::u_spline);
        interp::set(atmo::rho_spline);  interp::set(atmo::v_spline);
    
        geoac::set_limits();
        topo::set_bndlyr();

        if(rank == 0){
            cout << '\t' << "Propagation region limits:" << '\n';
            cout << '\t' << '\t' << "latitude = " << geoac::lat_min * (180.0 / Pi) << ", " << geoac::lat_max * (180.0 / Pi) << '\n';
            cout << '\t' << '\t' << "longitude = " << geoac::lon_min * (180.0 / Pi) << ", " << geoac::lon_max  * (180.0 / Pi) << '\n';
            cout << '\t' << '\t' << "altitude = " << topo::z0 << ", " << geoac::alt_max << '\n' << '\n';
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    return result;
}


int set_region(char* atmo_file, char* topo_file, char* atmo_format, bool invert_winds, int rank){
    int result = 1;
    if(rank == 0){
        cout << "Interpolating atmosphere data in '" << atmo_file << "' and topography data in '" << topo_file << "'..." << '\n';
        if(invert_winds){
            cout << '\t' << "Inverting wind fields for back projection analysis..." << '\n';
        }
        
        int n1, n2;
        ifstream file_in;
    
        file_2d_dims(topo_file, n1, n2);
        interp::prep(topo::spline, n1, n2);
    
        file_in.open(topo_file);
        if(!file_in.is_open()){
            cout << '\t' << "ERROR: Invalid terrain file (" << topo_file << ")" << '\n' << '\n';
            result = 0;
        } else {
            for (int n1 = 0; n1 < topo::spline.length_x; n1++){
                for (int n2 = 0; n2 < topo::spline.length_y; n2++){
                    file_in >> topo::spline.x_vals[n1];
                    file_in >> topo::spline.y_vals[n2];
                    file_in >> topo::spline.f_vals[n1][n2];
                }
            }
        }
        file_in.close();
    
        for (int n = 0; n < topo::spline.length_x; n++){ topo::spline.x_vals[n] *= (Pi / 180.0);}
        for (int n = 0; n < topo::spline.length_y; n++){ topo::spline.y_vals[n] *= (Pi / 180.0);}

        interp::prep(atmo::c_spline,    file_length(atmo_file));
        interp::prep(atmo::u_spline,    atmo::c_spline.length);
        interp::prep(atmo::v_spline,    atmo::c_spline.length);
        interp::prep(atmo::rho_spline,  atmo::c_spline.length);

        double temp;
        string line;
        n1 = 0;

        file_in.open(atmo_file);
        if(!file_in.is_open()){
            cout << '\t' << "ERROR: Invalid atmospheric specification (" << atmo_file << ")" << '\n' << '\n';
            result = 0;
        } else {
            while(!file_in.eof() && n1 < atmo::c_spline.length){
                getline (file_in, line);
                if(line.find("#") != 0){
                    stringstream ss(line);
                    if(strncmp(atmo_format, "zTuvdp", 6) == 0){
                        ss >> atmo::c_spline.x_vals[n1];    // Extract z_i value
                        ss >> temp;                        // Extract T(z_i) but don't store it
                        ss >> atmo::u_spline.f_vals[n1];    // Extract u(z_i)
                        ss >> atmo::v_spline.f_vals[n1];    // Extract v(z_i)
                        ss >> atmo::rho_spline.f_vals[n1];  // Extract rho(z_i)
                        ss >> atmo::c_spline.f_vals[n1];    // Extract p(z_i) into c(z_i) and convert below
                    } else if (strncmp(atmo_format, "zuvwTdp", 7) == 0){
                        ss >> atmo::c_spline.x_vals[n1];    // Extract z_i value
                        ss >> atmo::u_spline.f_vals[n1];    // Extract u(z_i)
                        ss >> atmo::v_spline.f_vals[n1];    // Extract v(z_i)
                        ss >> temp;                        // Extract w(z_i) but don't store it
                        ss >> temp;                        // Extract T(z_i) but don't store it
                        ss >> atmo::rho_spline.f_vals[n1];  // Extract rho(z_i)
                        ss >> atmo::c_spline.f_vals[n1];    // Extract p(z_i) into c(z_i) and convert below
                    } else if (strncmp(atmo_format, "zcuvd", 5) == 0){
                        ss >> atmo::c_spline.x_vals[n1];    // Extract z_i value
                        ss >> atmo::c_spline.f_vals[n1];    // Extract c(z_i)
                        ss >> atmo::u_spline.f_vals[n1];    // Extract u(z_i)
                        ss >> atmo::v_spline.f_vals[n1];    // Extract v(z_i)
                        ss >> atmo::rho_spline.f_vals[n1];  // Extract rho(z_i)
                    } else {
                        cout << "Unrecognized profile option: " << atmo_format << ".  Valid options are: zTuvdp, zuvwTdp, or zcuvd" << '\n';
                        break;
                    }
        
                    // Copy altitude values to other interpolations
                    atmo::u_spline.x_vals[n1] = atmo::c_spline.x_vals[n1];
                    atmo::v_spline.x_vals[n1] = atmo::c_spline.x_vals[n1];
                    atmo::rho_spline.x_vals[n1] = atmo::c_spline.x_vals[n1];
        
                    // Convert pressure and density to adiabatic sound speed unless c is specified and scale winds from m/s to km/s
                    if (strncmp(atmo_format, "zTuvdp", 6) == 0 || strncmp(atmo_format, "zuvwTdp", 7) == 0){
                        atmo::c_spline.f_vals[n1] = sqrt(0.1 * atmo::gam * atmo::c_spline.f_vals[n1] / atmo::rho_spline.f_vals[n1]) / 1000.0;
                    } else {
                        atmo::c_spline.f_vals[n1] /= 1000.0;
                    }
        
                    if(invert_winds){
                        atmo::u_spline.f_vals[n1] /= -1000.0;
                        atmo::v_spline.f_vals[n1] /= -1000.0;
                    } else {
                        atmo::u_spline.f_vals[n1] /= 1000.0;
                        atmo::v_spline.f_vals[n1] /= 1000.0;
                    }
                    n1++;
                }
            }
        }
        file_in.close();
    }
    MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(result == 1){    
        int spline_len [2];
        double spline_vals [5];
    
        // Broadcast topography spline
        if(rank == 0){
            spline_len[0] = topo::spline.length_x;
            spline_len[1] = topo::spline.length_y;
        }
        MPI_Bcast(&spline_len, 2, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank != 0){
            interp::prep(topo::spline, spline_len[0], spline_len[1]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    
        for (int nx = 0; nx < topo::spline.length_x; nx++){
            for (int ny = 0; ny < topo::spline.length_y; ny++){
                if(rank == 0){
                    spline_vals[0] = topo::spline.x_vals[nx];
                    spline_vals[1] = topo::spline.y_vals[ny];
                    spline_vals[2] = topo::spline.f_vals[nx][ny];
                }
                MPI_Bcast(&spline_vals, 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);
                MPI_Barrier(MPI_COMM_WORLD);
                if(rank != 0){
                    topo::spline.x_vals[nx] = spline_vals[0];
                    topo::spline.y_vals[ny] = spline_vals[1];
                    topo::spline.f_vals[nx][ny] = spline_vals[2];
                }
                MPI_Barrier(MPI_COMM_WORLD);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    
        // Broadcast atmosphere spline
        if(rank == 0) spline_len[0] = atmo::c_spline.length;
        MPI_Bcast(&spline_len, 2, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank != 0){
            interp::prep(atmo::c_spline, spline_len[0]);
            interp::prep(atmo::u_spline, spline_len[0]);
            interp::prep(atmo::v_spline, spline_len[0]);
            interp::prep(atmo::rho_spline, spline_len[0]);
        }
    
        for(int n = 0; n < spline_len[0]; n++){
            if(rank == 0){
                spline_vals[0] = atmo::c_spline.x_vals[n];
                spline_vals[1] = atmo::c_spline.f_vals[n];
                spline_vals[2] = atmo::u_spline.f_vals[n];
                spline_vals[3] = atmo::v_spline.f_vals[n];
                spline_vals[4] = atmo::rho_spline.f_vals[n];
            }
            MPI_Bcast(&spline_vals, 5, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
        
            if(rank != 0){
                atmo::c_spline.x_vals[n] = spline_vals[0];      atmo::c_spline.f_vals[n] = spline_vals[1];
                atmo::u_spline.x_vals[n] = spline_vals[0];      atmo::u_spline.f_vals[n] = spline_vals[2];
                atmo::v_spline.x_vals[n] = spline_vals[0];      atmo::v_spline.f_vals[n] = spline_vals[3];
                atmo::rho_spline.x_vals[n] = spline_vals[0];    atmo::rho_spline.f_vals[n] = spline_vals[4];
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    
        interp::set(topo::spline);
        interp::set(atmo::c_spline);    interp::set(atmo::u_spline);
        interp::set(atmo::rho_spline);  interp::set(atmo::v_spline);
    
        geoac::set_limits();
        topo::set_bndlyr();
    
        if(rank == 0){
            cout << '\t' << "Propagation region limits:" << '\n';
            cout << '\t' << '\t' << "latitude = " << geoac::lat_min * (180.0 / Pi) << ", " << geoac::lat_max * (180.0 / Pi) << '\n';
            cout << '\t' << '\t' << "longitude = " << geoac::lon_min * (180.0 / Pi) << ", " << geoac::lon_max  * (180.0 / Pi) << '\n';
            cout << '\t' << '\t' << "altitude = " << topo::z0 << ", " << geoac::alt_max << '\n' << '\n';
    
            cout << '\t' << "Maximum topography height: " << topo::z_max << '\n';
            cout << '\t' << "Boundary layer height: " << topo::z_bndlyr << '\n';
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    return result;
}


void clear_region(){
    if(geoac::is_topo){
        interp::clear(topo::spline);
    }
    
    interp::clear(atmo::c_spline);      interp::clear(atmo::u_spline);
    interp::clear(atmo::rho_spline);    interp::clear(atmo::v_spline);
}

#endif /* _ATMO_IO_SPH_STRAT_MPI_CPP_ */
