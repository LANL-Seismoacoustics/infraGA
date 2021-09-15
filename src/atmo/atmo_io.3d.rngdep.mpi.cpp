# ifndef _ATMO_IO_3D_RNGDEP_MPI_CPP_
# define _ATMO_IO_3D_RNGDEP_MPI_CPP_

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <sstream>
#include <math.h>
#include <mpi.h>

#include "atmo_state.h"
#include "atmo_io.3d.rngdep.h"
#include "../util/interpolation.h"
#include "../util/fileIO.h"
#include "../geoac/geoac.params.h"

using namespace std;

//-----------------------------------------------//
//---------Define the Propagation Region---------//
//------------For 3D Planar Propgation-----------//
//-----------------------------------------------//
void geoac::set_limits(){
    double buffer = 0.001;
    x_min = atmo::c_spline.x_vals[0] + buffer;  x_max = atmo::c_spline.x_vals[atmo::c_spline.length_x - 1] - buffer;
    y_min = atmo::c_spline.y_vals[0] + buffer;  y_max = atmo::c_spline.y_vals[atmo::c_spline.length_y - 1] - buffer;
    
    if(is_topo){
        x_min = max(x_min, topo::spline.x_vals[0] + buffer); x_max = min(x_max, topo::spline.x_vals[topo::spline.length_x - 1] - buffer);
        y_min = max(y_min, topo::spline.y_vals[0] + buffer); y_max = min(y_max, topo::spline.y_vals[topo::spline.length_y - 1] - buffer);
        if(x_max < x_min || y_max < y_min){
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
int set_region(char* atmo_prefix, char* atmo_locs_x, char* atmo_locs_y, char* atmo_format, bool invert_winds, int rank){
    if(rank == 0){
        cout << "Interpolating atmosphere data in '" << atmo_prefix << "'* using format '" << atmo_format << "'..." << '\n';
        int x_cnt, y_cnt, z_cnt;
        double temp;
        char output_buffer [512];
        string line;
        ifstream file_in;
    
        x_cnt = file_length(atmo_locs_x);
        y_cnt = file_length(atmo_locs_y);
    
        sprintf(output_buffer, "%s%i.met", atmo_prefix, 0);
        z_cnt = file_length(output_buffer);
    
        interp::prep(atmo::c_spline, x_cnt, y_cnt, z_cnt);
        interp::prep(atmo::u_spline, x_cnt, y_cnt, z_cnt);
        interp::prep(atmo::v_spline, x_cnt, y_cnt, z_cnt);
        interp::prep(atmo::rho_spline, x_cnt, y_cnt, z_cnt);
    
        file_in.open(atmo_locs_x);
        for (int nx = 0; nx < x_cnt; nx++){
            file_in >> atmo::c_spline.x_vals[nx];
            atmo::u_spline.x_vals[nx] = atmo::c_spline.x_vals[nx];
            atmo::v_spline.x_vals[nx] = atmo::c_spline.x_vals[nx];
            atmo::rho_spline.x_vals[nx] = atmo::c_spline.x_vals[nx];
        }
        file_in.close();
    
        file_in.open(atmo_locs_y);
        for (int ny = 0; ny < y_cnt; ny++){
            file_in >> atmo::c_spline.y_vals[ny];
            atmo::u_spline.y_vals[ny] = atmo::c_spline.y_vals[ny];
            atmo::v_spline.y_vals[ny] = atmo::c_spline.y_vals[ny];
            atmo::rho_spline.y_vals[ny] = atmo::c_spline.y_vals[ny];
        }
        file_in.close();
    
        for(int nx = 0; nx < x_cnt; nx++){
        for(int ny = 0; ny < y_cnt; ny++){
            sprintf(output_buffer, "%s%i.met", atmo_prefix, nx * y_cnt + ny);
            if((nx < 3) && (ny < 3)){
                cout << '\t' << "Setting grid node at (" << atmo::c_spline.x_vals[nx] << ", " << atmo::c_spline.y_vals[ny];
                cout << ") with profile " << output_buffer << '\n';
            } else if((nx < 3) && (ny == 3)) {
                cout << '\t' << "..." << '\n';
            }
        
            int nz = 0;
            file_in.open(output_buffer);
            while(!file_in.eof() && nz < z_cnt){
                getline (file_in, line);
                if(line.find("#") != 0){
                    stringstream ss(line);
                    if(strncmp(atmo_format, "zTuvdp", 6) == 0){
                        ss >> atmo::c_spline.z_vals[nz];                // Extract z_i value
                        ss >> temp;                                     // Extract T(z_i) but don't store it
                        ss >> atmo::u_spline.f_vals[nx][ny][nz];        // Extract u(z_i)
                        ss >> atmo::v_spline.f_vals[nx][ny][nz];        // Extract v(z_i)
                        ss >> atmo::rho_spline.f_vals[nx][ny][nz];      // Extract rho(z_i)
                        ss >> atmo::c_spline.f_vals[nx][ny][nz];        // Extract p(z_i) into c(z_i) and convert below
                    } else if (strncmp(atmo_format, "zuvwTdp", 7) == 0){
                        ss >> atmo::c_spline.z_vals[nz];                // Extract z_i value
                        ss >> atmo::u_spline.f_vals[nx][ny][nz];        // Extract u(z_i)
                        ss >> atmo::v_spline.f_vals[nx][ny][nz];        // Extract v(z_i)
                        ss >> temp;                                     // Extract w(z_i) but don't store it
                        ss >> temp;                                     // Extract T(z_i) but don't store it
                        ss >> atmo::rho_spline.f_vals[nx][ny][nz];      // Extract rho(z_i)
                        ss >> atmo::c_spline.f_vals[nx][ny][nz];        // Extract p(z_i) into c(z_i) and convert below
                    } else if (strncmp(atmo_format, "zcuvd", 5) == 0){
                        ss >> atmo::c_spline.z_vals[nz];                // Extract z_i value
                        ss >> atmo::c_spline.f_vals[nx][ny][nz];        // Extract c(z_i)
                        ss >> atmo::u_spline.f_vals[nx][ny][nz];        // Extract u(z_i)
                        ss >> atmo::v_spline.f_vals[nx][ny][nz];        // Extract v(z_i)
                        ss >> atmo::rho_spline.f_vals[nx][ny][nz];      // Extract rho(z_i)
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
                        atmo::c_spline.f_vals[nx][ny][nz] = sqrt(0.1 * atmo::gam * atmo::c_spline.f_vals[nx][ny][nz] / atmo::rho_spline.f_vals[nx][ny][nz]) / 1000.0;
                    } else {
                        atmo::c_spline.f_vals[nx][ny][nz] /= 1000.0;
                    }
            
                    if(invert_winds){
                        atmo::u_spline.f_vals[nx][ny][nz] /= -1000.0;
                        atmo::v_spline.f_vals[nx][ny][nz] /= -1000.0;
                    } else {
                        atmo::u_spline.f_vals[nx][ny][nz] /= 1000.0;
                        atmo::v_spline.f_vals[nx][ny][nz] /= 1000.0;
                    }
                    nz++;
                }
            }
            file_in.close();
        }}
        cout << '\n';
    }
    MPI_Barrier(MPI_COMM_WORLD);

    int spline_len [3];
    double spline_vals [7];
    
    if(rank == 0){
        spline_len[0] = atmo::c_spline.length_x;
        spline_len[1] = atmo::c_spline.length_y;
        spline_len[2] = atmo::c_spline.length_z;
    }
    MPI_Bcast(&spline_len, 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank != 0){
        interp::prep(atmo::c_spline, spline_len[0], spline_len[1], spline_len[2]);
        interp::prep(atmo::u_spline, spline_len[0], spline_len[1], spline_len[2]);
        interp::prep(atmo::v_spline, spline_len[0], spline_len[1], spline_len[2]);
        interp::prep(atmo::rho_spline, spline_len[0], spline_len[1], spline_len[2]);
    }

    for(int nx = 0; nx < spline_len[0]; nx++){
    for(int ny = 0; ny < spline_len[1]; ny++){
        for(int nz = 0; nz < spline_len[2]; nz++){
            if(rank == 0){
                spline_vals[0] = atmo::c_spline.x_vals[nx];
                spline_vals[1] = atmo::c_spline.y_vals[ny];
                spline_vals[2] = atmo::c_spline.z_vals[nz];
                spline_vals[3] = atmo::c_spline.f_vals[nx][ny][nz];
                spline_vals[4] = atmo::u_spline.f_vals[nx][ny][nz];
                spline_vals[5] = atmo::v_spline.f_vals[nx][ny][nz];
                spline_vals[6] = atmo::rho_spline.f_vals[nx][ny][nz];
            }
            MPI_Bcast(&spline_vals, 7, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            
            if(rank != 0){
                atmo::c_spline.x_vals[nx] = spline_vals[0];     atmo::c_spline.y_vals[ny] = spline_vals[1];     atmo::c_spline.z_vals[nz] = spline_vals[2];     atmo::c_spline.f_vals[nx][ny][nz] = spline_vals[3];
                atmo::u_spline.x_vals[nx] = spline_vals[0];     atmo::u_spline.y_vals[ny] = spline_vals[1];     atmo::u_spline.z_vals[nz] = spline_vals[2];     atmo::u_spline.f_vals[nx][ny][nz] = spline_vals[4];
                atmo::v_spline.x_vals[nx] = spline_vals[0];     atmo::v_spline.y_vals[ny] = spline_vals[1];     atmo::v_spline.z_vals[nz] = spline_vals[2];     atmo::v_spline.f_vals[nx][ny][nz] = spline_vals[5];
                atmo::rho_spline.x_vals[nx] = spline_vals[0];   atmo::rho_spline.y_vals[ny] = spline_vals[1];   atmo::rho_spline.z_vals[nz] = spline_vals[2];   atmo::rho_spline.f_vals[nx][ny][nz] = spline_vals[6];
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }}
    
    interp::set(atmo::c_spline);    interp::set(atmo::u_spline);
    interp::set(atmo::rho_spline);  interp::set(atmo::v_spline);
    
    geoac::set_limits();
    topo::set_bndlyr();

    if(rank == 0){
        cout << '\t' << "Propagation region limits:" << '\n';
        cout << '\t' << '\t' << "x = " << geoac::x_min << ", " << geoac::x_max << '\n';
        cout << '\t' << '\t' << "y = " << geoac::y_min << ", " << geoac::y_max << '\n';
        cout << '\t' << '\t' << "z = " << topo::z0 << ", " << geoac::alt_max << '\n';
    }
}


int set_region(char* atmo_prefix, char* atmo_locs_x, char* atmo_locs_y, char* topo_file, char* atmo_format, bool invert_winds, int rank){
    if(rank == 0){
        cout << "Interpolating atmosphere data in '" << atmo_prefix << "#' and topography data in '" << topo_file << "'..." << '\n';
    
        int x_cnt, y_cnt, z_cnt;
        double temp;
        char output_buffer [512];
        string line;
        ifstream file_in;
    
        // Set up topography interpolation
        file_2d_dims(topo_file, x_cnt, y_cnt);
        interp::prep(topo::spline, x_cnt, y_cnt);
    
        file_in.open(topo_file);
        for (int nx = 0; nx < topo::spline.length_x; nx++){
        for (int ny = 0; ny < topo::spline.length_y; ny++){
            file_in >> topo::spline.x_vals[nx];
            file_in >> topo::spline.y_vals[ny];
            file_in >> topo::spline.f_vals[nx][ny];
        }}
        file_in.close();
    
        x_cnt = file_length(atmo_locs_x);
        y_cnt = file_length(atmo_locs_y);
    
        sprintf(output_buffer, "%s%i.met", atmo_prefix, 0);
        z_cnt = file_length(output_buffer);
    
        interp::prep(atmo::c_spline, x_cnt, y_cnt, z_cnt);
        interp::prep(atmo::u_spline, x_cnt, y_cnt, z_cnt);
        interp::prep(atmo::v_spline, x_cnt, y_cnt, z_cnt);
        interp::prep(atmo::rho_spline, x_cnt, y_cnt, z_cnt);
    
        file_in.open(atmo_locs_x);
        for (int nx = 0; nx < x_cnt; nx++){
            file_in >> atmo::c_spline.x_vals[nx];
            atmo::u_spline.x_vals[nx] = atmo::c_spline.x_vals[nx];
            atmo::v_spline.x_vals[nx] = atmo::c_spline.x_vals[nx];
            atmo::rho_spline.x_vals[nx] = atmo::c_spline.x_vals[nx];
        }
        file_in.close();
    
        file_in.open(atmo_locs_y);
        for (int ny = 0; ny < y_cnt; ny++){
            file_in >> atmo::c_spline.y_vals[ny];
            atmo::u_spline.y_vals[ny] = atmo::c_spline.y_vals[ny];
            atmo::v_spline.y_vals[ny] = atmo::c_spline.y_vals[ny];
            atmo::rho_spline.y_vals[ny] = atmo::c_spline.y_vals[ny];
        }
        file_in.close();
    
        for(int nx = 0; nx < x_cnt; nx++){
        for(int ny = 0; ny < y_cnt; ny++){
            sprintf(output_buffer, "%s%i.met", atmo_prefix, nx * y_cnt + ny);
            if((nx < 3) && (ny < 3)){
                cout << '\t' << "Setting grid node at (" << atmo::c_spline.x_vals[nx] << ", " << atmo::c_spline.y_vals[ny];
                cout << ") with profile " << output_buffer << '\n';
            } else if((nx < 3) && (ny == 3)) {
                cout << '\t' << "..." << '\n';
            }
            
            int nz = 0;
            file_in.open(output_buffer);
            while(!file_in.eof() && nz < z_cnt){
                getline (file_in, line);
                if(line.find("#") != 0){
                    stringstream ss(line);
                    if(strncmp(atmo_format, "zTuvdp", 6) == 0){
                        ss >> atmo::c_spline.z_vals[nz];        // Extract z_i value
                        ss >> temp;                                     // Extract T(z_i) but don't store it
                        ss >> atmo::u_spline.f_vals[nx][ny][nz];        // Extract u(z_i)
                        ss >> atmo::v_spline.f_vals[nx][ny][nz];        // Extract v(z_i)
                        ss >> atmo::rho_spline.f_vals[nx][ny][nz];      // Extract rho(z_i)
                        ss >> atmo::c_spline.f_vals[nx][ny][nz];        // Extract p(z_i) into c(z_i) and convert below
                    } else if (strncmp(atmo_format, "zuvwTdp", 7) == 0){
                        ss >> atmo::c_spline.z_vals[nz];        // Extract z_i value
                        ss >> atmo::u_spline.f_vals[nx][ny][nz];        // Extract u(z_i)
                        ss >> atmo::v_spline.f_vals[nx][ny][nz];        // Extract v(z_i)
                        ss >> temp;                                     // Extract w(z_i) but don't store it
                        ss >> temp;                                     // Extract T(z_i) but don't store it
                        ss >> atmo::rho_spline.f_vals[nx][ny][nz];      // Extract rho(z_i)
                        ss >> atmo::c_spline.f_vals[nx][ny][nz];        // Extract p(z_i) into c(z_i) and convert below
                    } else if (strncmp(atmo_format, "zcuvd", 5) == 0){
                        ss >> atmo::c_spline.z_vals[nz];        // Extract z_i value
                        ss >> atmo::c_spline.f_vals[nx][ny][nz];        // Extract c(z_i)
                        ss >> atmo::u_spline.f_vals[nx][ny][nz];        // Extract u(z_i)
                        ss >> atmo::v_spline.f_vals[nx][ny][nz];        // Extract v(z_i)
                        ss >> atmo::rho_spline.f_vals[nx][ny][nz];      // Extract rho(z_i)
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
                        atmo::c_spline.f_vals[nx][ny][nz] = sqrt(0.1 * atmo::gam * atmo::c_spline.f_vals[nx][ny][nz] / atmo::rho_spline.f_vals[nx][ny][nz]) / 1000.0;
                    } else {
                        atmo::c_spline.f_vals[nx][ny][nz] /= 1000.0;
                    }
                
                    if(invert_winds){
                        atmo::u_spline.f_vals[nx][ny][nz] /= -1000.0;
                        atmo::v_spline.f_vals[nx][ny][nz] /= -1000.0;
                    } else {
                        atmo::u_spline.f_vals[nx][ny][nz] /= 1000.0;
                        atmo::v_spline.f_vals[nx][ny][nz] /= 1000.0;
                    }
                    nz++;
                }
            }
            file_in.close();
        }}
        cout << '\n';
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    int spline_len [3];
    double spline_vals [7];
    
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
    }}
    MPI_Barrier(MPI_COMM_WORLD);
    
    // Broadcast atmosphere spline
    if(rank == 0){
        spline_len[0] = atmo::c_spline.length_x;
        spline_len[1] = atmo::c_spline.length_y;
        spline_len[2] = atmo::c_spline.length_z;
    }
    MPI_Bcast(&spline_len, 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank != 0){
        interp::prep(atmo::c_spline, spline_len[0], spline_len[1], spline_len[2]);
        interp::prep(atmo::u_spline, spline_len[0], spline_len[1], spline_len[2]);
        interp::prep(atmo::v_spline, spline_len[0], spline_len[1], spline_len[2]);
        interp::prep(atmo::rho_spline, spline_len[0], spline_len[1], spline_len[2]);
    }
    
    if(rank == 0){
        spline_len[0] = atmo::c_spline.length_x;
        spline_len[1] = atmo::c_spline.length_y;
        spline_len[2] = atmo::c_spline.length_z;
    }
    MPI_Bcast(&spline_len, 3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank != 0){
        interp::prep(atmo::c_spline, spline_len[0], spline_len[1], spline_len[2]);
        interp::prep(atmo::u_spline, spline_len[0], spline_len[1], spline_len[2]);
        interp::prep(atmo::v_spline, spline_len[0], spline_len[1], spline_len[2]);
        interp::prep(atmo::rho_spline, spline_len[0], spline_len[1], spline_len[2]);
    }
    
    for(int nx = 0; nx < spline_len[0]; nx++){
    for(int ny = 0; ny < spline_len[1]; ny++){
        for(int nz = 0; nz < spline_len[2]; nz++){
            if(rank == 0){
                spline_vals[0] = atmo::c_spline.x_vals[nx];
                spline_vals[1] = atmo::c_spline.y_vals[ny];
                spline_vals[2] = atmo::c_spline.z_vals[nz];
                spline_vals[3] = atmo::c_spline.f_vals[nx][ny][nz];
                spline_vals[4] = atmo::u_spline.f_vals[nx][ny][nz];
                spline_vals[5] = atmo::v_spline.f_vals[nx][ny][nz];
                spline_vals[6] = atmo::rho_spline.f_vals[nx][ny][nz];
            }
            MPI_Bcast(&spline_vals, 7, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
                
            if(rank != 0){
                atmo::c_spline.x_vals[nx] = spline_vals[0];     atmo::c_spline.y_vals[ny] = spline_vals[1];     atmo::c_spline.z_vals[nz] = spline_vals[2];     atmo::c_spline.f_vals[nx][ny][nz] = spline_vals[3];
                atmo::u_spline.x_vals[nx] = spline_vals[0];     atmo::u_spline.y_vals[ny] = spline_vals[1];     atmo::u_spline.z_vals[nz] = spline_vals[2];     atmo::u_spline.f_vals[nx][ny][nz] = spline_vals[4];
                atmo::v_spline.x_vals[nx] = spline_vals[0];     atmo::v_spline.y_vals[ny] = spline_vals[1];     atmo::v_spline.z_vals[nz] = spline_vals[2];     atmo::v_spline.f_vals[nx][ny][nz] = spline_vals[5];
                atmo::rho_spline.x_vals[nx] = spline_vals[0];   atmo::rho_spline.y_vals[ny] = spline_vals[1];   atmo::rho_spline.z_vals[nz] = spline_vals[2];   atmo::rho_spline.f_vals[nx][ny][nz] = spline_vals[6];
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }}
    
    interp::set(topo::spline);
    interp::set(atmo::c_spline);    interp::set(atmo::u_spline);
    interp::set(atmo::rho_spline);  interp::set(atmo::v_spline);
    
    geoac::set_limits();
    topo::set_bndlyr();

    if(rank == 0){
        cout << '\t' << "Propagation region limits:" << '\n';
        cout << '\t' << '\t' << "x = " << geoac::x_min << ", " << geoac::x_max << '\n';
        cout << '\t' << '\t' << "y = " << geoac::y_min << ", " << geoac::y_max << '\n';
        cout << '\t' << '\t' << "z = " << topo::z0  << ", " << geoac::alt_max << '\n' << '\n';
    
        cout << '\t' << "Maximum topography height: " << topo::z_max << '\n';
        cout << '\t' << "Boundary layer height: " << topo::z_bndlyr << '\n';
    }
    MPI_Barrier(MPI_COMM_WORLD);

}


void clear_region(){
    if(geoac::is_topo){
        interp::clear(topo::spline);
    }
    
    interp::clear(atmo::c_spline);      interp::clear(atmo::u_spline);
    interp::clear(atmo::rho_spline);    interp::clear(atmo::v_spline);
}


#endif /*  _ATMO_IO_3D_RNGDEP_MPI_CPP_ */
