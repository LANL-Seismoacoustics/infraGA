#ifndef _GEOAC_INTERFACE_3D_CPP_
#define _GEOAC_INTERFACE_3D_CPP_

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "geoac.eqset.h"
#include "geoac.interface.h"
#include "geoac.params.h"

#include "../atmo/atmo_state.h"
#include "../util/rk4solver.h"

using namespace std;

//--------------------------------------//
//----Functions To Set Up The Solver----//
//--------------------------------------//
void geoac::configure(){
    if(calc_amp){   eq_cnt = 18;}   // Use auxiliary parameters to compute amplitude
    else {          eq_cnt = 6;}    // Don't use auxiliary parameters
}

//----------------------------------------------------------------//
//----Functions To Build, Clear, and Delete the Solution Array----//
//----------------------------------------------------------------//
void geoac::build_solution(double ** & solution, int length){
	solution = new double* [length];
	for(int n = 0; n < length; n++){ solution[n] = new double [eq_cnt];}
}

void geoac::clear_solution(double ** & solution, int index){
	for(int n1 = 0; n1 < index; n1++){ for(int n2 = 0; n2 < eq_cnt; n2++){ solution[n1][n2] = 0.0; }}
}

void geoac::delete_solution(double ** & solution, int length){
	for(int n = 0; n < length; n++){ delete [] solution[n];}
	delete [] solution;
}

//------------------------------------//
//---------Output The Profile---------//
//------------------------------------//
void geoac::write_prof(char* file_name, double x0, double y0, double azimuth){
	ofstream file_out;	file_out.open(file_name);
    
	if(!file_out.is_open()){
        cout << "Error opening file, check file name." << '\n';
    } else {
        for(int m = 0; m < alt_max * 10.; m++){
            double z0 = m / 10.0;

            file_out << z0 << '\t';
            file_out << atmo::c(x0, y0, z0) * 1000.0 << '\t';
            file_out << atmo::u(x0, y0, z0) * 1000.0 << '\t';
            file_out << atmo::v(x0, y0, z0) * 1000.0 << '\t';
            file_out << atmo::rho(x0, y0, z0) << '\t';
            file_out << (atmo::c(x0, y0, z0) + cos(azimuth) * atmo::u(x0, y0, z0) + sin(azimuth) * atmo::v(x0, y0, z0)) * 1000.0 << '\t';
            file_out << '\n';
        }
        file_out.close();
    }
}

#endif /* _GEOAC_INTERFACE_3D_CPP_ */
