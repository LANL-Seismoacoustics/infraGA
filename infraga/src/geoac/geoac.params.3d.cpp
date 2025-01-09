#ifndef _GEOAC_PARAMS_3D_CPP_
#define _GEOAC_PARAMS_3D_CPP_

#include "geoac.params.h"

//----------------------------//
//---Mathematical Constants---//
//----------------------------//
double Pi =  3.141592653589793238462643;

//----------------------------------------//
//------Initial angles at the source------//
//----------------------------------------//
double  geoac::theta;    // Initial inclination angle of ray path
double  geoac::phi;      // Initial azimuth angle of ray path

//-----------------------------//
//------Configuring GeoAc------//
//-----------------------------//
int     geoac::eq_cnt;      // Number of equations to solve
int     geoac::dim;         // Number of dimensions (2 or 3)

bool    geoac::calc_amp;    // Is amplitude to be calculated?
bool    geoac::is_strat;    // Is the medium stratified?
bool    geoac::is_topo;     // Is there topography?

double  geoac::s_max = 1000.0;  // Ray length at which to stop ray tracing
double  geoac::ds_min = 0.001;  // Smallest possible step size
double  geoac::ds_max = 0.05;   // Largest possible step size
double  geoac::ds_wvfrm = 0.5;  // Step size in weakly non-linear waveform solver

//-------------------------------------//
//----Set Propagation Region Limits----//
//-------------------------------------//
double  geoac::alt_max = 140.0;     // Altitude at which to stop ray tracing
double  geoac::rng_max = 1000.0;    // Range at which to stop ray tracing

double 	geoac::x_min = -1000.0;     // E-W range at which to stop ray tracing
double 	geoac::x_max =  1000.0;     // E-W range at which to stop ray tracing
double 	geoac::y_min = -1000.0;     // N-S range at which to stop ray tracing
double 	geoac::y_max =  1000.0;     // N-S range at which to stop ray tracing

#endif /* _GEOAC_PARAMS_3D_CPP_ */
