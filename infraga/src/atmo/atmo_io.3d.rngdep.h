# ifndef _ATMO_STATE_3D_RNGDEP_H_
# define _ATMO_STATE_3D_RNGDEP_H_

#include "../util/interpolation.h"

//-----------------------------------//
//-----Topography and Atmosphere-----//
//-----------Interpolations----------//
//-----------------------------------//
namespace topo {
    extern struct interp::natural_cubic_spline_2D spline;   // Topography interpolation
}

namespace atmo {
    extern struct interp::hybrid_spline_3D c_spline;       // Sound speed interpolation
    extern struct interp::hybrid_spline_3D u_spline;       // E-W wind interpolation
    extern struct interp::hybrid_spline_3D v_spline;       // N-S wind interpolation
    extern struct interp::hybrid_spline_3D rho_spline;     // Density interpolation
}

//-------------------------------------------//
//--------------Set Up or Clear--------------//
//---------Topography and Atmosphere---------//
//--------------Interpolations---------------//
//-------------------------------------------//
int set_region(char*, char*, char*, char*, bool);              // Input the atmosphere file prefix, limits, and atmosphere file format
int set_region(char*, char*, char*, char*, char*, bool);       // Input the atmosphere file prefix, limits, topography file, and atmosphere file format

int set_region(char*, char*, char*, char*, bool, int);         // Input the atmosphere file prefix, limits, and atmosphere file format using OpenMPI
int set_region(char*, char*, char*, char*, char*, bool, int);  // Input the atmosphere file prefix, limits, topography file, and atmosphere file format using OpenMPI

void clear_region();                                            // Clear the interpolations

#endif /* _ATMO_STATE_3D_RNGDEP_H_ */
