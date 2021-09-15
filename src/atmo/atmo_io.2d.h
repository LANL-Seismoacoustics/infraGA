# ifndef _ATMO_STATE_2D_H_
# define _ATMO_STATE_2D_H_

#include "../util/interpolation.h"

//-----------------------------------//
//-----Topography and Atmosphere-----//
//-----------Interpolations----------//
//-----------------------------------//
namespace topo {
    extern struct interp::natural_cubic_spline_1D spline;       // Topography interpolation
}

namespace atmo {
    extern struct interp::natural_cubic_spline_1D rho_spline;   // Density interpolation
    extern struct interp::natural_cubic_spline_1D c_spline;     // Temperature interpolation
    extern struct interp::natural_cubic_spline_1D u_spline;     // E-W wind interpolation
    extern struct interp::natural_cubic_spline_1D v_spline;     // N-S wind interpolation
}

//-----------------------------------------------//
//----------Set Up or Clear Topography-----------//
//---------and Atmosphere Interpolations---------//
//-----------------------------------------------//
int set_region(char*, char*, bool);        // Read in the atmosphere file and set up interpolations
int set_region(char*, char*, char*, bool); // Read in the atmosphere and topography files and set up interpolations
void clear_region();                        // Clear the interpolations

#endif /* _ATMO_STATE_2D_H_ */
