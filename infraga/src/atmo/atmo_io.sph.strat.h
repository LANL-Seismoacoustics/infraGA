# ifndef _ATMO_STATE_SPH_STRAT_H_
# define _ATMO_STATE_SPH_STRAT_H_

#include "../util/interpolation.h"

//-----------------------------------//
//-----Topography and Atmosphere-----//
//-----------Interpolations----------//
//-----------------------------------//
namespace topo {
    void set_bndlyr();                                      // Set the boundary layer altitude
    
    extern double bndlyr_param;                             // Parameter controlling boundary layer thickness
    extern struct interp::natural_cubic_spline_2D spline;   // Topography interpolation

    double bndlyr_sc(double, double, double);               // Scaling function to force a no-slip BC at the ground
    double bndlyr_dsc(double, double, double, int);         // Scaling function first derivative
    double bndlyr_ddsc(double, double, double, int, int);   // Scaling function second derivative

}

namespace atmo {
    extern struct interp::natural_cubic_spline_1D c_spline;       // Sound speed interpolation
    extern struct interp::natural_cubic_spline_1D u_spline;       // E-W wind interpolation
    extern struct interp::natural_cubic_spline_1D v_spline;       // N-S wind interpolation
    extern struct interp::natural_cubic_spline_1D rho_spline;     // Density interpolation
}

//-------------------------------------------//
//--------------Set Up or Clear--------------//
//---------Topography and Atmosphere---------//
//--------------Interpolations---------------//
//-------------------------------------------//
int set_region(char*, char*, bool);                // Read in the atmosphere file and set up interpolations
int set_region(char*, char*, char*, bool);         // Read in the atmosphere and topography files and set up interpolations

int set_region(char*, char*, bool, int);           // Read in the atmosphere file and set up interpolations using OpenMPI
int set_region(char*, char*, char*, bool, int);    // Read in the atmosphere and topography files and set up interpolations using OpenMPI

void clear_region();                                // Clear the interpolations


#endif /* _ATMO_STATE_SPH_STRAT_H_ */
