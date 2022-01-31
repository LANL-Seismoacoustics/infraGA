#ifndef _GEOAC_PARAMS_H_
#define _GEOAC_PARAMS_H_

//----------------------------//
//---Mathematical Constants---//
//----------------------------//
extern double   Pi;

namespace  geoac{

    //----------------------------------------//
    //------Initial angles at the source------//
    //----------------------------------------//
    extern double   theta;      // Initial inclination angle of ray path
    extern double   phi;        // Initial azimuth angle of ray path

    extern int      eq_cnt;     // Number of equations to solve
    extern int      dim;        // Number of dimensions (2 or 3)

    //-----------------------------//
    //------Configuring GeoAc------//
    //-----------------------------//
    extern bool     calc_amp;   // Is amplitude to be calculated?
    extern bool     is_strat;   // Is the medium stratified?
    extern bool     is_topo;    // Is there topography?

    extern double   s_max;      // Ray length at which to stop ray tracing
    extern double   ds_min;     // Smallest possible step size
    extern double   ds_max;     // Largest possible step size
    extern double   ds_wvfrm;   // Step size for weakly non-linear waveform solver

    //-------------------------------------//
    //----Set Propagation Region Limits----//
    //-------------------------------------//
    extern double   alt_max;    // Altitude at which to stop ray tracing
    extern double 	rng_max;    // Range at which to stop ray tracing
    extern double   time_max;   // Maximum propagation time (for global scale simulations)

    extern double 	x_min;      // E-W range at which to stop ray tracing in 3D.RngDep
    extern double 	x_max;      // E-W range at which to stop ray tracing in 3D.RngDep
    extern double 	y_min;      // N-S range at which to stop ray tracing in 3D.RngDep
    extern double 	y_max;      // N-S range at which to stop ray tracing in 3D.RngDep

    extern double 	lat_min;    // Minimum latitude of propagation region in Global.RngDep
    extern double 	lat_max;    // Maximum latitude of propagation region in Global.RngDep
    extern double 	lon_min;    // Minimum longitude of propagation region in Global.RngDep
    extern double 	lon_max;    // Maximum longitude of propagation region in Global.RngDep
}
#endif /* _GEOAC_PARAMS_H_ */
