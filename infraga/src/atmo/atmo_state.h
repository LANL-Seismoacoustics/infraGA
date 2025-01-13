# ifndef _ATMO_STATE_H_
# define _ATMO_STATE_H_

//--------------------------------------------//
//---------Set the Propagation Medium---------//
//--------------------------------------------//
namespace geoac{
    void set_limits();
}

//---------------------------------------//
//---------Define the Topography---------//
//---------------------------------------//
namespace topo{
    extern double z0;               // Difference between ground level and sea level (z = 0 in G2S file)
    extern double z_max;            // Maximum topographical altitude
    extern double z_bndlyr;         // Altitude of boundary layer (above this the topography is neglected)

    extern double bndlyr_param;     // Parameter controlling boundary layer thickness
    extern double vert_wind_grad;   // Boundary layer wind gradient for topography inclusion
    extern bool use_BLw;            // Control to use boundary layer vertical winds

    
    void set_bndlyr();      // Define the topography maximum and boundary layer height
    
    double bndlyr_sc(double, double, double);               // Scaling function to force a no-slip BC at the ground
    double bndlyr_dsc(double, double, double, int);         // Scaling function first derivative
    double bndlyr_ddsc(double, double, double, int, int);   // Scaling function second derivative
    
    double z(double);       // Function defining topography (range dependent)
    double dz(double);      // Function defining topography gradient (range dependent)
    double ddz(double);     // Function defining topography curvature (range dependent)

    double z(double, double);                       // Function defining topography (x,y or lat,lon dependent)
    double dz(double, double, int);                 // Function defining topography gradients (x,y or lat,lon dependent)
    double ddz(double, double, int, int);           // Function defining topography curvatures (x,y or lat,lon dependent)
    double dddz(double, double, int, int, int);     // Function defining topography curvature derivative(x,y or lat,lon dependent)
}
    
//---------------------------------------------------//
//---------Functions Defining the Atmosphere---------//
//---------------------------------------------------//
namespace atmo{
    extern double   gam;
    extern double   R;

    extern double   z_reflect;                  // Altitude to impose partial reflection

    double rho(double, double, double);         // Atmospheric density
    
    double c(double, double, double);           // Adiabatic sound speed and derivatives
    double dc(double, double, double, int);
    double ddc(double, double, double, int, int);

    double u(double, double, double);           // E-W component of the winds and derivatives
    double du(double, double, double, int);
    double ddu(double, double, double, int, int);

    double v(double, double, double);           // N-S component of the winds and derivatives
    double dv(double, double, double, int);
    double ddv(double, double, double, int, int);

    double w(double, double, double);           // Vertical component of the winds and derivatives
    double dw(double, double, double, int);
    double ddw(double, double, double, int, int);

    void calc_uvw(double, double, double, double &, double &, double &, double [], double [], double []);
    void calc_uvw(double, double, double, double &, double &, double &, double [], double [], double [], double [], double [], double []);
    
    //----------------------------------------------------//
    //---------Functions Defining the Attenuation---------//
    //----------------------------------------------------//
    extern double tweak_abs;                                // Coefficient to scale Sutherland Bass attenuation
    double SB_alpha(double, double, double, double);  // Attenuation model from Sutherland Bass
}

#endif  /* _ATMO_STATE_H_ */
