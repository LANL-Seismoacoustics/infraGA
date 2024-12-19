#ifndef _GEOAC_EQSET_H_
#define _GEOAC_EQSET_H_

namespace  geoac{
    void    set_system();                                       // Set dimensions and atmosphere type
    
    void    set_initial(double**&, double, double);             // Set initial conditions for a source at (r_0, z_0)
    void    set_initial(double**&, double, double, double);     // Set initial conditions for a source at (x_0, y_0, z_0) or (r0, lat0, lon0)
    
    void    approx_intercept(double**, int, double*&);          // Uses linear interpolation to estimate ground intercept values
    void    set_refl(double**&, int);                           // Set reflection conditions

    void    update_refs(double, double*);                       // Update values referenced in the equation evaluation
    double  set_ds(double*);                                    // Modify the solver step size
    double  eval_src_eq(double, double*, int);                  // Evaluate the source equations

    double  eval_eikonal(double, double*);                      // Evaluate the eikonal residual at s, current_values
    double  eval_eikonal(double**, int);                        // Evaluate the eikonal residual at index of solution**
    double  eval_eikonal_deriv(double**, int);                  // Evaluate the launch angle derivatives of the eikonal
    
    bool    break_check(double ** &, int);                      // Check for ray leaving the propagation region
    bool    ground_check(double **, int);                       // Check for ray returning to ground
    bool    reflect_check(double **, int);                      // Check for ray reaching partial reflection altitude

    double  travel_time(double **, int);                                        // Calculate travel time from ray origin to s = ds * index
    void    travel_time(double&, double**, int, int);                           // Increment travel time from ds*k_1 to d2*k_2

    void    travel_time_var(double **, int, double&, double&, double&);         // Calculate travel time with variances from origin to k^{th} step
    void    travel_time_var(double&, double&, double&, double**, int, int);     // Increment travel time and variances from step k1 to step k2

    double  atten(double **, int, double);                                      // Calculate atmospheric attenuation through a ray path
    void    atten(double&, double **, int , int, double);                       // Increment atmospheric attenuation through a ray path
    double  jacobian(double **, int);                                           // Calculate the Jacobian determinant
    double  amp(double **, int);                                                // Calculate the transport equation coefficient relative to spherical spreading
    double  est_dev(double **, int);                                            // Estimate the non-planar deviation in the back aizmuth for ESS.

    double    wnl_wvfrm(double **, double ** &, int, int, double, double, double, double, double, double);                   // Compute weakly non-linear prediction for waveform along ray path in 2D
    double    wnl_wvfrm(double **, double ** &, int, int, double, double, double, double, double, double, double, double);   // Compute weakly non-linear prediction for waveform along ray path in 3D/sph
}
#endif /* _GEOAC_EQSET_H_ */
