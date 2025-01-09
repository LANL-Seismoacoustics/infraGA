#ifndef _GEOAC_EIGENRAY_MPI_H_
#define _GEOAC_EIGENRAY_MPI_H_

namespace geoac {

    extern bool verbose;
    extern int verbose_opt;

    extern int eigenray_cnt;
    
    extern double damping;
    extern double tolerance;
    
    extern double dth_big;
    extern double dth_sml;

    extern std::ofstream eig_results;

    double mod_dth(double, double);                                                                                                 // Function to modify theta_step in est_eigenray function
    bool est_eigenray(double [3], double [2], double, double, double &, double &, double &, int, double, MPI_Comm, int, int, int);  // Identify inclination at fixed azimuth for specified arrival range
    bool find_eigenray(double [3], double [2], double &, double &, double, int, int, char [], int);                                 // Identify exact eigenray

}

#endif /* _GEOAC_EIGENRAY_MPI_H_ */
