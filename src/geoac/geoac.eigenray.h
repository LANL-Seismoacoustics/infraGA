#ifndef _GEOAC_EIGENRAY_H_
#define _GEOAC_EIGENRAY_H_

namespace geoac {

extern bool verbose;
extern int eigenray_cnt;
    
extern double damping;
extern double tolerance;
    
extern double dth_big;
extern double dth_sml;

extern std::ofstream eig_results;

double mod_dth(double, double);                      // Function to modify theta_step in est_eigenray function

bool est_eigenray(double [3], double [2], double, double, double &, double &, double &, int, double);               // Identify inclination at fixed azimuth for specified arrival range
void find_eigenray(double [3], double [2], double, double, double, int, int, char []);                              // Identify exact eigenray

}

#endif /* GEOAC_EIGENRAY_H_ */
