# ifndef GLOBE_CPP_
# define GLOBE_CPP_

#include <math.h>
#include "globe.h"

//---------------------------//
//----Radius of the Earth----//
//---------------------------//
double globe::r0 = 6370.0;

//----------------------------------------//
//----------Functions to Compute----------//
//-----Great Circle Bearing and Range-----//
//----------------------------------------//
double globe::bearing(double lat1, double lon1, double lat2, double lon2){
    double term1 = sin(lon2 - lon1);
    double term2 = cos(lat1) * tan(lat2) - sin(lat1) * cos(lon2 - lon1);
    
    return atan2(term1, term2);
}

double globe::gc_dist(double lat1, double lon1, double lat2, double lon2){
    double term1 = pow(sin((lat2 - lat1) / 2.0), 2);
    double term2 = cos(lat1) * cos(lat2) * pow(sin((lon2 - lon1) / 2.0), 2);
    
    return 2.0 * r0 * asin(sqrt(term1 + term2));
}

#endif /* GLOBE_CPP_ */
