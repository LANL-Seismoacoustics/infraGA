# ifndef GLOBE_H_
# define GLOBE_H_

//----------------------------------------//
//----------Functions to Compute----------//
//-----Great Circle Bearing and Range-----//
//----------------------------------------//
namespace  globe{
    extern double r0;                                   // Radius for spherical earth model
    
    double bearing(double, double, double, double);     // Bearing between two lat/lon locations
    double gc_dist(double, double, double, double);     // Great circle distance between two lat/lon locations
}

#endif /* GLOBE_H_ */
