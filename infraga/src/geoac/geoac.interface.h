#ifndef _GEOAC_INTERFACE_H_
#define _GEOAC_INTERFACE_H_

namespace geoac{

//-------------------------//
//----Set Up The Solver----//
//-------------------------//
void configure();

//--------------------------------------------------//
//----Build, Clear, or Delete the Solution Array----//
//--------------------------------------------------//
void build_solution(double ** &, int);
void clear_solution(double ** &, int);
void delete_solution(double ** &, int);

//----------------------------//
//-------Output Profile-------//
//----------------------------//
void write_prof(char*, double);
void write_prof(char*, double, double, double);
void write_prof2d(char*, double, double, double);

}
#endif /* _GEOAC_INTERFACE_H_ */
