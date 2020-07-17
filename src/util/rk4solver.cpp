#ifndef RK4_SOLVER_CPP_
#define RK4_SOLVER_CPP_

#include <iostream>
#include <math.h>

#include "rk4solver.h"
#include "../atmo/atmo_state.h"
#include "../geoac/geoac.params.h"
#include "../geoac/geoac.eqset.h"

using namespace std;

int geoac::prop_rk4(double ** & solution, bool & check){
	int k = 0;
    int step_limit = s_max * int(1.0 / (ds_min * 10));
    double s = 0, ds = ds_min;
    
	double *temp0 = new double [eq_cnt];    double *partial1 = new double [eq_cnt];
    double *temp1 = new double [eq_cnt];    double *partial2 = new double [eq_cnt];
    double *temp2 = new double [eq_cnt];	double *partial3 = new double [eq_cnt];
    double *temp3 = new double [eq_cnt];
	double *temp4 = new double [eq_cnt];

	check = false;
    
    for(k = 0; k < (step_limit - 1); k++){
		for (int i = 0; i < eq_cnt; i++){
			temp0[i] = solution[k][i];
		}
        update_refs(s, temp0);
		ds = set_ds(temp0);

		s += ds;
		for (int i = 0; i < eq_cnt; i++){
			temp1[i] = ds * eval_src_eq(s, temp0, i);
			partial1[i] = solution[k][i] + temp1[i]/2.0;
		}		

		update_refs(s + ds/2, partial1);
		for (int i = 0; i < eq_cnt; i++){
			temp2[i] = ds * eval_src_eq(s + ds/2.0, partial1, i);
			partial2[i] = solution[k][i] + temp2[i]/2.0;
		}

		update_refs(s + ds/2, partial2);
		for (int i = 0; i < eq_cnt; i++){
			temp3[i] = ds * eval_src_eq(s + ds/2.0, partial2,i);
			partial3[i] = solution[k][i] + temp3[i];
		}

		update_refs(s + ds, partial3);
		for (int i = 0; i < eq_cnt; i++){
			temp4[i] = ds * eval_src_eq(s+ds, partial3, i);
			solution[k + 1][i] = solution[k][i] + temp1[i] / 6.0 + temp2[i] / 3.0 + temp3[i] / 3.0 + temp4[i] / 6.0;
		}

        if(break_check(solution, k + 1)){   check = true;	break;}
		if(ground_check(solution, k + 1)){  check = false;  break;}
	}
	
    delete [] temp0;   delete [] partial1;
    delete [] temp1;   delete [] partial2;
    delete [] temp2;   delete [] partial3;
    delete [] temp3;
    delete [] temp4;
	
	return k + 1;
}

#endif /* RK4_SOLVER_CPP_ */
