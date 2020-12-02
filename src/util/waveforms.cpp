#ifndef WAVEFORMS_CPP_
#define WAVEFORMS_CPP_

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>
#include <sstream>
#include <fftw3.h>

#include "fileIO.h"
#include "waveforms.h"

using namespace std;

//-----------------------------//
//-----Waveform Parameters-----//
//-----------------------------//
double wvfrm::p0 = 10.0;
double wvfrm::t0 = 0.1;
double wvfrm::alpha = 4.0;

long int wvfrm::len = pow(2, 13);
int wvfrm::N_rise = 10;
int wvfrm::N_period = 40;

double wvfrm::filt1_g = 0.5;
double wvfrm::filt1_n1 = 4.0;
double wvfrm::filt1_n2 = 0.5;

double wvfrm::filt2_g = 0.75;
double wvfrm::filt2_n1 = 2.0;
double wvfrm::filt2_n2 = 1.0;

//------------------------------//
//-------Built-in Waveforms-----//
//------------------------------//
double wvfrm::impulse(double t){
    double x, x_pk, result;
    
    if(t >= 0.0){
        x = t / t0;
        x_pk = (alpha + 1.0) - sqrt(alpha + 1.0);
        
        result = pow(x, alpha) * (1.0 - x / (1.0 + alpha)) * exp(-x);
        result /= pow(x_pk, alpha) * (1.0 - x_pk / (1.0 + alpha)) * exp(-x_pk);
        
        return result;
    } else {
        return 0.0;
    }
}

double wvfrm::n_wave(double t){
    if(fabs(t) <= t0){  return t / t0;}
    else {              return 0.0;}
}

//-----------------------------//
//-------Methods to Define-----//
//----------the Waveform-------//
//-----------------------------//
void wvfrm::load_wvfrm(double** & wvfrm_array, char* wvfrm_file){
    cout << '\t' << "Loading waveform from file: " << wvfrm_file << "." << '\n';

    len = file_length(wvfrm_file);
    
    wvfrm_array = new double * [len];
    for(int n = 0; n < len; n++){
        wvfrm_array[n] = new double [2];
    }

    string line;
    ifstream file_in;        
    file_in.open(wvfrm_file);

    int n = 0;
    while(!file_in.eof() && n < len){
        getline(file_in, line);
        if(line.find("#") != 0){
            stringstream ss(line);
            ss >> wvfrm_array[n][0];
            ss >> wvfrm_array[n][1];
            n++;
        }
    }
    file_in.close();

}

void wvfrm::build_wvfrm(double** & wvfrm_array, char* option){
    double duration;
    
    wvfrm_array = new double * [len];
    if (strncmp(option, "impulse", 7) == 0){
        alpha = max(alpha, 0.01);
        duration = t0 * len * max(((1.0 + alpha) - sqrt(1.0 + alpha)) / N_rise, 2.0 * 3.14159 * sqrt(1 + alpha) / N_period);

        cout << '\t' << "Defining waveform from built in impulse with shaping parameter = " << alpha << "." << '\n';
        for(int n = 0; n < len; n++){
            wvfrm_array[n] = new double [2];

            double t = -duration / 2.0 + duration * n / len;
            wvfrm_array[n][0] = t;
            wvfrm_array[n][1] = p0 * impulse(t); 
        }
    } else if ((strncmp(option, "Uwave", 5) == 0) || (strncmp(option, "u-wave", 6) == 0)){
        alpha = max(alpha, 2.0);
        duration = t0 * len * max(((1.0 + alpha) - sqrt(1.0 + alpha)) / N_rise, 2.0 * 3.14159 * sqrt(1 + alpha) / N_period);

        cout << '\t' << "Defining waveform from built in U-wave with shaping parameter = " << alpha << "." << '\n';
        double eps=1.0e-2;
        for(int n = 0; n < len; n++){
            wvfrm_array[n] = new double [2];

            double t = -duration / 2.0 + duration * n / len;
            wvfrm_array[n][0] = t;
            wvfrm_array[n][1] = p0 * (impulse(t + eps) - impulse(t - eps)) / (2.0 * eps);
            wvfrm_array[n][1] /= 1.126 * (impulse(eps) - impulse(-eps)) / (2.0 * eps);
        }
    } else if ((strncmp(option, "Nwave", 5) == 0) || (strncmp(option, "n-wave", 6) == 0)){
        duration = (t0 / N_period) * len;

        cout << '\t' << "Defining waveform from built in N-wave." << '\n';        
        for(int n = 0; n < len; n++){
            wvfrm_array[n] = new double [2];

            double t = -duration / 2.0 + duration * n / len;
            wvfrm_array[n][0] = t;
            wvfrm_array[n][1] = - p0 * n_wave(t - t0 / 2.0);
        }
    } else {
        alpha = max(alpha, 0.01);
        duration = t0 * len * min(((1.0 + alpha) - sqrt(1.0 + alpha)) / N_rise, 2.0 * 3.14159 * sqrt(1 + alpha) / N_period);

        cout << '\t' << "Warning!!  Unrecognized option in waveform.  Options are 'impulse', 'Uwave', and 'Nwave'  Using impulse with shaping parameter " << alpha << "." << '\n';
        for(int n = 0; n < len; n++){

            wvfrm_array[n] = new double [2];

            double t = -duration / 2.0 + duration * n / len;
            wvfrm_array[n][0] = t;
            wvfrm_array[n][1] = p0 * impulse(t);
        }
    }
}

void wvfrm::delete_wvfrm(double** & wvfrm_array){
    for(int n = 0; n < len; n++){
        delete wvfrm_array[n];
    }
    delete wvfrm_array;
}



#endif /* WAVEFORMS_CPP_ */
