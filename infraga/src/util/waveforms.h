#ifndef WAVEFORMS_H_
#define WAVEFORMS_H_

namespace wvfrm {

    extern double p0;           // Amplitude of the defined waveform
    extern double t0;           // Time scale of the defined waveform
    extern double alpha;        // Shaping parameter of the defined waveform (only used for impulse)
    
    extern long int len;        // Integer length of the waveform array
    extern int N_rise;          // Points over the rise time for a weak shock solution
    extern int N_period;        // Points per period used to define dt

    extern double filt1_g;      // Nyquist frequency scalar for low-pass filter (filter 1 scales non-linear contribution)
    extern double filt1_n1;     // First exponential coefficient for low-pass filter (filter 1 scales non-linear contribution)
    extern double filt1_n2;     // Second exponential coefficient for low-pass filter (filter 1 scales non-linear contribution)
    
    extern double filt2_g;      // Nyquist frequency scalar for low-pass filter (filter 2 scales full waveform at each step)
    extern double filt2_n1;     // First exponential coefficient for low-pass filter (filter 2 scales full waveform at each step)
    extern double filt2_n2;     // Second exponential coefficient for low-pass filter (filter 2 scales full waveform at each step)
    
    double impulse(double);                 // A generalized acoustic impulse definition (asymptotes to blast wave) developed by Waxler
    double n_wave(double);                  // An N-wave waveform deffinition
    
    void load_wvfrm(double** &, char*);      // Load waveform from a file
    void build_wvfrm(double** &, char*);     // Build waveform from built in options
    void delete_wvfrm(double** &);     // Delete the waveform array
}

#endif /* WAVEFORMS_H_ */
