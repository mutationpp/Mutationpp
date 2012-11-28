#include "mutation++.h"

#include <time.h>

using namespace std;

double nanoseconds_per_call(
    const timespec& ts1, const timespec& ts2, const int iters)
{
    timespec ts;
    if ((ts2.tv_nsec-ts1.tv_nsec)<0) {
		ts.tv_sec = ts2.tv_sec-ts1.tv_sec-1;
		ts.tv_nsec = 1E9+ts2.tv_nsec-ts1.tv_nsec;
	} else {
		ts.tv_sec = ts2.tv_sec-ts1.tv_sec;
		ts.tv_nsec = ts2.tv_nsec-ts1.tv_nsec;
	}

    return (double)(ts.tv_sec * 1E9 + ts.tv_nsec) / (double) iters;
}

int main(int argc, char** argv)
{
    if (argc <= 1) exit(1);
    
    string mix_name(argv[argc-1]);
    cout << mix_name << endl;
    
    Mixture mix(mix_name);
        
    const int ne = mix.nElements();
    const int ns = mix.nSpecies();
    const int nr = mix.nReactions();
    
    double *p_c  = new double [ne];
    double *p_x  = new double [ns];
    double *p_cp = new double [ns];
    double *p_h  = new double [ns];
    double *p_s  = new double [ns];
    
    if (mix_name == "air5") {
        p_c[mix.elementIndex("N")]  = 0.79;
        p_c[mix.elementIndex("O")]  = 0.21;
    } else if (mix_name == "CO28") {
        p_c[mix.elementIndex("e-")] = 0.0;
        p_c[mix.elementIndex("C")]  = 1.0/3.0;
        p_c[mix.elementIndex("O")]  = 2.0/3.0;
    } else if (mix_name == "air11") {
        p_c[mix.elementIndex("e-")] = 0.0;
        p_c[mix.elementIndex("N")]  = 0.79;
        p_c[mix.elementIndex("O")]  = 0.21;
    } else if (mix_name == "air13") {
        p_c[mix.elementIndex("e-")] = 0.0;
        p_c[mix.elementIndex("N")]  = 0.78;
        p_c[mix.elementIndex("O")]  = 0.20;
        p_c[mix.elementIndex("Ar")] = 0.02;
    } else if (mix_name == "phenol40") {
        p_c[mix.elementIndex("e-")] = 0.0;
        p_c[mix.elementIndex("C")]  = 0.549;
        p_c[mix.elementIndex("H")]  = 0.124;
        p_c[mix.elementIndex("O")]  = 0.327;
    } else {
        exit(1);
    }
    
    
    
    timespec ts1, ts2;
    double T;    
    const int n = 100000;
    
    // Random seed
    srand((unsigned)time(NULL));
    
    // Time the setup functions needed in order to call the desired function
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts1);
    for (int i = 0; i < n; ++i) {
        T = 300.0 + 14700.0 * ((double)rand()/(double)RAND_MAX);
        mix.equilibrate(T, 101325.0, p_c, p_x);
        nd = mix.numberDensity(T, P);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts2);    
    cout << "setup:   " << 
        nanoseconds_per_call(ts1, ts2, n) << endl;
       
    // Now re-time with the function included
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts1);
    for (int i = 0; i < n; ++i) {
        T = 300.0 + 14700.0 / (double)n * (double)i;
        mix.equilibrate(T, 101325.0, p_c, p_x);
        nd = mix.numberDensity(T, P);
        eta = mix.eta(T, nd, p_x);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &ts2);    
    cout << "setup + function:  " << 
        nanoseconds_per_call(ts1, ts2, n) << endl;
}
