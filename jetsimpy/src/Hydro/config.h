#ifndef CONFIG
#define CONFIG

#include "../environment.h"

// simulation configuration
struct JetConfig {
    // initial conditions
    Array1D theta_edge;             // cell faces
    Array1D Eb;                     // Eb
    Array1D Ht;                     // Hb * beta_th
    Array1D Msw;                    // Msw
    Array1D Mej;                    // Mej
    Array1D R;                      // R
    double nwind;                   // wind density scale
    double nism;                    // interstellar medium density scale
    
    // configurations
    double tmin;                    // PDE starting time
    double tmax;                    // PDE ending time
    double rtol;                    // velocity solver tolerance
    double cfl;                     // Courant number
    bool spread;                    // Whether enabling spreading
    int calib_level;                // calibration level: 
                                    //   0: no calibration. 
                                    //   1: calibrate with Blandford-McKee all time.
                                    //   2: calibrate with Blandford-McKee in ultra-relativistic phase and Sedov-Taylor in Newtonian phase.
};

#endif