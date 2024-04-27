#ifndef TOOLS
#define TOOLS

#include "../environment.h"
#include "config.h"
#include "../Math/root.h"

// the tool to store useful functions
class Tool {
public:
    Tool(const JetConfig& jet_config);
    double solveDensity(double r);                                      // solve ambient density
    double solveBetaGammaSq(double msw_eb, double mej_eb, double r);    // solve (beta * gamma) by Msw/Eb and Mej/Eb
    double minmod(double x1, double x2);                                // Minmod function
    double solveS(double r, double beta_gamma_sq);                      // solve calibration coefficient
    void findIndex(const Array1D& x_array, const double x, int& index1, int& index2); // binary search for index range
    double linear(double x, double x1, double x2, double y1, double y2);  // linear interpolation

private:
    const double factor = 2.0;      // a constant factor to interpolate BM & ST
    double nwind;                   // wind density scale
    double nism;                    // ISM density scale
    double rtol;                    // relative tolerance of solving (beta * gamma)
    int calib_level;                // calibration level.
};

#endif