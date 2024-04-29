#ifndef JETSIMPY
#define JETSIMPY

#include "environment.h"
#include "Hydro/config.h"
#include "Hydro/tools.h"
#include "Hydro/sim_box.h"
#include "Hydro/interpolate.h"
#include "Afterglow/eats.h"
#include "Afterglow/afterglow.h"

// The object accessable from Python
class Jet {
public:
    Jet(const JetConfig& jet_config);

    // ---------- solve hydro ---------- //
    void solveJet();                     // solve hydro
    py::array_t<double> getY();          // get pde ys
    py::array_t<double> getT();          // get pde ts
    py::array_t<double> getTheta();      // get cell centers

    // ---------- hydro interpolation ---------- //
    double interpolateMsw(double t, double theta);
    double interpolateMej(double t, double theta);
    double interpolateBetaGamma(double t, double theta);
    double interpolateBetaTh(double t, double theta);
    double interpolateR(double t, double theta);
    double interpolateE0(double t, double theta);

    // ---------- afterglow calculation ---------- //
    void configParameters(const Dict& param);              // configurate parameter dictionary
    void configEmissivity(const std::string& model_name);  // configurate emissivity model
    void configAvgModel(const std::string& model_name);    // configure average models
    double calculateEATS(double Tobs, double theta, double phi, double theta_v, double z);    // calculate t of EATS
    double calculateIntensity(double Tobs, double nu, double theta, double phi);    // intensity in Jet coordinate
    double calculateLuminosity(double Tobs, double nu, double rtol);    // integrate luminosity
    double calculateAvgModel(double Tobs, double nu, double rtol);      // integrate average model
    double WeightedAverage(double Tobs, double nu, double rtol);
    double IntensityOfPixel(const double Tobs, const double nu, const double x_tilde, const double y_tilde);

private:
    JetConfig jet_config;            // configuration data
    Tool tool;                       // the tool containing many functions
    SimBox sim_box;                  // the simulation box
    Interpolator interpolator;       // the interpolation tool
    EATS eats;                       // equal arrival time surface solver
    Afterglow afterglow;             // afterglow algorithms
};

#endif