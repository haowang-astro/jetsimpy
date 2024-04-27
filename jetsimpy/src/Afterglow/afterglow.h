#ifndef JET_INTEGRATE
#define JET_INTEGRATE

#include "../environment.h"
#include "../Hydro/sim_box.h"
#include "../Hydro/tools.h"
#include "../Math/integral.h"
#include "../Math/optimize.h"
#include "blast.h"
#include "eats.h"
#include "models.h"

// This class provides a routine to calculate afterglow
class Afterglow {
public:
    Afterglow() {}
    void initialize(SimBox& sim_box, Tool& tool);          // initialize object
    void configParameters(const Dict& param);              // configurate parameter dictionary
    void configEmissivity(const std::string& model_name);  // configure emissivity model
    void configAvgModel(const std::string& model_name);    // configure average models
    
    // calculate the specific intensity
    double Intensity(const double Tobs, const double nu, const double theta, const double phi);

    // total luminosity
    double Luminosity(const double Tobs, const double nu, const double rtol);

    // integrate average model with dL_dOmega as the weight
    double integrateModel(const double Tobs, const double nu, const double rtol);

    // intensity of pixel (useful for sky map)
    double IntensityOfPixel(const double Tobs, const double nu, const double x_tilde, const double y_tilde);

private:
    Array1D* theta_data;
    Blast blast;    // blast object
    Dict param;     // parameter dictionary from python side

    // model function (pointer) to be called in integration
    double (*emissivity_model)(const double, const Dict&, const Blast&);
    double (*avg_model)(const double, const Dict&, const Blast&);

    // model object to store model functions
    Models models;

    // parameters we must have
    double theta_v;                 // observing angle
    double z;                       // redshift
    double d;                       // luminosity distance

    // equal arrival time surface solver
    EATS eats;

    // calculate dL/dOmega
    double dL_dOmega(const double Tobs_z, const double nu_z, const double theta, const double phi);

    // find peak
    double findPeak(const double Tobs_z, const double nu_z);
};

#endif