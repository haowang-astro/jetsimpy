#ifndef EQUAL_ARRIVAL_TIME_SURFACE
#define EQUAL_ARRIVAL_TIME_SURFACE

#include "../environment.h"
#include "../Hydro/sim_box.h"
#include "../Hydro/tools.h"
#include "blast.h"

// equal arrival time surface solver
class EATS {
public:
    EATS() {}
    void feedData(SimBox& sim_box, Tool& tool);    // feed PDE data

    // solve EATS & blast properties
    void solveBlast(double Tobs_z, double theta, double phi, double theta_v, Blast& blast);
    
private:
    // PDE data
    const Array3D* y_data;
    const Array1D* t_data;
    const Array1D* theta_data;

    Tool* tool;               // tool
    int ntheta;               // number of cells
    int nt;                   // number of time data
    double tmin;              // pde start time
    double tmax;              // pde end time
    double theta_min;         // minimum cell center
    double theta_max;         // maximum cell center

    // find index by binary search
    void findThetaIndex(double theta, int& theta_index1, int& theta_index2);
    void findTimeIndex(double mu, double Tobs_z, int theta_index, int& t_index1, int& t_index2);

    // solve t for EATS at theta_index
    double solveT(double mu, double Tobs_z, int theta_index, int t_index1, int t_index2);
    
    // solve primitive variables and return primitive & t
    Array1D solvePrimitive(double mu, double Tobs_z, int theta_index);

    // derive blast properties given the primitive variables
    void deriveBlast(double theta, double phi, double theta_v, const Array1D& val, Blast& blast);
};

#endif