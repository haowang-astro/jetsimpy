#ifndef INTERPOLATE
#define INTERPOLATE

#include "../environment.h"
#include "tools.h"
#include "sim_box.h"

class Interpolator {
public:
    Interpolator() {}
    void feedData(SimBox& sim_box, Tool& tool);     // feed PDE data

    // interpolate for thet y_index'th value
    double interpolateY(double t, double theta, int y_index);
    
private:
    const Array3D* y_data;          // pde data y
    const Array1D* t_data;          // pde data t
    const Array1D* theta_data;      // pde cell center
    Tool* tool;               // tool
    int ntheta;               // number of cells
    double tmin;              // pde start time
    double tmax;              // pde end time
    double theta_min;         // minimum cell center
    double theta_max;         // maximum cell center

    // find index by binary search
    void findThetaIndex(double theta, int& theta_index1, int& theta_index2);
    void findTimeIndex(double t, int& t_index1, int& t_index2);
};

#endif