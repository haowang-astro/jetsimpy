#include "interpolate.h"

void Interpolator::feedData(SimBox& sim_box, Tool& tool) {
    y_data = &(sim_box.getY());
    t_data = &(sim_box.getT());
    theta_data = &(sim_box.getTheta());
    this->tool = &tool;
    ntheta = theta_data->size();
    tmin = t_data->front();
    tmax = t_data->back();
    theta_min = theta_data->front();
    theta_max = theta_data->back();
}

void Interpolator::findThetaIndex(double theta, int& theta_index1, int& theta_index2) {
    if (theta < 0.0 || theta > PI) {
        // bound check error
        throw std::runtime_error("Interpolation: theta outside bounds.\n");
    }
    else if (theta < theta_min) {
        // north pole
        theta_index1 = 0;
        theta_index2 = 0;
    }
    else if (theta > theta_max) {
        // south pole
        theta_index1 = ntheta - 1;
        theta_index2 = ntheta - 1;
    }
    else {
        // middle
        tool->findIndex(*theta_data, theta, theta_index1, theta_index2);
    }
}

void Interpolator::findTimeIndex(double t, int& t_index1, int& t_index2) {
    if (t < 0.0 || t > tmax) {
        // bound check error
        throw std::runtime_error("Interpolation: t outside bounds.\n");
    }
    else if (t < tmin) {
        // before initial condition
        t_index1 = 0;
        t_index2 = 0;
    }
    else {
        // in the solution domain
        tool->findIndex(*t_data, t, t_index1, t_index2);
    }
}

// ----- interpolate over t at specific theta index ----- //

double Interpolator::interpolateY(double t, double theta, int y_index) {
    // find index
    int theta_index1, theta_index2;
    int t_index1, t_index2;
    findThetaIndex(theta, theta_index1, theta_index2);
    findTimeIndex(t, t_index1, t_index2);

    // coordinate values at the corners [t, theta]
    double t1 = (*t_data)[t_index1];
    double t2 = (*t_data)[t_index2];
    double theta1 = (*theta_data)[theta_index1];
    double theta2 = (*theta_data)[theta_index2];    
    
    // y values at the corners
    double y11 = (*y_data)[y_index][theta_index1][t_index1];
    double y12 = (*y_data)[y_index][theta_index2][t_index1];
    double y21 = (*y_data)[y_index][theta_index1][t_index2];
    double y22 = (*y_data)[y_index][theta_index2][t_index2];
    
    // interpolate over theta
    double y1 = (theta_index1 == theta_index2) ?  // in the poles?
                y11 // constant in the pole
                : 
                tool->linear(theta, theta1, theta2, y11, y12);
    double y2 = (theta_index1 == theta_index2) ?  // in the poles?
                y21 // constant in the pole
                : 
                tool->linear(theta, theta1, theta2, y21, y22);
    
    // interpolate over t
    double y = (t_index1 == t_index2) ?  // before initial condition?
                y1 // constant extrapolation before initial condition
                :
                tool->linear(t, t1, t2, y1, y2);
    return y;
}
