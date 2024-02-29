#ifndef EQUAL_ARRIVAL_TIME_SURFACE
#define EQUAL_ARRIVAL_TIME_SURFACE

// this header file provides a faster way of solving the equal-arrival-time surface. It performs binary search before doing interpolation.
#include "Definations.hpp"

class EATS_Solver {
public:
    const Array1D* t_data;
    const Array1D* theta_data;
    const Array3D* y_data;

    void config(const Array1D& t_data, const Array1D& theta_data, const Array3D& y_data, const double& nwind, const double& nism);

    // solve both t and primitive variables
    void solve_primitive(const double& Tobs_z, const double& theta, const double& phi, const double& theta_v, double& t, Array1D& y_pri);

private:
    double theta_min, theta_max;
    double t_min, t_max;
    int ntheta, nt;
    double nwind, nism;

    // find t_index at a specific theta_index
    void find_index(const double& mu, const double& Tobs_z, const int& theta_index, int& t_index1, int& t_index2);

    // solve t (EATS) at a given theta_index
    double solve_t(const double& mu, const double& Tobs_z, const int& theta_index, const int& t_index1, const int& t_index2);

    // solve primitive variables at a given theta_index
    void solve_primitive_at_a_cell(const double& t, const int& t_index1, const int& t_index2, const int& theta_index, Array1D& y_pri);

    // solve both t and primitive variables at a given theta_index
    void solve_a_cell(const double& mu, const double Tobs_z, const double& theta_index, double& t, Array1D& y_pri);
};

void EATS_Solver::config(const Array1D& t_data, const Array1D& theta_data, const Array3D& y_data, const double& nwind, const double& nism) {
    this->t_data = &t_data;
    this->theta_data = &theta_data;
    this->y_data = &y_data;
    this->nwind = nwind;
    this->nism = nism;
    theta_min = theta_data.front();
    theta_max = theta_data.back();
    t_min = t_data.front();
    t_max = t_data.front();
    ntheta = theta_data.size();
    nt = t_data.size();
}

void EATS_Solver::solve_primitive(const double& Tobs_z, const double& theta, const double& phi, const double& theta_v, double& t, Array1D& y_pri) {
    // ----- find theta index ----- //
    int theta_index1, theta_index2;

    // find theta index (binary search)
    if (theta <= theta_min) {
        theta_index1 = 0;
        theta_index2 = 0;
    }
    else if (theta >= theta_max) {
        theta_index1 = ntheta - 1;
        theta_index2 = ntheta - 1;
    }
    else {
        Find_Index(*theta_data, theta, theta_index1, theta_index2);
    }

    // ----- interpolate over theta ----- //
    if (theta <= theta_min) { // north pole
        // calculate cos angle
        double mu = std::cos(theta) * std::cos(theta_v) + std::sin(theta) * std::cos(phi) * std::sin(theta_v);

        // solve t and y_pri at theta_index
        solve_a_cell(mu, Tobs_z, theta_index1, t, y_pri);
    }
    else if (theta >= theta_max) { // south pole
        // calculate cos angle
        double mu = std::cos(theta) * std::cos(theta_v) + std::sin(theta) * std::cos(phi) * std::sin(theta_v);

        // solve t and y_pri at theta_index
        solve_a_cell(mu, Tobs_z, theta_index2, t, y_pri);
    }
    else { // middle
        double t1, t2;
        Array1D y_pri1 = Array1D(5);
        Array1D y_pri2 = Array1D(5);

        // calculate cos angle
        double mu1 = std::cos((*theta_data)[theta_index1]) * std::cos(theta_v) + std::sin((*theta_data)[theta_index1]) * std::cos(phi) * std::sin(theta_v);
        double mu2 = std::cos((*theta_data)[theta_index2]) * std::cos(theta_v) + std::sin((*theta_data)[theta_index2]) * std::cos(phi) * std::sin(theta_v);

        // solve t and y_pri at theta_index
        solve_a_cell(mu1, Tobs_z, theta_index1, t1, y_pri1);
        solve_a_cell(mu2, Tobs_z, theta_index2, t2, y_pri2);

        // interpolate time
        t = Linear((*theta_data)[theta_index1], (*theta_data)[theta_index2], t1, t2, theta);

        // interpolate primitive variables
        for (int i = 0; i < 5; ++i) {
            y_pri[i] = Linear((*theta_data)[theta_index1], (*theta_data)[theta_index2], y_pri1[i], y_pri2[i], theta);
        }
    }
}

void EATS_Solver::find_index(const double& mu, const double& Tobs_z, const int& theta_index, int& t_index1, int& t_index2) {
    // ----- find t index ----- //
    // equal arrival time surface function at t_index
    auto f = [&](const int& t_index) {
        // r at theta_index, t_index
        double r = (*y_data)[t_index][4][theta_index];

        // t at theta_index, r_index
        double t = (*t_data)[t_index];

        // equal arrival time surface function
        return t - r * mu / CONST.c - Tobs_z;
    };

    // find t index (binary search)
    if (f(0) > 0) { // smaller than tmin
        t_index1 = 0;
        t_index2 = 0;
    }
    else if (f(nt - 1) < 0) { // larger than tmax
        // throw error
        throw std::runtime_error("Observing time exceeds PDE maximum evolution time!\n");
    }
    else {
        t_index1 = 0;
        t_index2 = nt - 1;
        int index_mid;
        while (t_index2 - t_index1 > 1) {
            index_mid = (t_index1 + t_index2) / 2;
            if (f(index_mid) > 0) {
                t_index2 = index_mid;
            }
            else {
                t_index1 = index_mid;
            }
        }
    }
}

double EATS_Solver::solve_t(const double& mu, const double& Tobs_z, const int& theta_index, const int& t_index1, const int& t_index2) {
    // t at two ends
    double t1 = (*t_data)[t_index1];
    double t2 = (*t_data)[t_index2];

    // r values at two ends
    double r1 = (*y_data)[t_index1][4][theta_index];
    double r2 = (*y_data)[t_index2][4][theta_index];

    // if t < tmin
    if (t_index2 == 0) {
        return Tobs_z / (1 - r1 / (CONST.c * t_min) * mu);
    }
    else {
        // linear interpolation
        double slope = (r2 - r1) / (t2 - t1);
        double t = (Tobs_z + (r1 - slope * t1) * mu / CONST.c) / (1 - slope * mu / CONST.c);
        return t;
    }
}

void EATS_Solver::solve_primitive_at_a_cell(const double& t, const int& t_index1, const int& t_index2, const int& theta_index, Array1D& y_pri) {
    if (t < t_min) {
        double r = (*y_data)[0][4][theta_index] * t / t_min;
        y_pri[0] = nwind * CONST.mp * r / 1e17 * 1e51 + nism * CONST.mp * r * r * r / 3.0;
        y_pri[1] = (*y_data)[0][1][theta_index];
        y_pri[2] = (*y_data)[0][2][theta_index];
        y_pri[3] = (*y_data)[0][3][theta_index];
        y_pri[4] = r;
    }
    else {
        for (int i = 0; i < 5; ++i) {
            y_pri[i] = Linear((*t_data)[t_index1], (*t_data)[t_index2], (*y_data)[t_index1][i][theta_index], (*y_data)[t_index2][i][theta_index], t);
        }
    }
}

void EATS_Solver::solve_a_cell(const double& mu, const double Tobs_z, const double& theta_index, double& t, Array1D& y_pri) {
    // find t index
    int t_index1, t_index2;
    find_index(mu, Tobs_z, theta_index, t_index1, t_index2);

    // find t
    t = solve_t(mu, Tobs_z, theta_index, t_index1, t_index2);

    // find primitive
    solve_primitive_at_a_cell(t, t_index1, t_index2, theta_index, y_pri);
}

#endif