#ifndef DEFINATION
#define DEFINATION

#include <iostream>
#include <string>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <vector>
#include <list>
#include <exception>
#include <map>
#include <tuple>
#include <iterator>
#include <limits>
#include <chrono>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "math.hpp"

namespace py = pybind11;

// ------------------ //
//       Arrays       //
// ------------------ //

using Array1D = std::vector<double>;
using Array2D = std::vector<Array1D>;
using Array3D = std::vector<Array2D>;

Array1D Array(int n1, double fill = 0.0) {
    return std::vector<double>(n1, fill);
}
Array2D Array(int n1, int n2, double fill = 0.0) {
    Array1D array1d = Array(n2, fill);
    return std::vector<std::vector<double>>(n1, array1d);
}
Array3D Array(int n1, int n2, int n3, double fill = 0.0) {
    Array2D array2d = Array(n2, n3, fill);
    return std::vector<std::vector<std::vector<double>>>(n1, array2d);
}

// ------------------ //
//       structs      //
// ------------------ //

struct Constants {
    double c = 	29979245800.0;
    double mp = 1.672622e-24;
    double me = 9.109384e-28;
    double sigma_t = 6.6524587e-25;
    double e = 4.803204673e-10;
    double mpc = 3.09e24;
    double pi = M_PI;
    double mas = 1.0 / 206264806.24709466;
}CONST;

struct Blast {
    double t;
    double Tobs;
    double theta;
    double phi;
    double R, dR;
    double beta, gamma, beta_th, beta_r;
    double beta_f, gamma_f;
    double e;
    double n_blast, n_ambient;
    double p;
    double doppler;
    double cos_theta_r;
    double cos_theta_beta;
    double s;
};

struct Stokes {
    double Pi_lin = 0.0;    // total linear polarization
    double Pi_circ = 0.0;   // circular polarization
    double emissivity = 0.0; // emissivity

    void zero() {
        Pi_lin = 0.0;
        Pi_circ = 0.0;
        emissivity = 0.0;
    }
};

// ------------------ //
//        alias       //
// ------------------ //

using Dict = std::map<std::string, double>;
using RadDict = std::map<std::string, void(*)(const double&, const Dict&, const Blast&, Stokes&)>;

#endif