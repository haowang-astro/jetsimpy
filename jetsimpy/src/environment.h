#ifndef ENVIRONMENT
#define ENVIRONMENT

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
#include <functional>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

// global constants in cgs unit
const double CSpeed = 29979245800.0;
const double MassP = 1.672622e-24;
const double MassE = 9.109384e-28;
const double SigmaT = 6.6524587e-25;
const double ECharge = 4.803204673e-10;
const double MPC = 3.09e24;
const double PI = 3.14159265358979323846;
const double MAS = 1.0 / 206264806.24709466;

// type alias
namespace py = pybind11;
using Array1D = std::vector<double>;
using Array2D = std::vector<Array1D>;
using Array3D = std::vector<Array2D>;
using Dict = std::map<std::string, double>;

// array constructions
inline Array1D Array(int n1, double fill = 0.0) {
    return std::vector<double>(n1, fill);
}
inline Array2D Array(int n1, int n2, double fill = 0.0) {
    Array1D array1d = Array(n2, fill);
    return std::vector<std::vector<double>>(n1, array1d);
}
inline Array3D Array(int n1, int n2, int n3, double fill = 0.0) {
    Array2D array2d = Array(n2, n3, fill);
    return std::vector<std::vector<std::vector<double>>>(n1, array2d);
}

// print function
template <typename T>
void print(T var1) {
    std::cout << var1 << "\n";
}

template <typename T, typename... Types>
void print(T var1, Types... var2) {
    std::cout << var1 << " ";
    print(var2...);
}



#endif