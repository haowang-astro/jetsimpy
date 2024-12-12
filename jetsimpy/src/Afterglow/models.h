#ifndef MODELS
#define MODELS

#include "../environment.h"
#include "blast.h"

using Dict = std::map<std::string, double>;
using ModelDict = std::map<std::string, std::function<double(const double, const Dict&, const Blast&)>>;

struct Models {
    // constants (define your constants here)
    const double CSpeed = 29979245800.0;
    const double MassP = 1.672622e-24;
    const double MassE = 9.109384e-28;
    const double SigmaT = 6.6524587e-25;
    const double ECharge = 4.803204673e-10;
    const double MPC = 3.09e24;
    const double PI = 3.14159265358979323846;
    const double MAS = 1.0 / 206264806.24709466;

    // models to be registered in "models.cpp"
    ModelDict radiation_models;
    ModelDict avg_models;

    // register radiation models
    void registerIntensity();

    // register average models for weighted average calculation
    void registerAvgModels();
};

#endif