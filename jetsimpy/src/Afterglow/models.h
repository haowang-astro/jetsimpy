#ifndef MODELS
#define MODELS

#include "../environment.h"
#include "blast.h"

using ModelDict = std::map<std::string, std::function<double(const double, const Dict&, const Blast&)>>;

struct Models {
    // models to be registered in "models.cpp"
    ModelDict emissivity_models;
    ModelDict avg_models;

    // register emissivity models
    void registerEmissivity();

    // register average models for weighted average calculation
    void registerAvgModels();
};

#endif