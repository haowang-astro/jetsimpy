// this file do what

#ifndef GRBAFTERGLOW
#define GRBAFTERGLOW

#include "Definations.hpp"
#include "Jet.hpp"
#include "JetIntegrator.hpp"
#include "Radiation.hpp"
#include "JetImaging.hpp"

class GRBafterglow {
public:
    // jet module
    void config_jet(const Array1D& theta_edge, const Array2D& yini, const bool& spreading, const double& tmin, const double& tmax, const double& rtol, const double& cfl, const double& nwind, const double& nism);
    void calibrate_jet();
    void solve_jet();
    Array1D& get_t_pde();
    Array3D& get_y_pde();
    Interpolator& get_interpolator();

    // integrator module
    void config_integrator(const std::string& model, const Dict& para_rad);
    double FluxDensity(const double& Tobs, const double& nu, const double& xtol, const double& rtol);
    double Offset(const double& Tobs, const double& nu, const double& xtol, const double& rtol);
    double Xscale(const double& Tobs, const double& nu, const double& xtol, const double& rtol);
    double Yscale(const double& Tobs, const double& nu, const double& xtol, const double& rtol);
    double Pi_lin(const double& Tobs, const double& nu, const double& xtol, const double& rtol);
    double intensity(const double& Tobs, const double& nu, const double& theta, const double& phi);
    double test_function(const double& Tobs, const double& nu, const double& theta, const double& phi);

    // imaging module
    void config_image();
    void solve_image(const double& Tobs, const double& nu, const int& inum, const int& pnum);
    AfterglowImage& image();

private:
    Jet1D jet;
    JetIntegrator jet_integrator;
    JetImaging jet_imaging;
};

// ---------- implementations ---------- //
void GRBafterglow::config_jet(const Array1D& theta_edge, const Array2D& yini, const bool& spreading, const double& tmin, const double& tmax, const double& rtol, const double& cfl, const double& nwind, const double& nism) {
    jet.config(theta_edge, yini, spreading, tmin, tmax, rtol, cfl, nwind, nism);
}

void GRBafterglow::calibrate_jet() {
    jet.pde_solver.calibration = true;
}

void GRBafterglow::solve_jet() {
    if (jet.spread) {
        jet.solve();
    }
    else {
        jet.solve_no_spreading();
    }
}

Array1D& GRBafterglow::get_t_pde() {
    return jet.ts;
}

Array3D& GRBafterglow::get_y_pde() {
    return jet.ys;
}

Interpolator& GRBafterglow::get_interpolator() {
    return jet.interpolator;
}

void GRBafterglow::config_integrator(const std::string& model, const Dict& para_rad) {
    jet_integrator.config(model, para_rad, &jet);
}

double GRBafterglow::FluxDensity(const double& Tobs, const double& nu, const double& xtol, const double& rtol) {
    return jet_integrator.FluxDensity(Tobs, nu, xtol, rtol);
}

double GRBafterglow::Offset(const double& Tobs, const double& nu, const double& xtol, const double& rtol) {
    return jet_integrator.Offset(Tobs, nu, xtol, rtol);
}

double GRBafterglow::Xscale(const double& Tobs, const double& nu, const double& xtol, const double& rtol) {
    return jet_integrator.Xscale(Tobs, nu, xtol, rtol);
}

double GRBafterglow::Yscale(const double& Tobs, const double& nu, const double& xtol, const double& rtol) {
    return jet_integrator.Yscale(Tobs, nu, xtol, rtol);
}

double GRBafterglow::Pi_lin(const double& Tobs, const double& nu, const double& xtol, const double& rtol) {
    return jet_integrator.Pi_lin(Tobs, nu, xtol, rtol);
}

double GRBafterglow::intensity(const double& Tobs, const double& nu, const double& theta, const double& phi) {
    double Tobs_z = Tobs / (1.0 + jet_integrator.para_rad.at("z"));
    double nu_z = nu * (1.0 + jet_integrator.para_rad.at("z"));
    return jet_integrator.Intensity(Tobs_z, nu_z, theta, phi);
}

double GRBafterglow::test_function(const double& Tobs, const double& nu, const double& theta, const double& phi) {
    double Tobs_z = Tobs / (1.0 + jet_integrator.para_rad.at("z"));
    double nu_z = nu * (1.0 + jet_integrator.para_rad.at("z"));
    return jet_integrator.test_function(Tobs_z, nu_z, theta, phi);
}

void GRBafterglow::config_image() {
    jet_imaging.config(&jet_integrator);
}

void GRBafterglow::solve_image(const double& Tobs, const double& nu, const int& inum, const int& pnum) {
    jet_imaging.solve_image(Tobs, nu, inum, pnum);
}

AfterglowImage& GRBafterglow::image() {
    return jet_imaging.image;
}

#endif