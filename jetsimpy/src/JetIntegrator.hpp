// this file do what
#ifndef JETINTEGRATOR
#define JETINTEGRATOR

#include "Definations.hpp"
#include "RootSolver.hpp"
#include "Integral.hpp"
#include "Optimization.hpp"
#include "Jet.hpp"

class JetIntegrator {
public:
    // configuration
    void config(const std::string& model, const Dict& para_rad, Jet1D* jet);
    Jet1D* jet;

    RadDict rad_dict;           // Radiation models
    Dict para_rad;              // radiation parameters
    Blast blast;                // fluid element properties
    void(*radiation)(const double&, const Dict&, const Blast&, Stokes&); // radiation function pointer
    Stokes S;         // polarization properties

    // solve fluid element property and save it in .blast
    void get_blast(const double& Tobs_z, const double& theta, const double& phi);

    // integrate and get flux, polarization, and superluminal motion
    double Intensity(const double& Tobs_z, const double& nu_z, const double& theta, const double& phi);
    double dL_dOmega(const double& Tobs_z, const double& nu_z, const double& theta, const double& phi);
    double FluxDensity(const double& Tobs, const double& nu, const double& xtol, const double& rtol);
    double Offset(const double& Tobs, const double& nu, const double& xtol, const double& rtol);
    double Pi_lin(const double& Tobs, const double& nu, const double& xtol, const double& rtol);
    double Pi_circ(const double& Tobs, const double& nu, const double& xtol, const double& rtol);
    double Xscale(const double& Tobs, const double& nu, const double& xtol, const double& rtol);
    double Yscale(const double& Tobs, const double& nu, const double& xtol, const double& rtol);
    double SigmaX(const double& Tobs, const double& nu, const double& xtol, const double& rtol);
    double SigmaY(const double& Tobs, const double& nu, const double& xtol, const double& rtol);
    double test_function(const double& Tobs_z, const double& nu_z, const double& theta, const double& phi);

    // imaging
    void solve_image(const double& Tobs, const double& nu, const int& inum, const int& pnum);

private:
    Array1D cos_theta_samples;  // initial integration samples of cos(theta)
    Array1D phi_samples;        // initial integration samples of phi
    int max_iter = 50;          // maximum iteration

    // solve initial integration samples
    void cos_theta_phi_samples(const double& Tobs_z, const double& nu_z, double& theta_peak);

    // register radiation model
    void register_models();
};

// implementation
void JetIntegrator::config(const std::string& model, const Dict& para_rad, Jet1D* jet) {
    register_models();

    // ----- setup radiation models
    try {
        this->radiation = rad_dict.at(model);
    }
    catch(const std::exception& e) {
        std::string text = "Radiation model '" + model + "' doesn't exist!";
        throw std::runtime_error(text);
    }

    this->para_rad = para_rad;
    this->jet = jet;
}

// solve fluid element
void JetIntegrator::get_blast(const double& Tobs_z, const double& theta, const double& phi) {
    using namespace std;
    double& theta_v = para_rad.at("theta_v");

    // solve eats and primitive variables
    double t;
    Array1D y_pri = Array(5);
    jet->eats_solver.solve_primitive(Tobs_z, theta, phi, theta_v, t, y_pri);

    // solve derived variables
    double Msw, u, beta_t, R, s; 
    double dR, beta, gamma, beta_r, e, n_blast, p, doppler, n;
    double nr, nth, mu_alpha, mu_beta, mu_r;
    
    // PDE variables
    Msw = y_pri[0];
    u = std::sqrt(y_pri[2]);
    beta_t = y_pri[3];
    R = y_pri[4];
    if (jet->pde_solver.calibration) {
        double sbm = calib_BM(jet->nwind, jet->nism, R);
        double sst = calib_ST(jet->nwind, jet->nism, R);
        s = (sst + sbm * 2 * y_pri[2]) / (1.0 + 2 * y_pri[2]);
    }
    else {
        s = 1.0;
    }
    
    // derived variables
    n = density0(jet->nwind, jet->nism, R);
    gamma = sqrt(u * u + 1.0);
    beta = u / gamma;
    beta_r = (beta * beta - beta_t * beta_t > 0) ? sqrt(beta * beta - beta_t * beta_t) : 0.0; // prevent beta_th > beta
    n_blast = 4.0 * gamma * n;
    e = (gamma - 1.0) * n_blast * CONST.mp * CONST.c * CONST.c * s;
    p = 4.0 / 3.0 * u * u * n * CONST.mp * CONST.c * CONST.c * s;
    dR = Msw / R / R / n_blast / CONST.mp;
    
    // doppler factor
    nr = beta_r / beta;
    nth = beta_t / beta;
    mu_beta = (nr * sin(theta) * cos(phi) + nth * cos(theta) * cos(phi)) * sin(theta_v)
            + (nr * cos(theta) - nth * sin(theta)) * cos(theta_v);
    mu_r = std::cos(theta) * std::cos(theta_v) + std::sin(theta) * std::cos(phi) * std::sin(theta_v);
    doppler = 1.0 / gamma / (1.0 - beta * mu_beta);
    //doppler = 1.0 / gamma / (1.0 - beta * mu_r);

    // insert to blast
    blast.t = t;
    blast.Tobs = Tobs_z;
    blast.theta = theta;
    blast.phi = phi;
    blast.R = R;
    blast.dR = dR;
    blast.beta = beta;
    blast.beta_f = 4.0 * beta * gamma * gamma / (4.0 * gamma * gamma - 2.0);
    blast.gamma = gamma;
    blast.gamma_f = (4.0 * gamma * gamma - 1.0) / std::sqrt(8 * gamma * gamma + 1.0);
    blast.beta_th = beta_t;
    blast.beta_r = beta_r;
    blast.e = e;
    blast.n_blast = n_blast;
    blast.n_ambient = n;
    blast.p = p;
    blast.doppler = doppler;
    blast.cos_theta_beta = (mu_beta - beta) / (1.0 - beta * mu_beta);
    blast.s = s;
}

double JetIntegrator::Intensity(const double& Tobs_z, const double& nu_z, const double& theta, const double& phi) {
    double& theta_v = para_rad.at("theta_v");

    // solve blast properties
    this->get_blast(Tobs_z, theta, phi);
    
    // nu in comoving frame
    double nu_src = nu_z / blast.doppler;

    // zero the polarization
    S.zero();
    
    // solve radiation
    (*this->radiation)(nu_src, para_rad, blast, S);

    // intensity
    return S.emissivity / 4.0 / CONST.pi * blast.dR * blast.doppler * blast.doppler * blast.doppler;
}

double JetIntegrator::test_function(const double& Tobs_z, const double& nu_z, const double& theta, const double& phi) {
    double& theta_v = para_rad.at("theta_v");

    // solve blast properties
    this->get_blast(Tobs_z, theta, phi);
    
    // nu in comoving frame
    double nu_src = nu_z / blast.doppler;

    // zero the polarization
    S.zero();
    
    // solve radiation
    (*this->radiation)(nu_src, para_rad, blast, S);

    return blast.t;
}

double JetIntegrator::dL_dOmega(const double& Tobs_z, const double& nu_z, const double& theta, const double& phi) {
    double intensity = Intensity(Tobs_z, nu_z, theta, phi);
    return 4.0 * CONST.pi * intensity * blast.R * blast.R;
}

void JetIntegrator::cos_theta_phi_samples(const double& Tobs_z, const double& nu_z, double& theta_peak) {
    double& theta_v = para_rad.at("theta_v");

    // ------------------ find peak --------------------- //
    // define function
    auto f = [&](const double& theta) {
        return - std::log(dL_dOmega(Tobs_z, nu_z, theta, 0.0));
    };

    // solve theta_peak
    try {
        theta_peak = minimization(f, jet->theta_edge, 1e-6);
    }
    catch(const std::exception& e) {
        std::string text = "\nFailed to find angle of peak intensity. \n↓\n";
        throw std::runtime_error(text + e.what());
    }

    // solve beaming angle
    get_blast(Tobs_z, theta_peak, 0.0);
    double beaming_angle = 1.0 / blast.gamma;

    cos_theta_samples = {-1, std::cos(beaming_angle), std::cos(beaming_angle / 2.0), 1};
    phi_samples = {0.0, CONST.pi};
}

// slove flux
double JetIntegrator::FluxDensity(const double& Tobs, const double& nu, const double& xtol, const double& rtol) {
    double z = para_rad.at("z");
    double dl = para_rad.at("d");
    double Tobs_z = Tobs / (1.0 + z);
    double nu_z = nu * (1.0 + z);

    // initial integration samples
    double theta_peak;
    cos_theta_phi_samples(Tobs_z, nu_z, theta_peak);
    
    // function to integrate
    auto f = [&](const double& cos_theta_rot, const double& phi_rot) -> Array1D {
        double theta_rot = std::acos(cos_theta_rot);

        // interruption detection
        if (PyErr_CheckSignals() != 0) {
            throw std::runtime_error("Interrupted at Flux calculation.");
        }
        
        // transform from "peak" coordinate to jet coordinate
        double x = std::sin(theta_rot) * std::cos(phi_rot) * std::cos(theta_peak) + std::cos(theta_rot) * std::sin(theta_peak);
        double y = std::sin(theta_rot) * std::sin(phi_rot);
        double z = - std::sin(theta_rot) * std::cos(phi_rot) * std::sin(theta_peak) + std::cos(theta_rot) * std::cos(theta_peak);

        double cos_theta = std::min(1.0, std::fabs(z)) * z / std::fabs(z);
        double theta = std::acos(cos_theta);

        double phi = std::atan2(y, x);
        phi = (phi < 0.0) ? phi + CONST.pi * 2.0 : phi;

        // luminosity
        double dl_domega = this->dL_dOmega(Tobs_z, nu_z, theta, phi);
        
        // integral body
        return {dl_domega};
    };

    double result;
    try {
        result = Adaptive_2D(f, cos_theta_samples, phi_samples, xtol, rtol, max_iter)[0];
    }
    catch(const std::exception& e) {
        std::string text = "\nFailed to perform adaptive integration. \n↓\n";
        throw std::runtime_error(text + e.what());
    }

    return 2.0 * result * (1.0 + z) / 4.0 / CONST.pi / dl / dl / CONST.mpc / CONST.mpc;
}

double JetIntegrator::Offset(const double& Tobs, const double& nu, const double& xtol, const double& rtol) {
    double theta_v = para_rad.at("theta_v");
    double z = para_rad.at("z");
    double dl = para_rad.at("d");
    double Tobs_z = Tobs / (1.0 + z);
    double nu_z = nu * (1.0 + z);

    // avoid the hard calculation of theta_v = 0
    if (theta_v == 0.0) {
        return 0.0;
    }

    // avoid Tobs = 0
    if (Tobs == 0.0) {
        return 0.0;
    }

    // solve integrating samples
    double theta_peak;
    cos_theta_phi_samples(Tobs_z, nu_z, theta_peak);

    // function to integrate
    auto f = [&](const double& cos_theta_rot, const double& phi_rot) -> Array1D {
        double theta_rot = std::acos(cos_theta_rot);

        // interruption detection
        if (PyErr_CheckSignals() != 0) {
            throw std::runtime_error("Interrupted at flux centroid calculation.");
        }

        // transform from "peak" coordinate to jet coordinate
        double x = std::sin(theta_rot) * std::cos(phi_rot) * std::cos(theta_peak) + std::cos(theta_rot) * std::sin(theta_peak);
        double y = std::sin(theta_rot) * std::sin(phi_rot);
        double z = - std::sin(theta_rot) * std::cos(phi_rot) * std::sin(theta_peak) + std::cos(theta_rot) * std::cos(theta_peak);

        double cos_theta = std::min(1.0, std::fabs(z)) * z / std::fabs(z);
        double theta = std::acos(cos_theta);

        double phi = std::atan2(y, x);
        phi = (phi < 0.0) ? phi + CONST.pi * 2.0 : phi;

        // luminosity
        double dl_domega = this->dL_dOmega(Tobs_z, nu_z, theta, phi); // get_blast run

        // LOS coordinate
        double x_tilde = std::sin(theta) * std::cos(phi) * std::cos(theta_v) - std::cos(theta) * std::sin(theta_v);

        // integral body
        double dX_domega = - x_tilde * blast.R * dl_domega;

        return {dl_domega, dX_domega};
    };

    // integrate over phi
    Array1D result;
    try {
        result = Adaptive_2D(f, cos_theta_samples, phi_samples, xtol, rtol, max_iter);
    }
    catch(const std::exception& e) {
        std::string text = "\nFailed to perform adaptive integration. \n↓\n";
        throw std::runtime_error(text + e.what());
    }

    return result[1] / result[0] / dl / CONST.mpc / (1.0 + z) / (1.0 + z) / CONST.mas;
}

double JetIntegrator::Xscale(const double& Tobs, const double& nu, const double& xtol, const double& rtol) {
    double theta_v = para_rad.at("theta_v");
    double z = para_rad.at("z");
    double dl = para_rad.at("d");
    double Tobs_z = Tobs / (1.0 + z);
    double nu_z = nu * (1.0 + z);

    // centroid
    double xc_offset = Offset(Tobs, nu, xtol, rtol);
    double xc = xc_offset * dl * CONST.mpc * (1.0 + z) * (1.0 + z) * CONST.mas;

    // solve integrating samples
    double theta_peak;
    cos_theta_phi_samples(Tobs_z, nu_z, theta_peak);

    // function to integrate
    auto f = [&](const double& cos_theta_rot, const double& phi_rot) -> Array1D {
        double theta_rot = std::acos(cos_theta_rot);

        // interruption detection
        if (PyErr_CheckSignals() != 0) {
            throw std::runtime_error("Interrupted at flux centroid calculation.");
        }

        // transform from "peak" coordinate to jet coordinate
        double x = std::sin(theta_rot) * std::cos(phi_rot) * std::cos(theta_peak) + std::cos(theta_rot) * std::sin(theta_peak);
        double y = std::sin(theta_rot) * std::sin(phi_rot);
        double z = - std::sin(theta_rot) * std::cos(phi_rot) * std::sin(theta_peak) + std::cos(theta_rot) * std::cos(theta_peak);

        double cos_theta = std::min(1.0, std::fabs(z)) * z / std::fabs(z);
        double theta = std::acos(cos_theta);

        double phi = std::atan2(y, x);
        phi = (phi < 0.0) ? phi + CONST.pi * 2.0 : phi;

        // luminosity
        double dl_domega = this->dL_dOmega(Tobs_z, nu_z, theta, phi); // get_blast run

        // LOS coordinate
        double x_tilde = std::sin(theta) * std::cos(phi) * std::cos(theta_v) - std::cos(theta) * std::sin(theta_v);

        // integral body
        double dX_domega = std::abs(- x_tilde * blast.R - xc) * dl_domega;
        return {dl_domega, dX_domega};
    };

    // integrate over phi
    Array1D result;
    try {
        result = Adaptive_2D(f, cos_theta_samples, phi_samples, xtol, rtol, max_iter);
    }
    catch(const std::exception& e) {
        std::string text = "\nFailed to perform adaptive integration. \n↓\n";
        throw std::runtime_error(text + e.what());
    }

    return 4.0 * result[1] / result[0] / dl / CONST.mpc / (1.0 + z) / (1.0 + z) / CONST.mas;
}

double JetIntegrator::Yscale(const double& Tobs, const double& nu, const double& xtol, const double& rtol) {
    double z = para_rad.at("z");
    double dl = para_rad.at("d");
    double Tobs_z = Tobs / (1.0 + z);
    double nu_z = nu * (1.0 + z);

    // solve integrating samples
    double theta_peak;
    cos_theta_phi_samples(Tobs_z, nu_z, theta_peak);

    // function to integrate
    auto f = [&](const double& cos_theta_rot, const double& phi_rot) -> Array1D {
        double theta_rot = std::acos(cos_theta_rot);

        // interruption detection
        if (PyErr_CheckSignals() != 0) {
            throw std::runtime_error("Interrupted at flux centroid calculation.");
        }

        // transform from "peak" coordinate to jet coordinate
        double x = std::sin(theta_rot) * std::cos(phi_rot) * std::cos(theta_peak) + std::cos(theta_rot) * std::sin(theta_peak);
        double y = std::sin(theta_rot) * std::sin(phi_rot);
        double z = - std::sin(theta_rot) * std::cos(phi_rot) * std::sin(theta_peak) + std::cos(theta_rot) * std::cos(theta_peak);

        double cos_theta = std::min(1.0, std::fabs(z)) * z / std::fabs(z);
        double theta = std::acos(cos_theta);

        double phi = std::atan2(y, x);
        phi = (phi < 0.0) ? phi + CONST.pi * 2.0 : phi;

        // luminosity
        double dl_domega = this->dL_dOmega(Tobs_z, nu_z, theta, phi); // get_blast run

        // LOS coordinate
        //double y_tilde = y;

        // integral body
        double dY_domega = std::abs(y) * blast.R * dl_domega;
        return {dl_domega, dY_domega};
    };

    // integrate over phi
    Array1D result;
    try {
        result = Adaptive_2D(f, cos_theta_samples, phi_samples, xtol, rtol, max_iter);
    }
    catch(const std::exception& e) {
        std::string text = "\nFailed to perform adaptive integration. \n↓\n";
        throw std::runtime_error(text + e.what());
    }

    return 4.0 * result[1] / result[0] / dl / CONST.mpc / (1.0 + z) / (1.0 + z) / CONST.mas;
}

double JetIntegrator::SigmaX(const double& Tobs, const double& nu, const double& xtol, const double& rtol) {
    double theta_v = para_rad.at("theta_v");
    double z = para_rad.at("z");
    double dl = para_rad.at("d");
    double Tobs_z = Tobs / (1.0 + z);
    double nu_z = nu * (1.0 + z);

    // centroid
    double xc_offset = Offset(Tobs, nu, xtol, rtol);
    double xc = xc_offset * dl * CONST.mpc * (1.0 + z) * (1.0 + z) * CONST.mas;

    // solve integrating samples
    double theta_peak;
    cos_theta_phi_samples(Tobs_z, nu_z, theta_peak);

    // function to integrate
    auto f = [&](const double& cos_theta_rot, const double& phi_rot) -> Array1D {
        double theta_rot = std::acos(cos_theta_rot);

        // interruption detection
        if (PyErr_CheckSignals() != 0) {
            throw std::runtime_error("Interrupted at flux centroid calculation.");
        }

        // transform from "peak" coordinate to jet coordinate
        double x = std::sin(theta_rot) * std::cos(phi_rot) * std::cos(theta_peak) + std::cos(theta_rot) * std::sin(theta_peak);
        double y = std::sin(theta_rot) * std::sin(phi_rot);
        double z = - std::sin(theta_rot) * std::cos(phi_rot) * std::sin(theta_peak) + std::cos(theta_rot) * std::cos(theta_peak);

        double cos_theta = std::min(1.0, std::fabs(z)) * z / std::fabs(z);
        double theta = std::acos(cos_theta);

        double phi = std::atan2(y, x);
        phi = (phi < 0.0) ? phi + CONST.pi * 2.0 : phi;

        // luminosity
        double dl_domega = this->dL_dOmega(Tobs_z, nu_z, theta, phi); // get_blast run

        // LOS coordinate
        double x_tilde = - std::sin(theta) * std::cos(phi) * std::cos(theta_v) + std::cos(theta) * std::sin(theta_v);

        // integral body
        double dX_domega = (x_tilde * blast.R - xc) * (x_tilde * blast.R - xc) * dl_domega;
        return {dl_domega, dX_domega};
    };

    // integrate over phi
    Array1D result;
    try {
        result = Adaptive_2D(f, cos_theta_samples, phi_samples, xtol, rtol, max_iter);
    }
    catch(const std::exception& e) {
        std::string text = "\nFailed to perform adaptive integration. \n↓\n";
        throw std::runtime_error(text + e.what());
    }

    double sigma_x_sq = result[1] / result[0];
    return std::sqrt(sigma_x_sq) / dl / CONST.mpc / (1.0 + z) / (1.0 + z) / CONST.mas;
}

double JetIntegrator::SigmaY(const double& Tobs, const double& nu, const double& xtol, const double& rtol) {
    double z = para_rad.at("z");
    double dl = para_rad.at("d");
    double Tobs_z = Tobs / (1.0 + z);
    double nu_z = nu * (1.0 + z);

    // solve integrating samples
    double theta_peak;
    cos_theta_phi_samples(Tobs_z, nu_z, theta_peak);

    // function to integrate
    auto f = [&](const double& cos_theta_rot, const double& phi_rot) -> Array1D {
        double theta_rot = std::acos(cos_theta_rot);

        // interruption detection
        if (PyErr_CheckSignals() != 0) {
            throw std::runtime_error("Interrupted at flux centroid calculation.");
        }

        // transform from "peak" coordinate to jet coordinate
        double x = std::sin(theta_rot) * std::cos(phi_rot) * std::cos(theta_peak) + std::cos(theta_rot) * std::sin(theta_peak);
        double y = std::sin(theta_rot) * std::sin(phi_rot);
        double z = - std::sin(theta_rot) * std::cos(phi_rot) * std::sin(theta_peak) + std::cos(theta_rot) * std::cos(theta_peak);

        double cos_theta = std::min(1.0, std::fabs(z)) * z / std::fabs(z);
        double theta = std::acos(cos_theta);

        double phi = std::atan2(y, x);
        phi = (phi < 0.0) ? phi + CONST.pi * 2.0 : phi;

        // luminosity
        double dl_domega = this->dL_dOmega(Tobs_z, nu_z, theta, phi); // get_blast run

        // LOS coordinate
        //double y_tilde = y;

        // integral body
        double dY_domega = y * blast.R * y * blast.R * dl_domega;
        return {dl_domega, dY_domega};
    };

    // integrate over phi
    Array1D result;
    try {
        result = Adaptive_2D(f, cos_theta_samples, phi_samples, xtol, rtol, max_iter);
    }
    catch(const std::exception& e) {
        std::string text = "\nFailed to perform adaptive integration. \n↓\n";
        throw std::runtime_error(text + e.what());
    }

    double sigma_y_sq = result[1] / result[0];
    return std::sqrt(sigma_y_sq) / dl / CONST.mpc / (1.0 + z) / (1.0 + z) / CONST.mas;
}

double JetIntegrator::Pi_lin(const double& Tobs, const double& nu, const double& xtol, const double& rtol) {
    double theta_v = para_rad.at("theta_v");
    double z = para_rad.at("z");
    double Tobs_z = Tobs / (1.0 + z);
    double nu_z = nu * (1.0 + z);

    // avoid the hard calculation of theta_v = 0
    if (theta_v == 0.0) {
        return 0.0;
    }

    // initial integration samples
    double theta_peak;
    cos_theta_phi_samples(Tobs_z, nu_z, theta_peak);

    // function to integrate
    auto f = [&](const double& cos_theta_rot, const double& phi_rot) -> Array1D {
        double theta_rot = std::acos(cos_theta_rot);

        // interruption detection
        if (PyErr_CheckSignals() != 0) {
            throw std::runtime_error("Interrupted at Linear Polarization calculation.");
        }
        
        // transform from "peak" coordinate to jet coordinate
        double x = std::sin(theta_rot) * std::cos(phi_rot) * std::cos(theta_peak) + std::cos(theta_rot) * std::sin(theta_peak);
        double y = std::sin(theta_rot) * std::sin(phi_rot);
        double z = - std::sin(theta_rot) * std::cos(phi_rot) * std::sin(theta_peak) + std::cos(theta_rot) * std::cos(theta_peak);

        double cos_theta = std::min(1.0, std::fabs(z)) * z / std::fabs(z);
        double theta = std::acos(cos_theta);

        double phi = std::atan2(y, x);
        phi = (phi < 0.0) ? phi + CONST.pi * 2.0 : phi;
        
        // luminosity
        double dl_domega = this->dL_dOmega(Tobs_z, nu_z, theta, phi); // get_blast run

        // LOS coordinate
        double x_tilde = std::sin(theta) * std::cos(phi) * std::cos(theta_v) - std::cos(theta) * std::sin(theta_v);
        double y_tilde = std::sin(theta) * std::sin(phi);
        double tan_phi_tilde = (x_tilde == 0.0) ? 0.0 : y_tilde / x_tilde;
        double cos_2phi_tilde = (1.0 - tan_phi_tilde * tan_phi_tilde) / (1.0 + tan_phi_tilde * tan_phi_tilde);

        // integral body
        double dLin_domega = S.Pi_lin * cos_2phi_tilde * dl_domega;
        return {dl_domega, dLin_domega};
    };

    // integrate over phi
    Array1D result;
    try {
        result = Adaptive_2D(f, cos_theta_samples, phi_samples, xtol, rtol, max_iter);
    }
    catch(const std::exception& e) {
        std::string text = "\nFailed to perform adaptive integration. \n↓\n";
        throw std::runtime_error(text + e.what());
    }

    return result[1] / result[0];
}

#endif