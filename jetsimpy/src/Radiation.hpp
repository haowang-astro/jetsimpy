#ifndef RADIATION
#define RADIATION

#include "Definations.hpp"

// ------------------------------------------------------------------------ //
//            Define your radiation models inside this namespace            //
// ------------------------------------------------------------------------ //
namespace Rad_Model {

    // The template of a radiation model:
    //
    // void Model_Name(const double& nu,        // frequency
    //                 const Dict& P,           // same as the parameter disctionary from Python side
    //                 const Blast& blast,      // the object to access all fluid element properties in its coming frame
    //                 Stokes& S)               // the return values
    // {
    //     ...
    //
    //     S.emissivity = xxx;           // isotropic emissivity
    //     S.Pi_lin = xxx;               // local linear polarization degree, i.e., Stokes √(Q^2 + U^2) / I
    //     S.Pi_circ = xxx;              // (not working yet) local circular polarization degree, i.e., Stokes V / I
    //     
    // }

    // ~~~~~~~~~~ Built-in Radiation model ~~~~~~~~~~ //
    // standard synchrotron radiation spectrum + linear polarization (Gruzinov 1999; Sari 1999)
    void Synchrotron(const double& nu, const Dict& P, const Blast& blast, Stokes& S) {
        double eps_e = P.at("eps_e");
        double eps_b = P.at("eps_b");
        double p = P.at("p");
        double b = P.at("b");
        
        double n_blast = blast.n_blast;
        double t = blast.t;
        double gamma = blast.gamma;
        double e = blast.e;
        double sin_theta_beta_sq = 1.0 - blast.cos_theta_beta * blast.cos_theta_beta;
    
        double emissivity;
        double gamma_m, gamma_c, B, nu_m, nu_c, e_p;
        
        gamma_m = (p - 2.0) / (p - 1.0) * (eps_e * CONST.mp / CONST.me * (gamma - 1.0));
        B = std::sqrt(8.0 * CONST.pi * eps_b * e);
        gamma_c = 6.0 * CONST.pi * CONST.me * gamma * CONST.c / CONST.sigma_t / B / B / t;
        nu_m = 3.0 * CONST.e * B * gamma_m * gamma_m / 4.0 / CONST.pi / CONST.c / CONST.me;
        nu_c = 3.0 * CONST.e * B * gamma_c * gamma_c / 4.0 / CONST.pi / CONST.c / CONST.me;
        e_p = std::sqrt(3.0) * CONST.e * CONST.e * CONST.e * B * n_blast / CONST.me / CONST.c / CONST.c;

        if (nu_m < nu_c) {
            if (nu < nu_m) {
                emissivity = e_p * std::cbrt(nu / nu_m);
            }
            else if (nu < nu_c) {
                emissivity = e_p * std::pow(nu / nu_m, - (p - 1) / 2.0);
            }
            else {
                emissivity = e_p * std::pow(nu_c / nu_m, - (p - 1) / 2.0) * std::pow(nu / nu_c, - p / 2);
            }
        }
        else {
            if (nu < nu_c) {
                emissivity = e_p * std::cbrt(nu / nu_c);
            }
            else if (nu < nu_m) {
                emissivity = e_p / std::sqrt(nu / nu_c);
            }
            else {
                emissivity = e_p / std::sqrt(nu_m / nu_c) * std::pow(nu / nu_m, - p / 2);
            }
        }

        // emissivity and polarization in the [comoving] frame.
        S.emissivity = emissivity;
        S.Pi_lin = (p + 1.0) / (p + 7.0 / 3.0) * sin_theta_beta_sq * (b - 1.0) / (2.0 + sin_theta_beta_sq * (b - 1.0));
        S.Pi_circ = 0.0;
    }

    // ---------- define your emissivity models below ---------- //

    // Deep Newtonian Phase
    void Synchrotron_DNP(const double& nu, const Dict& P, const Blast& blast, Stokes& S) {
        double eps_e = P.at("eps_e");
        double eps_b = P.at("eps_b");
        double p = P.at("p");
        double b = P.at("b");
        
        double n_blast = blast.n_blast;
        double t = blast.t;
        double gamma = blast.gamma;
        double e = blast.e;
        double sin_theta_beta_sq = 1.0 - blast.cos_theta_beta * blast.cos_theta_beta;
    
        double emissivity;
        double gamma_m, gamma_c, B, nu_m, nu_c, e_p;
        
        gamma_m = (p - 2.0) / (p - 1.0) * eps_e * CONST.mp / CONST.me * (gamma - 1.0);
        double f = 1.0;
        if (gamma_m <= 1) {
            gamma_m = 1.0;
            f = (p - 2.0) / (p - 1.0) * eps_e * CONST.mp / CONST.me * (gamma - 1.0) / gamma_m;
        }
        B = std::sqrt(8.0 * CONST.pi * eps_b * e);
        gamma_c = 6.0 * CONST.pi * CONST.me * gamma * CONST.c / CONST.sigma_t / B / B / t;
        nu_m = 3.0 * CONST.e * B * gamma_m * gamma_m / 4.0 / CONST.pi / CONST.c / CONST.me;
        nu_c = 3.0 * CONST.e * B * gamma_c * gamma_c / 4.0 / CONST.pi / CONST.c / CONST.me;
        e_p = std::sqrt(3.0) * CONST.e * CONST.e * CONST.e * B * f * n_blast / CONST.me / CONST.c / CONST.c;

        if (nu_m < nu_c) {
            if (nu < nu_m) {
                emissivity = e_p * std::cbrt(nu / nu_m);
            }
            else if (nu < nu_c) {
                emissivity = e_p * std::pow(nu / nu_m, - (p - 1) / 2.0);
            }
            else {
                emissivity = e_p * std::pow(nu_c / nu_m, - (p - 1) / 2.0) * std::pow(nu / nu_c, - p / 2);
            }
        }
        else {
            if (nu < nu_c) {
                emissivity = e_p * std::cbrt(nu / nu_c);
            }
            else if (nu < nu_m) {
                emissivity = e_p / std::sqrt(nu / nu_c);
            }
            else {
                emissivity = e_p / std::sqrt(nu_m / nu_c) * std::pow(nu / nu_m, - p / 2);
            }
        }

        // emissivity and polarization in the [comoving] frame.
        S.emissivity = emissivity;
        S.Pi_lin = (p + 1.0) / (p + 7.0 / 3.0) * sin_theta_beta_sq * (b - 1.0) / (2.0 + sin_theta_beta_sq * (b - 1.0));
        S.Pi_circ = 0.0;
    }
}

// ----- register radiation models ----- //
// this function must be put at the end of this file //
void JetIntegrator::register_models() {
    // built-in emissivity
    rad_dict["sync"] = &Rad_Model::Synchrotron;

    // custom emissivity
    rad_dict["dnp"] = &Rad_Model::Synchrotron_DNP;
}

#endif