#include "models.h"

void Models::registerEmissivity() {
    // The template of a emissivity model:
    //
    // emissivity_models["model_name"] = [](    // "model_name" is the keyword used from Python side
    //    const double nu,               // frequency
    //    const Dict& P,                 // parameter dictionary from Python side
    //    const Blast& blast,            // an object that contains all fluid element properties
    // ) {
    //
    //     ...
    //
    //     return emissivity;            // return value must be a "double"
    // };

    // ---------- default radiation model (Sari 1998) ---------- //
    emissivity_models["sync"] = [&](const double nu, const Dict& P, const Blast& blast) {
        double eps_e = P.at("eps_e");
        double eps_b = P.at("eps_b");
        double p = P.at("p");
        
        double n_blast = blast.n_blast;
        double t = blast.t;
        double gamma = blast.gamma;
        double e = blast.e_density;
    
        double emissivity;
        double gamma_m, gamma_c, B, nu_m, nu_c, e_p;
        
        gamma_m = (p - 2.0) / (p - 1.0) * (eps_e * MassP / MassE * (gamma - 1.0));
        B = std::sqrt(8.0 * PI * eps_b * e);
        gamma_c = 6.0 * PI * MassE * gamma * CSpeed / SigmaT / B / B / t;
        nu_m = 3.0 * ECharge * B * gamma_m * gamma_m / 4.0 / PI / CSpeed / MassE;
        nu_c = 3.0 * ECharge * B * gamma_c * gamma_c / 4.0 / PI / CSpeed / MassE;
        e_p = std::sqrt(3.0) * ECharge * ECharge * ECharge * B * n_blast / MassE / CSpeed / CSpeed;

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

        return emissivity;
    };

    // Bonus! deep newtonian phase correction
    emissivity_models["sync_dnp"] = [&](const double nu, const Dict& P, const Blast& blast) {
        double eps_e = P.at("eps_e");
        double eps_b = P.at("eps_b");
        double p = P.at("p");
        
        double n_blast = blast.n_blast;
        double t = blast.t;
        double gamma = blast.gamma;
        double e = blast.e_density;
    
        double emissivity;
        double gamma_m, gamma_c, B, nu_m, nu_c, e_p;
        
        gamma_m = (p - 2.0) / (p - 1.0) * eps_e * MassP / MassE * (gamma - 1.0);
        double f = 1.0;
        if (gamma_m <= 1) {
            gamma_m = 1.0;
            f = (p - 2.0) / (p - 1.0) * eps_e * MassP / MassE * (gamma - 1.0) / gamma_m;
        }
        B = std::sqrt(8.0 * PI * eps_b * e);
        gamma_c = 6.0 * PI * MassE * gamma * CSpeed / SigmaT / B / B / t;
        nu_m = 3.0 * ECharge * B * gamma_m * gamma_m / 4.0 / PI / CSpeed / MassE;
        nu_c = 3.0 * ECharge * B * gamma_c * gamma_c / 4.0 / PI / CSpeed / MassE;
        e_p = std::sqrt(3.0) * ECharge * ECharge * ECharge * B * f * n_blast / MassE / CSpeed / CSpeed;

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

        return emissivity;
    };

    // ---------- define your own model below ---------- //
    // emissivity_models["model_name"] = [](const double nu, const Dict& P, const Blast& blast) {
    //     ...
    //     return emissivity;
    // }
}

// weighted average models
void Models::registerAvgModels() {
    // for offset
    avg_models["offset"] = [](const double nu, const Dict& P, const Blast& blast) {
        double theta_v = P.at("theta_v");
        double x_tilde = - std::sin(blast.theta) * std::cos(blast.phi) * std::cos(theta_v) + std::cos(blast.theta) * std::sin(theta_v);
        return x_tilde * blast.R;
    };

    // for sigma_x
    avg_models["sigma_x"] = [](const double nu, const Dict& P, const Blast& blast) {
        double theta_v = P.at("theta_v");
        double x_tilde = - std::sin(blast.theta) * std::cos(blast.phi) * std::cos(theta_v) + std::cos(blast.theta) * std::sin(theta_v);
        return x_tilde * blast.R * x_tilde * blast.R;
    };

    // for sigma_y
    avg_models["sigma_y"] = [](const double nu, const Dict& P, const Blast& blast) {
        double y = std::sin(blast.theta) * std::sin(blast.phi);
        return y * blast.R * y * blast.R;
    };

    // ---------- define your own model below ---------- //
    // avg_models["model_name"] = [](const double nu, const Dict& P, const Blast& blast) {
    //     ...
    //     return xxx;
    // }
}