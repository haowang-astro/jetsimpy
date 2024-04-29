#include "afterglow.h"

void Afterglow::initialize(SimBox& sim_box, EATS& eats) {
    this->eats = &eats;
    //eats.feedData(sim_box, tool);
    theta_data = &(sim_box.getTheta());
    models.registerEmissivity();
    models.registerAvgModels();
}

void Afterglow::configParameters(const Dict& param) {
    // save important parameters and perform parameter checking
    theta_v = param.at("theta_v");
    d = param.at("d");
    z = param.at("z");
    this->param = param;
}

void Afterglow::configEmissivity(const std::string& model_name) {
    try {
        emissivity_model = models.emissivity_models.at(model_name);
    }
    catch (const std::exception& e) {
        throw std::runtime_error("Emissivity: Model name not found!");
    }
}

void Afterglow::configAvgModel(const std::string& model_name) {
    try {
        avg_model = models.avg_models[model_name];
    }
    catch (const std::exception& e) {
        throw std::runtime_error("Weighted average: Model name not found!");
    }
}

double Afterglow::Intensity(const double Tobs, const double nu, const double theta, const double phi) {
    double Tobs_z = Tobs / (1 + z);
    double nu_z = nu * (1 + z);
    
    // solve equal-arrival-time-surface and get blast properties
    eats->solveBlast(Tobs_z, theta, phi, theta_v, blast);
    
    // nu in comoving frame
    double nu_src = nu_z / blast.doppler;
    
    // solve emissivity at comoving frame
    double emissivity = (*emissivity_model)(nu_src, param, blast);

    // convert to intensity at observer's frame
    return emissivity / 4.0 / PI * blast.dR * blast.doppler * blast.doppler * blast.doppler;
}

double Afterglow::dL_dOmega(const double Tobs_z, const double nu_z, const double theta, const double phi) {
    // nu in comoving frame
    double nu_src = nu_z / blast.doppler;
    
    // solve emissivity at comoving frame
    double emissivity = (*emissivity_model)(nu_src, param, blast);

    return emissivity * blast.dR * blast.R * blast.R * blast.doppler * blast.doppler * blast.doppler;
}

double Afterglow::Luminosity(const double Tobs, const double nu, const double rtol) {
    double Tobs_z = Tobs / (1 + z);
    double nu_z = nu * (1 + z);
    double theta_peak = findPeak(Tobs_z, nu_z);

    auto f = [&](double cos_theta_rot, double phi_rot) {
        // interruption detection
        if (PyErr_CheckSignals() != 0) {
            throw std::runtime_error("Luminosity: Keyboard interruption.");
        }

        // transform from "peak" coordinate to jet coordinate
        double theta_rot = std::acos(cos_theta_rot);

        double x = std::sin(theta_rot) * std::cos(phi_rot) * std::cos(theta_peak) + std::cos(theta_rot) * std::sin(theta_peak);
        double y = std::sin(theta_rot) * std::sin(phi_rot);
        double z = - std::sin(theta_rot) * std::cos(phi_rot) * std::sin(theta_peak) + std::cos(theta_rot) * std::cos(theta_peak);

        double cos_theta = std::min(1.0, std::fabs(z)) * z / std::fabs(z);
        double theta = std::acos(cos_theta);

        double phi = std::atan2(y, x);
        phi = (phi < 0.0) ? phi + PI * 2.0 : phi;

        // solve equal-arrival-time-surface and get blast properties
        eats->solveBlast(Tobs_z, theta, phi, theta_v, blast);

        return dL_dOmega(Tobs_z, nu_z, theta, phi);
    };

    // beaming angle
    eats->solveBlast(Tobs_z, theta_peak, 0.0, theta_v, blast);
    double beaming_angle = 1.0 / blast.gamma;

    // initial inegral samples
    Array1D cos_theta_samples = {-1.0, std::cos(beaming_angle), std::cos(beaming_angle / 2.0), 1.0};
    Array1D phi_samples = {0.0, PI};

    // multiply by 2 because the integral domain for phi is [0, pi].
    double luminosity = Adaptive_2D(f, cos_theta_samples, phi_samples, 0.0, rtol) * 2.0;

    return luminosity;
}

double Afterglow::integrateModel(const double Tobs, const double nu, const double rtol) {
    double Tobs_z = Tobs / (1 + z);
    double nu_z = nu * (1 + z);
    double theta_peak = findPeak(Tobs_z, nu_z);
    double nu_src = nu_z / blast.doppler;

    auto f = [&](double cos_theta_rot, double phi_rot) {
        // interruption detection
        if (PyErr_CheckSignals() != 0) {
            throw std::runtime_error("Average Model: Keyboard interruption.");
        }

        // transform from "peak" coordinate to jet coordinate
        double theta_rot = std::acos(cos_theta_rot);

        double x = std::sin(theta_rot) * std::cos(phi_rot) * std::cos(theta_peak) + std::cos(theta_rot) * std::sin(theta_peak);
        double y = std::sin(theta_rot) * std::sin(phi_rot);
        double z = - std::sin(theta_rot) * std::cos(phi_rot) * std::sin(theta_peak) + std::cos(theta_rot) * std::cos(theta_peak);

        double cos_theta = std::min(1.0, std::fabs(z)) * z / std::fabs(z);
        double theta = std::acos(cos_theta);

        double phi = std::atan2(y, x);
        phi = (phi < 0.0) ? phi + PI * 2.0 : phi;

        // solve equal-arrival-time-surface and get blast properties
        eats->solveBlast(Tobs_z, theta, phi, theta_v, blast);

        return (*avg_model)(nu_src, param, blast) * dL_dOmega(Tobs_z, nu_z, theta, phi);
    };

    // beaming angle
    eats->solveBlast(Tobs_z, theta_peak, 0.0, theta_v, blast);
    double beaming_angle = 1.0 / blast.gamma;

    // initial inegral samples (phi from 0 to 2 * pi)
    Array1D cos_theta_samples = {-1.0, std::cos(beaming_angle), std::cos(beaming_angle / 2.0), 1.0};
    Array1D phi_samples = {0.0, PI, 2.0 * PI};

    double integral = Adaptive_2D(f, cos_theta_samples, phi_samples, 0.0, rtol);

    return integral;
}

double Afterglow::findPeak(const double Tobs_z, const double nu_z) {
    // define function
    auto f = [&](const double& theta) {
        return - dL_dOmega(Tobs_z, nu_z, theta, 0.0);
    };

    // solve theta_peak
    double theta_peak = minimization(f, *theta_data, 1e-6);
    return theta_peak;
}

double Afterglow::IntensityOfPixel(const double Tobs, const double nu, const double x_tilde, const double y_tilde) {
    // function to solve intensity from LOS spherical coordinate
    auto f_intensity = [&](const double theta_tilde, const double phi_tilde) {
        // convert to source coordinate (cartisan)
        double x = std::sin(theta_tilde) * std::cos(phi_tilde) * std::cos(theta_v) + std::cos(theta_tilde) * std::sin(theta_v);
        double y = std::sin(theta_tilde) * std::sin(phi_tilde);
        double z = - std::sin(theta_tilde) * std::cos(phi_tilde) * std::sin(theta_v) + std::cos(theta_tilde) * std::cos(theta_v);

        // convert to spherical coordinate
        double theta = std::acos(z);
        double phi = std::atan2(y, x);
        phi = (phi < 0.0) ? phi + PI * 2.0 : phi;

        double intensity = Intensity(Tobs, nu, theta, phi);
        return intensity;
    };

    // projection to the axis of LOS coordinate
    double projection = std::sqrt(x_tilde * x_tilde + y_tilde * y_tilde) * d * MPC * (1.0 + z) * (1.0 + z) * MAS;

    // solve azimuthal angle
    double phi_tilde = std::atan2(y_tilde, x_tilde);
    phi_tilde = (phi_tilde < 0.0) ? phi_tilde + PI * 2.0 : phi_tilde;

    // function to solve root and optimize
    auto f_root = [&](const double theta_tilde) {
        // convert to source coordinate (cartisan)
        double x = std::sin(theta_tilde) * std::cos(phi_tilde) * std::cos(theta_v) + std::cos(theta_tilde) * std::sin(theta_v);
        double y = std::sin(theta_tilde) * std::sin(phi_tilde);
        double z = - std::sin(theta_tilde) * std::cos(phi_tilde) * std::sin(theta_v) + std::cos(theta_tilde) * std::cos(theta_v);

        // convert to spherical coordinate
        double theta = std::acos(z);
        double phi = std::atan2(y, x);
        phi = (phi < 0.0) ? phi + PI * 2.0 : phi;

        // solve EATS
        eats->solveBlast(Tobs / (1 + this->z), theta, phi, theta_v, blast);

        // why do I use arcsinh? I don't even remember!
        return projection - blast.R * std::sin(theta_tilde);
    };

    // perform minimization (because there might be two intersections on two sides of theta_tilde_peak)
    double theta_tilde_peak = minimization(f_root, 0.0, PI, 10, 1e-4, "linear");

    // initialize intensity
    double intensity_tot = 0.0;

    // root function value at three critical points
    double f1 = f_root(0.0);
    double f2 = f_root(theta_tilde_peak);
    double f3 = f_root(PI);

    // line of sight and 2D surface may have 0 or 2 intersections
    if (f1 * f2 <= 0.0) {
        double theta_tilde = brentq(f_root, 0.0, theta_tilde_peak, 1e-6, 1e-6);
        double intensity = f_intensity(theta_tilde, phi_tilde);
        intensity_tot += intensity;
    }
    if (f2 * f3 <= 0.0) {
        double theta_tilde = brentq(f_root, theta_tilde_peak, PI, 1e-6, 1e-6);
        double intensity = f_intensity(theta_tilde, phi_tilde);
        intensity_tot += intensity;
    }

    return intensity_tot;
}
