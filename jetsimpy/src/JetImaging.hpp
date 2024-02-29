// This file do what

#ifndef JETIMAGING
#define JETIMAGING

#include "Definations.hpp"
#include "Jet.hpp"
#include "JetIntegrator.hpp"

struct AfterglowImage {
    double half_width;            // [mas]
    double xscale;                // [mas]
    double yscale;                // [mas]
    double offset;                // [mas]
    Array2D intensity_map;        // [inum, inum] (Inu)
    Array2D pi_scale_map;         // [pnum, pnum] (%)
    Array2D pi_angle_map;         // [pnum, pnum] (radian)
};

class JetImaging {
public:
    JetIntegrator* jet_integrator;
    AfterglowImage image;

    JetImaging() {}
    void config(JetIntegrator* jet_integrator);
    void solve_image(const double& Tobs, const double& nu, const int& inum, const int& pnum);
private:

};

// ---------- implementations ---------- //

void JetImaging::config(JetIntegrator* jet_integrator) {
    this->jet_integrator = jet_integrator;
}

void JetImaging::solve_image(const double& Tobs, const double& nu, const int& inum, const int& pnum) {
    // create an alias
    Dict& para_rad = jet_integrator->para_rad;

    double redshift = para_rad.at("z");
    double theta_v = para_rad.at("theta_v");
    double dl = para_rad.at("d");
    double Tobs_z = Tobs / (1.0 + redshift);
    double nu_z = nu * (1.0 + redshift);

    double dxy, dxy_pi, x_tilde, y_tilde, theta_tilde, phi_tilde, projection;

    // initialize the containers
    image.intensity_map = Array(inum, inum);
    image.pi_scale_map = Array(pnum, pnum);
    image.pi_angle_map = Array(pnum, pnum);
    
    // find scales
    //image.xscale = jet_integrator->Xscale(Tobs, nu, 0.0, 1e-2);
    //image.yscale = jet_integrator->Yscale(Tobs, nu, 0.0, 1e-2);
    image.xscale = jet_integrator->SigmaX(Tobs, nu, 0.0, 1e-2);
    image.yscale = jet_integrator->SigmaY(Tobs, nu, 0.0, 1e-2);
    image.offset = jet_integrator->Offset(Tobs, nu, 0.0, 1e-2);

    // setup image size
    double factor = 1.1;
    image.half_width = (image.offset + image.xscale * 2) * factor;
    dxy = image.half_width * 2.0 / (inum - 1);
    dxy_pi = image.half_width * 2.0 / (pnum - 1);
    
    // function to solve intensity & polarization
    Array1D polar_array = Array(3);
    auto f_intensity = [&](const double& theta_tilde, const double& phi_tilde) -> Array1D {
        // convert to source coordinate (cartisan)
        double x = std::sin(theta_tilde) * std::cos(phi_tilde) * std::cos(theta_v) + std::cos(theta_tilde) * std::sin(theta_v);
        double y = std::sin(theta_tilde) * std::sin(phi_tilde);
        double z = - std::sin(theta_tilde) * std::cos(phi_tilde) * std::sin(theta_v) + std::cos(theta_tilde) * std::cos(theta_v);

        // convert to spherical coordinate
        double theta = std::acos(z);
        double phi = std::atan2(y, x);
        phi = (phi < 0.0) ? phi + CONST.pi * 2.0 : phi;

        double intensity = jet_integrator->Intensity(Tobs_z, nu_z, theta, phi);
        double Pi_lin = jet_integrator->S.Pi_lin;
        double Pi_circ = jet_integrator->S.Pi_circ;

        return {intensity, Pi_lin, Pi_circ};
    };
    
    // solve intensity map
    for (int j = 0; j < inum; ++j) {
        // set x grid points
        x_tilde = (- image.half_width + j * dxy);
        for (int k = 0; k < inum; ++k) {
            // interruption detection
            if (PyErr_CheckSignals() != 0) {
                throw std::runtime_error("Interrupted at image calculation.");
            }

            // set y grid points
            y_tilde = (- image.half_width + k * dxy);

            // solve azimuthal angle
            phi_tilde = std::atan2(y_tilde, x_tilde);
            phi_tilde = (phi_tilde < 0.0) ? phi_tilde + CONST.pi * 2.0 : phi_tilde;

            // offset of this point
            projection = std::sqrt(x_tilde * x_tilde + y_tilde * y_tilde) * dl * CONST.mpc * (1.0 + redshift) * (1.0 + redshift) * CONST.mas;

            // function to solve root and optimize
            auto f_root = [&](const double& theta_tilde) {
                // convert to source coordinate (cartisan)
                double x = std::sin(theta_tilde) * std::cos(phi_tilde) * std::cos(theta_v) + std::cos(theta_tilde) * std::sin(theta_v);
                double y = std::sin(theta_tilde) * std::sin(phi_tilde);
                double z = - std::sin(theta_tilde) * std::cos(phi_tilde) * std::sin(theta_v) + std::cos(theta_tilde) * std::cos(theta_v);

                // convert to spherical coordinate
                double theta = std::acos(z);
                double phi = std::atan2(y, x);
                phi = (phi < 0.0) ? phi + CONST.pi * 2.0 : phi;

                // solve EATS
                jet_integrator->get_blast(Tobs_z, theta, phi);

                return std::asinh(projection) - std::asinh(jet_integrator->blast.R * std::sin(theta_tilde));
            };

            // perform minimization
            double theta_tilde_peak = minimization(f_root, 0.0, CONST.pi, 10, 1e-4, "linear");
            
            // initialize intensity
            image.intensity_map[j][k] = 0.0;
            
            // solve intensity
            double f1 = f_root(0.0);
            double f2 = f_root(theta_tilde_peak);
            double f3 = f_root(CONST.pi);
            // line of sight and 2D surface have 0 or 2 intersections
            if (f1 * f2 <= 0.0) {
                theta_tilde = brentq(f_root, 0.0, theta_tilde_peak, 1e-6, 1e-6);
                polar_array = f_intensity(theta_tilde, phi_tilde);
                image.intensity_map[j][k] += polar_array[0];
            }
            if (f2 * f3 <= 0.0) {
                theta_tilde = brentq(f_root, theta_tilde_peak, CONST.pi, 1e-6, 1e-6);
                polar_array = f_intensity(theta_tilde, phi_tilde);
                image.intensity_map[j][k] += polar_array[0];
            }
        }
    }

    // solve polarization map
    for (int j = 0; j < pnum; ++j) {
        // set x grid points
        x_tilde = (- image.half_width + j * dxy_pi);
        for (int k = 0; k < pnum; ++k) {
            // interruption detection
            if (PyErr_CheckSignals() != 0) {
                throw std::runtime_error("Interrupted at image calculation.");
            }

            // set y grid points
            y_tilde = (- image.half_width + k * dxy_pi);

            // solve azimuthal angle
            phi_tilde = std::atan2(y_tilde, x_tilde);
            phi_tilde = (phi_tilde < 0.0) ? phi_tilde + CONST.pi * 2.0 : phi_tilde;

            // offset of this point
            projection = std::sqrt(x_tilde * x_tilde + y_tilde * y_tilde) * dl * CONST.mpc * (1.0 + redshift) * (1.0 + redshift) * CONST.mas;

            // function to solve root and optimize
            auto f_root = [&](const double& theta_tilde) {
                // convert to source coordinate (cartisan)
                double x = std::sin(theta_tilde) * std::cos(phi_tilde) * std::cos(theta_v) + std::cos(theta_tilde) * std::sin(theta_v);
                double y = std::sin(theta_tilde) * std::sin(phi_tilde);
                double z = - std::sin(theta_tilde) * std::cos(phi_tilde) * std::sin(theta_v) + std::cos(theta_tilde) * std::cos(theta_v);

                // convert to spherical coordinate
                double theta = std::acos(z);
                double phi = std::atan2(y, x);
                phi = (phi < 0.0) ? phi + CONST.pi * 2.0 : phi;

                // solve EATS
                jet_integrator->get_blast(Tobs_z, theta, phi);

                return std::asinh(projection) - std::asinh(jet_integrator->blast.R * std::sin(theta_tilde));
            };

            // perform minimization
            double theta_tilde_peak = minimization(f_root, 0.0, CONST.pi, 10, 1e-4, "linear");
            
            // initialize polarization
            image.pi_angle_map[j][k] = 0.0;
            image.pi_scale_map[j][k] = 0.0;
            double qu = 0.0;
            double intensity = 0.0;
            
            // solve polarization
            double f1 = f_root(0.0);
            double f2 = f_root(theta_tilde_peak);
            double f3 = f_root(CONST.pi);
            // line of sight and 2D surface have 0 or 2 intersections
            if (f1 * f2 <= 0.0) {
                theta_tilde = brentq(f_root, 0.0, theta_tilde_peak, 1e-6, 1e-6);
                polar_array = f_intensity(theta_tilde, phi_tilde);
                qu += polar_array[0] * polar_array[1];
                intensity += polar_array[0];
            }
            if (f2 * f3 <= 0.0) {
                theta_tilde = brentq(f_root, theta_tilde_peak, CONST.pi, 1e-6, 1e-6);
                polar_array = f_intensity(theta_tilde, phi_tilde);
                qu += polar_array[0] * polar_array[1];
                intensity += polar_array[0];
            }
            image.pi_scale_map[j][k] = (intensity == 0.0) ? 0.0 : qu / intensity;

            // determine the rotation angle of polarization
            if (image.pi_scale_map[j][k] > 0.0) {
                image.pi_angle_map[j][k] = phi_tilde + CONST.pi / 2.0;
                if (image.pi_angle_map[j][k] > CONST.pi * 2.0) {
                    image.pi_angle_map[j][k] -= CONST.pi * 2.0;
                }
            }
            else if (image.pi_scale_map[j][k] < 0.0) {
                image.pi_angle_map[j][k] = phi_tilde;
            }
            else {
                image.pi_angle_map[j][k] = 0.0;
            }
        }
    }
}

#endif