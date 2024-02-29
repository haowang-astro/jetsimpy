#ifndef INTERPOLATION
#define INTERPOLATION

#include "Definations.hpp"

class Interpolator {
public:
    Interpolator(const Array1D& t_data, const Array1D& theta_data, const Array3D& y_pri_data, const double& nwind, const double& nism);
    Interpolator() {}

    double Msw(const double& t, const double& theta);
    double Mej(const double& t, const double& theta);
    double beta_gamma(const double& t, const double& theta);
    double beta_th(const double& t, const double& theta);
    double R(const double& t, const double& theta);
    double s(const double& t, const double& theta);
    double Eb0(const double& t, const double& theta);

private:
    const Array1D* t_data;
    const Array1D* theta_data;
    const Array3D* y_pri_data;
    double nwind, nism;
    double tmin, tmax;
    double theta_min, theta_max;
    int ntheta;
    double interpolate(const double& t, const double& theta, const int& pde_index);
    double interpolate(const double& theta, const int& pde_index);
};

// -------------------------------- Jet1D interpolator ----------------------------------- //
Interpolator::Interpolator(const Array1D& t_data, const Array1D& theta_data, const Array3D& y_pri_data, const double& nwind, const double& nism) {
    this->t_data = &t_data;
    this->theta_data = &theta_data;
    this->y_pri_data = &y_pri_data;
    this->nwind = nwind;
    this->nism = nism;
    this->tmin = t_data.front();
    this->tmax = t_data.back() - 1.0;
    this->theta_min = theta_data.front();
    this->theta_max = theta_data.back();
    this->ntheta = theta_data.size();
}

double Interpolator::interpolate(const double& t, const double& theta, const int& pde_index) {
    // find t, theta index
    int t_index1, t_index2;
    int theta_index1, theta_index2;
    
    const Array1D& ts = *t_data;
    const Array1D& thetas = *theta_data;
    const Array3D& ys = *y_pri_data;

    // find t index
    Find_Index(ts, t, t_index1, t_index2);
    
    // interpolate
    if (theta <= theta_min) {
        double result = Linear(ts[t_index1], ts[t_index2],
                      ys[t_index1][pde_index][0], ys[t_index2][pde_index][0],
                      t);
        return result;
    }
    else if (theta >= theta_max) {
        double result = Linear(ts[t_index1], ts[t_index2],
                      ys[t_index1][pde_index][ntheta - 1], ys[t_index2][pde_index][ntheta - 1],
                      t);
        return result;
    }
    else {
        // find theta index
        Find_Index(thetas, theta, theta_index1, theta_index2);

        double result = Bilinear(ts[t_index1], ts[t_index2],
                      thetas[theta_index1], thetas[theta_index2],
                      ys[t_index1][pde_index][theta_index1],
                      ys[t_index1][pde_index][theta_index2],
                      ys[t_index2][pde_index][theta_index1],
                      ys[t_index2][pde_index][theta_index2],
                      t, theta);

        return result;
    }
}

double Interpolator::interpolate(const double& theta, const int& pde_index) {
    int theta_index1, theta_index2;

    const Array1D& thetas = *theta_data;
    const Array3D& ys = *y_pri_data;
    
    // interpolate
    if (theta <= theta_min) {
        return ys[0][pde_index][0];
    }
    else if (theta>= theta_max) {
        return ys[0][pde_index][ntheta - 1];
    }
    else {
        // find theta index
        Find_Index(thetas, theta, theta_index1, theta_index2);
        double result = Linear(thetas[theta_index1], thetas[theta_index2],
                      ys[0][pde_index][theta_index1], ys[0][pde_index][theta_index2],
                      theta);
        return result;
    }
}

double Interpolator::Msw(const double& t, const double& theta) {
    int pde_index = 0;
    double r, msw;
    r = R(t, theta);

    // bound check
    if (t < 0.0 || t > tmax || theta < 0.0 || theta > CONST.pi) {
        // throw error
        throw std::runtime_error("a value is out of the interpolation range.");
    }
    // find if t < tmin (constant)
    else if (t < tmin) {
        // calculate extrapolation
        msw = nwind * CONST.mp * r / 1e17 * 1e51 + nism * CONST.mp * r * r * r / 3.0;
    }
    else {
        msw = interpolate(t, theta, pde_index);
    }

    return msw;
}

double Interpolator::Mej(const double& t, const double& theta) {
    int pde_index = 1;
    double mej;

    // bound check
    if (t < 0.0 || t > tmax || theta < 0.0 || theta > CONST.pi) {
        // throw error
        throw std::runtime_error("a value is out of the interpolation range.");
    }
    // find if t < tmin (constant)
    else if (t < tmin) {
        mej = interpolate(theta, pde_index);
    }
    else {
        mej = interpolate(t, theta, pde_index);
    }

    return mej;
}

// The spatial interpolation of u must be cubic, in order to increase accuracy
double Interpolator::beta_gamma(const double& t, const double& theta) {
    int pde_index = 2;
    double beta_gamma_sq;

    // bound check
    if (t < 0.0 || t > tmax || theta < 0.0 || theta > CONST.pi) {
        // throw error
        throw std::runtime_error("a value is out of the interpolation range.");
    }
    // find if t < tmin (constant)
    else if (t < tmin) {
        beta_gamma_sq = interpolate(theta, pde_index);
    }
    else {
       beta_gamma_sq = interpolate(t, theta, pde_index);
    }

    return std::sqrt(beta_gamma_sq);
}

double Interpolator::beta_th(const double& t, const double& theta) {
    int pde_index = 3;
    double beta_theta;

    // bound check
    if (t < 0.0 || t > tmax || theta < 0.0 || theta > CONST.pi) {
        // throw error
        throw std::runtime_error("a value is out of the interpolation range.");
    }
    // find if t < tmin (constant)
    else if (t < tmin) {
        beta_theta = 0.0;
    }
    else {
        beta_theta = interpolate(t, theta, pde_index);
    }

    return beta_theta;
}

// The spatial interpolation of u must be cubic, in order to increase accuracy
double Interpolator::R(const double& t, const double& theta) {
    int pde_index = 4;
    double r;

    // bound check
    if (t < 0.0 || t > tmax || theta < 0.0 || theta > CONST.pi) {
        // throw error
        throw std::runtime_error("a value is out of the interpolation range.");
    }
    // find if t < tmin (constant)
    else if (t < tmin) {
        r = interpolate(theta, pde_index);
        r = t / tmin * r;
    }
    else {
        r = interpolate(t, theta, pde_index);
    }

    return r;
}

double Interpolator::s(const double& t, const double& theta) {
    int pde_index = 5;
    double ss;

    // bound check
    if (t < 0.0 || t > tmax || theta < 0.0 || theta > CONST.pi) {
        // throw error
        throw std::runtime_error("a value is out of the interpolation range.");
    }
    // find if t < tmin (constant)
    else if (t < tmin) {
        ss = interpolate(theta, pde_index);
    }
    else {
        ss = interpolate(t, theta, pde_index);
    }

    return ss;
}

double Interpolator::Eb0(const double& t, const double& theta) {
    double msw = Msw(t, theta);
    double mej = Mej(t, theta);
    double r = R(t, theta);
    double beta_gamma_sq = std::pow(beta_gamma(t, theta), 2);
    double beta_sq = beta_gamma_sq / (beta_gamma_sq + 1);
    double gamma = std::sqrt(beta_gamma_sq + 1.0);
    double s = this->s(t, theta);
    return s * gamma * gamma * (1.0 + beta_sq * beta_sq / 3.0) * msw + (1.0 - s) * gamma * msw + gamma * mej - msw - mej;
}

#endif