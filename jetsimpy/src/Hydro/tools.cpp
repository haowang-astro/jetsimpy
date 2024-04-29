#include "tools.h"

Tool::Tool(const JetConfig& jet_config) {
    nwind = jet_config.nwind;
    nism = jet_config.nism;
    rtol = jet_config.rtol;
    calib_level = jet_config.calib_level;
}

double Tool::solveDensity(double r) {
    return nwind / (r / 1e17) / (r / 1e17) + nism;
}

double Tool::solveS(double r, double beta_gamma_sq) {
    // density slope
    double k = 2.0 * nwind / (nwind + nism * (r / 1e17) * (r / 1e17));

    // calibration coefficient
    double s;
    if (calib_level == 0) {
        s = 1.0;
    }
    else if (calib_level == 1) {
        s = 0.52935729 - 0.05698377 * k - 0.00158176 * k * k - 0.00939548 * k * k * k;
    }
    else if (calib_level == 2) {
        double sBM = 0.52935729 - 0.05698377 * k - 0.00158176 * k * k - 0.00939548 * k * k * k;
        double sST = 1.635 - 0.651 * k;
        s = (sST + sBM * factor * beta_gamma_sq) / (1.0 + factor * beta_gamma_sq);
    }
    else {
        throw std::runtime_error("Hydro: invalid calibration level!");
    }

    return s;
}

double Tool::solveBetaGammaSq(double msw_eb, double mej_eb, double r) {
    double beta_gamma_min_sq = 0.0; //(0.75 / (msw_eb + mej_eb) - 1) * 0.99;
    double gamma_max = 1.0 / (msw_eb + mej_eb);
    double beta_gamma_max_sq = (gamma_max * gamma_max - 1.0) * 1.01;

    // solve root
    auto f = [&](const double& u_sq) {
        double beta_sq = u_sq / (u_sq + 1.0);
        double gamma = std::sqrt(u_sq + 1.0);
        double s = solveS(r, u_sq);
        return s * gamma * gamma * (1.0 + beta_sq * beta_sq / 3.0) * msw_eb + gamma * ((1.0 - s) * msw_eb + mej_eb) - 1.0;
    };

    double u_sq = brentq(f, beta_gamma_min_sq, beta_gamma_max_sq, 0.0, rtol);
    return u_sq;
}

double Tool::minmod(double x1, double x2) {
    if (x1 * x2 > 0) {
        if (std::fabs(x1) < std::fabs(x2)) {
            return x1;
        }
        else {
            return x2;
        }
    }
    else {
        return 0;
    }
}

void Tool::findIndex(const Array1D& x_array, const double x, int& index1, int& index2) {
    index1 = 0;
    index2 = x_array.size() - 1;

    int index_mid;
    while (index2 - index1 > 1) {
        index_mid = (index1 + index2) / 2;
        if (x > x_array[index_mid]) {
            index1 = index_mid;
        }
        else {
            index2 = index_mid;
        }
    }
}

double Tool::linear(double x, double x1, double x2, double y1, double y2) {
    return (x1 == x2) ? y1 : (y2 - y1) / (x2 - x1) * (x - x1) + y1;
}