#include "sim_box.h"

// ---------- public functions ---------- //
SimBox::SimBox(const JetConfig& jet_config, Tool& tool) {
    // member setup
    this->tool = &tool;
    cfl = jet_config.cfl;
    tmin = jet_config.tmin;
    tmax = jet_config.tmax;
    spread = jet_config.spread;

    // mesh
    ntheta = jet_config.Eb.size();
    theta_edge = jet_config.theta_edge;
    theta = Array(ntheta);
    for (int i = 0; i < ntheta; ++i) {
        theta[i] = (theta_edge[i] + theta_edge[i + 1]) / 2.0;
    }

    // conserved variables
    Eb = jet_config.Eb;
    Ht = jet_config.Ht;
    Msw = jet_config.Msw;
    Mej = jet_config.Mej;
    R = jet_config.R;

    // primitive variables
    beta_gamma_sq = Array(ntheta);
    beta_th = Array(ntheta);

    // variables for convinience
    beta = Array(ntheta);
    gamma = Array(ntheta);
    Psw = Array(ntheta);
    Hb = Array(ntheta);
    s = Array(ntheta);

    // eigenvalues
    eigenvalues = Array(ntheta);
    alpha_R = Array(ntheta);

    // solpe
    slope = Array(5, ntheta);
    R_slope_l = Array(ntheta);
    R_slope_r = Array(ntheta);

    // numerical flux
    numerical_flux = Array(4, ntheta + 1);
    dR_dt = Array(ntheta);

    // dy/dt
    dy_dt = Array(5, ntheta);

    // solve initial condition
    solvePrimitive();
    solveEigen();
}

void SimBox::solvePDE() {
    if (spread) {
        solveSpread();
    }
    else {
        solveNoSpread();
    }
}

Array3D& SimBox::getY() {
    return ys;
}

Array1D& SimBox::getT() {
    return ts;
}

Array1D& SimBox::getTheta() {
    return theta;
}

// ---------- private functions ---------- //
void SimBox::solvePrimitive() {
    for (int i = 0; i < ntheta; ++i) {
        // velocity
        try {
            beta_gamma_sq[i] = tool->solveBetaGammaSq(Msw[i] / Eb[i], Mej[i] / Eb[i], R[i]);
        }
        catch (const std::exception& e) {
            std::string text = "Hydro Primitive solver: ";
            throw std::runtime_error(text + e.what());
        }

        // convenient variables
        s[i] = tool->solveS(R[i], beta_gamma_sq[i]);
        gamma[i] = std::sqrt(beta_gamma_sq[i] + 1.0);
        beta[i] = std::sqrt(beta_gamma_sq[i] / (beta_gamma_sq[i] + 1.0));
        Psw[i] = s[i] * beta[i] * beta[i] * Msw[i] / 3.0;
        Hb[i] = Eb[i] + Psw[i];

        // tangent velocity
        beta_th[i] =  Ht[i] / Hb[i];
    }
}

void SimBox::solveEigen() {
    for (int i = 0; i < ntheta; ++i) {
        double A = 2.0 * s[i] / 3.0 * Msw[i] * (4.0 * gamma[i] * gamma[i] * gamma[i] * gamma[i] - 1.0) + ((1.0 - s[i]) * Msw[i] + Mej[i]) * gamma[i] * gamma[i] * gamma[i];
        double dPsw_dEb = 2.0 * s[i] / 3.0 * Msw[i] / A;
        double dPsw_dMsw = s[i] * beta[i] * beta[i] / 3.0 - 2.0 * s[i] / 3.0 * (Eb[i] - gamma[i] * Mej[i]) / A;
        double dPsw_dMej = - 2.0 * s[i] / 3.0 * gamma[i] * Msw[i] / A;

        double B = Mej[i] / Hb[i] * dPsw_dMej + Msw[i] / Hb[i] * dPsw_dMsw;
        double C = std::sqrt((1.0 - beta_th[i] * beta_th[i]) * (dPsw_dEb + B) + beta_th[i] * beta_th[i] / 4.0 * B * B);
        double alpha1 = beta_th[i];
        double alpha2 = beta_th[i] * (1.0 - B / 2.0) + C;
        double alpha3 = beta_th[i] * (1.0 - B / 2.0) - C;

        alpha1 = std::abs(alpha1);
        alpha2 = std::abs(alpha2);
        alpha3 = std::abs(alpha3);

        //eigenvalues[i] = std::max(std::max(alpha1, alpha2), alpha3) * CSpeed / R[i];
        eigenvalues[i] = std::max(std::max(alpha1, alpha2), alpha3) * CSpeed;
        alpha_R[i] = std::abs(beta_th[i]) / R[i];
    }
}

void SimBox::solveSlope() {
    // array for convinience
    std::vector<Array1D*> vars_ptr = {&Msw, &Mej, &beta_gamma_sq, &beta_th, &R};
    
    // variables to use
    int index1, index2;
    double diff1, diff2;
    double slope1, slope2;

    // msw, mej, beta_gamma_sq, beta_th
    for (int j = 0; j < 4; ++j) {
        Array1D& var = *(vars_ptr[j]);
        for (int i = 0; i < ntheta; ++i) {
            index1 = std::max(i - 1, 0);           // left cell index
            index2 = std::min(i + 1, ntheta - 1);  // right cell index
            diff1 = var[i] - var[index1];          // difference to left cell
            diff2 = var[index2] - var[i];          // difference to right cell
            slope1 = (i == index1) ? 0.0 : diff1 / (theta[i] - theta[index1]);  // left biased slope
            slope2 = (i == index2) ? 0.0 : diff2 / (theta[index2] - theta[i]);  // right biased slope
            slope[j][i] = tool->minmod(slope1, slope2);  // slope limiter
        }
    }

    // r
    for (int i = 0; i < ntheta; ++i) {
        index1 = std::max(i - 1, 0);           // left cell index (reflective boundary condition)
        index2 = std::min(i + 1, ntheta - 1);  // right cell index (reflective boundary condition)
        diff1 = R[i] - R[index1];          // difference to left cell
        diff2 = R[index2] - R[i];          // difference to right cell
        R_slope_l[i] = (i == index1) ? 0.0 : diff1 / (theta[i] - theta[index1]);  // left biased slope
        R_slope_r[i] = (i == index2) ? 0.0 : diff2 / (theta[index2] - theta[i]);  // right biased slope
        //slope[4][i] = (R_slope_l[i] + R_slope_r[i]) / 2.0;  // total slope
        slope[4][i] = tool->minmod(R_slope_l[i], R_slope_r[i]);
    }
}

void SimBox::solveNumericalFlux() {
    Array1D var_l = Array(5);        // left-biased reconstructed variable 
    Array1D var_r = Array(5);        // right-biased reconstructed variable
    Array1D F_l = Array(5);          // left-biased physical flux
    Array1D F_r = Array(5);          // right-biased physical flux
    double alpha;                    // numerical viscocity
    double Eb_l, Psw_l, Ht_l;        // left-biased conserved
    double Eb_r, Psw_r, Ht_r;        // right-biased conserved
    double s_l;                      // left biased calibration coefficient
    double s_r;                      // right biased calibration coefficient

    // alias
    double& Msw_l = var_l[0];
    double& Mej_l = var_l[1];
    double& beta_gamma_sq_l = var_l[2];
    double& beta_th_l = var_l[3];
    double& R_l = var_l[4];

    double& Msw_r = var_r[0];
    double& Mej_r = var_r[1];
    double& beta_gamma_sq_r = var_r[2];
    double& beta_th_r = var_r[3];
    double& R_r = var_r[4];

    // array for convinience
    std::vector<Array1D*> vars_ptr = {&Msw, &Mej, &beta_gamma_sq, &beta_th, &R};

    // only need to loop over middle faces. Polar faces have flux=0
    for (int i = 1; i < ntheta; ++i) {
        // reconstruct primitive variables
        for (int j = 0; j < 5; ++j) {
            Array1D& var = *(vars_ptr)[j];    // variable to reconstruct
            
            // left biased reconstruction
            var_l[j] = var[i - 1] + slope[j][i - 1] * (theta_edge[i] - theta[i - 1]);

            // right biased reconstruction
            var_r[j] = var[i] + slope[j][i] * (theta_edge[i] - theta[i]);
        }
        
        // solve calibration coefficients
        s_l = tool->solveS(R_l, beta_gamma_sq_l);
        s_r = tool->solveS(R_r, beta_gamma_sq_r);

        // solve left-biased conserved variables
        Eb_l = s_l * (1.0 + beta_gamma_sq_l * beta_gamma_sq_l / (beta_gamma_sq_l + 1) / (beta_gamma_sq_l + 1) / 3.0) * (beta_gamma_sq_l + 1) * Msw_l
             + (1.0 - s_l) * std::sqrt(beta_gamma_sq_l + 1) * Msw_l
             + std::sqrt(beta_gamma_sq_l + 1) * Mej_l;
        Psw_l = s_l * beta_gamma_sq_l / (beta_gamma_sq_l + 1) * Msw_l / 3.0;
        Ht_l = (Eb_l + Psw_l) * beta_th_l;

        // solve right-biased conserved variables
        Eb_r = s_r * (1.0 + beta_gamma_sq_r * beta_gamma_sq_r / (beta_gamma_sq_r + 1) / (beta_gamma_sq_r + 1) / 3.0) * (beta_gamma_sq_r + 1) * Msw_r 
             + (1.0 - s_r) * std::sqrt(beta_gamma_sq_r + 1) * Msw_r
             + std::sqrt(beta_gamma_sq_r + 1) * Mej_r;
        Psw_r = s_r * beta_gamma_sq_r / (beta_gamma_sq_r + 1) * Msw_r / 3.0;
        Ht_r = (Eb_r + Psw_r) * beta_th_r;

        // physical flux (left)
        F_l[0] = Ht_l / R_l * CSpeed;
        F_l[1] = (Ht_l * beta_th_l + Psw_l) / R_l * CSpeed;
        F_l[2] = Msw_l * beta_th_l / R_l * CSpeed;
        F_l[3] = Mej_l * beta_th_l / R_l * CSpeed;

        // physical flux (right)
        F_r[0] = Ht_r / R_r * CSpeed;
        F_r[1] = (Ht_r * beta_th_r + Psw_r) / R_r * CSpeed;
        F_r[2] = Msw_r * beta_th_r / R_r * CSpeed;
        F_r[3] = Mej_r * beta_th_r / R_r * CSpeed;

        // viscosity (maximum eigenvalue over 4 neighbor cells)
        int index_l = std::max(i - 2, 0);
        int index_r = std::min(i + 1, ntheta - 1);
        alpha = *std::max_element(eigenvalues.begin() + index_l, eigenvalues.begin() + index_r);

        // numerical flux
        numerical_flux[0][i] = 0.5 * (F_l[0] + F_r[0] - alpha * (Eb_r / R_r - Eb_l / R_l));
        numerical_flux[1][i] = 0.5 * (F_l[1] + F_r[1] - alpha * (Ht_r / R_r - Ht_l / R_l));
        numerical_flux[2][i] = 0.5 * (F_l[2] + F_r[2] - alpha * (Msw_r / R_r - Msw_l / R_l));
        numerical_flux[3][i] = 0.5 * (F_l[3] + F_r[3] - alpha * (Mej_r / R_r - Mej_l / R_l));
        
        //numerical_flux[0][i] = 0.5 * (F_l[0] + F_r[0] - alpha * (Eb_r - Eb_l));
        //numerical_flux[1][i] = 0.5 * (F_l[1] + F_r[1] - alpha * (Ht_r - Ht_l));
        //numerical_flux[2][i] = 0.5 * (F_l[2] + F_r[2] - alpha * (Msw_r - Msw_l));
        //numerical_flux[3][i] = 0.5 * (F_l[3] + F_r[3] - alpha * (Mej_r - Mej_l));
        for (int j = 0 ; j < 4; ++j) {
            numerical_flux[j][i] *= std::sin(theta_edge[i]);
        }
    }
}

double SimBox::solveDeltaT() {
    double omega, omega_all;
    double delta_t_min, delta_t;
    
    // find minimum value
    delta_t_min = std::numeric_limits<double>::max();
    for (int i = 0; i < ntheta; ++i) {
        // maximum possible omega
        omega = beta[i] * CSpeed / R[i];

        // combined signal speed and maximum omega
        //omega_all = eigenvalues[i] + 0.05 * omega;
        omega_all = eigenvalues[i] / R[i] + 0.05 * omega;
        //omega_all = omega;

        // delta_t of a cell
        delta_t = cfl * (theta_edge[i + 1] - theta_edge[i]) / omega_all;

        // find minimum
        delta_t_min = std::min(delta_t_min, delta_t);
    }

    return delta_t_min;
}

void SimBox::solveDyDt() {
    for (int i = 0; i < ntheta; ++i) {
        // solve needed variables
        double beta_f = 4.0 * beta[i] * (beta_gamma_sq[i] + 1) / (4.0 * beta_gamma_sq[i] + 3.0);
        double beta_r = (beta_th[i] <= beta[i]) ? std::sqrt(beta[i] * beta[i] - beta_th[i] * beta_th[i]) : 0.0;
        double vol = std::cos(theta_edge[i]) - std::cos(theta_edge[i + 1]);

        // solve dR_dt
        int index_l = std::max(i - 1, 0);
        int index_r = std::min(i + 1, ntheta - 1);
        double alpha = std::max(std::max(alpha_R[index_l], alpha_R[i]), alpha_R[index_r]);
        dy_dt[4][i] = beta_f - slope[4][i] * beta_th[i] / R[i] + 0.5 * alpha * (R_slope_r[i] - R_slope_l[i]);
        dy_dt[4][i] *= CSpeed;

        // conserved variables
        double rho = tool->solveDensity(R[i]) * MassP;
        dy_dt[0][i] = (numerical_flux[0][i] - numerical_flux[0][i + 1]) / vol
                    + dy_dt[4][i] * rho * R[i] * R[i];
        dy_dt[1][i] = (numerical_flux[1][i] - numerical_flux[1][i + 1]) / vol
                    + (std::cos(theta[i]) / std::sin(theta[i]) * Psw[i] - Ht[i] * beta_r) * CSpeed / R[i];
        dy_dt[2][i] = (numerical_flux[2][i] - numerical_flux[2][i + 1]) / vol
                    + dy_dt[4][i] * rho * R[i] * R[i];
        dy_dt[3][i] = (numerical_flux[3][i] - numerical_flux[3][i + 1]) / vol
                    + 0.0;
    }
}

void SimBox::solveDyDt_no_spread() {
    for (int i = 0; i < ntheta; ++i) {
        // solve needed variables
        double beta_f = 4.0 * beta[i] * (beta_gamma_sq[i] + 1) / (4.0 * beta_gamma_sq[i] + 3.0);

        // solve dR_dt
        dy_dt[4][i] = beta_f * CSpeed;

        // conserved variables
        double rho = tool->solveDensity(R[i]) * MassP;
        dy_dt[0][i] = dy_dt[4][i] * rho * R[i] * R[i];
        dy_dt[1][i] = 0.0;
        dy_dt[2][i] = dy_dt[4][i] * rho * R[i] * R[i];
        dy_dt[3][i] = 0.0;
    }
}

void SimBox::oneStepRK2(double dt) {
    // alias
    std::vector<Array1D*> conserved = {&Eb, &Ht, &Msw, &Mej, &R};

    // copy initial conserved variables
    Array2D conserved_ini = {Eb, Ht, Msw, Mej, R};

    // ---------- step 1 ---------- //
    // solve dydy
    solveSlope();
    solveNumericalFlux();
    solveDyDt();

    // update variables
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            (*conserved[i])[j] += dt * dy_dt[i][j];
        }
    }

    // update primitive & eigenvalues
    solvePrimitive();
    solveEigen();

    // ---------- step 2 ---------- //
    // solve dydy
    solveSlope();
    solveNumericalFlux();
    solveDyDt();

    // update variables
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            (*conserved[i])[j] = 0.5 * conserved_ini[i][j] + 0.5 * (*conserved[i])[j] + 0.5 * dt * dy_dt[i][j];
        }
    }

    // update primitive & eigenvalues
    solvePrimitive();
    solveEigen();
}

void SimBox::oneStepRK45(double& dt, const double rtol, bool& succeeded) {
    // alias of conserved variables
    std::vector<Array1D*> conserved = {&Eb, &Ht, &Msw, &Mej, &R};

    // save initial conserved variables
    Array2D conserved_ini = {Eb, Ht, Msw, Mej, R};

    // middle steps
    Array2D k1, k2, k3, k4, k5, k6;
    k1 = k2 = k3 = k4 = k5 = k6 = Array(5, ntheta);

    // ---------- step 1 ---------- //
    // solve dydt
    solveDyDt_no_spread();

    // update variables
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            k1[i][j] = dt * dy_dt[i][j];
            (*conserved[i])[j] = conserved_ini[i][j] + k1[i][j] * 2.0 / 9.0;
        }
    }

    // update primitives
    solvePrimitive();

    // ---------- step 2 ---------- //
    // solve dydt
    solveDyDt_no_spread();

    // update variables
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            k2[i][j] = dt * dy_dt[i][j];
            (*conserved[i])[j] = conserved_ini[i][j] + k1[i][j] / 12.0 + k2[i][j] / 4.0;
        }
    }

    // update primitives
    solvePrimitive();

    // ---------- step 3 ---------- //
    // solve dydt
    solveDyDt_no_spread();

    // update variables
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            k3[i][j] = dt * dy_dt[i][j];
            (*conserved[i])[j] = conserved_ini[i][j] + k1[i][j] * 69.0 / 128.0 - k2[i][j] * 243.0 / 128.0 + k3[i][j] * 135.0 / 64.0;
        }
    }

    // update primitives
    solvePrimitive();

    // ---------- step 4 ---------- //
    // solve dydt
    solveDyDt_no_spread();

    // update variables
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            k4[i][j] = dt * dy_dt[i][j];
            (*conserved[i])[j] = conserved_ini[i][j] - k1[i][j] * 17.0 / 12.0 + k2[i][j] * 27.0 / 4.0 - k3[i][j] * 27.0 / 5.0 + k4[i][j] * 16.0 / 15.0;
        }
    }

    // update primitives
    solvePrimitive();

    // ---------- step 5 ---------- //
    // solve dydt
    solveDyDt_no_spread();

    // update variables
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            k5[i][j] = dt * dy_dt[i][j];
            (*conserved[i])[j] = conserved_ini[i][j] + k1[i][j] * 65.0 / 432.0 - k2[i][j] * 5.0 / 16.0 + k3[i][j] * 13.0 / 16.0 + k4[i][j] * 4.0 / 27.0 + k5[i][j] * 5.0 / 144.0;
        }
    }

    // update primitives
    solvePrimitive();

    // ---------- step 6 ---------- //
    // solve dydt
    solveDyDt_no_spread();

    // update variables
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            k6[i][j] = dt * dy_dt[i][j];
            (*conserved[i])[j] = conserved_ini[i][j] + k1[i][j] * 47.0 / 450.0 + k3[i][j] * 12.0 / 25.0 + k4[i][j] * 32.0 / 225.0 + k5[i][j] / 30.0 + k6[i][j] * 6.0 / 25.0;
        }
    }

    // ---------- error estimate ---------- //
    double error;
    double rerror = 0.0;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            error = std::abs(k1[i][j] / 150.0 - k3[i][j] * 3.0 / 100.0 + k4[i][j] * 16.0 / 75.0 + k5[i][j] / 20.0 - k6[i][j] * 6.0 / 25.0);
            rerror = std::max(rerror, error / std::abs((*conserved[i])[j]));
        }
    }

    if (rerror < rtol) {
        // good! update primitives
        solvePrimitive();

        // mark
        succeeded = true;
    }
    else {
        // roll back to the initial conserved variables
        for (int i = 0; i < 5; ++i) {
            *conserved[i] = conserved_ini[i];
        }

        // update primitives
        solvePrimitive();

        // mark
        succeeded = false;
    }

    // update dt
    double boost_factor = 0.9 * std::pow(rtol / rerror, 0.2);
    boost_factor = std::min(1.5, boost_factor);
    dt *= boost_factor;
}

void SimBox::solveSpread() {
    // alias
    std::vector<Array1D*> primitives = {&Msw, &Mej, &beta_gamma_sq, &beta_th, &R};

    // record initial condition
    ts.push_back(tmin);
    ys = Array(5, ntheta, 1);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            ys[i][j][0] = (*primitives[i])[j];
        }
    }

    // solve PDE
    double dt;
    double t = tmin;
    while (t < tmax) {
        // interruption detection
        if (PyErr_CheckSignals() != 0) {
            throw std::runtime_error("Hydro: Keyboard interruption.");
        }

        // solve one step
        dt = std::min(solveDeltaT(), tmax - t + 1e-6);
        oneStepRK2(dt);
        t += dt;

        // save variables
        ts.push_back(t);
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < ntheta; ++j) {
                ys[i][j].push_back((*primitives[i])[j]);
            }
        }
    }
}

void SimBox::solveNoSpread() {
    // alias
    std::vector<Array1D*> primitives = {&Msw, &Mej, &beta_gamma_sq, &beta_th, &R};

    // record initial condition
    ts.push_back(tmin);
    ys = Array(5, ntheta, 1);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            ys[i][j][0] = (*primitives[i])[j];
        }
    }

    // solve PDE
    double delta_t = 1.0;
    double t = tmin;
    while (t < tmax) {
        // interruption detection
        if (PyErr_CheckSignals() != 0) {
            throw std::runtime_error("Hydro: Keyboard interruption.");
        }

        // solve one step
        bool succeeded;
        double dt = delta_t;
        oneStepRK45(dt, 1e-6, succeeded);
        
        if (succeeded) {
            t += delta_t;
            delta_t = dt;
            
            // save variables
            ts.push_back(t);
            for (int i = 0; i < 5; ++i) {
                for (int j = 0; j < ntheta; ++j) {
                    ys[i][j].push_back((*primitives[i])[j]);
                }
            }
        }
        else {
            delta_t = dt;
        }
    }
}