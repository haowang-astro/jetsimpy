// this file do what

#ifndef JET
#define JET

#include "Definations.hpp"
#include "Interpolation.hpp"
#include "RootSolver.hpp"
#include "TimeIntegration.hpp"
#include "EqualArrivalTimeSurface.hpp"

struct PDE_Solver {
    // external medium
    double nwind, nism;

    // velocity solver tolerance
    double rtol;

    // cells
    Array1D theta;
    int ntheta;
    Array1D theta_edge;
    Array1D voltheta;

    // conserved
    Array1D Eb;
    Array1D Ht;
    Array1D Msw;
    Array1D Mej;
    Array1D R;
    
    // primitive
    Array1D Psw;
    Array1D Hb;
    Array1D beta_gamma_sq;
    Array1D beta;
    Array1D gamma;
    Array1D beta_th;
    Array1D beta_r;
    Array1D beta_f;

    // eigenvalues
    Array1D eigenvalues;
    Array1D alpha_R;

    // calibration
    bool calibration = false;
    Array1D sST;
    Array1D sBM;
    Array1D s;
    double factor;

    // delta_t
    double dt;

    // slope
    Array1D Msw_slope;
    Array1D Mej_slope;
    Array1D beta_gamma_sq_slope;
    Array1D beta_th_slope;
    Array1D R_slope;
    Array1D R_slope1;
    Array1D R_slope2;

    // numerical flux
    Array2D numerical_flux;

    // performance test
    double time1 = 0.0;
    double time2 = 0.0;

    // configuration
    void config(const double& nwind, const double& nism, const int& ntheta, const Array1D& theta, const Array1D& theta_edge, const Array1D& voltheta, const double& rtol);
    void solve_delta_t(const double& cfl);
    void solve_dy_dt(const Array2D& y, Array2D& dy_dt);
    void solve_dy_dt_no_spreading(const Array2D& y, Array2D& dy_dt);

private:
    void update_conserved(const Array2D& y);
    double solve_beta_gamma_sq(const double& msw_eb, const double& mej_eb, const double& sbm, const double& sst);
    void solve_state();
    void solve_eigen();
    void solve_slope();
    void solve_numerical_flux();
};

void PDE_Solver::config(const double& nwind, const double& nism, const int& ntheta, const Array1D& theta, const Array1D& theta_edge, const Array1D& voltheta, const double& rtol) {
    // external medium
    this->nwind = nwind;
    this->nism = nism;

    // velocity solver tolerance
    this->rtol = rtol;

    // cells
    this->ntheta = ntheta;
    this->theta = theta;
    this->theta_edge = theta_edge;
    this->voltheta = voltheta;

    // conserved
    Eb = Array(ntheta);
    Msw = Array(ntheta);
    Mej = Array(ntheta);
    Ht = Array(ntheta);
    R = Array(ntheta);

    // primitive
    Hb = Array(ntheta);
    Psw = Array(ntheta);
    beta_gamma_sq = Array(ntheta);
    beta = Array(ntheta);
    gamma = Array(ntheta);
    beta_th = Array(ntheta);
    beta_r = Array(ntheta);
    beta_f = Array(ntheta);

    // eigenvalues
    eigenvalues = Array(ntheta);
    alpha_R = Array(ntheta);

    // calibration
    sST = Array(ntheta, 1.0);
    sBM = Array(ntheta, 1.0);
    s = Array(ntheta, 1.0);
    factor = 2.0;

    // slope
    Msw_slope = Array(ntheta);
    Mej_slope = Array(ntheta);
    beta_gamma_sq_slope = Array(ntheta);
    beta_th_slope = Array(ntheta);
    R_slope = Array(ntheta);
    R_slope1 = Array(ntheta);
    R_slope2 = Array(ntheta);

    // numerical flux
    numerical_flux = Array(4, ntheta + 1);
}

void PDE_Solver::update_conserved(const Array2D& y) {
    Eb = y[0];
    Ht = y[1];
    Msw = y[2];
    Mej = y[3];
    R = y[4];
}

double PDE_Solver::solve_beta_gamma_sq(const double& msw_eb, const double& mej_eb, const double& sbm, const double& sst) {
    // newtonian
    //double beta_gamma_sq_newton = (1.0 - msw_eb - mej_eb) / (msw_eb + 0.5 * mej_eb);
    //if (beta_gamma_sq_newton < 0.1) {
    //    return beta_gamma_sq_newton;
    //}

    double beta_gamma_min_sq = 0.0;//(0.75 / (msw_eb + mej_eb) - 1) * 0.99;
    double gamma_max = 1.0 / (msw_eb + mej_eb);
    double beta_gamma_max_sq = (gamma_max * gamma_max - 1.0) * 1.01;

    // solve root
    auto f = [&](const double& u_sq) {
        double beta_sq = u_sq / (u_sq + 1.0);
        double gamma = std::sqrt(u_sq + 1.0);
        double s = (sst + sbm * factor * u_sq) / (1.0 + factor * u_sq);
        return s * gamma * gamma * (1.0 + beta_sq * beta_sq / 3.0) * msw_eb + gamma * ((1.0 - s) * msw_eb + mej_eb) - 1.0;
    };

    double u_sq = brentq(f, beta_gamma_min_sq, beta_gamma_max_sq, 0.0, rtol);
    return u_sq;
}

void PDE_Solver::solve_state() {
    // calibration factors
    if (calibration) {
        for (int i = 0; i < ntheta; ++i) {
            sBM[i] = calib_BM(nwind, nism, R[i]);
            sST[i] = calib_ST(nwind, nism, R[i]);
            //sST[i] = calib_BM(nwind, nism, R[i]);
            //sST[i] = 1.0;
        }
    }

    // solve velocity and calibration factor
    for (int i = 0; i < ntheta; ++i) {
        // solve velocity
        beta_gamma_sq[i] = solve_beta_gamma_sq(Msw[i] / Eb[i], Mej[i] / Eb[i], sBM[i], sST[i]);
        s[i] = (sST[i] + sBM[i] * factor * beta_gamma_sq[i]) / (1.0 + factor * beta_gamma_sq[i]);
    }
    
    // solve primitive variables
    double beta_sq;
    for (int i = 0; i < ntheta; ++i) {
        // velocity
        gamma[i] = std::sqrt(beta_gamma_sq[i] + 1.0);
        beta_sq = beta_gamma_sq[i] / (beta_gamma_sq[i] + 1.0);
        beta[i] = std::sqrt(beta_sq);

        // pressure and enthalpy
        Psw[i] = s[i] * beta_sq * Msw[i] / 3.0;
        Hb[i] = Eb[i] + Psw[i];

        // velocity components
        beta_th[i] = Ht[i] / Hb[i];
        beta_r[i] = (beta_th[i] <= beta[i]) ? std::sqrt(beta_sq - beta_th[i] * beta_th[i]) : 0.0;
        
        // forward shock
        beta_f[i] = 4.0 * beta[i] * (beta_gamma_sq[i] + 1) / (4.0 * beta_gamma_sq[i] + 3.0);   
    }
}

void PDE_Solver::solve_eigen() {
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

        eigenvalues[i] = std::max(std::max(alpha1, alpha2), alpha3) * CONST.c / R[i];
        alpha_R[i] = std::abs(beta_th[i]) / R[i];
    }
}

void PDE_Solver::solve_delta_t(const double& cfl) {
    double omega, omega_all;
    double delta_t_min, delta_t;
    
    // find minimum value
    delta_t_min = std::numeric_limits<double>::max();
    for (int i = 0; i < ntheta; ++i) {
        // maximum possible omega
        omega = beta[i] * CONST.c / R[i];

        // combined signal speed and maximum omega
        omega_all = eigenvalues[i] + 0.05 * omega;
        //omega_all = omega;

        // delta_t of a cell
        delta_t = cfl * (theta_edge[i + 1] - theta_edge[i]) / omega_all;

        // find minimum
        delta_t_min = std::min(delta_t_min, delta_t);
    }

    dt = delta_t_min;
}

void PDE_Solver::solve_slope() {
    // slope limit
    int index1, index2;
    double diff1, diff2;
    double slope1, slope2;

    for (int i = 0; i < ntheta; ++i) {
        index1 = std::max(i - 1, 0);
        index2 = std::min(i + 1, ntheta - 1);
        
        // Msw
        diff1 = Msw[i] - Msw[index1];
        diff2 = Msw[index2] - Msw[i];
        slope1 = (diff1 == 0.0) ? 0.0 : diff1 / (theta[i] - theta[index1]);
        slope2 = (diff2 == 0.0) ? 0.0 : diff2 / (theta[index2] - theta[i]);
        Msw_slope[i] = Minmod(slope1, slope2);

        // Mej
        diff1 = Mej[i] - Mej[index1];
        diff2 = Mej[index2] - Mej[i];
        slope1 = (diff1 == 0.0) ? 0.0 : diff1 / (theta[i] - theta[index1]);
        slope2 = (diff2 == 0.0) ? 0.0 : diff2 / (theta[index2] - theta[i]);
        Mej_slope[i] = Minmod(slope1, slope2);

        // beta_gamma_sq
        diff1 = beta_gamma_sq[i] - beta_gamma_sq[index1];
        diff2 = beta_gamma_sq[index2] - beta_gamma_sq[i];
        slope1 = (diff1 == 0.0) ? 0.0 : diff1 / (theta[i] - theta[index1]);
        slope2 = (diff2 == 0.0) ? 0.0 : diff2 / (theta[index2] - theta[i]);
        beta_gamma_sq_slope[i] = Minmod(slope1, slope2);

        // beta_th
        diff1 = beta_th[i] - beta_th[index1];
        diff2 = beta_th[index2] - beta_th[i];
        slope1 = (diff1 == 0.0) ? 0.0 : diff1 / (theta[i] - theta[index1]);
        slope2 = (diff2 == 0.0) ? 0.0 : diff2 / (theta[index2] - theta[i]);
        beta_th_slope[i] = Minmod(slope1, slope2);

        // R
        diff1 = R[i] - R[index1];
        diff2 = R[index2] - R[i];
        slope1 = (diff1 == 0.0) ? 0.0 : diff1 / (theta[i] - theta[index1]);
        slope2 = (diff2 == 0.0) ? 0.0 : diff2 / (theta[index2] - theta[i]);
        R_slope[i] = (slope1 + slope2) / 2.0;
        R_slope1[i] = slope1;
        R_slope2[i] = slope2;
    }
}

void PDE_Solver::solve_numerical_flux() {
    // reconstructed variable from the "left" and "right" cells
    double Msw_l, Mej_l, beta_gamma_sq_l, beta_th_l, R_l, Eb_l, Ht_l, Psw_l;
    double Msw_r, Mej_r, beta_gamma_sq_r, beta_th_r, R_r, Eb_r, Ht_r, Psw_r;

    // physical flux from the left and right cells
    Array1D F_l = Array(4);
    Array1D F_r = Array(4);

    // viscosity
    double alpha;

    for (int i = 1; i < ntheta; ++i) {
        // reconstruct primitive (left)
        Msw_l = Msw[i - 1] + Msw_slope[i - 1] * (theta_edge[i] - theta[i - 1]);
        Mej_l = Mej[i - 1] + Mej_slope[i - 1] * (theta_edge[i] - theta[i - 1]);
        beta_gamma_sq_l = beta_gamma_sq[i - 1] + beta_gamma_sq_slope[i - 1] * (theta_edge[i] - theta[i - 1]);
        beta_th_l = beta_th[i - 1] + beta_th_slope[i - 1] * (theta_edge[i] - theta[i - 1]);
        R_l = R[i - 1] + R_slope[i - 1] * (theta_edge[i] - theta[i - 1]);

        // reconstruct primitive (right)
        Msw_r = Msw[i] + Msw_slope[i] * (theta_edge[i] - theta[i]);
        Mej_r = Mej[i] + Mej_slope[i] * (theta_edge[i] - theta[i]);
        beta_gamma_sq_r = beta_gamma_sq[i] + beta_gamma_sq_slope[i] * (theta_edge[i] - theta[i]);
        beta_th_r = beta_th[i] + beta_th_slope[i] * (theta_edge[i] - theta[i]);
        R_r = R[i] + R_slope[i] * (theta_edge[i] - theta[i]);

        // reconstruct conserved (left)
        Eb_l = s[i] * (1.0 + beta_gamma_sq_l * beta_gamma_sq_l / (beta_gamma_sq_l + 1) / (beta_gamma_sq_l + 1) / 3.0) * (beta_gamma_sq_l + 1) * Msw_l
             + (1.0 - s[i]) * std::sqrt(beta_gamma_sq_l + 1) * Msw_l
             + std::sqrt(beta_gamma_sq_l + 1) * Mej_l;
        Psw_l = s[i] * beta_gamma_sq_l / (beta_gamma_sq_l + 1) * Msw_l / 3.0;
        Ht_l = (Eb_l + Psw_l) * beta_th_l;

        // reconstruct conserved (right)
        Eb_r = s[i] * (1.0 + beta_gamma_sq_r * beta_gamma_sq_r / (beta_gamma_sq_r + 1) / (beta_gamma_sq_r + 1) / 3.0) * (beta_gamma_sq_r + 1) * Msw_r 
             + (1.0 - s[i]) * std::sqrt(beta_gamma_sq_r + 1) * Msw_r
             + std::sqrt(beta_gamma_sq_r + 1) * Mej_r;
        Psw_r = s[i] * beta_gamma_sq_r / (beta_gamma_sq_r + 1) * Msw_r / 3.0;
        Ht_r = (Eb_r + Psw_r) * beta_th_r;

        // physical flux (left)
        F_l[0] = Ht_l / R_l * CONST.c;
        F_l[1] = (Ht_l * beta_th_l + Psw_l) / R_l * CONST.c;
        F_l[2] = Msw_l * beta_th_l / R_l * CONST.c;
        F_l[3] = Mej_l * beta_th_l / R_l * CONST.c;

        // physical flux (right)
        F_r[0] = Ht_r / R_r * CONST.c;
        F_r[1] = (Ht_r * beta_th_r + Psw_r) / R_r * CONST.c;
        F_r[2] = Msw_r * beta_th_r / R_r * CONST.c;
        F_r[3] = Mej_r * beta_th_r / R_r * CONST.c;

        // viscosity
        int index_l = std::max(i - 2, 0);
        int index_r = std::min(i + 1, ntheta - 1);
        alpha = *std::max_element(eigenvalues.begin() + index_l, eigenvalues.begin() + index_r);

        // numerical flux
        numerical_flux[0][i] = 0.5 * (F_l[0] + F_r[0] - alpha * (Eb_r - Eb_l));
        numerical_flux[1][i] = 0.5 * (F_l[1] + F_r[1] - alpha * (Ht_r - Ht_l));
        numerical_flux[2][i] = 0.5 * (F_l[2] + F_r[2] - alpha * (Msw_r - Msw_l));
        numerical_flux[3][i] = 0.5 * (F_l[3] + F_r[3] - alpha * (Mej_r - Mej_l));
        for (int j = 0 ; j < 4; ++j) {
            numerical_flux[j][i] *= std::sin(theta_edge[i]);
        }
    }
}

void PDE_Solver::solve_dy_dt(const Array2D& y, Array2D& dy_dt) {
    auto t1 = std::chrono::high_resolution_clock::now();
    // update conserved
    update_conserved(y);

    // solve pde_state
    solve_state();
    auto t2 = std::chrono::high_resolution_clock::now();

    // solve slope
    solve_slope();

    // solve eigenvalues
    solve_eigen();

    // solve numerical flux
    solve_numerical_flux();

    // solve dy_dt
    for (int i = 0; i < ntheta; ++i) {
        int index_l = std::max(i - 1, 0);
        int index_r = std::min(i + 1, ntheta - 1);
        double alpha = std::max(std::max(alpha_R[index_l], alpha_R[i]), alpha_R[index_r]);

        // dR_dt
        dy_dt[4][i] = beta_f[i] - R_slope[i] * beta_th[i] / R[i] + 0.5 * alpha * (R_slope2[i] - R_slope1[i]);
        dy_dt[4][i] *= CONST.c;

        // conserved variables
        double rho = density0(nwind, nism, R[i]) * CONST.mp;
        dy_dt[0][i] = (numerical_flux[0][i] - numerical_flux[0][i + 1]) / voltheta[i]
                    + dy_dt[4][i] * rho * R[i] * R[i];
        dy_dt[1][i] = (numerical_flux[1][i] - numerical_flux[1][i + 1]) / voltheta[i]
                    + (std::cos(theta[i]) / std::sin(theta[i]) * Psw[i] - Ht[i] * beta_r[i]) * CONST.c / R[i];
        dy_dt[2][i] = (numerical_flux[2][i] - numerical_flux[2][i + 1]) / voltheta[i]
                    + dy_dt[4][i] * rho * R[i] * R[i];
        dy_dt[3][i] = (numerical_flux[3][i] - numerical_flux[3][i + 1]) / voltheta[i]
                    + 0.0;
    }
    auto t3 = std::chrono::high_resolution_clock::now();
    time1 += std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() * 1e-9;
    time2 += std::chrono::duration_cast<std::chrono::nanoseconds>(t3 - t2).count() * 1e-9;
}

void PDE_Solver::solve_dy_dt_no_spreading(const Array2D& y, Array2D& dy_dt) {
    // update conserved
    update_conserved(y);

    // solve pde_state
    solve_state();

    // solve dy_dt
    for (int i = 0; i < ntheta; ++i) {
        // dR_dt
        dy_dt[4][i] = beta_f[i] * CONST.c;

        // conserved variables
        double rho = density0(nwind, nism, R[i]) * CONST.mp;
        dy_dt[0][i] = dy_dt[4][i] * rho * R[i] * R[i];
        dy_dt[1][i] = 0.0;
        dy_dt[2][i] = dy_dt[4][i] * rho * R[i] * R[i];
        dy_dt[3][i] = 0.0;
    }
}

// ******************** 1D Jet solver ******************** //

class Jet1D {
public:
    // theta
    int ntheta;                 // cell number
    Array1D theta;              // (ntheta) cell center
    Array1D theta_edge;         // (ntheta + 1) cell edge
    Array1D voltheta;           // (ntheta) cell volume
    
    // configuration
    bool spread = true;      // spreading
    double tmin, tmax;          // PDE maximum time
    double xtol = 0; 
    double rtol = 1e-6;         // root solver accuracy
    double cfl;                 // cfl number
    double nwind, nism;             // define ambient density
    Array2D y0;

    // PDE intermediate
    Array1D ts;                 // (nt)
    Array3D ys;                 // (nt, 5, ntheta) y: Eb, HtR, Msw, Mej, R
    Array3D y_pris;             // (nt, 5, ntheta) y: Msw, Mej, beta_gamma_sq, beta_th, R
    PDE_Solver pde_solver;

    // PDE
    Interpolator interpolator;  // interpolator based on t_interp, y_interp
    EATS_Solver eats_solver;

    // max Tobs for EATS
    double Tobs_z_max;          // Maximum observing time
    
    // configuration
    void config(const Array1D& theta_edge, const Array2D& yini, const bool& spread, const double& tmin, const double& tmax, const double& rtol, const double& cfl, const double& nwind, const double& nism);

    // solve PDE
    void solve();
    void solve_no_spreading();

private:

};

// 1d Jet setup //
// y0: Eb, HtR, Msw, Mej, R
void Jet1D::config(const Array1D& theta_edge, const Array2D& y0, const bool& spread, const double& tmin, const double& tmax, const double& rtol, const double& cfl, const double& nwind, const double& nism) {
    // configuration
    this->theta_edge = theta_edge;
    this->spread = spread;
    this->tmin = tmin;
    this->tmax = tmax;
    this->rtol = rtol;
    this->cfl = cfl * 2;
    this->nwind = nwind;
    this->nism = nism;
    this->y0 = y0;

    // setup theta
    ntheta = theta_edge.size() - 1;
    theta = Array(ntheta);
    voltheta = Array(ntheta);
    for (int i = 0; i < ntheta; ++i) {
        theta[i] = (theta_edge[i] + theta_edge[i + 1]) / 2.0;
        voltheta[i] = std::cos(theta_edge[i]) - std::cos(theta_edge[i + 1]);
    }

    // PDE data initialization
    ts = Array(1);
    ts[0] = tmin;
    ys = Array(1, 5, ntheta);
    ys[0] = y0;

    // pde solver
    pde_solver.config(nwind, nism, ntheta, theta, theta_edge, voltheta, rtol);
}

void Jet1D::solve() {
    // initial conditions
    double t = ts[0];
    Array2D y = ys[0];
    Array2D y_pri = Array(6, ntheta);
    
    // dy_dt function for SSPRK 
    auto f = [&](Array2D& y, Array2D& dy_dt, Array2D& y_pri, double& delta_t) {
        // solve dy_dt
        pde_solver.solve_dy_dt(y, dy_dt);

        // get primitive variables
        y_pri[0] = pde_solver.Msw;
        y_pri[1] = pde_solver.Mej;
        y_pri[2] = pde_solver.beta_gamma_sq;
        y_pri[3] = pde_solver.beta_th;
        y_pri[4] = pde_solver.R;
        y_pri[5] = pde_solver.s;

        // get delta_t
        pde_solver.solve_delta_t(cfl);
        delta_t = pde_solver.dt;
    };

    // update
    double delta_t;
    int n = 0;
    while (t <= tmax) {
        // update t, y
        try {
            //SSP_RK2(f, y, y_pri, delta_t); // 2 stage SSP-RK2 with cfl < 1 (cfl_eff < 0.5)
            SSP_RK3_4(f, y, y_pri, delta_t); // 4 stage SSP-RK3 with cfl < 2 (cfl_eff < 0.5)
        }
        catch (const std::exception& e) {
            std::string text = "\nFailed to perform SSP_RK2 time integration. \n↓";
            throw std::runtime_error(text + e.what());
        }

        // update time
        t += delta_t;

        //std::cout << t << "\n";

        // save result
        ts.push_back(t);
        y_pris.push_back(y_pri);
        //ys.push_back(y);
    }

    // solve last primitive
    Array2D dy_dt_tmp = Array(5, ntheta);
    f(y, dy_dt_tmp, y_pri, delta_t);
    y_pris.push_back(y_pri);

    //std::cout << pde_solver.time1 << " " << pde_solver.time2 << "\n";

    // interpolate primitive variables
    interpolator = Interpolator(ts, theta, y_pris, nwind, nism);
    eats_solver.config(ts, theta, y_pris, nwind, nism);

    // maximum time
    Tobs_z_max = ts.back() - *std::max_element(y_pris.back()[4].begin(), y_pris.back()[4].end()) / CONST.c;
}

void Jet1D::solve_no_spreading() {
    double t = ts[0];
    Array2D y = ys[0];
    Array2D y_pri = Array(6, ntheta);
    double delta_t = 1.0;
    
    // function
    auto f = [&](Array2D& y, Array2D& dy_dt, Array2D& y_pri) {
        // solve dy_dt
        pde_solver.solve_dy_dt_no_spreading(y, dy_dt);

        // get primitive variables
        y_pri[0] = pde_solver.Msw;
        y_pri[1] = pde_solver.Mej;
        y_pri[2] = pde_solver.beta_gamma_sq;
        y_pri[3] = pde_solver.beta_th;
        y_pri[4] = pde_solver.R;
        y_pri[5] = pde_solver.s;
    };

    // update
    while (t <= tmax) {
        // update t, y
        bool succeeded;
        double dt = delta_t;
        try {
            RK45(f, y, y_pri, dt, 1e-6, succeeded);
        }
        catch (const std::exception& e) {
            std::string text = "\nFailed to perform RK45 time integration. \n↓";
            throw std::runtime_error(text + e.what());
        }

        if (succeeded) {
            t += delta_t;
            delta_t = dt;

            // save primitive variable
            y_pris.push_back(y_pri);
            
            // save result
            ts.push_back(t);
            //ys.push_back(y);
        }
        else {
            delta_t = dt;
        }
    }

    // solve last primitive
    Array2D dy_dt_tmp = Array(5, ntheta);
    f(y, dy_dt_tmp, y_pri);
    y_pris.push_back(y_pri);

    // interpolate primitive variables
    interpolator = Interpolator(ts, theta, y_pris, nwind, nism);
    eats_solver.config(ts, theta, y_pris, nwind, nism);

    // maximum time
    Tobs_z_max = ts.back() - *std::max_element(y_pris.back()[4].begin(), y_pris.back()[4].end()) / CONST.c;
}

// ---------- other functions ---------- //



#endif