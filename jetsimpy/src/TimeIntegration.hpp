#ifndef TIME_INTEGRATION
#define TIME_INTEGRATION

#include "Definations.hpp"

template<typename F1>
void SSP_RK2(F1& f, Array2D& y, Array2D& y_pri, double& dt) {
    // size
    int ntheta = y[0].size();
    Array2D dy_dt = Array(5, ntheta);
    Array2D ymid = Array(5, ntheta);
    double dt_tmp;
    Array2D y_pri_tmp = y_pri;

    // step 1
    f(y, dy_dt, y_pri, dt);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            ymid[i][j] = y[i][j] + dt * dy_dt[i][j];
        }
    }

    // step 2
    f(ymid, dy_dt, y_pri_tmp, dt_tmp);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            y[i][j] = 0.5 * y[i][j] + 0.5 * ymid[i][j] + 0.5 * dt * dy_dt[i][j];
        }
    }
}

template<typename F1>
void SSP_RK2_10(F1& f, Array2D& y, Array2D& y_pri, double& dt) {
    int ntheta = y[0].size();
    double dt_tmp;
    Array2D y_pri_tmp = y_pri;

    int nstage = 10;

    Array3D ys = Array(nstage, 5, ntheta);
    Array3D dy_dts = Array(nstage, 5, ntheta);

    // step 0
    ys[0] = y;

    // step 1
    f(ys[0], dy_dts[0], y_pri, dt);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            ys[1][i][j] = ys[0][i][j] + dt * dy_dts[0][i][j] / (nstage - 1);
        }
    }

    // step 2 - (nstage - 1)
    for (int i = 1; i < (nstage - 1); ++i) {
        f(ys[i], dy_dts[i], y_pri_tmp, dt_tmp);
        ys[i + 1] = y;
        for (int j = 0; j < 5; ++j) {
            for (int k = 0; k < ntheta; ++k) {
                for (int l = 0; l < i + 1; ++l) {
                    ys[i + 1][j][k] += dt * dy_dts[l][j][k] / (nstage - 1);
                }
            }
        }
    }

    // step nstage
    f(ys[nstage - 1], dy_dts[nstage - 1], y_pri_tmp, dt_tmp);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            for (int k = 0; k < nstage; ++k) {
                y[i][j] += dt * dy_dts[k][i][j] / nstage;
            }
        }
    }
}

template<typename F1>
void SSP_RK3_4(F1& f, Array2D& y, Array2D& y_pri, double& dt) {
    int ntheta = y[0].size();
    double dt_tmp;
    Array2D y_pri_tmp = y_pri;
    Array3D ys = Array(4, 5, ntheta);
    Array3D dy_dts = Array(4, 5, ntheta);

    // step 0
    ys[0] = y;

    // step 1
    f(ys[0], dy_dts[0], y_pri, dt);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            ys[1][i][j] = y[i][j] + dt * dy_dts[0][i][j] / 2.0;
        }
    }

    // step 2
    f(ys[1], dy_dts[1], y_pri_tmp, dt_tmp);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            ys[2][i][j] = y[i][j] + dt * dy_dts[0][i][j] / 2.0 + dt * dy_dts[1][i][j] / 2.0;
        }
    }

    // step 3
    f(ys[2], dy_dts[2], y_pri_tmp, dt_tmp);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            ys[3][i][j] = y[i][j] + dt * dy_dts[0][i][j] / 6.0 + dt * dy_dts[1][i][j] / 6.0 + dt * dy_dts[2][i][j] / 6.0;
        }
    }

    // step 4
    f(ys[3], dy_dts[3], y_pri_tmp, dt_tmp);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            y[i][j] += dt * dy_dts[0][i][j] / 6.0 + dt * dy_dts[1][i][j] / 6.0 + dt * dy_dts[2][i][j] / 6.0 + dt * dy_dts[3][i][j] / 2.0;
        }
    }
}

template<typename F1>
void RK45(F1& f, Array2D& y, Array2D& y_pri, double& dt, const double& rtol, bool& succeeded) {
    int ntheta = y[0].size();
    Array2D k1, k2, k3, k4, k5, k6;
    k1 = k2 = k3 = k4 = k5 = k6 = Array(5, ntheta);
    Array2D yy = Array(5, ntheta);
    Array2D dy_dt = Array(5, ntheta);
    Array2D ynew = Array(5, ntheta);
    Array2D y_pri_tmp = y_pri;

    // step 1
    f(y, dy_dt, y_pri);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            k1[i][j] = dt * dy_dt[i][j];
            yy[i][j] = y[i][j] + k1[i][j] * 2.0 / 9.0;
        }
    }

    // step 2
    f(yy, dy_dt, y_pri_tmp);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            k2[i][j] = dt * dy_dt[i][j];
            yy[i][j] = y[i][j] + k1[i][j] / 12.0 + k2[i][j] / 4.0;
        }
    }

    // step 3
    f(yy, dy_dt, y_pri_tmp);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            k3[i][j] = dt * dy_dt[i][j];
            yy[i][j] = y[i][j] + k1[i][j] * 69.0 / 128.0 - k2[i][j] * 243.0 / 128.0 + k3[i][j] * 135.0 / 64.0;
        }
    }

    // step 4
    f(yy, dy_dt, y_pri_tmp);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            k4[i][j] = dt * dy_dt[i][j];
            yy[i][j] = y[i][j] - k1[i][j] * 17.0 / 12.0 + k2[i][j] * 27.0 / 4.0 - k3[i][j] * 27.0 / 5.0 + k4[i][j] * 16.0 / 15.0;
        }
    }

    // step 5
    f(yy, dy_dt, y_pri_tmp);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            k5[i][j] = dt * dy_dt[i][j];
            yy[i][j] = y[i][j] + k1[i][j] * 65.0 / 432.0 - k2[i][j] * 5.0 / 16.0 + k3[i][j] * 13.0 / 16.0 + k4[i][j] * 4.0 / 27.0 + k5[i][j] * 5.0 / 144.0;
        }
    }

    // step 6
    f(yy, dy_dt, y_pri_tmp);
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            k6[i][j] = dt * dy_dt[i][j];
        }
    }

    // error estimate
    double error;
    double rerror = 0.0;
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j < ntheta; ++j) {
            ynew[i][j] = y[i][j] + k1[i][j] * 47.0 / 450.0 + k3[i][j] * 12.0 / 25.0 + k4[i][j] * 32.0 / 225.0 + k5[i][j] / 30.0 + k6[i][j] * 6.0 / 25.0;
            error = std::abs(k1[i][j] / 150.0 - k3[i][j] * 3.0 / 100.0 + k4[i][j] * 16.0 / 75.0 + k5[i][j] / 20.0 - k6[i][j] * 6.0 / 25.0);
            rerror = std::max(rerror, error / std::abs(ynew[i][j]));
        }
    }

    if (rerror < rtol) {
        // update y
        for (int i = 0; i < 5; ++i) {
            for (int j = 0; j < ntheta; ++j) {
                y[i][j] = ynew[i][j];
            }
        }
        succeeded = true;
    }
    else {
        succeeded = false;
    }

    // update dt
    double boost_factor = 0.9 * std::pow(rtol / rerror, 0.2);
    boost_factor = std::min(1.5, boost_factor);
    dt *= boost_factor;
}

#endif