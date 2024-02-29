#ifndef REFINEMENT
#define REFINEMENT

#include "Definations.hpp"
#include "Interpolation.hpp"
#include "Integral.hpp"

double diff_high_order(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4) {
    double x = (x2 + x3) / 2.0;
    double p1 = ((x - x2) * (x - x3) + (x - x2) * (x - x4) + (x - x3) * (x - x4)) / (x1 - x2) / (x1 - x3) / (x1 - x4) * y1;
    double p2 = ((x - x1) * (x - x3) + (x - x1) * (x - x4) + (x - x3) * (x - x4)) / (x2 - x1) / (x2 - x3) / (x2 - x4) * y2;
    double p3 = ((x - x1) * (x - x2) + (x - x1) * (x - x4) + (x - x2) * (x - x4)) / (x3 - x1) / (x3 - x2) / (x3 - x4) * y3;
    double p4 = ((x - x1) * (x - x2) + (x - x1) * (x - x3) + (x - x2) * (x - x3)) / (x4 - x1) / (x4 - x2) / (x4 - x3) * y4;
    return p1 + p2 + p3 + p4;
}

double Lagrange_4pts(double x1, double x2, double x3, double x4, double y1, double y2, double y3, double y4) {
    double x = (x2 + x3) / 2.0;
    double l1 = (x - x2) / (x1 - x2) * (x - x3) / (x1 - x3) * (x - x4) / (x1 - x4);
    double l2 = (x - x1) / (x2 - x1) * (x - x3) / (x2 - x3) * (x - x4) / (x2 - x4);
    double l3 = (x - x1) / (x3 - x1) * (x - x2) / (x3 - x2) * (x - x4) / (x3 - x4);
    double l4 = (x - x1) / (x4 - x1) * (x - x2) / (x4 - x2) * (x - x3) / (x4 - x3);
    return y1 * l1 + y2 * l2 + y3 * l3 + y4 * l4;
}

Array1D Smoothen(const Array1D& x_data, const Array1D& y_data, const double& threshold) {
    int n = x_data.size();
    Array1D ynew = Array(n);

    // extrapolate y
    Array1D y = Array(n + 6);
    y[0] = y_data[3];
    y[1] = y_data[2];
    y[2] = y_data[1];
    for ( int i = 0; i < n; ++i) {
        y[i + 3] = y_data[i];
    }
    y[n + 3] = y_data[n - 2];
    y[n + 4] = y_data[n - 3];
    y[n + 5] = y_data[n - 4];

    for (int i = 0; i < n; ++i) {
        if (ynew[i] <= threshold) {
            ynew[i] = (
                      0.0625 * y[i]
                    + 0.1250 * y[i + 1]
                    + 0.1875 * y[i + 2]
                    + 0.2500 * y[i + 3]
                    + 0.1875 * y[i + 4]
                    + 0.1250 * y[i + 5]
                    + 0.0625 * y[i + 6]);
        }
    }

    return ynew;
}

Array1D Refine1(const Array1D& x, const Array1D& y, const double& xtol) {
    int n = x.size();
    Array1D xx = x;
    Array1D yy = y;
    Array1D x_new = x;

    double ymax = *std::max_element(yy.begin(), yy.end());

    // extrapolate
    xx.insert(xx.begin(), - xx[1]);
    xx.insert(xx.end(), xx[xx.size() - 2]);
    yy.insert(yy.begin(), yy[1]);
    yy.insert(yy.end(), y[y.size() - 2]);
    
    for (int i = 0; i < n - 3; ++i) {
        // estimate value
        double y_l = (yy[i + 1] + yy[i + 2]) / 2.0;
        double y_h = Lagrange_4pts(xx[i], xx[i + 1], xx[i + 2], xx[i + 3], yy[i], yy[i + 1], yy[i + 2], yy[i + 3]);

        // estimate derivative
        double slope_l = (yy[i + 2] - yy[i + 1]) / (xx[i + 2] - xx[i + 1]);
        double slope_h = diff_high_order(xx[i], xx[i + 1], xx[i + 2], xx[i + 3], yy[i], yy[i + 1], yy[i + 2], yy[i + 3]);

        // refine
        if (ymax - std::max(y_l, y_h) < 6.0 && (std::abs(y_l - y_h) > xtol * y_h || std::abs(slope_l - slope_h) > xtol)) {
            x_new.push_back((xx[i + 1] + xx[i + 2]) / 2.0);
        }
    }

    std::sort(x_new.begin(), x_new.end());

    return x_new;
}

Array1D Refine2(const Array1D& x, const Array1D& y1, const Array1D& y2, const double& xtol) {
    int n = x.size();
    Array1D xx = x;
    Array1D yy1 = y1;
    Array1D yy2 = y2;
    Array1D x_new = x;

    double y1max = *std::max_element(yy1.begin(), yy1.end());
    double y2max = *std::max_element(yy2.begin(), yy2.end());

    // extrapolate
    xx.insert(xx.begin(), - xx[1]);
    xx.insert(xx.end(), xx[xx.size() - 2]);
    yy1.insert(yy1.begin(), yy1[1]);
    yy1.insert(yy1.end(), yy1[yy1.size() - 2]);
    yy2.insert(yy2.begin(), yy2[1]);
    yy2.insert(yy2.end(), yy2[yy2.size() - 2]);
    
    for (int i = 0; i < n - 1; ++i) {
        // estimate value
        double y1_l = (yy1[i + 1] + yy1[i + 2]) / 2.0;
        double y1_h = Lagrange_4pts(xx[i], xx[i + 1], xx[i + 2], xx[i + 3], yy1[i], yy1[i + 1], yy1[i + 2], yy1[i + 3]);
        double y2_l = (yy2[i + 1] + yy2[i + 2]) / 2.0;
        double y2_h = Lagrange_4pts(xx[i], xx[i + 1], xx[i + 2], xx[i + 3], yy2[i], yy2[i + 1], yy2[i + 2], yy2[i + 3]);

        // estimate derivative
        double slope1_l = (yy1[i + 2] - yy1[i + 1]) / (xx[i + 2] - xx[i + 1]);
        double slope1_h = diff_high_order(xx[i], xx[i + 1], xx[i + 2], xx[i + 3], yy1[i], yy1[i + 1], yy1[i + 2], yy1[i + 3]);
        double slope2_l = (yy2[i + 2] - yy2[i + 1]) / (xx[i + 2] - xx[i + 1]);
        double slope2_h = diff_high_order(xx[i], xx[i + 1], xx[i + 2], xx[i + 3], yy2[i], yy2[i + 1], yy2[i + 2], yy2[i + 3]);

        // refine
        if ((y1max - std::max(y1_l, y1_h) < 6.0) && (std::abs(y1_l - y1_h) > xtol * y1_h || std::abs(y2_l - y2_h) > xtol * y2_h || std::abs(slope1_l - slope1_h) > xtol || std::abs(slope2_l - slope2_h) > xtol)) {
            x_new.push_back((xx[i + 1] + xx[i + 2]) / 2.0);
        }
    }

    std::sort(x_new.begin(), x_new.end());

    return x_new;
}

#endif