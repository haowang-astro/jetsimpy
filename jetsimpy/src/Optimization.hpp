#ifndef OPTIMIZATION
#define OPTIMIZATION

#include "Definations.hpp"

// this algorithms is only valid to the setup of Flux calculation
// this algorithm does minimization
template<typename F>
double golden_section(F& f, const double& xx1, const double& xx2, const double& atol, const int& max_iter = 100) {
    double r = 0.618;

    double x1 = xx1;
    double x2 = xx2;
    double xm1 = (x1 + r * x2) / (1.0 + r);
    double f1 = f(x1);
    double fm1 = f(xm1);
    double f2 = f(x2);

    // braket the range
    int n = 0;
    while (n < 10) {
        if (fm1 > f1 && fm1 <= f2) {
            x2 = xm1;
            f2 = fm1;

            xm1 = (x1 + r * x2) / (1.0 + r);
            fm1 = f(xm1);
            n += 1;
        }
        else if (fm1 <= f1 && fm1 > f2) {
            x1 = xm1;
            f1 = fm1;

            xm1 = (x1 + r * x2) / (1.0 + r);
            fm1 = f(xm1);
            n += 1;
        }
        else if (fm1 <= f1 && fm1 <= f2) {
            break;
        }
        else {
            if (f1 < f2) {
                x2 = xm1;
                f2 = fm1;
                
                xm1 = (x1 + r * x2) / (1.0 + r);
                fm1 = f(xm1);
                n += 1;
            }
            else {
                x1 = xm1;
                f1 = fm1;
                
                xm1 = (x1 + r * x2) / (1.0 + r);
                fm1 = f(xm1);
                n += 1;
            }
        }
    }

    // if the braketing failed
    if (n >= 10) {
        if (fm1 > f1) {
            return x1;
        }
        else {
            return x2;
        }
    }

    // another side
    double xm2 = (r * x1 + x2) / (1.0 + r);
    double fm2 = f(xm2);

    // loop
    n = 0;
    while(std::fabs(xm1 - xm2) > atol && n < max_iter) {
        if (fm1 < fm2) {
            x2 = xm2;
            f2 = fm2;

            xm2 = xm1;
            fm2 = fm1;

            xm1 = (x1 + r * x2) / (1.0 + r);
            fm1 = f(xm1);
        }
        else {
            x1 = xm1;
            f1 = fm1;

            xm1 = xm2;
            fm1 = fm2;

            xm2 = (r * x1 + x2) / (1.0 + r);
            fm2 = f(xm2);
        }
        n += 1;
    }

    if (n >= max_iter) {
        std::string text = "Optimizer doesn't converge after " + std::to_string(max_iter) + " interations!\n";
        throw std::runtime_error(text);
    }

    return (xm1 + xm2) / 2.0;
}

// minimization givn an initial sampling
template<typename F>
double minimization(F& f, const double& x1, const double& x2, const int& n_opt, const double& atol, const std::string& type_, const int& max_iter = 100) {
    Array1D x_ini = Array(n_opt);
    Array1D f_ini = Array(n_opt);

    // evaluate initial sampling
    if (type_ == "linear") { //linear spacing
        for (int i = 0; i < n_opt; ++i) {
            x_ini[i] = x1 + i * (x2 - x1) / (n_opt - 1);
            f_ini[i] = f(x_ini[i]);
        }
    }
    else if (type_ == "log") { // spacing by arcsinh
        double asinh_x1 = std::asinh(x1);
        double asinh_x2 = std::asinh(x2);
        for (int i = 0; i < n_opt; ++i) {
            x_ini[i] = asinh_x1 + i * (asinh_x2 - asinh_x1) / (n_opt - 1);
            x_ini[i] = std::sinh(x_ini[i]);
            f_ini[i] = f(x_ini[i]);
        }
    }
    else {
        std::string text = "No such optimization method: '" + type_ + "'\n";
        throw std::runtime_error(text);
    }

    // find minimum place
    double x_peak;
    int argmin = std::distance(f_ini.begin(), std::min_element(f_ini.begin(), f_ini.end()));
    if (argmin == 0) {
        x_peak = golden_section(f, x_ini[0], x_ini[1], atol, max_iter);
    }
    else if (argmin == n_opt - 1) {
        x_peak = golden_section(f, x_ini[n_opt - 2], x_ini[n_opt - 1], atol, max_iter);
    }
    else {
        x_peak = golden_section(f, x_ini[argmin - 1], x_ini[argmin + 1], atol, max_iter);
    }

    return x_peak;
}

template<typename F>
double minimization(F& f, const Array1D& x_ini, const double& atol, const int& max_iter = 100) {
    int n = x_ini.size();
    Array1D f_ini = Array(n);

    // evaluate initial sampling
    for (int i = 0; i < n; ++i) {
        f_ini[i] = f(x_ini[i]);
    }

    // find minimum place
    double x_peak;
    int argmin = std::distance(f_ini.begin(), std::min_element(f_ini.begin(), f_ini.end()));
    if (argmin == 0) {
        x_peak = golden_section(f, x_ini[0], x_ini[1], atol, max_iter);
    }
    else if (argmin == n - 1) {
        x_peak = golden_section(f, x_ini[n - 2], x_ini[n - 1], atol, max_iter);
    }
    else {
        x_peak = golden_section(f, x_ini[argmin - 1], x_ini[argmin + 1], atol, max_iter);
    }

    return x_peak;
}

#endif