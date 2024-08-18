#ifndef INTEGRAL
#define INTEGRAL

#include "../environment.h"

// representation of a 5-point interval for integration
struct Interval {
    double x1;                   // start point x value of the interval
    std::array<double, 5> y;     // array of function values
    double h;                    // interval length
    double lorder;               // low order integral estimation
    double horder;               // high order integral estimation

    // compute the low order and high order estimations
    void integrate() {
        // Simpson's rule
        lorder = (y[0] + 4.0 * y[2] + y[4]) * h / 6.0;

        // Composite Simpson's rule
        //horder = (y[0] + 4.0 * y[1] + 2.0 * y[2] + 4.0 * y[3] + y[4]) * h / 12.0;

        // Boole's rule (1 step Romberg extrapolation)
        horder = (7.0 * y[0] + 32.0 * y[1] + 12.0 * y[2] + 32.0 * y[3] + 7.0 * y[4]) * h / 90.0;
    }
};

// Adaptive 1D integral with initial x samples
template<class F>
double Adaptive_1D(F& f, const Array1D xini, const double xtol, const double rtol, const int max_iter = 50) {
    // initialize interval
    Interval new_itv;
    std::vector<Interval> intervals = {};
    double integral_tot = 0.0;
    int xsize = xini.size();
    for (int i = 0; i < xsize - 1; ++i) {
        // assign values to a new interval
        new_itv.x1 = xini[i];
        new_itv.h = xini[i + 1] - xini[i];
        new_itv.y = {f(xini[i]), f(xini[i] + new_itv.h / 4.0), f(xini[i] + 2.0 * new_itv.h / 4.0), f(xini[i] + 3.0 * new_itv.h / 4.0), f(xini[i + 1])};
        new_itv.integrate();

        // push_back to interval vectors
        intervals.push_back(new_itv);
        integral_tot += new_itv.horder;
    }
    double h_tot = xini.back() - xini.front();

    // control variables
    bool need_refine;          // if a single interval needs refinement
    bool any_refined;          // if any of the interval is refined
    double integral_new;       // a variable to update total integral

    // adptive refinement
    for (int i = 0; i < max_iter; ++i) {
        // initialize some variables
        any_refined = false;
        int intervals_size = intervals.size();
        integral_new = integral_tot;

        // loop over all intervals and make necessary refinement
        for (int i = 0; i < intervals_size; ++i) {
            // alias
            Interval& itv = intervals[i];

            // test if refinement is needed
            need_refine = (
                // local relative error
                std::abs(itv.lorder - itv.horder) > xtol * itv.h / h_tot + rtol * std::abs(itv.horder)
                &&
                // globally relative local error (Gander & Gautschi 2001 ?)
                std::abs(itv.lorder - itv.horder) * intervals_size > xtol + rtol * std::abs(integral_tot)
            );
            
            // refine
            if (need_refine) {
                // subtract the contribution from this interval
                integral_new -= itv.horder;

                // refine right half part by a new interval
                new_itv.x1 = itv.x1 + itv.h / 2.0;
                new_itv.h = itv.h / 2.0;
                new_itv.y = {itv.y[2], f(new_itv.x1 + new_itv.h / 4.0), itv.y[3], f(new_itv.x1 + 3.0 * new_itv.h / 4.0), itv.y[4]};
                new_itv.integrate();

                // refine left half part by changing original interval
                itv.h = itv.h / 2.0;
                itv.y = {itv.y[0], f(itv.x1 + itv.h / 4.0), itv.y[1], f(itv.x1 + 3.0 * itv.h / 4.0), itv.y[2]};
                itv.integrate();

                // update total integral
                integral_new += (new_itv.horder + itv.horder);

                // mark
                any_refined = true;

                // add new interval to the last
                intervals.push_back(new_itv);
            }
            else {
                continue;
            }
        }

        // update total integral
        integral_tot = integral_new;

        // if any of the intervals are refined?
        if (any_refined) {
            // do nothing, keep refinement
            continue;
        }
        else {
            // if no refinement, return the value
            return integral_tot;
        }
    }

    // if integral does not converge in max_iter loops, throw error
    throw std::runtime_error("Integral: Not converged after " + std::to_string(max_iter) + " iterations!");
}

// Adaptive 2D integral with initial x and y samples
template<class F>
double Adaptive_2D(F& f, const Array1D& xini, const Array1D& yini, const double xtol, const double rtol, const int max_iter = 50) {
    auto g = [&](const double y) {
        auto h = [&](const double x) {
            return f(x, y);
        };
        double result = Adaptive_1D(h, xini, xtol, rtol, max_iter);
        return result;
    };
    return Adaptive_1D(g, yini, xtol, rtol, max_iter);
}

#endif