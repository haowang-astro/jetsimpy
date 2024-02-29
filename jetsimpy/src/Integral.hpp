#ifndef INTEGRAL
#define INTEGRAL

#include "Definations.hpp"

struct Point {
    double x;
    Array1D y;
};

struct Interval {
    std::list<Point>::iterator pt1, pt2, pt3, pt4, pt5;
    double h;
    Array1D simps_1, simps_2;
    Array1D simps_sum, simps_tot;
};

template<typename F>
Array1D Adaptive_1D(F& f, const Array1D& xini, const double& xtol, const double& rtol, const int max_iter = 20) {
    // point and interval list
    std::list<Point> points;
    std::list<Interval> intervals;

    // initialize points
    double h_tot = xini.back() - xini.front();
    int xini_size = xini.size();
    Point pt;
    for (int i = 0; i < xini_size - 1; ++i) {
        pt.x = xini[i];
        pt.y = f(pt.x);
        points.push_back(pt);

        pt.x = 0.75 * xini[i] + 0.25 * xini[i + 1];
        pt.y = f(pt.x);
        points.push_back(pt);
        
        pt.x = 0.5 * xini[i] + 0.5 * xini[i + 1];
        pt.y = f(pt.x);
        points.push_back(pt);

        pt.x = 0.25 * xini[i] + 0.75 * xini[i + 1];
        pt.y = f(pt.x);
        points.push_back(pt);
    }
    pt.x = xini[xini_size - 1];
    pt.y = f(pt.x);
    points.push_back(pt);

    // initialize lists
    int n = pt.y.size();
    Interval itv;
    itv.simps_1 = itv.simps_2 = itv.simps_sum = itv.simps_tot = Array(n);
    Array1D integral = Array(n);
    for (auto it = ++points.begin(); it != points.end(); ++++++++it) {
        itv.pt1 = std::prev(it);
        itv.pt2 = it;
        itv.pt3 = std::next(it);
        itv.pt4 = std::next(itv.pt3);
        itv.pt5 = std::next(itv.pt4);
        itv.h = itv.pt5->x - itv.pt1->x;
        for (int i = 0; i < n; ++i) {
            itv.simps_1[i] = (itv.pt1->y[i] / 6.0 + itv.pt2->y[i] * 2.0 / 3.0 + itv.pt3->y[i] / 6.0) * itv.h / 2.0;
            itv.simps_2[i] = (itv.pt3->y[i] / 6.0 + itv.pt4->y[i] * 2.0 / 3.0 + itv.pt5->y[i] / 6.0) * itv.h / 2.0;
            itv.simps_sum[i] = itv.simps_1[i] + itv.simps_2[i];
            itv.simps_tot[i] = (itv.pt1->y[i] / 6.0 + itv.pt3->y[i] * 2.0 / 3.0 + itv.pt5->y[i] / 6.0) * itv.h;
            integral[i] += itv.simps_sum[i];
        }
        intervals.push_back(itv);
    }

    // adaptive refinement
    bool should_refine = false;
    bool refined = false;
    Array1D integral_new = Array(n);
    for (int i = 0; i < max_iter; ++i) {
        integral_new = integral;
        refined = false;
        for (auto it = intervals.begin(); it != intervals.end(); ++it) {
            // test if refine is needed
            should_refine = false;
            for (int j = 0; j < n; ++j) {
                if (std::abs(it->simps_tot[j] - it->simps_sum[j]) > xtol * it->h / h_tot + rtol * std::abs(it->simps_sum[j])
                    && 
                    std::abs(it->simps_sum[j]) > xtol + rtol * std::abs(integral[j])) {
                    should_refine = true;
                    break;
                }
            }

            // refine
            if (should_refine) {
                refined = true;

                // refine left
                pt.x = (it->pt1->x + it->pt2->x) / 2.0;
                pt.y = f(pt.x);
                points.insert(it->pt2, pt);
                
                pt.x = (it->pt2->x + it->pt3->x) / 2.0;
                pt.y = f(pt.x);
                points.insert(it->pt3, pt);
                
                itv.pt1 = it->pt1;
                itv.pt2 = std::next(itv.pt1);
                itv.pt3 = std::next(itv.pt2);
                itv.pt4 = std::next(itv.pt3);
                itv.pt5 = std::next(itv.pt4);
                itv.h = it->h / 2.0;
                for (int j = 0; j < n; ++j) {
                    itv.simps_1[j] = (itv.pt1->y[j] / 6.0 + itv.pt2->y[j] * 2.0 / 3.0 + itv.pt3->y[j] / 6.0) * itv.h / 2.0;
                    itv.simps_2[j] = (itv.pt3->y[j] / 6.0 + itv.pt4->y[j] * 2.0 / 3.0 + itv.pt5->y[j] / 6.0) * itv.h / 2.0;
                    itv.simps_sum[j] = itv.simps_1[j] + itv.simps_2[j];
                    itv.simps_tot[j] = it->simps_1[j];
                    integral_new[j] -= it->simps_sum[j];
                    integral_new[j] += itv.simps_sum[j];
                }
                intervals.insert(it, itv);

                // refine right
                pt.x = (it->pt3->x + it->pt4->x) / 2.0;
                pt.y = f(pt.x);
                points.insert(it->pt4, pt);

                pt.x = (it->pt4->x + it->pt5->x) / 2.0;
                pt.y = f(pt.x);
                points.insert(it->pt5, pt);

                it->pt1 = it->pt3;
                it->pt2 = std::next(it->pt1);
                it->pt3 = std::next(it->pt2);
                it->pt4 = std::next(it->pt3);
                it->pt5 = std::next(it->pt4);
                it->h /= 2;
                for (int j = 0; j < n; ++j) {
                    it->simps_tot[j] = it->simps_2[j];
                    it->simps_1[j] = (it->pt1->y[j] / 6.0 + it->pt2->y[j] * 2.0 / 3.0 + it->pt3->y[j] / 6.0) * it->h / 2.0;
                    it->simps_2[j] = (it->pt3->y[j] / 6.0 + it->pt4->y[j] * 2.0 / 3.0 + it->pt5->y[j] / 6.0) * it->h / 2.0;
                    it->simps_sum[j] = it->simps_1[j] + it->simps_2[j];
                    integral_new[j] += it->simps_sum[j];
                }
            }
        }
        // if no refinement, return result
        if (!refined) {
            return integral;
        }
        integral = integral_new;
    }

    // solution doesn't converge in max_iter iterations
    //std::cout << "integral may not be accurate: maximum depth = " << max_iter << " is achieved. \n";
    throw std::runtime_error("Integration doesn't converge after " + std::to_string(max_iter) + " iterations!");
    return integral;
}

template<typename F>
Array1D Adaptive_2D(F& func, const Array1D& xini, const Array1D& yini, const double& xtol, const double& rtol, const int max_iter = 20) {
    auto f = [&](const double& y) {
        auto g = [&](const double& x) {
            return func(x, y);
        };
        // integrate over theta
        Array1D result = Adaptive_1D(g, xini, xtol, rtol, max_iter);
        return result;
    };
    return Adaptive_1D(f, yini, xtol, rtol, max_iter);
}

#endif