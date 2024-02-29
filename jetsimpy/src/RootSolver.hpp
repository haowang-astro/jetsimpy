#ifndef ROOTSOLVER
#define ROOTSOLVER

#include "Definations.hpp"

// Brent's method. Code from Scipy.
template<typename F>
double brentq(F& f, double xa, double xb, double xtol, double rtol = 1e-8, int iter = 100) {
    using namespace std;
    double xpre = xa, xcur = xb;
    double xblk = 0., fpre, fcur, fblk = 0., spre = 0., scur = 0., sbis;
    // the tolerance is 2*delta
    double delta;
    double stry, dpre, dblk;
    int i;
    
    fpre = f(xpre);
    fcur = f(xcur);
    
    // check if the root is at boundary
    if (fabs(fpre) == 0) {
        return xpre;
    }
    else if (fabs(fcur) == 0) {
        return xcur;
    }

    // check if f(a) and f(b) have different signs
    //try {
    if (fpre * fcur > 0) {
        string text = "f(a) and f(b) must have different signs!\n";
        throw std::runtime_error(text);
    }

    // iteration
    for (i = 0; i < iter; i++) {
        
        if (fpre != 0 && fcur != 0 &&
	    (signbit(fpre) != signbit(fcur))) {
            xblk = xpre;
            fblk = fpre;
            spre = scur = xcur - xpre;
        }
        if (fabs(fblk) < fabs(fcur)) {
            xpre = xcur;
            xcur = xblk;
            xblk = xpre;

            fpre = fcur;
            fcur = fblk;
            fblk = fpre;
        }

        delta = (xtol + rtol * fabs(xcur)) / 2;
        sbis = (xblk - xcur) / 2;
        if (fcur == 0 || fabs(sbis) < delta) {
            return xcur;
        }

        if (fabs(spre) > delta && fabs(fcur) < fabs(fpre)) {
            if (xpre == xblk) {
                // interpolate
                stry = -fcur * (xcur - xpre) / (fcur - fpre);
            }
            else {
                // extrapolate
                dpre = (fpre - fcur) / (xpre - xcur);
                dblk = (fblk - fcur) / (xblk - xcur);
                stry = - fcur * (fblk * dblk - fpre * dpre) 
                       / (dblk * dpre * (fblk - fpre));
            }
            if (2 * fabs(stry) < min(fabs(spre), 3 * fabs(sbis) - delta)) {
                // good short step
                spre = scur;
                scur = stry;
            } else {
                // bisect
                spre = sbis;
                scur = sbis;
            }
        }
        else {
            // bisect
            spre = sbis;
            scur = sbis;
        }

        xpre = xcur; fpre = fcur;
        if (fabs(scur) > delta) {
            xcur += scur;
        }
        else {
            xcur += (sbis > 0 ? delta : -delta);
        }

        fcur = f(xcur);
    }

    // the solver does not converge after xxx iterations!
    string text = "Solver doesn't converge after " + to_string(iter) + " interations!\n";
    throw std::runtime_error(text);

    return xcur;
}

template<typename F1, typename F2>
double newton_raphson(F1& f, F2& df_dx, double x0, double xtol = 1e-8, int iter = 100) {
    double xnew;

    for (int i = 0; i < iter; ++i) {
        xnew = x0 - f(x0) / df_dx(x0);

        if (std::abs(xnew - x0) < xtol) {
            return xnew;
        }
        else {
            x0 = xnew;
        }
    }

    // the solver does not converge after xxx iterations!
    std::string text = "Solver doesn't converge after " + std::to_string(iter) + " interations!\n";
    throw std::runtime_error(text);

    return xnew;
}
#endif
