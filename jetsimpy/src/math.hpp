#ifndef MATH
#define MATH

#include <vector>

// ------------------ //
//      functions     //
// ------------------ //

double Minmod(const double& x1, const double& x2) {
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

// Binary search (without bound check)
void Find_Index(const std::vector<double>& x_array, const double& x, int& index1, int& index2) {
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

// ---------- density ---------- //
double density0(const double& nwind, const double& nism, const double& r) {
    return nwind / (r / 1e17) / (r / 1e17) + nism;
}

double density_k(const double& nwind, const double& nism, const double& r) {
    return 2.0 * nwind / (nwind + nism * (r / 1e17) * (r / 1e17));
}

double calib_BM(const double& nwind, const double& nism, const double& r) {
    double k = density_k(nwind, nism, r);
    return 0.52935729 - 0.05698377 * k - 0.00158176 * k * k - 0.00939548 * k * k * k;
}

double calib_ST(const double& nwind, const double& nism, const double& r) {
    double k = density_k(nwind, nism, r);
    return 1.635 - 0.651 * k;
}

// -------------------------------- Linear interpolation ----------------------------------- //
double Linear(const double& x1, const double& x2, const double& f1, const double& f2, const double& x) {
    return (f1 * (x2 - x) + f2 * (x - x1)) / (x2 - x1);
}

double Bilinear(const double& x1, const double& x2, const double& y1, const double& y2, 
              const double& f11, const double& f12, const double& f21, const double& f22,
              const double& x, const double& y) {
    return (f11 * (x2 - x) * (y2 - y) + f21 * (x - x1) * (y2 - y) + f12 * (x2 - x) * (y - y1) + f22 * (x - x1) * (y - y1))
           / (x2 - x1) / (y2 - y1);
}

#endif