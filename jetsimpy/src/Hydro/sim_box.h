#ifndef SIM_BOX
#define SIM_BOX

#include "../environment.h"
#include "config.h"
#include "tools.h"

class SimBox {
public:
    SimBox(const JetConfig& jet_config, Tool& tool);

    // ---------- interfaces ---------- //
    void solvePDE();             // solve PDE
    Array3D& getY();             // get PDE ys
    Array1D& getT();             // get PDE ts
    Array1D& getTheta();         // get cell centers
private:
    // tools
    Tool& tool;                  // useful functions

    // configuration
    double tmin;                 // minimum pde time
    double tmax;                 // maximum pde time
    double cfl;                  // cfl number
    bool spread;                 // spread or not

    // mesh
    int ntheta;                  // number of cells
    Array1D theta;               // cell center position
    Array1D theta_edge;          // cell edge position

    // conserved variables
    Array1D Eb;                  // Eb
    Array1D Ht;                  // Hb * beta_th
    Array1D Msw;                 // Msw
    Array1D Mej;                 // Mej
    Array1D R;                   // R

    // primitive variables
    Array1D beta_gamma_sq;       // (beta * gamma) ^ 2
    Array1D beta_th;             // beta_theta
    Array1D Psw;                 // Psw
    Array1D Hb;                  // Hb
    Array1D s;                   // s

    // variables for convinience
    Array1D beta;                // beta
    Array1D gamma;               // gamma

    // eigenvalues
    Array1D eigenvalues;         // maximum eigenvalues of the four
    Array1D alpha_R;             // viscosity for R equation

    // solpe
    Array2D slope;                       // (5, ntheta) [Msw, Mej, u_sq, beta_th, R] slope
    Array1D R_slope_l;                   // left biased R slope
    Array1D R_slope_r;                   // right biased R slope

    // numerical flux
    Array2D numerical_flux;              // (4, ntheta + 1) numerical flux
    Array1D dR_dt;                       // dR / dt

    // dy / dt
    Array2D dy_dt;                       // (5, ntheta) dy/dt

    // PDE solution
    Array3D ys;                          // (5, ntheta, nt)
    Array1D ts;                          // (nt)

    // ---------- functions ---------- //
    void solvePrimitive();               // solve primitive and convenient variables
    void solveEigen();                   // solve eigenvalues
    void solveSlope();                   // reconstruct
    void solveNumericalFlux();           // solve riemann problem
    double solveDeltaT();                // solve delta_t
    void solveDyDt();                    // solve dy_dt (with spread)
    void solveDyDt_no_spread();          // solve dy_dt (without spread)
    void oneStepRK2(double dt);          // one step forward (with spread)
    void oneStepRK45(double& dt, const double rtol, bool& succeeded);        // one step forward (without spread)
    void solveSpread();                  // solve PDE with spreading
    void solveNoSpread();                // solve PDE without spreading
};

#endif