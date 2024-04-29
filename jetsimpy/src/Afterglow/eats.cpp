#include "eats.h"

void EATS::feedData(SimBox& sim_box, Tool& tool) {
    y_data = &(sim_box.getY());
    t_data = &(sim_box.getT());
    theta_data = &(sim_box.getTheta());
    this->tool = &tool;
    ntheta = theta_data->size();
    nt = t_data->size();
    tmin = t_data->front();
    tmax = t_data->back();
    theta_min = theta_data->front();
    theta_max = theta_data->back();
}

void EATS::findThetaIndex(double theta, int& theta_index1, int& theta_index2) {
    if (theta < 0.0 || theta > PI) {
        // bound check error
        throw std::runtime_error("EATS: theta outside bounds.\n");
    }
    else if (theta < theta_min) {
        // north pole
        theta_index1 = 0;
        theta_index2 = 0;
    }
    else if (theta > theta_max) {
        // south pole
        theta_index1 = ntheta - 1;
        theta_index2 = ntheta - 1;
    }
    else {
        // middle
        tool->findIndex(*theta_data, theta, theta_index1, theta_index2);
    }
}

void EATS::findTimeIndex(double mu, double Tobs_z, int theta_index, int& t_index1, int& t_index2) {
    // equal arrival time surface function at t_index
    auto f = [&](const int& t_index) {
        // r at theta_index, t_index
        double r = (*y_data)[4][theta_index][t_index];

        // t at theta_index, r_index
        double t = (*t_data)[t_index];

        // equal arrival time surface function
        return t - r * mu / CSpeed - Tobs_z;
    };

    // find t index (binary search)
    if (f(0) > 0) { // smaller than tmin
        t_index1 = 0;
        t_index2 = 0;
    }
    else if (f(nt - 1) < 0) { // larger than tmax
        // throw error
        throw std::runtime_error("EATS: Observing time exceeds PDE maximum evolution time!\n");
    }
    else {
        t_index1 = 0;
        t_index2 = nt - 1;
        int index_mid;
        while (t_index2 - t_index1 > 1) {
            index_mid = (t_index1 + t_index2) / 2;
            if (f(index_mid) > 0) {
                t_index2 = index_mid;
            }
            else {
                t_index1 = index_mid;
            }
        }
    }
}

double EATS::solveT(double mu, double Tobs_z, int theta_index, int t_index1, int t_index2) {
    if (t_index1 == t_index2) { // before tmin
        double r = (*y_data)[4][theta_index][t_index1];
        double t = Tobs_z / (1 - r / (CSpeed * tmin) * mu);
        return t;
    }
    else {
        // t at two ends
        double t1 = (*t_data)[t_index1];
        double t2 = (*t_data)[t_index2];

        // r values at two ends
        double r1 = (*y_data)[4][theta_index][t_index1];
        double r2 = (*y_data)[4][theta_index][t_index2];

        // linear interpolation
        double slope = (r2 - r1) / (t2 - t1);
        double t = (Tobs_z + (r1 - slope * t1) * mu / CSpeed) / (1 - slope * mu / CSpeed);
        return t;
    }
}

Array1D EATS::solvePrimitive(double mu, double Tobs_z, int theta_index) {
    // find t_index
    int t_index1, t_index2;
    findTimeIndex(mu, Tobs_z, theta_index, t_index1, t_index2);

    // solve t
    double t = solveT(mu, Tobs_z, theta_index, t_index1, t_index2);

    // solve primitive
    Array1D val = Array(6);      // [Msw, Mej, beta_gamma_sq, beta_th, R, t]
    for (int i = 0; i < 5; ++i) {
        val[i] = tool->linear(
            t, (*t_data)[t_index1], (*t_data)[t_index2], 
            (*y_data)[i][theta_index][t_index1], (*y_data)[i][theta_index][t_index2]
        );
    }
    val[5] = t;

    return val;
}

void EATS::deriveBlast(double theta, double phi, double theta_v, const Array1D& val, Blast& blast) {
    // PDE variables
    double Msw = val[0];
    double beta_gamma = std::sqrt(val[2]);
    double beta_th = val[3];
    double R = val[4];

    // coordinate values (burster frame)
    blast.t = val[5];
    blast.theta = theta;
    blast.phi = phi;
    blast.R = R;

    // blast velocity (burster frame)
    blast.gamma = std::sqrt(val[2] + 1.0);
    blast.beta = beta_gamma / blast.gamma;
    blast.beta_th = beta_th;
    blast.beta_r = (beta_th <= blast.beta) ? std::sqrt(blast.beta * blast.beta - beta_th * beta_th) : 0.0;
    blast.beta_f = 4.0 * blast.beta * (val[2] + 1) / (4.0 * val[2] + 3.0);
    blast.gamma_f = (4.0 * val[2] + 3.0) / std::sqrt(8 * val[2] + 9.0);
    blast.s = tool->solveS(R, val[2]);

    // angles
    double nr = blast.beta_r / blast.beta;
    double nth = beta_th / blast.beta;
    double mu_beta = (nr * sin(theta) * cos(phi) + nth * cos(theta) * cos(phi)) * sin(theta_v)
                   + (nr * cos(theta) - nth * sin(theta)) * cos(theta_v);
    double mu_r = std::cos(theta) * std::cos(theta_v) + std::sin(theta) * std::cos(phi) * std::sin(theta_v);
    blast.doppler = 1.0 / blast.gamma / (1.0 - blast.beta * mu_beta);
    blast.cos_theta_beta = (mu_beta - blast.beta) / (1.0 - blast.beta * mu_beta);

    // thermaldynamic properties (comoving frame)
    blast.n_ambient = tool->solveDensity(R);
    blast.n_blast = 4.0 * blast.gamma * blast.n_ambient;
    blast.e_density = (blast.gamma - 1.0) * blast.n_blast * MassP * CSpeed * CSpeed * blast.s;
    blast.pressure = 4.0 / 3.0 * val[2] * blast.n_ambient * MassP * CSpeed * CSpeed * blast.s;
    blast.dR = Msw / R / R / blast.n_blast / MassP;
}

void EATS::solveBlast(double Tobs_z, double theta, double phi, double theta_v, Blast& blast) {
    // find theta index
    int theta_index1;
    int theta_index2;
    findThetaIndex(theta, theta_index1, theta_index2);

    Array1D val;
    if (theta_index1 == theta_index2) { // near poles
        // compute mu
        double mu = std::cos((*theta_data)[theta_index1]) * std::cos(theta_v) + std::sin((*theta_data)[theta_index1]) * std::cos(phi) * std::sin(theta_v);

        // solve primitive
        val = solvePrimitive(mu, Tobs_z, theta_index1);

        // construct blast object
        deriveBlast(theta, phi, theta_v, val, blast);
    }
    else {
        Array1D val_l, val_r;
        // solve val_l
        {
            // compute mu
            double mu = std::cos((*theta_data)[theta_index1]) * std::cos(theta_v) + std::sin((*theta_data)[theta_index1]) * std::cos(phi) * std::sin(theta_v);

            // solve primitive
            val_l = solvePrimitive(mu, Tobs_z, theta_index1);
        }

        // solve val_r
        {
            // compute mu
            double mu = std::cos((*theta_data)[theta_index2]) * std::cos(theta_v) + std::sin((*theta_data)[theta_index2]) * std::cos(phi) * std::sin(theta_v);

            // solve primitive
            val_r = solvePrimitive(mu, Tobs_z, theta_index2);
        }

        // interpolate val over theta
        val = Array(6);
        for (int i = 0; i < 6; ++i) {
            val[i] = tool->linear(
                theta,
                (*theta_data)[theta_index1],
                (*theta_data)[theta_index2],
                val_l[i],
                val_r[i]
            );
        }

        // construct blast object
        deriveBlast(theta, phi, theta_v, val, blast);
    }
}

double EATS::solveEATS(double Tobs_z, double theta, double phi, double theta_v) {
    // find theta index
    int theta_index1;
    int theta_index2;
    findThetaIndex(theta, theta_index1, theta_index2);

    Array1D val;
    if (theta_index1 == theta_index2) { // near poles
        // compute mu
        double mu = std::cos((*theta_data)[theta_index1]) * std::cos(theta_v) + std::sin((*theta_data)[theta_index1]) * std::cos(phi) * std::sin(theta_v);

        // solve primitive
        val = solvePrimitive(mu, Tobs_z, theta_index1);
    }
    else {
        Array1D val_l, val_r;
        // solve val_l
        {
            // compute mu
            double mu = std::cos((*theta_data)[theta_index1]) * std::cos(theta_v) + std::sin((*theta_data)[theta_index1]) * std::cos(phi) * std::sin(theta_v);

            // solve primitive
            val_l = solvePrimitive(mu, Tobs_z, theta_index1);
        }

        // solve val_r
        {
            // compute mu
            double mu = std::cos((*theta_data)[theta_index2]) * std::cos(theta_v) + std::sin((*theta_data)[theta_index2]) * std::cos(phi) * std::sin(theta_v);

            // solve primitive
            val_r = solvePrimitive(mu, Tobs_z, theta_index2);
        }

        // interpolate val over theta
        val = Array(6);
        for (int i = 0; i < 6; ++i) {
            val[i] = tool->linear(
                theta,
                (*theta_data)[theta_index1],
                (*theta_data)[theta_index2],
                val_l[i],
                val_r[i]
            );
        }
    }

    return val[5];
}