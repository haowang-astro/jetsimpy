#ifndef BLAST
#define BLAST

struct Blast {
public:
    // observer's parameters
    //double Tobs_z;          // Tobs / (1 + z)

    // coordinate values (burster frame)
    double t;
    double theta;
    double phi;
    double R;

    // blast velocity (burster frame)
    double beta;
    double gamma;
    double beta_th;
    double beta_r;
    double beta_f;
    double gamma_f;
    double s;

    // angles (burster frame)
    double doppler;
    //double cos_theta_r;
    double cos_theta_beta;

    // thermodynamic values (comoving frame)
    double n_blast;
    double e_density;
    double pressure;
    double n_ambient;
    double dR;
};

#endif