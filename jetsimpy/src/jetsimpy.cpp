#include "jetsimpy.h"

Jet::Jet(const JetConfig& jet_config)
  : jet_config (jet_config),
    tool (jet_config),
    sim_box (jet_config, tool)
{

}

// ---------- solve hydro ---------- //

void Jet::solveJet() {
    // solve PDE
    sim_box.solvePDE();

    // feed data to interpolator
    interpolator.feedData(sim_box, tool);

    // feed data to eats solver
    eats.feedData(sim_box, tool);

    // initialize afterglow object
    afterglow.initialize(sim_box, eats);
}

py::array_t<double> Jet::getY() {
    // get data
    Array3D& y_data = sim_box.getY();

    // create np.ndarray object
    size_t nt = sim_box.getT().size();
    size_t ntheta = sim_box.getTheta().size();
    size_t shape[3] = {5, ntheta, nt};
    size_t strides[3] = {ntheta * nt * sizeof(double), nt * sizeof(double), sizeof(double)};
    auto array = py::array_t<double>(
        shape, 
        strides
    );

    // register memory layout
    auto view = array.mutable_unchecked<3>();
    for (size_t i = 0; i < 5; ++i) {
        for (size_t j = 0; j < ntheta; ++j) {
            for (size_t k = 0; k < nt; ++k) {
                view(i, j, k) = y_data[i][j][k];
            }
        }
    }

    return array;
}

py::array_t<double> Jet::getT() {
    // get data
    Array1D& t_data = sim_box.getT();

    // create np.ndarray object
    size_t nt = sim_box.getT().size();
    size_t shape[1] = {nt};
    size_t strides[1] = {sizeof(double)};
    auto array = py::array_t<double>(
        shape, 
        strides
    );

    // register memory layout
    auto view = array.mutable_unchecked<1>();
    for (size_t i = 0; i < array.shape(0); ++i) {
        view(i) = t_data[i];
    }

    return array;
}

py::array_t<double> Jet::getTheta() {
    // get data
    Array1D& theta_data = sim_box.getTheta();

    // create np.ndarray object
    auto array = py::array_t<double>(
        {theta_data.size()}, 
        {sizeof(double)}
    );

    // register memory layout
    auto view = array.mutable_unchecked<1>();
    for (size_t i = 0; i < array.shape(0); ++i) {
        view(i) = theta_data[i];
    }

    return array;
}

// ---------- hydro interpolation ---------- //

double Jet::interpolateMsw(double t, double theta) {
    return interpolator.interpolateY(t, theta, 0);
}

double Jet::interpolateMej(double t, double theta) {
    return interpolator.interpolateY(t, theta, 1);
}

double Jet::interpolateBetaGamma(double t, double theta) {
    // remember what we record and interpolate is (beta * gamma)^2
    return std::sqrt(interpolator.interpolateY(t, theta, 2));
}

double Jet::interpolateBetaTh(double t, double theta) {
    return interpolator.interpolateY(t, theta, 3);
}

double Jet::interpolateR(double t, double theta) {
    return interpolator.interpolateY(t, theta, 4);
}

double Jet::interpolateE0(double t, double theta) {
    double msw = interpolateMsw(t, theta);
    double mej = interpolateMej(t, theta);
    double beta_gamma_sq = std::pow(interpolateBetaGamma(t, theta), 2.0);
    double R = interpolateR(t, theta);
    double s = tool.solveS(R, beta_gamma_sq);

    double E0 = s * (1.0 + beta_gamma_sq * beta_gamma_sq / (beta_gamma_sq + 1.0) / (beta_gamma_sq + 1.0) / 3.0) * (beta_gamma_sq + 1.0) * msw
              + (1.0 - s) * std::sqrt(beta_gamma_sq + 1.0) * msw
              + std::sqrt(beta_gamma_sq + 1.0) * mej
              - msw
              - mej;
    return E0 * CSpeed * CSpeed;
}

// ---------- afterglow calculation ---------- //

double Jet::calculateEATS(double Tobs, double theta, double phi, double theta_v, double z) {
    double Tobs_z = Tobs / (1.0 + z);
    return eats.solveEATS(Tobs_z, theta, phi, theta_v);
}

double Jet::calculateIntensity(double Tobs, double nu, double theta, double phi) {
    return afterglow.Intensity(Tobs, nu, theta, phi);
}

void Jet::configParameters(const Dict& param) {
    afterglow.configParameters(param);
}

void Jet::configEmissivity(const std::string& model_name) {
    afterglow.configEmissivity(model_name);
}

void Jet::configAvgModel(const std::string& model_name) {
    afterglow.configAvgModel(model_name);
}

void Jet::configEmissivityPy(py::function py_f) {
    afterglow.configEmissivityPy(py_f);
}

void Jet::configAvgModelPy(py::function py_f) {
    afterglow.configAvgModelPy(py_f);
}

double Jet::calculateLuminosity(double Tobs, double nu, double rtol) {
    return afterglow.Luminosity(Tobs, nu, rtol);
}

double Jet::calculateAvgModel(double Tobs, double nu, double rtol) {
    return afterglow.integrateModel(Tobs, nu, rtol);
}

double Jet::WeightedAverage(double Tobs, double nu, double rtol) {
    // save time
    if (Tobs == 0.0) return 0.0;

    double luminosity = calculateLuminosity(Tobs, nu, rtol);
    double integral = calculateAvgModel(Tobs, nu, rtol);

    return integral / luminosity;
}

double Jet::IntensityOfPixel(const double Tobs, const double nu, const double x_tilde, const double y_tilde) {
    return afterglow.IntensityOfPixel(Tobs, nu, x_tilde, y_tilde);
}
