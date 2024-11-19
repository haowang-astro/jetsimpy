#include "jetsimpy.h"

PYBIND11_MODULE(jetsimpy_extension, m) {
    // parameter object
    py::class_<JetConfig>(m, "JetConfig")
        .def(py::init<>())
        .def_readwrite("theta_edge", &JetConfig::theta_edge)
        .def_readwrite("Eb", &JetConfig::Eb)
        .def_readwrite("Ht", &JetConfig::Ht)
        .def_readwrite("Msw", &JetConfig::Msw)
        .def_readwrite("Mej", &JetConfig::Mej)
        .def_readwrite("R", &JetConfig::R)
        .def_readwrite("nwind", &JetConfig::nwind)
        .def_readwrite("nism", &JetConfig::nism)
        .def_readwrite("tmin", &JetConfig::tmin)
        .def_readwrite("tmax", &JetConfig::tmax)
        .def_readwrite("rtol", &JetConfig::rtol)
        .def_readwrite("cfl", &JetConfig::cfl)
        .def_readwrite("spread", &JetConfig::spread)
        .def_readwrite("cal_level", &JetConfig::cal_level)
    ;

    // jet & afterglow calculation class
    py::class_<Jet>(m, "Jet")
        .def(py::init<const JetConfig&>())
        // ---------- solve hydro ---------- //
        .def("solveJet", &Jet::solveJet)
        .def("getY", &Jet::getY)
        .def("getT", &Jet::getT)
        .def("getTheta", &Jet::getTheta)

        // ---------- hydro interpolation ---------- //
        .def("interpolateMsw", py::vectorize(&Jet::interpolateMsw))
        .def("interpolateMej", py::vectorize(&Jet::interpolateMej))
        .def("interpolateBetaGamma", py::vectorize(&Jet::interpolateBetaGamma))
        .def("interpolateBetaTh", py::vectorize(&Jet::interpolateBetaTh))
        .def("interpolateR", py::vectorize(&Jet::interpolateR))
        .def("interpolateE0", py::vectorize(&Jet::interpolateE0))

        // ---------- afterglow calculation ---------- //
        .def("configParameters", &Jet::configParameters)
        .def("configIntensity", &Jet::configIntensity)
        .def("configAvgModel", &Jet::configAvgModel)
        .def("configIntensityPy", &Jet::configIntensityPy)
        .def("configAvgModelPy", &Jet::configAvgModelPy)
        .def("calculateEATS", py::vectorize(&Jet::calculateEATS))
        .def("calculateIntensity", py::vectorize(&Jet::calculateIntensity))
        .def("calculateLuminosity", py::vectorize(&Jet::calculateLuminosity))
        .def("WeightedAverage", py::vectorize(&Jet::WeightedAverage))
        .def("IntensityOfPixel", py::vectorize(&Jet::IntensityOfPixel))
    ;

    // bind blast object
    py::class_<Blast>(m, "Blast")
        .def_readonly("t", &Blast::t)
        .def_readonly("theta", &Blast::theta)
        .def_readonly("phi", &Blast::phi)
        .def_readonly("R", &Blast::R)
        .def_readonly("beta", &Blast::beta)
        .def_readonly("gamma", &Blast::gamma)
        .def_readonly("beta_th", &Blast::beta_th)
        .def_readonly("beta_r", &Blast::beta_r)
        .def_readonly("beta_f", &Blast::beta_f)
        .def_readonly("gamma_f", &Blast::gamma_f)
        .def_readonly("s", &Blast::s)
        .def_readonly("doppler", &Blast::doppler)
        .def_readonly("n_blast", &Blast::n_blast)
        .def_readonly("e_density", &Blast::e_density)
        .def_readonly("pressure", &Blast::pressure)
        .def_readonly("n_ambient", &Blast::n_ambient)
        .def_readonly("dR", &Blast::dR)
    ;
}
