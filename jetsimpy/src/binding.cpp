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
        .def_readwrite("calib_level", &JetConfig::calib_level)
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
        .def("configEmissivity", &Jet::configEmissivity)
        .def("configAvgModel", &Jet::configAvgModel)
        .def("calculateEATS", py::vectorize(&Jet::calculateEATS))
        .def("calculateIntensity", py::vectorize(&Jet::calculateIntensity))
        .def("calculateLuminosity", py::vectorize(&Jet::calculateLuminosity))
        .def("WeightedAverage", py::vectorize(&Jet::WeightedAverage))
        //.def("calculateAvgModel", py::vectorize(&Jet::calculateAvgModel))
        .def("IntensityOfPixel", py::vectorize(&Jet::IntensityOfPixel))
    ;

}
