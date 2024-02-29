#include "Definations.hpp"
#include "Jet.hpp"
#include "Grbafterglow.hpp"
#include "Refinement.hpp"

PYBIND11_MODULE(jetsimpy_extension, m) {
    // load module data
    m.def("Refine1", &Refine1);
    m.def("Refine2", &Refine2);

    // interpolator class
    py::class_<Interpolator>(m, "Interpolator")
        .def("Eb0", py::vectorize(&Interpolator::Eb0))
        .def("Msw", py::vectorize(&Interpolator::Msw))
        .def("Mej", py::vectorize(&Interpolator::Mej))
        .def("beta_gamma", py::vectorize(&Interpolator::beta_gamma))
        .def("beta_th", py::vectorize(&Interpolator::beta_th))
        .def("R", py::vectorize(&Interpolator::R))
        ;

    // imaging class
    py::class_<AfterglowImage>(m, "afterglow_image")
        .def_readonly("half_width", &AfterglowImage::half_width)
        .def_readonly("xscale", &AfterglowImage::xscale)
        .def_readonly("yscale", &AfterglowImage::yscale)
        .def_readonly("offset", &AfterglowImage::offset)
        .def_readonly("intensity_image", &AfterglowImage::intensity_map)
        .def_readonly("polarization_scale_image", &AfterglowImage::pi_scale_map)
        .def_readonly("polarization_angle_image", &AfterglowImage::pi_angle_map)
        ;

    // jet & afterglow calculation class
    py::class_<GRBafterglow>(m, "afterglow")
        .def(py::init<>())

        // jet module
        .def("config_jet", &GRBafterglow::config_jet)
        .def("calibrate_jet", &GRBafterglow::calibrate_jet)
        .def("solve_jet", &GRBafterglow::solve_jet)
        .def("get_t_pde", &GRBafterglow::get_t_pde)
        .def("get_y_pde", &GRBafterglow::get_y_pde)
        .def("get_interpolator", &GRBafterglow::get_interpolator)

        // integrator module
        .def("config_integrator", &GRBafterglow::config_integrator)
        .def("FluxDensity", py::vectorize(&GRBafterglow::FluxDensity))
        .def("Offset", py::vectorize(&GRBafterglow::Offset))
        .def("Xscale", py::vectorize(&GRBafterglow::Xscale))
        .def("Yscale", py::vectorize(&GRBafterglow::Yscale))
        .def("Pi_lin", py::vectorize(&GRBafterglow::Pi_lin))
        .def("intensity", py::vectorize(&GRBafterglow::intensity))
        .def("test_function", py::vectorize(&GRBafterglow::test_function))

        // image module
        .def("config_image", &GRBafterglow::config_image)
        .def("solve_image", &GRBafterglow::solve_image)
        .def("image", &GRBafterglow::image)
        ;
}
