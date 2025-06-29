#include "solver_pcg.hpp"
#include "model.hpp"
#include "lmm.hpp"

PYBIND11_MODULE(evolm, m)
{
    m.doc() = "some description ..."; // Optional module docstring

    pybind11::class_<evolm::Model>(m, "Model")

        .def(pybind11::init<>())
#ifdef UTEST
        .def("size_of", &evolm::Model::size_of)
        .def("shape_of", &evolm::Model::shape_of)
        .def("print", &evolm::Model::print)
#endif
        .def("clear_residuals", &evolm::Model::clear_residuals)
        .def("clear_observations", &evolm::Model::clear_observations)
        .def("clear_effects", &evolm::Model::clear_effects)
        .def("clear_corrstruct", &evolm::Model::clear_corrstruct)
        .def("clear_traitstruct", &evolm::Model::clear_traitstruct)
        .def("clear", &evolm::Model::clear)        

        .def( "append_residual", static_cast< int (evolm::Model::*)(pybind11::array_t<float>, size_t) >(&evolm::Model::append_residual) )
        .def( "append_residual", static_cast< int (evolm::Model::*)(const std::string &) >(&evolm::Model::append_residual) )        
        .def( "append_observation", static_cast< int (evolm::Model::*)(pybind11::array_t<float>, size_t) >(&evolm::Model::append_observation) )
        .def( "append_observation", static_cast< int (evolm::Model::*)(const std::string &) >(&evolm::Model::append_observation) )
        
        .def( "append_effect", static_cast< int (evolm::Model::*)(pybind11::array_t<float>, size_t, size_t) >(&evolm::Model::append_effect) )

        .def( "append_effect", static_cast< int (evolm::Model::*)(const std::string &, size_t) >(&evolm::Model::append_effect) )
        .def( "append_effect", static_cast< int (evolm::Model::*)(const std::string &, bool) >(&evolm::Model::append_effect) )
        .def( "append_effect", static_cast< int (evolm::Model::*)(const std::string &, float) >(&evolm::Model::append_effect) )
        .def( "append_effect", static_cast< int (evolm::Model::*)(const std::string &, double) >(&evolm::Model::append_effect) )

        .def( "append_corrstruct", static_cast< int (evolm::Model::*)(pybind11::array_t<float>, size_t, pybind11::array_t<float>, size_t, pybind11::array_t<int>) >(&evolm::Model::append_corrstruct) )
        .def( "append_corrstruct", static_cast< int (evolm::Model::*)(pybind11::array_t<float>, size_t, const std::string &, pybind11::array_t<int>) >(&evolm::Model::append_corrstruct) )
        .def( "append_corrstruct", static_cast< int (evolm::Model::*)(const std::string &, const std::string &, pybind11::array_t<int>) >(&evolm::Model::append_corrstruct) )

        .def( "append_corrstruct", static_cast< int (evolm::Model::*)(pybind11::array_t<float>, size_t, std::string &, size_t, pybind11::array_t<int>) >(&evolm::Model::append_corrstruct) )
        .def( "append_corrstruct", static_cast< int (evolm::Model::*)(const std::string &, std::string &, size_t, pybind11::array_t<int>) >(&evolm::Model::append_corrstruct) )
        
        .def( "append_traitstruct", static_cast< int (evolm::Model::*)(int, pybind11::array_t<int>) >(&evolm::Model::append_traitstruct) );
    
    pybind11::class_<evolm::Pcg>(m, "Pcg")
        .def(pybind11::init<>())
        .def("append_model", &evolm::Pcg::append_model)
        .def("solve", pybind11::overload_cast<>(&evolm::Pcg::solve))
        .def("solve", pybind11::overload_cast<int>(&evolm::Pcg::solve))
        .def( "get_solution", static_cast< int (evolm::Pcg::*)(const std::string &) >(&evolm::Pcg::get_solution) )
        .def( "get_solution", static_cast< pybind11::array_t<float> (evolm::Pcg::*)() >(&evolm::Pcg::get_solution) );

    pybind11::class_<evolm::lmm>(m, "lmm")
        .def(pybind11::init<>())
        .def("clear", &evolm::lmm::clear)
        .def( "define", &evolm::lmm::define )
        .def( "define_infile", &evolm::lmm::define_infile )
        .def( "snp_to_bv", &evolm::lmm::snp_to_bv )
        .def( "solve", pybind11::overload_cast<const std::string &, int, int, const std::string &, const std::string &>(&evolm::lmm::solve) )
        .def( "solve", pybind11::overload_cast<const std::string &, int, int, const std::string &, const std::string &, double>(&evolm::lmm::solve) )
        .def( "solve", pybind11::overload_cast<const std::string &, int, int, const std::string &, const std::string &, double, size_t>(&evolm::lmm::solve) );
}
