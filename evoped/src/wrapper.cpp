#include "Amat.hpp"
#include "Gmat.hpp"
#include "Hmat.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#endif

PYBIND11_MODULE(evoped, m)
{
    m.doc() = "some description ...";

    pybind11::class_< evoped::Amat<double> >( m, "Amatd" )
        .def(pybind11::init<>())
        .def(pybind11::init<double>())
        .def("make_matrix", pybind11::overload_cast<const std::string &, bool>(&evoped::Amat<double>::make_matrix))
        .def("make_matrix", pybind11::overload_cast<const std::string &, const std::string &, bool>(&evoped::Amat<double>::make_matrix))
        .def("make_matrix_forgenotyped", pybind11::overload_cast<const std::string &, const std::string &, bool>(&evoped::Amat<double>::make_matrix_forgenotyped))
        .def("make_all", pybind11::overload_cast<const std::string &, const std::string &>(&evoped::Amat<double>::make_all))
        .def("get_inbreeding", pybind11::overload_cast<const std::string &>(&evoped::Amat<double>::get_inbreeding))
        .def("clear", &evoped::Amat<double>::clear)
        .def("set_sparsiity_threshold", &evoped::Amat<double>::set_sparsiity_threshold)
        .def("save_matrix", &evoped::Amat<double>::save_matrix)
        .def("save_ids", &evoped::Amat<double>::save_ids)
        ;
    
    pybind11::class_<evoped::Amat<float> >( m, "Amat" )
        .def(pybind11::init<>())
        .def(pybind11::init<double>())
        .def("make_matrix", pybind11::overload_cast<const std::string &, bool>(&evoped::Amat<float>::make_matrix))
        .def("make_matrix", pybind11::overload_cast<const std::string &, const std::string &, bool>(&evoped::Amat<float>::make_matrix))
        .def("make_matrix_forgenotyped", pybind11::overload_cast<const std::string &, const std::string &, bool>(&evoped::Amat<float>::make_matrix_forgenotyped))
        .def("make_all", pybind11::overload_cast<const std::string &, const std::string &>(&evoped::Amat<float>::make_all))
        .def("get_inbreeding", pybind11::overload_cast<const std::string &>(&evoped::Amat<float>::get_inbreeding))
        .def("clear", &evoped::Amat<float>::clear)
        .def("set_sparsiity_threshold", &evoped::Amat<float>::set_sparsiity_threshold)
        .def("save_matrix", &evoped::Amat<float>::save_matrix)
        .def("save_ids", &evoped::Amat<float>::save_ids)
        ;
    
    pybind11::class_< evoped::Gmat<double> >( m, "Gmatd" )
        .def(pybind11::init<>())
        .def("clear", &evoped::Gmat<double>::clear)
        .def("read_matrix", pybind11::overload_cast<const std::string &>(&evoped::Gmat<double>::read_matrix))
        .def("make_matrix", pybind11::overload_cast<const std::string &>(&evoped::Gmat<double>::make_matrix))
        .def("make_matrix", pybind11::overload_cast<const std::string &, const std::string &>(&evoped::Gmat<double>::make_matrix))
        .def("scale_genotypes", pybind11::overload_cast<const std::string &>(&evoped::Gmat<double>::scale_genotypes))
        .def("scale_genotypes", pybind11::overload_cast<const std::string &, const std::string &>(&evoped::Gmat<double>::scale_genotypes))
        .def("scale_genotypes", pybind11::overload_cast<const std::string &, const std::string &, const std::string &>(&evoped::Gmat<double>::scale_genotypes))
        .def("scale_genotypes", pybind11::overload_cast<const std::string &, const std::string &, const std::string &, const std::string &>(&evoped::Gmat<double>::scale_genotypes))
        //.def("impute_genotypes", pybind11::overload_cast<const std::string &, const std::string &, const std::string &>(&evoped::Gmat<double>::impute_genotypes))
        .def("impute_genotypes", pybind11::overload_cast<const std::string &, const std::string &, const std::string &, const std::string &, const std::string &>(&evoped::Gmat<double>::impute_genotypes))
        .def("scale_diag", &evoped::Gmat<double>::scale_diag)
        .def("scale_matrix", pybind11::overload_cast<const std::string &, double>(&evoped::Gmat<double>::scale_matrix))
        .def("invert_matrix", pybind11::overload_cast<bool>(&evoped::Gmat<double>::invert_matrix))
        .def("invert_matrix", pybind11::overload_cast<>(&evoped::Gmat<double>::invert_matrix))
        .def("save_matrix", pybind11::overload_cast<const std::string &>(&evoped::Gmat<double>::save_matrix))
        .def("save_matrix2", pybind11::overload_cast<const std::string &>(&evoped::Gmat<double>::save_matrix2))
        .def("save_ids", &evoped::Gmat<double>::save_ids)
        ;

    pybind11::class_< evoped::Gmat<float> >( m, "Gmat" )
        .def(pybind11::init<>())
        .def("clear", &evoped::Gmat<float>::clear)
        .def("read_matrix", pybind11::overload_cast<const std::string &>(&evoped::Gmat<float>::read_matrix))
        .def("make_matrix", pybind11::overload_cast<const std::string &>(&evoped::Gmat<float>::make_matrix))
        .def("make_matrix", pybind11::overload_cast<const std::string &, const std::string &>(&evoped::Gmat<float>::make_matrix))
        .def("scale_genotypes", pybind11::overload_cast<const std::string &>(&evoped::Gmat<float>::scale_genotypes))
        .def("scale_genotypes", pybind11::overload_cast<const std::string &, const std::string &>(&evoped::Gmat<float>::scale_genotypes))
        .def("scale_genotypes", pybind11::overload_cast<const std::string &, const std::string &, const std::string &>(&evoped::Gmat<float>::scale_genotypes))
        .def("scale_genotypes", pybind11::overload_cast<const std::string &, const std::string &, const std::string &, const std::string &>(&evoped::Gmat<float>::scale_genotypes))
        //.def("impute_genotypes", pybind11::overload_cast<const std::string &, const std::string &, const std::string &>(&evoped::Gmat<float>::impute_genotypes))
        .def("impute_genotypes", pybind11::overload_cast<const std::string &, const std::string &, const std::string &, const std::string &, const std::string &>(&evoped::Gmat<float>::impute_genotypes))
        .def("scale_diag", &evoped::Gmat<float>::scale_diag)
        .def("scale_matrix", pybind11::overload_cast<const std::string &, float>(&evoped::Gmat<float>::scale_matrix))
        .def("invert_matrix", pybind11::overload_cast<bool>(&evoped::Gmat<float>::invert_matrix))
        .def("invert_matrix", pybind11::overload_cast<>(&evoped::Gmat<float>::invert_matrix))
        .def("save_matrix", pybind11::overload_cast<const std::string &>(&evoped::Gmat<float>::save_matrix))
        .def("save_matrix2", pybind11::overload_cast<const std::string &>(&evoped::Gmat<float>::save_matrix2))
        .def("save_ids", &evoped::Gmat<float>::save_ids)
        ;

    pybind11::class_< evoped::Hmat<double> >( m, "Hmatd" )
        .def(pybind11::init<>())
        .def("clear", &evoped::Hmat<double>::clear)
        .def("make_matrix", pybind11::overload_cast<const std::string &, const std::string &, const std::string &, const std::string &, const std::string &, const std::string &>(&evoped::Hmat<double>::make_matrix))
        .def("make_matrix", pybind11::overload_cast<const std::string &, const std::string &, const std::string &, const std::string &>(&evoped::Hmat<double>::make_matrix))
        .def("make_matrix", pybind11::overload_cast<const std::string &, const std::string &, const std::string &>(&evoped::Hmat<double>::make_matrix))
        .def("save_matrix", pybind11::overload_cast<const std::string &>(&evoped::Hmat<double>::save_matrix))
        ;
    
    pybind11::class_< evoped::Hmat<float> >( m, "Hmat" )
        .def(pybind11::init<>())
        .def("clear", &evoped::Hmat<float>::clear)
        .def("make_matrix", pybind11::overload_cast<const std::string &, const std::string &, const std::string &, const std::string &, const std::string &, const std::string &>(&evoped::Hmat<float>::make_matrix))
        .def("make_matrix", pybind11::overload_cast<const std::string &, const std::string &, const std::string &, const std::string &>(&evoped::Hmat<float>::make_matrix))
        .def("make_matrix", pybind11::overload_cast<const std::string &, const std::string &, const std::string &>(&evoped::Hmat<float>::make_matrix))
        .def("save_matrix", pybind11::overload_cast<const std::string &>(&evoped::Hmat<float>::save_matrix))
        ;
}