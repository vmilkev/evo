#include "Population.hpp"
#include "Group.hpp"
#include "Trait.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#endif

PYBIND11_MODULE(evogen, m)
{
    pybind11::class_<evogen::Population>( m, "Population" )
        .def( pybind11::init<>() )
        .def( "size", &evogen::Population::size )
        .def( "capacity", &evogen::Population::capacity )
        .def( "clear", &evogen::Population::clear )
        .def( "reshape", &evogen::Population::reshape )
        .def( "aging", &evogen::Population::aging )        
        .def( "set_population", pybind11::overload_cast< size_t, const std::string &, float, int >(&evogen::Population::set_population) )
        .def( "set_population", pybind11::overload_cast< const std::string &, const std::string &, bool >(&evogen::Population::set_population) )
        .def( "id_at", pybind11::overload_cast< size_t, unsigned long >(&evogen::Population::id_at) )
        .def( "id_at", pybind11::overload_cast< size_t >(&evogen::Population::id_at) )
        .def( "sire_at", pybind11::overload_cast< size_t, unsigned long >(&evogen::Population::sire_at) )
        .def( "sire_at", pybind11::overload_cast< size_t >(&evogen::Population::sire_at) )
        .def( "dame_at", pybind11::overload_cast< size_t, unsigned long >(&evogen::Population::dame_at) )
        .def( "dame_at", pybind11::overload_cast< size_t >(&evogen::Population::dame_at) )
        .def( "age_at", pybind11::overload_cast< size_t, int >(&evogen::Population::age_at) )
        .def( "age_at", pybind11::overload_cast< size_t >(&evogen::Population::age_at) )
        .def( "alive_at", pybind11::overload_cast< size_t, bool >(&evogen::Population::alive_at) )
        .def( "alive_at", pybind11::overload_cast< size_t >(&evogen::Population::alive_at) )
        .def( "isgenotyped_at", pybind11::overload_cast< size_t, bool >(&evogen::Population::isgenotyped_at) )
        .def( "isgenotyped_at", pybind11::overload_cast< size_t >(&evogen::Population::isgenotyped_at) )
        .def( "sex_at", pybind11::overload_cast< size_t, int >(&evogen::Population::sex_at) )
        .def( "sex_at", pybind11::overload_cast< size_t >(&evogen::Population::sex_at) )
        .def( "phenotype_at", pybind11::overload_cast< size_t, pybind11::array_t<float> >(&evogen::Population::phenotype_at) )
        .def( "phenotype_at", pybind11::overload_cast< size_t >(&evogen::Population::phenotype_at) )
        .def( "breedingvalue_at", pybind11::overload_cast< size_t, pybind11::array_t<float> >(&evogen::Population::breedingvalue_at) )
        .def( "breedingvalue_at", pybind11::overload_cast< size_t >(&evogen::Population::breedingvalue_at) );

    pybind11::class_<evogen::Group>( m, "Group" )
        .def( pybind11::init<>() )
        .def( "size", &evogen::Group::size )
        .def( "size_at", &evogen::Group::size_at )
        .def( "clear", &evogen::Group::clear )
        .def( "remove", &evogen::Group::remove )
        .def( "move", &evogen::Group::move )
        .def( "add", pybind11::overload_cast< evogen::Population & >(&evogen::Group::add) )
        .def( "add", pybind11::overload_cast< evogen::Population &, size_t >(&evogen::Group::add) )
        .def( "add", pybind11::overload_cast< evogen::Group & >(&evogen::Group::add) )
        .def( "mate", pybind11::overload_cast<  >(&evogen::Group::mate) )
        .def( "mate", pybind11::overload_cast< bool, int, float >(&evogen::Group::mate) )
        .def( "regroup_newborn", &evogen::Group::regroup_newborn )
        .def( "aging", &evogen::Group::aging )
        .def( "genotype", &evogen::Group::genotype )
        .def( "kill", &evogen::Group::kill )
        .def( "make_observation", pybind11::overload_cast< evogen::Trait &, pybind11::array_t<float> >(&evogen::Group::make_observation) )
        .def( "make_observation", pybind11::overload_cast< evogen::Trait &, pybind11::array_t<float>, const std::string & >(&evogen::Group::make_observation) )
        .def( "make_observation", pybind11::overload_cast< evogen::Trait &, pybind11::array_t<float>, const std::string &, const std::string & >(&evogen::Group::make_observation) )
        .def( "make_observation", pybind11::overload_cast< evogen::Trait &, pybind11::array_t<float>, pybind11::array_t<float> >(&evogen::Group::make_observation) )
        .def( "make_observation", pybind11::overload_cast< evogen::Trait &, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<int> >(&evogen::Group::make_observation) );

    pybind11::class_<evogen::Trait>( m, "Trait" )
        .def( pybind11::init<>() )
        .def( pybind11::init< evogen::Population &, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>, size_t, pybind11::array_t<float> >() )
        //.def( "set_trait", &evogen::Trait::set_trait )
        .def( "set_trait", pybind11::overload_cast< evogen::Population &, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>, size_t, pybind11::array_t<float> >(&evogen::Trait::set_trait) )
        .def( "reset_trait", pybind11::overload_cast< evogen::Population &, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>, size_t, pybind11::array_t<float> >(&evogen::Trait::reset_trait) )
        //.def( "reset_trait", &evogen::Trait::reset_trait )
        .def( "clear", &evogen::Trait::clear )
        .def( "is_cleared", &evogen::Trait::is_cleared )
        .def( "get_observations", pybind11::overload_cast< evogen::Population &, pybind11::array_t<float> >(&evogen::Trait::get_observations) )
        .def( "get_observations", pybind11::overload_cast< evogen::Population &, pybind11::array_t<float>, const std::string & >(&evogen::Trait::get_observations) )
        .def( "get_observations", pybind11::overload_cast< evogen::Population &, pybind11::array_t<float>, const std::string &, const std::string & >(&evogen::Trait::get_observations) )
        .def( "get_observations", pybind11::overload_cast< evogen::Population &, pybind11::array_t<float>, pybind11::array_t<float> >(&evogen::Trait::get_observations) )
        .def( "get_observations", pybind11::overload_cast< evogen::Population &, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<int> >(&evogen::Trait::get_observations) );

}
