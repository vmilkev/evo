#ifndef sparse_pcg_hpp__
#define sparse_pcg_hpp__

#include <fstream>
#include <thread>
#include <ctime>

#include "sparse_solver.hpp"

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#endif

namespace evolm
{
        class sparse_pcg : public sparse_solver
        {

        public:
                void solve();
                void solve2();

                void set_tolerance(double tol);
                void set_maxiter(size_t iter);

                sparse_pcg();
                ~sparse_pcg();

#ifdef PYBIND
                pybind11::array_t<float> get_solution();
#else
                std::vector<float> get_solution();
#endif

                void get_solution(const std::string &fname);

#ifdef UTEST
                std::vector<double> test_dval();
#endif

        private:
                void jacobi_pcg();
                void construct_dval(std::vector<double> &inverted_diagonal);
                void update_vect(std::vector<double> &out_vect, std::vector<double> &in_vect);
                double v_dot_v(std::vector<double> &v1, std::vector<double> &v2);

                void jacobi_pcg2();
                void construct_dval(std::vector<long double> &inverted_diagonal);
                void update_vect(std::vector<long double> &out_vect, std::vector<long double> &in_vect);
                long double v_dot_v(std::vector<long double> &v1, std::vector<long double> &v2);

                std::vector<double> sol;
                std::vector<long double> sol2;

                double tolerance = 1e-6; // default;
                size_t iterations;
                size_t max_iterations;
                bool befault_max_iter = true;                
        };

} // end of namespace evolm

#endif // sparse_pcg_hpp__