#ifndef sparse_pcg_hpp__
#define sparse_pcg_hpp__

#include <fstream>
#include <thread>
#include <ctime>
#include <iomanip>

#include "sparse_solver.hpp"
//#include <omp.h>

#ifdef PYBIND
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#endif

namespace evolm
{
        class sparse_pcg : public sparse_solver
        {
                using sparse_solver::sparse_solver;

        public:
                void solve();

                void set_tolerance(double tol);
                void set_maxiter(size_t iter);

                sparse_pcg() = default;
                ~sparse_pcg();

#ifdef PYBIND
                pybind11::array_t<float> get_solution();
#else
                std::vector<float> get_solution();
#endif

                void get_solution(const std::string &fname);
                void get_solution(std::vector<double> &out_sol);

#ifdef UTEST
                std::vector<double> test_dval();
#endif

        private:
                void jacobi_pcg();
                void construct_dval(std::vector<double> &inverted_diagonal);
                void update_vect(std::vector<double> &out_vect, std::vector<double> &in_vect);
                double v_dot_v(std::vector<double> &v1, std::vector<double> &v2);
                double v_dot_v2(std::vector<double> &v1, std::vector<double> &v2);

                std::vector<double> sol;

                double tolerance = 1e-6; // default;
                size_t iterations;
                size_t max_iterations;
                bool befault_max_iter = true;                
        };

} // end of namespace evolm

#endif // sparse_pcg_hpp__
