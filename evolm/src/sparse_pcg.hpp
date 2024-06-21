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
                sparse_pcg()
                {
                        tolerance = 1e-6; // default
                        befault_max_iter = true;
                        amatrix_ondisk = false;
                        amatrix_onmem = false;

                        pipeline_val = 4; // default

                        std::random_device rd;
                        srand(rd());

                        int iNum = std::rand() % 100000;

                        binFilename = "amatrix_" + std::to_string(iNum);
                }

                ~sparse_pcg()
                {
                        if (fA.is_open())
                                fA.close();
                }

                void solve();
                void solve(int pipeline_index);

                void set_tolerance(double tol);
                void set_maxiter(size_t iter);

#ifdef PYBIND
                pybind11::array_t<float> get_solution();
#else
                std::vector<float> get_solution();
#endif

                int get_solution(const std::string &fname);

#ifdef UTEST
                std::vector<double> test_dval();
#endif

        private:

                template <typename T>
                T v_dot_v(const matrix<T> &v1, const matrix<T> &v2);

                void jacobi_pcg();

                matrix<float> construct_dval(std::vector<std::vector<size_t>> &cov_offsets,
                                             size_t num_levels,
                                             std::vector<size_t> &ordered_levels);

                void construct_dval2(std::vector<double> &inverted_diagonal);
                void update_vect2( std::vector<double> &out_vect, std::vector<double> &in_vect );
                double v_dot_v2(std::vector<double> &v1, std::vector<double> &v2);

                void set_amatr(std::vector<std::vector<size_t>> &cov_offsets,
                               size_t num_levels,
                               std::vector<size_t> &ordered_levels);

                void set_amatr(std::vector<std::vector<size_t>> &cov_offsets,
                               size_t num_levels,
                               std::vector<size_t> &ordered_levels, bool on_mem);

                void update_vect(std::vector<std::vector<size_t>> &cov_offsets,
                                 size_t num_levels,
                                 std::vector<size_t> &ordered_levels,
                                 matrix<double> &out_vect,
                                 matrix<double> &in_vect);

                matrix<float> fget_vect(size_t all_tr_levels, size_t row);

                void set_amatrix();
                void memload_amatr();
                void set_pipeline(int which_pipeline);

                //matrix<double> sol;
                std::vector<double> sol;
                //matrix<float> amatr;
                double tolerance;
                size_t iterations;
                size_t max_iterations;
                bool befault_max_iter;
                bool amatrix_ondisk;
                bool amatrix_onmem;

                int pipeline_val;

                std::string binFilename; /* Name of binary file to store A on disck. */
                std::fstream fA;
        };

} // end of namespace evolm

#endif // sparse_pcg_hpp__
