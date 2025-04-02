#include "sparse_pcg.hpp"

namespace evolm
{
    sparse_pcg::~sparse_pcg()
    {
        clear_model_matrix();
    }

    void sparse_pcg::solve()
    {
        try
        {
            log_message("Loading data ...");

            auto start = std::chrono::high_resolution_clock::now();

            memload_effects();
            memload_cor();
            memload_cor_effects();

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

            log_message("completed, elapsed time:", duration.count());
            log_message("Constructing RHS vector ...");

            start = std::chrono::high_resolution_clock::now();
            
            construct_rhs();
            
            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

            log_message("completed, elapsed time:", duration.count());

            double size_of_maps = r_map.size() * (sizeof(size_t) + r_map.begin()->second.size() * sizeof(float)) * 3.0;
            double size_of_data = (size_of_maps + (double)model->get_size_of_data()) / 1073741824.0;

            set_data_size(size_of_data);

            log_message("Constructing model matrix ...");

            start = std::chrono::high_resolution_clock::now();

            set_model_matrix();

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

            log_message("completed, elapsed time:", duration.count());
            log_message("Cleaning the model data ...");

            start = std::chrono::high_resolution_clock::now();

            remove_model();

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

            log_message("completed, elapsed time:", duration.count());
            log_message("Loading model matrix into memory ...");
            
            start = std::chrono::high_resolution_clock::now();

            if (first_row_ondisk < model_matrix.size() - 1) // after cleaning data, load as much as possible model_matrix to memory
            {
                size_t rows_per_memory = std::floor( ( ( get_memory_limit() - 2.0 * size_of_data) / row_size_upper_bound ) - 2.0 );
                
                size_t last_row = rows_per_memory - 1;
                if (last_row >= model_matrix.size())
                    last_row = model_matrix.size();

                std::fstream fA;
                fA.open(bin_fname, fA.binary | fA.in);
                fA.seekg(model_matrix[first_row_ondisk].bin_file_read_position, fA.beg);

                for (size_t i = first_row_ondisk; i < last_row; i++)
                    model_matrix[i].fread(fA);

                fA.close();
                if ( rows_per_memory < first_row_ondisk )
                    throw std::string("rows_per_memory < first_row_ondisk");
                
                first_row_ondisk = last_row;
            }

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

            log_message("completed, elapsed time:", duration.count());
            
            log_message(" ");
            log_message("Starting PCG iterations (stops when rate values lower or equal to 0.0) ...");
            log_message(" ");

            jacobi_pcg();

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
            
            log_message(" ");
            log_message("completed PCG iterations, elapsed time:", duration.count());
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::solve()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_pcg::solve()." << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::solve()." << '\n';
            throw;
        }
    }

#ifdef PYBIND

    pybind11::array_t<float> sparse_pcg::get_solution()
    {
        try
        {
            auto solution = pybind11::array_t<float>(sol.size());

            pybind11::buffer_info buf = solution.request();

            if (buf.ndim != 1)
                throw std::runtime_error("Number of dimensions must be one");

            float *ptr = static_cast<float *>(buf.ptr);

            for (pybind11::ssize_t i = 0; i < buf.shape[0]; i++)
            {
                ptr[i] = sol[i];
            }

            return solution;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_pcg::(): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::get_solution()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::get_solution()." << '\n';
            throw;
        }
    }

#else

    std::vector<float> sparse_pcg::get_solution()
    {
        try
        {
            std::vector<float> solution;

            for (size_t i = 0; i < sol.size(); i++)
                solution.push_back(sol[i]);

            return solution;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::get_solution()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::get_solution()." << '\n';
            throw;
        }
    }

#endif

    void sparse_pcg::get_solution(const std::string &fname)
    {
        try
        {
            std::ofstream solution(fname);

            if (solution.is_open())
            {
                for (size_t i = 0; i < sol.size(); i++)
                    solution << std::setprecision(16) << sol[i] << "\n";

                solution.close();
            }
            else
                throw "Unable to open solution file for output!";
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_pcg::get_solution(const std::string &): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::get_solution(const std::string &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::get_solution(const std::string &)." << '\n';
            throw;
        }
    }

    void sparse_pcg::get_solution(std::vector<double> &out_sol)
    {
        try
        {
            out_sol = sol;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_pcg::get_solution(std::vector<double> &): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::get_solution(std::vector<double> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::get_solution(std::vector<double> &)." << '\n';
            throw;
        }
    }

    void sparse_pcg::construct_dval(std::vector<double> &inverted_diagonal)
    {
        try
        {
            if (inverted_diagonal.size() != model_matrix.size())
                throw std::string("The size allocated for the vector of model matrix diagonals is not correct!");

#pragma omp parallel for num_threads(available_cpu)
            for (size_t j = 0; j < first_row_ondisk; j++) // process data which is already in memory
            {
                float d = model_matrix[j].value_at(0, j);

                if (d == 0.0f)
                    inverted_diagonal[j] = 1.0;
                else
                    inverted_diagonal[j] = 1.0 / static_cast<double>(d);
            }

            /*for (size_t i = 0; i < get_num_of_mem_blocks(); i++) // variant 1, on blocks with consecutive reading
            {
                size_t first_row = 0;
                size_t last_row = 0;

                load_model_matrix(i);
                get_mem_block_range(i, first_row, last_row);

#pragma omp parallel for num_threads( available_cpu )
                for (size_t j = first_row; j <= last_row; j++) // process in parallel
                {
                    float d = model_matrix[j].value_at(0,j);

                    if (d == 0.0f)
                        throw std::string("sparse_pcg::construct_dval2(std::vector<float> &) => The diagonal element of the model matrix is 0.0!");

                    inverted_diagonal[j] = 1.0 / static_cast<double>(d);
                }
                unload_model_matrix(i);
            }*/
            //std::vector<bool> reading_complete(model_matrix.size(), false);
            std::fstream fA;
            fA.open(bin_fname, fA.binary | fA.in);
            fA.seekg(model_matrix[first_row_ondisk].bin_file_read_position, fA.beg);

/*
#pragma omp parallel sections
            {
#pragma omp section
                {
                    for (size_t j = first_row_ondisk; j < model_matrix.size(); j++)
                    {
                        model_matrix[j].fread(fA);
                        reading_complete[j] = true;
                    }
                }
#pragma omp section
                {
                    size_t j = first_row_ondisk;
                    while (j < model_matrix.size())
                    {
                        if (reading_complete[j] == true)
                        {
                            float d = model_matrix[j].value_at(0, j);
                            if (d == 0.0f)
                                inverted_diagonal[j] = 1.0;
                            else
                                inverted_diagonal[j] = 1.0 / static_cast<double>(d);
                            model_matrix[j].remove_data();
                            j = j + 1;
                        }
                    }
                }
            }
            fA.close();*/

            for (size_t j = first_row_ondisk; j < model_matrix.size(); j++)
            {
                model_matrix[j].fread(fA);
                float d = model_matrix[j].value_at(0, j);
                if (d == 0.0f)
                    inverted_diagonal[j] = 1.0;
                else
                    inverted_diagonal[j] = 1.0 / static_cast<double>(d);
                model_matrix[j].remove_data();
            }

            fA.close();
        }
        catch (const std::string &err)
        {
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::construct_dval(std::vector<float> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::construct_dval(std::vector<float> &)." << '\n';
            throw;
        }
    }

    void sparse_pcg::update_vect(std::vector<double> &out_vect, std::vector<double> &in_vect)
    {
        try
        {
#pragma omp parallel for num_threads(available_cpu)
            for (size_t j = 0; j < first_row_ondisk; j++) // process data which is already in memory
            {
                double result = 0.0;
                model_matrix[j].vect_dot_vect(in_vect, result);
                out_vect[j] = result;
            }

            /*for (size_t i = 0; i < get_num_of_mem_blocks(); i++) // variant 1, on blocks with consecutive reading
            {
                size_t first_row = 0;
                size_t last_row = 0;

                load_model_matrix(i);
                get_mem_block_range(i, first_row, last_row);

#pragma omp parallel for num_threads( available_cpu )
                for (size_t j = first_row; j <= last_row; j++)
                {
                    double result = 0.0;
                    model_matrix[j].vect_dot_vect(in_vect, result);
                    out_vect[j] = result;
                }
                unload_model_matrix(i);
            }*/

            //std::vector<bool> reading_complete(model_matrix.size(), false);
            std::fstream fA;
            fA.open(bin_fname, fA.binary | fA.in);
            fA.seekg(model_matrix[first_row_ondisk].bin_file_read_position, fA.beg);
            /*
#pragma omp parallel sections
            {
#pragma omp section
                {
                    for (size_t j = first_row_ondisk; j < model_matrix.size(); j++)
                    {                        
                        model_matrix[j].fread(fA);
                        reading_complete[j] = true;
                    }
                }
#pragma omp section
                {
                    size_t j = first_row_ondisk;
                    while (j < model_matrix.size())
                    {
                        if (reading_complete[j] == true)
                        {
                            double result = 0.0;
                            model_matrix[j].vect_dot_vect(in_vect, result);
                            out_vect[j] = result;
                            model_matrix[j].remove_data();
                            j = j + 1;
                        }
                    }
                }
            }*/

            for (size_t j = first_row_ondisk; j < model_matrix.size(); j++)
            {                        
                model_matrix[j].fread(fA);
                double result = 0.0;
                model_matrix[j].vect_dot_vect(in_vect, result);
                out_vect[j] = result;
                model_matrix[j].remove_data();
            }

            fA.close();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::update_vect2(std::vector<double> &, std::vector<double> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::update_vect2(std::vector<double> &, std::vector<double> &)." << '\n';
            throw;
        }
    }

    void sparse_pcg::jacobi_pcg()
    {
        try
        {
            matrix<size_t> shape_rhs = rhs.shape();

            if (shape_rhs[0] == 0)
                throw std::string("The size of RHS is 0!");

            size_t unknowns = shape_rhs[0];

            if (befault_max_iter)
                max_iterations = unknowns;

            std::vector<double> Mi(unknowns, 0.0);     // inverse of diagonal of A
            std::vector<double> _rhs(unknowns, 0.0);   // rhs
            std::vector<double> tVect(unknowns, 0.0);  // vector to keep the result of operation: A*x
            std::vector<double> r_vect(unknowns, 0.0); // residual vector: r = rhs - A*x
            std::vector<double> d(unknowns, 0.0);
            std::vector<double> q(unknowns, 0.0);
            std::vector<double> s(unknowns, 0.0);

            sol.resize(unknowns, 0.0); // solution vector

            log_message("Getting the inverce of diagonal elements of the coefficient matrix ...");

            auto start = std::chrono::high_resolution_clock::now();

            construct_dval(Mi);

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

            log_message("completed, elapsed time:", duration.count());

            for (size_t i = 0; i < unknowns; i++)
                _rhs[i] = static_cast<double>(rhs[i]);

            for (size_t i = 0; i < unknowns; i++)
                sol[i] = _rhs[i] * Mi[i]; // initial solution

            log_message("Expected duration of one pcg iteration");

            start = std::chrono::high_resolution_clock::now();

            update_vect(tVect, sol); // A*x --> tVect

            for (size_t i = 0; i < unknowns; i++)
                r_vect[i] = _rhs[i] - tVect[i]; // residual = rhs - A*x

            tVect.clear();

            for (size_t i = 0; i < unknowns; i++)
                d[i] = Mi[i] * r_vect[i];

            double delta_new = v_dot_v(r_vect, d);

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

            log_message("(upper bound) ", (size_t)( 1.3 * (double)duration.count() ));
            log_message("Expected overall solution time");
            log_message("using the maximum number of iterations (the upper bound) ", (size_t)( 1.3 * (double)max_iterations * (double)duration.count() ));
            
            std::cout<<"stop when error < "<< delta_new * tolerance * tolerance<<'\n';
            std::cout<<"delta_zero = "<< delta_new <<'\n';
            std::cout<<"tolerance =  "<< tolerance<<'\n';

            double delta_zero = delta_new;

            iterations = 1;

            double i_value = 0.0;
            double alpha = 0.0;
            double delta_old = 0.0;
            double betha = 0.0;

            log_message("Iteration & Error:");
            log_message(" ");

            while ( iterations < max_iterations && delta_new > delta_zero * tolerance * tolerance )
            {
                update_vect(q, d); // here we casting vector from float to double every time by calling this function

                i_value = v_dot_v(d, q);

                if (i_value == 0.0)
                    throw std::string("sparse_pcg::jacobi_pcg() => Expected division by 0.0: v_dot_v(d, q) == 0.0!");

                alpha = delta_new / i_value;

                for (size_t i = 0; i < sol.size(); i++)
                    sol[i] = sol[i] + alpha * d[i];

                if (!(iterations % 50))
                {
                    update_vect(tVect, sol);

                    for (size_t i = 0; i < unknowns; i++)
                        r_vect[i] = _rhs[i] - tVect[i]; // r = b - A*x
                }
                else
                {
                    for (size_t i = 0; i < r_vect.size(); i++)
                        r_vect[i] = r_vect[i] - alpha * q[i];
                }

                for (size_t i = 0; i < unknowns; i++)
                    s[i] = Mi[i] * r_vect[i];

                delta_old = delta_new;

                delta_new = v_dot_v(r_vect, s);

                if (delta_old == 0.0)
                    throw std::string("sparse_pcg::jacobi_pcg() => Expected division by 0.0: delta_old == 0.0!");

                betha = delta_new / delta_old;

                for (size_t i = 0; i < unknowns; i++)
                    d[i] = s[i] + betha * d[i];

                //if ( !( iterations % logging_rate ) )
                    log_message(iterations, (delta_new/(delta_zero * tolerance * tolerance))-1.0);

                iterations = iterations + 1;
            }
            iterations = iterations - 1;

            std::cout<<"delta_new = "<<delta_new<<'\n';
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_pcg::jacobi_pcg(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::jacobi_pcg(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::jacobi_pcg(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            throw;
        }
    }

    double sparse_pcg::v_dot_v(std::vector<double> &v1, std::vector<double> &v2)
    {
        try
        {
            double res = 0.0;

            for (size_t i = 0; i < v1.size(); i++)
                res = res + v1[i] * v2[i];

            return res;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::v_dot_v(std::vector<double> &, std::vector<double> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::v_dot_v(std::vector<double> &, std::vector<double> &)." << '\n';
            throw;
        }
    }

    double sparse_pcg::v_dot_v2(std::vector<double> &v1, std::vector<double> &v2)
    {
        try
        {
            // it seems like this approach more acurate?
            // At least needs less iterations to converge!

            double *v11 = new double[v1.size()];
            double *v22 = new double[v2.size()];

            for (size_t i = 0; i < v1.size(); i++)
            {
                v11[i] = v1[i];
                v22[i] = v2[i];
            }

            double res = 0.0;
            for (size_t i = 0; i < v1.size(); i++)
            {
                res = res + v11[i] * v22[i];
            }

            delete[] v11;
            delete[] v22;

            return res;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::v_dot_v2(std::vector<double> &, std::vector<double> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::v_dot_v2(std::vector<double> &, std::vector<double> &)." << '\n';
            throw;
        }
    }

    void sparse_pcg::set_tolerance(double tol)
    {
        try
        {
            tolerance = tol;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::set_tolerance(double)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::set_tolerance(double)." << '\n';
            throw;
        }
    }

    void sparse_pcg::set_maxiter(size_t iter)
    {
        try
        {
            max_iterations = iter;
            befault_max_iter = false;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::set_maxiter(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::set_maxiter(size_t)." << '\n';
            throw;
        }
    }

#ifdef UTEST

    std::vector<double> sparse_pcg::test_dval()
    {
        try
        {
            memload_cor_effects();
            memload_cor();

            set_model_matrix();

            construct_rhs();

            std::vector<double> out(model_matrix.size(), 0.0);

            construct_dval(out);

            for (size_t i = 0; i < out.size(); i++)
                out[i] = 1.0 / out[i]; // because out is the inverse of diagonal values

            return out;
        }
        catch (const std::string &err)
        {
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::test_dval()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::test_dval()." << '\n';
            throw;
        }
    }
#endif
}
