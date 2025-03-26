#include "solver_pcg.hpp"

namespace evolm
{
    void Pcg::solve()
    {
        try
        {
            set_pipeline(pipeline_val); // defines a memory management

            construct_rhs(); // constructing RHS

            size_t n_all_levels = num_all_levels(); // number of all random effects in the model

            std::vector<size_t> ordered_random_levels = get_ordered_levels(); // for each trait stack the number of levels for each effect matrix, in the order as it appiars in the model

            std::vector<std::vector<size_t>> rcov_offsets = get_cov_offsets(ordered_random_levels);

            jacobi_pcg(rcov_offsets, n_all_levels, ordered_random_levels);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::solve()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Pcg::solve()." << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::solve()." << '\n';
            throw;
        }
    }

    void Pcg::solve(int pipeline_index)
    {
        try
        {
            pipeline_val = pipeline_index;

            set_pipeline(pipeline_val);

            construct_rhs();

            size_t n_all_levels = num_all_levels(); // number of all random effects in the model

            std::vector<size_t> ordered_random_levels = get_ordered_levels();

            std::vector<std::vector<size_t>> rcov_offsets = get_cov_offsets(ordered_random_levels);

            jacobi_pcg(rcov_offsets, n_all_levels, ordered_random_levels);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::solve(int)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::solve(int)." << '\n';
            throw;
        }
    }

    void Pcg::set_pipeline(int which_pipeline)
    {
        try
        {
            switch (which_pipeline)
            {
            case 1:
                /* everything on disk */

                break;

            case 2:
                /* only data on memory */

                memload_effects();
                memload_var();
                memload_cor();

                break;

            case 3:
                /* everything on memory */

                memload_effects();
                memload_var();
                memload_cor();
                memload_amatr();

                diskload_effects();
                diskload_var();
                diskload_cor();

                break;

            case 4:
                /* building coefficient matrix on disk */

                memload_effects();
                memload_var();
                memload_cor();
                set_amatrix();

                diskload_effects();
                diskload_var();
                diskload_cor();

                break;

            default:
                throw std::runtime_error("Wrong set-up of computational pipelne!");
                break;
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Pcg::set_pipeline(int): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::set_pipeline(int). " << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::set_pipeline(int)." << '\n';
            throw;
        }
    }

#ifdef PYBIND

    pybind11::array_t<float> Pcg::get_solution()
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
            std::cerr << "Exception in Pcg::(): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::get_solution()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::get_solution()." << '\n';
            throw;
        }
    }

#else

    std::vector<float> Pcg::get_solution()
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
            std::cerr << "Exception in Pcg::get_solution()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::get_solution()." << '\n';
            throw;
        }
    }

#endif

    int Pcg::get_solution(const std::string &fname)
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
            {
                throw "Unable to open solution file for output!";
            }

            return 0;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Pcg::get_solution(const std::string &): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::get_solution(const std::string &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::get_solution(const std::string &)." << '\n';
            throw;
        }
    }

    matrix<double> Pcg::construct_dval(std::vector<std::vector<size_t>> &cov_offsets, size_t num_levels, std::vector<size_t> &ordered_levels)
    {
        try
        {
            matrix<size_t> shape_rhs = rhs.shape();

            if (shape_rhs[0] == 0)
                throw std::string("The size of RHS is 0!");

            matrix<double> dval(shape_rhs[0], 1);

            size_t cmatr_row = 0;

            for (size_t i_trate = 0; i_trate < n_trait; i_trate++)
            {
                size_t tr_levels = get_levels(i_trate);

                for (size_t i_eff = 0; i_eff < tr_levels; i_eff++)
                {
                    matrix<float> vect;

                    if (amatrix_ondisk && !amatrix_onmem)
                    {
                        vect = fget_vect(shape_rhs[0], cmatr_row);
                    }
                    else if (amatrix_onmem)
                    {
                        size_t rows[] = {cmatr_row, cmatr_row};
                        size_t cols[] = {0, shape_rhs[0] - 1};

                        amatr.cast_fget(rows, cols, vect);
                    }
                    else if (!amatrix_ondisk && !amatrix_onmem)
                    {
                        vect = get_row_cmatr(shape_rhs[0], i_trate, i_eff, cov_offsets, num_levels, ordered_levels, cmatr_row);
                    }
                    else
                        throw std::string("Wrong computational pipeline. In construct_dval(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)");

                    if (vect(0, cmatr_row) == 0.0f)
                        throw std::string("Pcg::construct_dval(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &) => Expected division by 0.0!");

                    dval(cmatr_row, 0) = 1.0 / static_cast<double>(vect(0, cmatr_row));

                    cmatr_row = cmatr_row + 1;
                }
            }

            return dval;
        }
        catch (const std::string &err)
        {
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::construct_dval(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::construct_dval(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            throw;
        }
    }

    void Pcg::memload_amatr()
    {
        try
        {
            if (amatrix_ondisk && !amatrix_onmem)
            {
                size_t all_tr_levels = get_all_levels();

                if (all_tr_levels == 0)
                    throw std::string("The size of Unknowns is 0! In memload_amatr().");

                amatr.resize(all_tr_levels, all_tr_levels);

                std::ifstream f(binFilename.c_str());

                if (f.good())
                {
                    fA.open(binFilename, fA.binary | fA.in);

                    if (!fA.is_open())
                        throw std::string("Error while opening a binary file in Pcg::memload_amatr()");

                    fA.seekg(0 * all_tr_levels * sizeof(float));

                    float *A;
#ifdef intelmkl
                    A = (float *)mkl_malloc(all_tr_levels * all_tr_levels * sizeof(float), sizeof(float) * 8);
#else
                    //A = (float *)aligned_alloc(sizeof(float) * 8, all_tr_levels * all_tr_levels * sizeof(float));
                    A = (float *)malloc(all_tr_levels * all_tr_levels * sizeof(float));
#endif
                    if (A == NULL)
                    {
#ifdef intelmkl
                        mkl_free(A);
#else
                        free(A);
#endif
                        throw std::string("Memory allocation error in memload_amatr()");
                    }

                    fA.read((char *)A, all_tr_levels * all_tr_levels * sizeof(float));

                    amatr.insert_array(A);

                    amatrix_onmem = true;

#ifdef intelmkl
                    mkl_free(A);
#else
                    free(A);
#endif

                    fA.close();
                }
                else
                    throw std::string("Error while getting a binary file name in Pcg::memload_amatr()");
            }
            else if (!amatrix_ondisk && !amatrix_onmem)
            {
                // build amatr directly on memory

                size_t n_all_levels = num_all_levels();

                size_t all_levels = get_all_levels();

                std::vector<size_t> ordered_random_levels = get_ordered_levels();

                std::vector<std::vector<size_t>> rcov_offsets = get_cov_offsets(ordered_random_levels);

                amatr.resize(all_levels, all_levels);

//auto start = std::chrono::high_resolution_clock::now();

                set_amatr(rcov_offsets, n_all_levels, ordered_random_levels, true);

//auto stop = std::chrono::high_resolution_clock::now();
//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
//std::cout <<"set_amatr() (milliseconds): "<< duration.count() << std::endl;

                amatrix_onmem = true;
            }
            else
            {
                throw std::string("The matrix is already loaded to the memory in memload_amatr()");
            }
        }
        catch (const std::string &err)
        {
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::memload_amatr()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::memload_amatr()." << '\n';
            throw;
        }
    }

    void Pcg::set_amatrix()
    {
        try
        {
            size_t n_all_levels = num_all_levels();

            std::vector<size_t> ordered_random_levels = get_ordered_levels();

            std::vector<std::vector<size_t>> rcov_offsets = get_cov_offsets(ordered_random_levels);

            set_amatr(rcov_offsets, n_all_levels, ordered_random_levels);

            amatrix_ondisk = true;
        }
        catch (const std::string &err)
        {
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::set_amatrix()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::set_amatrix()." << '\n';
            throw;
        }
    }

    void Pcg::set_amatr(std::vector<std::vector<size_t>> &cov_offsets, size_t num_levels, std::vector<size_t> &ordered_levels)
    {
        // Build coefficient matrix on disk

        // Consider trait 1 effects: x1, z1;
        //          trait 2 effects: z2;
        // Model effects stacked consecutively as [ x1 z1 z2 ]';
        // Model observations stacked as [ y1 y2 ];
        // Than coefficient matrix calculated as (row-by-row):
        //                       [ x1 ]                   [ x1' * x1 * R11  x1' * z1 * R11  x1' * z2 * R12 ]
        //                       | z1 | * [ x1 z1 z2 ] =  | z1' * x1 * R11  z1' * z1 * R11  z1' * z2 * R12 |
        //                       [ z2 ]                   [ z2' * x1 * R21  z2' * z1 * R21  z2' * z2 * R22 ]

        try
        {
            size_t all_tr_levels = get_all_levels();

            if (all_tr_levels == 0)
                throw std::string("The size of Unknowns is 0!");

            fA.open(binFilename, fA.binary | fA.app | fA.out);

            if (!fA.is_open())
                throw std::string("Error while opening a binary file. Pcg::set_amatr()");

            size_t cmatr_row = 0;

            for (size_t i_trate = 0; i_trate < n_trait; i_trate++)
            {
                size_t tr_levels = get_levels(i_trate);

                for (size_t i_eff = 0; i_eff < tr_levels; i_eff++)
                {
                    matrix<float> vect = get_row_cmatr(all_tr_levels, i_trate, i_eff, cov_offsets, num_levels, ordered_levels, cmatr_row);

                    float *b;

                    b = vect.return_array();

                    fA.write((char *)b, all_tr_levels * sizeof(float));

                    cmatr_row = cmatr_row + 1;
                }
            }

            fA.close();
        }
        catch (const std::string &err)
        {
            std::cerr << err << "  In in Pcg::set_amatr(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::set_amatr(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::set_amatr(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            throw;
        }
    }

    void Pcg::set_amatr(std::vector<std::vector<size_t>> &cov_offsets, size_t num_levels, std::vector<size_t> &ordered_levels, bool on_mem)
    {
        // Build coefficient matrix on memory

        try
        {
            if (on_mem)
            {
                const auto processor_count = std::thread::hardware_concurrency();

                int n_threads = 1;

                if (var_onmem && cor_onmem && z_on_memory)
                {
                    n_threads = processor_count;
                }

                //std::cout << "n_threads: " << n_threads << "\n";

                for (size_t i_trate = 0; i_trate < n_trait; i_trate++)
                {
#pragma omp parallel for num_threads(n_threads)
                    for (size_t i_eff = 0; i_eff < get_levels(i_trate); i_eff++)
                    {
                        size_t private_raw = get_all_levels(i_trate) + i_eff;

                        size_t all_tr_levels = get_all_levels();

                        matrix<float> vect = get_row_cmatr(all_tr_levels, i_trate, i_eff, cov_offsets, num_levels, ordered_levels, private_raw);
                        for (size_t i = 0; i < all_tr_levels; i++)
                        {
                            amatr(private_raw, i) = vect(0, i);
                        }
                    }
                }
            }
        }
        catch (const std::string &err)
        {
            std::cerr << err << "  In in Pcg::set_amatr(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &, bool)." << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::set_amatr(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &, bool)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::set_amatr(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &, bool)." << '\n';
            throw;
        }
    }

    matrix<float> Pcg::fget_vect(size_t all_tr_levels, size_t row)
    {
        try
        {
            matrix<float> vect(1, all_tr_levels);

            std::ifstream f(binFilename.c_str());

            if (f.good())
            {
                fA.open(binFilename, fA.binary | fA.in);

                if (!fA.is_open())
                    throw std::string("Error while opening a binary file in Pcg::fget_vect(size_t row)");

                fA.seekg(row * all_tr_levels * sizeof(float));

                float *A;
#ifdef intelmkl
                A = (float *)mkl_malloc(1 * all_tr_levels * sizeof(float), sizeof(float) * 8);
#else
                //A = (float *)aligned_alloc(sizeof(float) * 8, 1 * all_tr_levels * sizeof(float));
                A = (float *)malloc(1 * all_tr_levels * sizeof(float));
#endif
                if (A == NULL)
                {
#ifdef intelmkl
                    mkl_free(A);
#else
                    free(A);
#endif
                    throw std::string("Memory allocation errorin fget_vect(size_t all_tr_levels, size_t row)");
                }

                fA.read((char *)A, all_tr_levels * sizeof(float));

                vect.insert_array(A);

#ifdef intelmkl
                mkl_free(A);
#else
                free(A);
#endif

                fA.close();
            }
            else
                throw std::string("Error while getting a binary file name in Pcg::fget_vect(size_t row)");

            return vect;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Pcg::fget_vect(size_t, size_t): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::fget_vect(size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::fget_vect(size_t, size_t)." << '\n';
            throw;
        }
    }

    void Pcg::update_vect(std::vector<std::vector<size_t>> &cov_offsets,
                          size_t num_levels,
                          std::vector<size_t> &ordered_levels,
                          matrix<double> &out_vect,
                          matrix<double> &in_vect)
    {
        try
        {
            matrix<size_t> shape = out_vect.shape();

            size_t unknowns = shape[0];

            size_t cmatr_row = 0;

            for (size_t i_trate = 0; i_trate < n_trait; i_trate++)
            {
                size_t tr_levels = get_levels(i_trate);

                for (size_t i_eff = 0; i_eff < tr_levels; i_eff++)
                {
                    matrix<float> vect;

                    if (amatrix_ondisk && !amatrix_onmem)
                    {
                        vect = fget_vect(unknowns, cmatr_row);
                    }
                    else if (amatrix_onmem)
                    {
                        size_t rows[] = {cmatr_row, cmatr_row};
                        size_t cols[] = {0, unknowns - 1};

                        amatr.cast_fget(rows, cols, vect);
                    }
                    else if (!amatrix_ondisk && !amatrix_onmem)
                    {
                        vect = get_row_cmatr(unknowns, i_trate, i_eff, cov_offsets, num_levels, ordered_levels, cmatr_row);
                    }

                    matrix<double> vect2 = vect._double(); // cast to double; may take time, so needs to be checked!

                    vect.clear();

                    matrix<double> res = vect2 * in_vect;
                    
                    out_vect(cmatr_row, 0) = res[0];
                    
                    cmatr_row = cmatr_row + 1;
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::update_vect(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &, matrix<double> &, matrix<double> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::update_vect(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &, matrix<double> &, matrix<double> &)." << '\n';
            throw;
        }
    }

    void Pcg::jacobi_pcg(std::vector<std::vector<size_t>> &cov_offsets, size_t num_levels, std::vector<size_t> &ordered_levels)
    {
        try
        {
            matrix<size_t> shape_rhs = rhs.shape();

            if (shape_rhs[0] == 0)
                throw std::string("The size of RHS is 0!");

            size_t unknowns = shape_rhs[0];

            if (befault_max_iter)
                max_iterations = unknowns * 10;

            sol.resize(unknowns, 1);

            matrix<double> Mi = construct_dval(cov_offsets, num_levels, ordered_levels); // actually, returning the inverse of dval (== Mi)

            //matrix<double> Mi = _Mi._double(); // cast to double
            matrix<double> _rhs = rhs._double(); // cast to double

            for (size_t i = 0; i < unknowns; i++)
                sol[i] = _rhs[i] * Mi[i]; // initial solution

            // -------------------------------------------------
            matrix<double> tVect(unknowns, 1); // vector to keep the result of operation: A*x

//auto start = std::chrono::high_resolution_clock::now();
            update_vect(cov_offsets, num_levels, ordered_levels, tVect, sol); // A*x(==sol)
//auto stop = std::chrono::high_resolution_clock::now();
//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
//std::cout <<"update_vect() (milliseconds): "<< duration.count() << std::endl;

            matrix<double> r_vect = _rhs - tVect; // r = b - A*x

            tVect.clear();
            // -------------------------------------------------

            matrix<double> d(unknowns, 1);

            for (size_t i = 0; i < unknowns; i++)
            {
                d[i] = Mi[i] * r_vect[i];
            }

//start = std::chrono::high_resolution_clock::now();
            double delta_new = v_dot_v(r_vect, d);
//stop = std::chrono::high_resolution_clock::now();
//duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
//std::cout <<"v_dot_v() (milliseconds): "<< duration.count() << std::endl;

            //double delta_zero = delta_new;

            iterations = 1;

//std::cout << "max_iterations: "<< max_iterations << "; iter: " << iterations << "; delta_new: " << delta_new << " condition: "<< /*delta_zero */ tolerance * tolerance << "\n";
            while (iterations < max_iterations && delta_new > /*delta_zero */ tolerance * tolerance)
            {
                matrix<double> q(unknowns, 1); // !!! can be declared outside the loop ?

                update_vect(cov_offsets, num_levels, ordered_levels, q, d); // here we casting vector from float to double every time by calling this function

                double i_value = v_dot_v(d, q);

                if (i_value == 0.0)
                    throw std::string("Pcg::jacobi_pcg() => Expected division by 0.0: v_dot_v(d, q) == 0.0!");

                double alpha = delta_new / i_value;

                sol = sol + alpha * d;

                if (!(iterations % 50))
                {
                    matrix<double> _tVect(unknowns, 1);

                    update_vect(cov_offsets, num_levels, ordered_levels, _tVect, sol);

                    r_vect = _rhs - _tVect;
                    _tVect.clear();
                }
                else
                {
                    r_vect = r_vect - alpha * q;
                }

                matrix<double> s(unknowns, 1); // !!! should be declared outside the loop ?

                for (size_t i = 0; i < unknowns; i++)
                    s[i] = Mi[i] * r_vect[i];

                double delta_old = delta_new;

                delta_new = v_dot_v(r_vect, s);

                if (delta_old == 0.0)
                    throw std::string("Pcg::jacobi_pcg() => Expected division by 0.0: delta_old == 0.0!");

                double betha = delta_new / delta_old;

                d = s + betha * d;

                // debugging
                //if (!(iterations % 10))
                //    std::cout << "iter: " << iterations << "; delta_new: " << delta_new <<" alpha: "<<alpha<<" i_value: "<<i_value<< "\n";

                iterations = iterations + 1;
            }
            iterations = iterations - 1;
            //std::cout << "no iterations: " << iterations << " delta_new: " << delta_new << " condition: "<< /*delta_zero */ tolerance * tolerance  <<"\n";
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Pcg::jacobi_pcg(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::jacobi_pcg(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::jacobi_pcg(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            throw;
        }
    }

    matrix<float>
    Pcg::get_row_cmatr(size_t rhs_size,
                        size_t i_trate,
                        size_t i_eff,
                        std::vector<std::vector<size_t>> &cov_offsets,
                        size_t num_levels,
                        std::vector<size_t> &ordered_levels,
                        size_t i_row)
    {
        try
        {
            matrix<float> vect_a(1, rhs_size);

            size_t correlations = corr_size();

            size_t last_random_level = 0;

            for (size_t j_trate = 0; j_trate < n_trait; j_trate++)
            {
                size_t tr_levels = get_levels(j_trate);

                size_t first_random_level = last_random_level + 1;
                last_random_level = last_random_level + tr_levels;

                size_t which_r = i_trate*(i_trate+1)/2 + j_trate;

                if ( j_trate > i_trate )
                    which_r = j_trate*(j_trate+1)/2 + i_trate;

                matrix<float> row_vect = z_dot_z(i_eff, tr_levels, i_trate, j_trate, which_r);
                //matrix<float> row_vect = z_dot_z(i_eff, tr_levels, i_trate, j_trate, r(i_trate,j_trate)); // if we use constant Rij(-1) for all observations

                for (size_t i = first_random_level - 1; i < last_random_level; i++)
                {
                    vect_a(0, i) = row_vect(0, i - (first_random_level - 1));
                }
            }

            // adding covariance structure:
            for (size_t i1 = 0; i1 < correlations; i1++)
            {
                matrix<int> which_effects = get_corr_effects(i1);
                // which_effects.fread();
                matrix<size_t> shape_eff = which_effects.shape();

                for (size_t i2 = 0; i2 < shape_eff[0]; i2++)
                {
                    for (size_t i3 = 0; i3 < shape_eff[0]; i3++)
                    {
                        size_t iblock_row = which_effects(i2, 0);
                        size_t iblock_col = which_effects(i3, 0);

                        std::vector<size_t> ioffset = cov_offsets[(iblock_row)*num_levels + iblock_col];
                        size_t first_row = ioffset[0];
                        size_t first_col = ioffset[1];

                        size_t last_row = first_row + ordered_levels[iblock_row] - 1;
                        size_t last_col = first_col + ordered_levels[iblock_col] - 1;

                        if (i_row >= first_row && i_row <= last_row)
                        {
                            size_t t_row = i_row - first_row;
                            size_t t_col1 = 0;
                            size_t t_col2 = ordered_levels[iblock_col] - 1;

                            bool identity = identity_correlation(i1);

                            matrix<float> var = get_variance(i1, i2, i3);

                            if (identity)
                            {
                                vect_a(0, first_col + t_row) = vect_a(0, first_col + t_row) + 1.0f * var[0];
                            }
                            else
                            {
                                matrix<float> corr = get_correlation(i1, t_row, t_col1, t_col2);

                                for (size_t l = first_col; l <= last_col; l++)
                                    vect_a(0, l) = vect_a(0, l) + corr(0, l - first_col) * var[0];
                            }
                        }
                    }
                }
            }

            return vect_a;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::get_row_cmatr(size_t, size_t, size_t, std::vector<std::vector<size_t>> &, size_t, std::vector<int> &, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::get_row_cmatr(size_t, size_t, size_t, std::vector<std::vector<size_t>> &, size_t, std::vector<int> &, size_t)." << '\n';
            throw;
        }   
    }

    matrix<float> Pcg::z_dot_z(size_t row, size_t vect_size, size_t i_matr, size_t j_matr, size_t r_index)
    {
        try
        {
            matrix<float> out_vect(1, vect_size); // (1,1002) => (1, n_levels)

            size_t vect_dim = n_obs[i_matr];

            // ---- Solution 1: -----------------------------
            // matrix<float> v1 = get_vect_z_uni(i_matr, row);
            // v1.transpose();
            // ----------------------------------------------

            // ---- Solution 2: -----------------------------
            // std::vector<std::vector<float>> v11(1, std::vector<float>(n_obs[0]));
            // std::vector<std::vector<float>> v22(1, std::vector<float>(n_obs[0]));
            // ----------------------------------------------

            // ---- Solution 3: -----------------------------
            float **v11 = new float *[1];
            float **v22 = new float *[1];

            v11[0] = new float[vect_dim];
            v22[0] = new float[vect_dim];
            // ----------------------------------------------

            get_vect_z_uni2(i_matr, row, v11);

            for (size_t i = 0; i < vect_dim; i++) // Multiply each element by correct R(-1) according to observations pattern
                v11[0][i] = v11[0][i] * r_map[ R_hash[i] ][ r_index ];

            for (size_t i = 0; i < vect_size; i++)
            {
                // matrix<float> v2 = get_vect_z_uni(j_matr, i); // !!! very expensive: x 115; (1, 489) => (1, n_obs)

                get_vect_z_uni2(j_matr, i, v22);

                float t = 0.0f;
                for (size_t j = 0; j < vect_dim; j++)
                    t = t + v11[0][j] * v22[0][j];

                out_vect(0, i) = t;
                // matrix<float> res = v2 * v1; // x 1.6
            }

            delete[] v11[0];
            delete[] v22[0];
            delete[] v11;
            delete[] v22;

            return out_vect;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::z_dot_z(size_t, size_t, size_t, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::z_dot_z(size_t, size_t, size_t, size_t, size_t)." << '\n';
            throw;
        }
    }

    matrix<float> Pcg::z_dot_z(size_t row, size_t vect_size, size_t i_matr, size_t j_matr, float r_val)
    {
        try
        {
            matrix<float> out_vect(1, vect_size); // (1,1002) => (1, n_levels)

            size_t vect_dim = n_obs[i_matr];

            // ---- Solution 1: -----------------------------
            // matrix<float> v1 = get_vect_z_uni(i_matr, row);
            // v1.transpose();
            // ----------------------------------------------

            // ---- Solution 2: -----------------------------
            // std::vector<std::vector<float>> v11(1, std::vector<float>(n_obs[0]));
            // std::vector<std::vector<float>> v22(1, std::vector<float>(n_obs[0]));
            // ----------------------------------------------

            // ---- Solution 3: -----------------------------
            float **v11 = new float *[1];
            float **v22 = new float *[1];

            v11[0] = new float[vect_dim];
            v22[0] = new float[vect_dim];
            // ----------------------------------------------

            get_vect_z_uni2(i_matr, row, v11);

            // size_t vect_dim = v11[0].size();

            for (size_t i = 0; i < vect_size; i++)
            {
                // matrix<float> v2 = get_vect_z_uni(j_matr, i); // !!! very expensive: x 115; (1, 489) => (1, n_obs)

                get_vect_z_uni2(j_matr, i, v22);

                float t = 0.0;
                for (size_t j = 0; j < vect_dim; j++)
                    t = t + v11[0][j] * v22[0][j];

                out_vect(0, i) = t * r_val;
                // matrix<float> res = v2 * v1; // x 1.6
                //  out_vect(0, i) = res[0] * r_val;
            }

            delete[] v11[0];
            delete[] v22[0];
            delete[] v11;
            delete[] v22;

            return out_vect;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::z_dot_z(size_t, size_t, size_t, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::z_dot_z(size_t, size_t, size_t, size_t, size_t)." << '\n';
            throw;
        }
    }

    void Pcg::construct_rhs()
    {
        // Consider trait 1 effects: x1, z1;
        //          trait 2 effects: z2;
        // Model effects stacked consecutively as [ x1 z1 z2 ]';
        // Model observations stacked as [ y1 y2 ];
        // Than RHS calculated as:
        //                       [ x1 ]                                            [ x1 * y1 * R11  x1 * y2 * R12 ]
        //                       | z1 | * [ y1 y2 ]; then summ by rows the matrix: | z1 * y1 * R11  z1 * y2 * R12 |
        //                       [ z2 ]                                            [ z2 * y1 * R21  z2 * y2 * R22 ]
        
        try
        {
            size_t all_tr_levels = get_all_levels(); // number of levels of all the traits

            rhs.resize(all_tr_levels, 1); // resize the dimension of RHS, it is the size of all levels in a model

            size_t levels_2 = 0;

            for (size_t i = 0; i < n_trait; i++)
            {
                size_t tr_levels = get_levels(i); // number of levels for a specific trait

                size_t levels_1 = levels_2 + 1; // lower index for the trait input in the RHS, starts at level_1

                levels_2 = levels_2 + tr_levels; // upper index for the trait input in the RHS, ends at level_2

                for (size_t j = 0; j < n_trait; j++)
                {
                    size_t which_r = i*(i+1)/2 + j;

                    if ( j > i )
                        which_r = j*(j+1)/2 + i;

                    matrix<float> vect = z_dot_y(tr_levels, i, j, which_r );
                    //matrix<float> vect = z_dot_y(tr_levels, i, j, r(i,j) ); // if we use constant Rij(-1) for all observations

                    size_t k2 = 0;

                    for (size_t k = levels_1 - 1; k < levels_2; k++)
                    {
                        rhs(k, 0) = rhs(k, 0) + vect(k2, 0); // tr(X1)*R11*y1 + tr(X1)*R12*y2 + ...
                        k2++;
                    }
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::construct_rhs()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::construct_rhs()." << '\n';
            throw;
        }
    }

    matrix<float> Pcg::z_dot_y(size_t vect_size, size_t i_trait, size_t j_trait, size_t r_index)
    {
        // Returns a column vector which will be inserted into a specific location of the RHS vector
        //
        // Here we do matrix-by-vector multiplication:
        // get_vect_z_uni(i_trait, i) return a row of a specific effect matrix which then
        // multiplied by the observation vector v2 and then by the corresponding component
        // of an inversed residual matrix r
        
        try
        {
            matrix<float> out_vect(vect_size, 1);

            matrix<float> v2 = y[j_trait].fget(); // observations for a specific trait

            for (size_t i = 0; i < v2.size(); i++) // multiply each v2 by relevant R(-1), corrected for missing observations
                v2(i,0) = v2(i,0) * r_map[ R_hash[i] ][r_index];

            for (size_t i = 0; i < vect_size; i++)
            {
                matrix<float> v1 = get_vect_z_uni(i_trait, i);
                matrix<float> res = v1 * v2;
                out_vect(i, 0) = res[0];
            }

            return out_vect;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::z_dot_y(size_t, size_t, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::z_dot_y(size_t, size_t, size_t, size_t)." << '\n';
            throw;
        }        
    }

    matrix<float> Pcg::z_dot_y(size_t vect_size, size_t i_trait, size_t j_trait, float r_val)
    {
        // giving a column vector
        // Here we do matrix-by-vector multiplication:
        // get_vect_z_uni(i_trait, i) return a row of a specific effect matrix which then
        // multiplied by the observation vector v2 and then by the corresponding component
        // of an inversed residual matrix r
        
        try
        {
            matrix<float> out_vect(vect_size, 1);

            matrix<float> v2 = y[j_trait].fget(); // observations for a specific trait

            for (size_t i = 0; i < vect_size; i++)
            {
                matrix<float> v1 = get_vect_z_uni(i_trait, i);
                matrix<float> res = v1 * v2;
                out_vect(i, 0) = res[0] * r_val;
            }

            return out_vect;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::z_dot_y(size_t, size_t, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::z_dot_y(size_t, size_t, size_t, size_t)." << '\n';
            throw;
        }        
    }

    size_t Pcg::get_levels(size_t which_trait)
    {
        try
        {
            size_t tr_levels = 0;

            for (auto const &e : n_lev[which_trait])
            {
                tr_levels = tr_levels + e;
            }

            return tr_levels;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::get_levels(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::get_levels(size_t)." << '\n';
            throw;
        }        
    }

    size_t Pcg::get_all_levels()
    {
        try
        {
            size_t tr_levels = 0;

            for (size_t i = 0; i < n_trait; i++)
            {
                for (auto const &e : n_lev[i])
                {
                    tr_levels = tr_levels + e;
                }
            }

            return tr_levels;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::get_all_levels()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::get_all_levels()." << '\n';
            throw;
        }
    }

    size_t Pcg::get_all_levels(size_t before_trait)
    {
        try
        {
            size_t tr_levels = 0;

            for (size_t i = 0; i < before_trait; i++)
            {
                for (auto const &e : n_lev[i])
                {
                    tr_levels = tr_levels + e;
                }
            }

            return tr_levels;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::get_all_levels(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::get_all_levels(size_t)." << '\n';
            throw;
        }        
    }

    size_t Pcg::num_all_levels()
    {
        try
        {
            size_t levels = 0;

            for (size_t i = 0; i < n_trait; i++)
            {
                levels = levels + n_lev[i].size();
            }

            return levels;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::num_all_levels()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::num_all_levels()." << '\n';
            throw;
        }        
    }

    std::vector<size_t> Pcg::get_ordered_levels()
    {
        try
        {
            std::vector<size_t> ordered_levels(num_all_levels(), 0);

            size_t olevels = 0;

            for (size_t i = 0; i < n_trait; i++)
            {
                for (auto const &e : n_lev[i])
                {
                    ordered_levels[olevels] = e;
                    olevels = olevels + 1;
                }
            }

            return ordered_levels;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::get_ordered_levels()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::get_ordered_levels()." << '\n';
            throw;
        }
    }

    std::vector<std::vector<size_t>> Pcg::get_cov_offsets(const std::vector<size_t> &ordered_levels)
    {
        /*
        %-------- Create coordinate vectors for random effect covar. blocks ------
        %
        % Here is the strategy:
        % create a matrix 'rcov_offsets' which holds [row col] values of a very first (top left corner)
        % element of each variance/covariance block of random effects in the order
        % these blocks appear in the coefficient matrix A. The matrix
        % 'rcov_offsets' is a vector representation of full matrix
        % rcov_offsets(l) = {[row col]}.
        % Indexing:
        % consider the part of the coefficient matrix A which corresponds to
        % the vsriance/covariance block of random effects B[ i = n_eff_random, j = n_eff_random ],
        % indexing is as follows B[ i, j ] => rcov_offsets( (i-1)*n_eff_random + j ).

        % n_eff_random := total number of random effects
        % nTrait       := number of traits
        % fixed_levels := number of all fixed effects' levels
        */
        try
        {
            size_t n_eff_random = ordered_levels.size();

            std::vector<std::vector<size_t>> rcov_offsets(n_eff_random * n_eff_random, std::vector<size_t>(2));

            size_t cont_index = 0;
            size_t left_shift = 0;

            for (size_t i = 0; i < n_eff_random; i++)
            {
                size_t right_shift = 0;

                for (size_t j = 0; j < n_eff_random; j++)
                {
                    rcov_offsets[cont_index][0] = left_shift;
                    rcov_offsets[cont_index][1] = right_shift;

                    cont_index = cont_index + 1;

                    right_shift = right_shift + ordered_levels[j];
                }
                left_shift = left_shift + ordered_levels[i];
            }

            return rcov_offsets;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::get_cov_offsets(const std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::get_cov_offsets(const std::vector<int> &)." << '\n';
            throw;
        }
    }

    template <typename T>
    T Pcg::v_dot_v(const matrix<T> &v1, const matrix<T> &v2)
    {
        try
        {
            matrix<T> res;

            res = v1;

            res.transpose();

            res = res * v2;

            return res[0];

            /*T res2 = (T)0.0;
            for (size_t i = 0; i < v2.size(); i++)
                res2 = res2 + v1[i]*v2[i];        
            return res2;*/
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::v_dot_v(const matrix<T> &, const matrix<T> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::v_dot_v(const matrix<T> &, const matrix<T> &)." << '\n';
            throw;
        }
    }

    template double Pcg::v_dot_v(const matrix<double> &v1, const matrix<double> &v2);
    template float Pcg::v_dot_v(const matrix<float> &v1, const matrix<float> &v2);

    void Pcg::set_tolerance(double tol)
    {
        try
        {
            tolerance = tol;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::set_tolerance(double)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::set_tolerance(double)." << '\n';
            throw;
        }
    }

    void Pcg::set_maxiter(size_t iter)
    {
        try
        {
            max_iterations = iter;
            befault_max_iter = false;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::set_maxiter(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::set_maxiter(size_t)." << '\n';
            throw;
        }
    }

#ifdef UTEST
    matrix<float> Pcg::test_z_dot_y(size_t vect_size, size_t i_trait, size_t j_trait, size_t r_index)
    {

        try
        {
            matrix<float> out_vect = z_dot_y(vect_size, i_trait, j_trait, r_index);
            return out_vect;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::test_z_dot_y(size_t, size_t, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::test_z_dot_y(size_t, size_t, size_t, size_t)." << '\n';
            throw;
        }
    }

    matrix<float> Pcg::test_rhs()
    {

        try
        {
            construct_rhs();

            matrix<float> rhs_out = rhs;

            rhs.clear();

            return rhs_out;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::test_rhs()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::test_rhs()." << '\n';
            throw;
        }
    }

    size_t Pcg::test_num_all_levels()
    {
        try
        {
            size_t n_all_levels = num_all_levels();

            return n_all_levels;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::test_num_all_levels()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::test_num_all_levels()." << '\n';
            throw;
        }
    }

    std::vector<size_t> Pcg::test_ordered_levels()
    {
        try
        {
            std::vector<size_t> ordered_random_levels = get_ordered_levels();

            return ordered_random_levels;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::test_ordered_levels()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::test_ordered_levels()." << '\n';
            throw;
        }
    }

    std::vector<std::vector<size_t>> Pcg::test_cov_offsets()
    {
        try
        {
            std::vector<size_t> ordered_random_levels = get_ordered_levels();

            std::vector<std::vector<size_t>> rcov_offsets = get_cov_offsets(ordered_random_levels);

            return rcov_offsets;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::test_cov_offsets()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::test_cov_offsets()." << '\n';
            throw;
        }
    }

    std::vector<float> Pcg::test_dval()
    {
        try
        {
            std::vector<float> out;

            construct_rhs();

            size_t n_all_levels = num_all_levels(); // number of all random effects in the model

            std::vector<size_t> ordered_random_levels = get_ordered_levels();

            std::vector<std::vector<size_t>> rcov_offsets = get_cov_offsets(ordered_random_levels);

            // set_pipeline(1);

            matrix<double> Mi = construct_dval(rcov_offsets, n_all_levels, ordered_random_levels);

            for (size_t i = 0; i < Mi.size(); i++)
            {
                if (Mi[i] == 0.0)
                    throw std::string("Pcg::test_dval() => Expected division by 0.0: Mi[i] == 0.0!");

                out.push_back( static_cast<float>(1.0 / Mi[i]) );
            }

            rhs.clear();

            return out;
        }
        catch (const std::string &err)
        {
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::test_dval()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::test_dval()." << '\n';
            throw;
        }
    }

    std::vector<std::vector<float>> Pcg::construct_A(std::vector<std::vector<size_t>> &cov_offsets, size_t num_levels, std::vector<size_t> &ordered_levels)
    {
        try
        {
            matrix<size_t> shape_rhs = rhs.shape();

            if (shape_rhs[0] == 0)
                throw std::string("The size of RHS is 0!");

            std::vector<std::vector<float>> A(shape_rhs[0], std::vector<float>(shape_rhs[0]));

            size_t cmatr_row = 0;

            for (size_t i_trate = 0; i_trate < n_trait; i_trate++)
            {
                size_t tr_levels = get_levels(i_trate);

                for (size_t i_eff = 0; i_eff < tr_levels; i_eff++)
                {
                    matrix<float> vect = get_row_cmatr(shape_rhs[0], i_trate, i_eff, cov_offsets, num_levels, ordered_levels, cmatr_row);
                    for (size_t i2 = 0; i2 < shape_rhs[0]; i2++)
                    {
                        A[cmatr_row][i2] = vect(0, i2);
                    }
                    cmatr_row = cmatr_row + 1;
                }
            }

            return A;
        }
        catch (const std::string &err)
        {
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::construct_A(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::construct_A(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            throw;
        }
    }

    std::vector<std::vector<float>> Pcg::test_A()
    {
        try
        {
            construct_rhs();

            size_t n_all_levels = num_all_levels(); // number of all random effects in the model

            std::vector<size_t> ordered_random_levels = get_ordered_levels();

            std::vector<std::vector<size_t>> rcov_offsets = get_cov_offsets(ordered_random_levels);

            std::vector<std::vector<float>> A = construct_A(rcov_offsets, n_all_levels, ordered_random_levels);

            rhs.clear();

            return A;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Pcg::test_A()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Pcg::test_A()." << '\n';
            throw;
        }
    }

#endif
}
