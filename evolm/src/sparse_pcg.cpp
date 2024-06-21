#include "sparse_pcg.hpp"

namespace evolm
{
    void sparse_pcg::solve()
    {
        try
        {
            //set_pipeline(pipeline_val); // defines a memory management

            memload_effects();
            //memload_var();
            memload_cor();
            memload_cor_effects();

            construct_rhs();

            set_model_matrix();

            diskload_effects();
            diskload_var();
            diskload_cor();
            diskload_cor_effects();
            
            jacobi_pcg();

            for (size_t i = 0; i < get_num_of_mem_blocks(); i++)
                unload_model_matrix(i);
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
            std::cerr <<"Reason: "<< e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::solve()." << '\n';
            throw;
        }
    }

    void sparse_pcg::solve(int pipeline_index)
    {
        try
        {
            pipeline_val = pipeline_index;

            set_pipeline(pipeline_val);

            construct_rhs();

            //size_t n_all_levels = num_all_levels(); // number of all random effects in the model

            //std::vector<size_t> ordered_random_levels = get_ordered_levels();

            //std::vector<std::vector<size_t>> rcov_offsets = get_cov_offsets(ordered_random_levels);

            jacobi_pcg();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::solve(int)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::solve(int)." << '\n';
            throw;
        }
    }

    void sparse_pcg::set_pipeline(int which_pipeline)
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
            std::cerr << "Exception in sparse_pcg::set_pipeline(int): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::set_pipeline(int). " << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::set_pipeline(int)." << '\n';
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

    int sparse_pcg::get_solution(const std::string &fname)
    {
        try
        {
            std::ofstream solution(fname);

            if (solution.is_open())
            {
                for (size_t i = 0; i < sol.size(); i++)
                    solution << sol[i] << "\n";

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

    matrix<float> sparse_pcg::construct_dval(std::vector<std::vector<size_t>> &cov_offsets, size_t num_levels, std::vector<size_t> &ordered_levels)
    {
        try
        {
            matrix<size_t> shape_rhs = rhs.shape();

            if (shape_rhs[0] == 0)
                throw std::string("The size of RHS is 0!");

            matrix<float> dval(shape_rhs[0], 1);

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
                        throw std::string("sparse_pcg::construct_dval(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &) => Expected division by 0.0!");

                    dval(cmatr_row, 0) = 1.0f / vect(0, cmatr_row);

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
            std::cerr << "Exception in sparse_pcg::construct_dval(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::construct_dval(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            throw;
        }
    }

    void sparse_pcg::construct_dval2(std::vector<double> &inverted_diagonal)
    {
        try
        {
            if (inverted_diagonal.size() != model_matrix.size())
                throw std::string("The size allocated for the vector of model matrix diagonals is not correct!");

            size_t first_row = 0;
            size_t last_row = 0;

            for (size_t i = 0; i < get_num_of_mem_blocks(); i++)
            {
                load_model_matrix(i);
                get_mem_block_range(i, first_row, last_row);

#pragma omp parallel for //num_threads(3)
                for (size_t j = first_row; j <= last_row; j++)
                {
                    float d = model_matrix[j].value_at(0,j);

                    if (d == 0.0f)
                        throw std::string("sparse_pcg::construct_dval2(std::vector<float> &) => The diagonal element of the model matrix is 0.0!");
                    
                    inverted_diagonal[j] = static_cast<double>(1.0 / d);
                }

                unload_model_matrix(i);
            }
        }
        catch (const std::string &err)
        {
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::construct_dval2(std::vector<float> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::construct_dval2(std::vector<float> &)." << '\n';
            throw;
        }
    }

    void sparse_pcg::memload_amatr()
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
                        throw std::string("Error while opening a binary file in sparse_pcg::memload_amatr()");

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
                    throw std::string("Error while getting a binary file name in sparse_pcg::memload_amatr()");
            }
            else if (!amatrix_ondisk && !amatrix_onmem)
            {
                // build amatr directly on memory

                size_t n_all_levels = num_all_levels();

                size_t all_levels = get_all_levels();

                std::vector<size_t> ordered_random_levels = get_ordered_levels();

                std::vector<std::vector<size_t>> rcov_offsets = get_cov_offsets(ordered_random_levels);

                amatr.resize(all_levels, all_levels);

                set_amatr(rcov_offsets, n_all_levels, ordered_random_levels, true);

                //---------------------------
                /*std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

                set_amatr(rcov_offsets, n_all_levels, ordered_random_levels, true);

                std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

                std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

                std::cout << "time to set amatrix: " << time_span.count() << " seconds."
                          << "\n";*/
                //-----------------------------

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
            std::cerr << "Exception in sparse_pcg::memload_amatr()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::memload_amatr()." << '\n';
            throw;
        }
    }

    void sparse_pcg::set_amatrix()
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
            std::cerr << "Exception in sparse_pcg::set_amatrix()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::set_amatrix()." << '\n';
            throw;
        }
    }

    void sparse_pcg::set_amatr(std::vector<std::vector<size_t>> &cov_offsets, size_t num_levels, std::vector<size_t> &ordered_levels)
    {
        // Build coefficient matrix on disk

        // Consider trait 1 effects: x1, z1;
        //          trait 2 effects: z2;
        // Model effects stacked consecutively as [ x1 z1 z2 ]';
        // Model observations stacked as [ y1 y2 ];
        // Than coefficient matrix calculated as (row-by-row):
        // [ x1 ]                   [ x1' * x1 * R11  x1' * z1 * R11  x1' * z2 * R12 ]
        // | z1 | * [ x1 z1 z2 ] =  | z1' * x1 * R11  z1' * z1 * R11  z1' * z2 * R12 |
        // [ z2 ]                   [ z2' * x1 * R21  z2' * z1 * R21  z2' * z2 * R22 ]

        try
        {
            size_t all_tr_levels = get_all_levels();

            if (all_tr_levels == 0)
                throw std::string("The size of Unknowns is 0!");

            fA.open(binFilename, fA.binary | fA.app | fA.out);

            if (!fA.is_open())
                throw std::string("Error while opening a binary file. sparse_pcg::set_amatr()");

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
            std::cerr << err << "  In in sparse_pcg::set_amatr(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::set_amatr(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::set_amatr(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            throw;
        }
    }

    void sparse_pcg::set_amatr(std::vector<std::vector<size_t>> &cov_offsets, size_t num_levels, std::vector<size_t> &ordered_levels, bool on_mem)
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

                std::cout << "n_threads: " << n_threads << "\n";

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
            std::cerr << err << "  In in sparse_pcg::set_amatr(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &, bool)." << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::set_amatr(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &, bool)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::set_amatr(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &, bool)." << '\n';
            throw;
        }
    }

    matrix<float> sparse_pcg::fget_vect(size_t all_tr_levels, size_t row)
    {
        try
        {
            matrix<float> vect(1, all_tr_levels);

            std::ifstream f(binFilename.c_str());

            if (f.good())
            {
                fA.open(binFilename, fA.binary | fA.in);

                if (!fA.is_open())
                    throw std::string("Error while opening a binary file in sparse_pcg::fget_vect(size_t row)");

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
                throw std::string("Error while getting a binary file name in sparse_pcg::fget_vect(size_t row)");

            return vect;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_pcg::fget_vect(size_t, size_t): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::fget_vect(size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::fget_vect(size_t, size_t)." << '\n';
            throw;
        }
    }

    void sparse_pcg::update_vect(std::vector<std::vector<size_t>> &cov_offsets,
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
            std::cerr << "Exception in sparse_pcg::update_vect(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &, matrix<double> &, matrix<double> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::update_vect(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &, matrix<double> &, matrix<double> &)." << '\n';
            throw;
        }
    }

    void sparse_pcg::update_vect2( std::vector<double> &out_vect, std::vector<double> &in_vect )
    {
        try
        {
            size_t first_row = 0;
            size_t last_row = 0;

            for (size_t i = 0; i < get_num_of_mem_blocks(); i++)
            {
                load_model_matrix(i);
                get_mem_block_range(i, first_row, last_row);

#pragma omp parallel for //num_threads(3)
                for (size_t j = first_row; j <= last_row; j++)
                {
                    double result = 0.0;
                    model_matrix[j].vect_dot_vect(in_vect, result);
                    out_vect[j] = result;
                }

                unload_model_matrix(i);
            }
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
                max_iterations = unknowns * 10;

            std::vector<double> Mi(unknowns, 0.0);

            construct_dval2(Mi);

            sol.resize(unknowns, 0.0);
            for (size_t i = 0; i < unknowns; i++)
                sol[i] = static_cast<double>(rhs[i]) * Mi[i]; // initial solution

            std::vector<double> _rhs(unknowns, 0.0);

            for (size_t i = 0; i < unknowns; i++)
                _rhs[i] = static_cast<double>(rhs[i]);

            std::vector<double> tVect(unknowns, 1); // vector to keep the result of operation: A*x

            update_vect2(tVect, sol); // A*x(==sol)

            std::vector<double> r_vect(unknowns, 0.0);

            for (size_t i = 0; i < unknowns; i++)
                r_vect[i] = _rhs[i] - tVect[i]; // r = b - A*x

            tVect.clear();

            std::vector<double> d(unknowns, 0.0);

            for (size_t i = 0; i < unknowns; i++)
                d[i] = Mi[i] * r_vect[i];

            double delta_new = v_dot_v2(r_vect, d);

            double delta_zero = delta_new;

            iterations = 1;

            std::vector<double> q(unknowns, 0.0);
            std::vector<double> s(unknowns, 0.0);
            double i_value = 0.0;
            double alpha = 0.0;
            double delta_old = 0.0;
            double betha = 0.0;

            while (iterations < max_iterations && delta_new > /*delta_zero */ tolerance * tolerance)
            {
                update_vect2(q, d); // here we casting vector from float to double every time by calling this function

                i_value = v_dot_v2(d, q);

                if (i_value == 0.0)
                    throw std::string("sparse_pcg::jacobi_pcg() => Expected division by 0.0: v_dot_v(d, q) == 0.0!");

                alpha = delta_new / i_value;

                for (size_t i = 0; i < sol.size(); i++)
                    sol[i] = sol[i] + alpha * d[i];

                if ( !(iterations % 50) )
                {
                    update_vect2(tVect, sol);

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

                delta_new = v_dot_v2(r_vect, s);

                if (delta_old == 0.0)
                    throw std::string("sparse_pcg::jacobi_pcg() => Expected division by 0.0: delta_old == 0.0!");

                betha = delta_new / delta_old;

                for (size_t i = 0; i < unknowns; i++)
                    d[i] = s[i] + betha * d[i];

                // debugging
                //if (!(iterations % 1))
                //    std::cout << "max_iterations: "<< max_iterations << "; iter: " << iterations << "; delta_new: " << delta_new << " condition: "<< /*delta_zero */ tolerance * tolerance << "\n";

                iterations = iterations + 1;
            }

            iterations = iterations - 1;
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

    template <typename T>
    T sparse_pcg::v_dot_v(const matrix<T> &v1, const matrix<T> &v2)
    {
        try
        {
            matrix<T> res;

            res = v1;

            res.transpose();

            res = res * v2;

            return res[0];
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_pcg::v_dot_v(const matrix<T> &, const matrix<T> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_pcg::v_dot_v(const matrix<T> &, const matrix<T> &)." << '\n';
            throw;
        }
    }

    template double sparse_pcg::v_dot_v(const matrix<double> &v1, const matrix<double> &v2);
    template float sparse_pcg::v_dot_v(const matrix<float> &v1, const matrix<float> &v2);

    double sparse_pcg::v_dot_v2(std::vector<double> &v1, std::vector<double> &v2)
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

            construct_dval2(out);

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
