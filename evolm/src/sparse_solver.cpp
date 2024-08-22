#include "sparse_solver.hpp"

namespace evolm
{
    sparse_solver::sparse_solver()
    {
        const auto processor_count = std::thread::hardware_concurrency();
        available_cpu = processor_count;
    }

    sparse_solver::~sparse_solver()
    {
        try
        {
            remove_model();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::~sparse_solver()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver destructor." << '\n';
        }
    }

    void sparse_solver::append_model(const model_sparse &m)
    {
        try
        {
            model = m;
            process_model();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::append_model(const model_sparse &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::append_model(const model_sparse &)." << '\n';
            throw;
        }
    }

    void sparse_solver::remove_model()
    {
        try
        {
            model.clear();

            std::vector<size_t>().swap(n_obs);

            std::vector<std::vector<size_t>>().swap(n_lev);

            for (auto &e : z_uni)
            {
                for (auto &p : e)
                {
                    p.clear();
                }
            }
            std::vector<std::vector<effects_storage>>().swap(z_uni);

            for (auto &e : y)
            {
                e.fclear();
                e.clear();
            }
            std::vector<matrix<float>>().swap(y);

            s.clear();
            s.shrink_to_fit();

            R_hash.clear();
            R_hash.shrink_to_fit();

            r_map.clear();

            z_on_memory = false;

            adj_effects_order.clear();

            /*clear_model_matrix();

            bin_fnames.clear();
            bin_fnames.shrink_to_fit();

            blocks_ranges.clear();
            blocks_ranges.shrink_to_fit();*/
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::remove_model()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::remove_model()." << '\n';
            throw;
        }
    }

    void sparse_solver::process_model()
    {
        try
        {
            n_trait = model.observation_trait.size(); // get the number of traits by the number of indicated associations
                                                      // between submitted observations and considered traits

            size_t i_adj_effects = 0;

            for (size_t i = 0; i < n_trait; i++)
            {
                int trait_obs_id = model.observation_trait[i]; // index of submitted observations vector in the observations std::vector

                if (((size_t)trait_obs_id > model.observations.size() - 1) || (trait_obs_id < 0))
                    throw std::string("Wrong index of observations vector associated to a trait!");

                matrix<int> trait_eff_ids = model.effects_trait[i]; // ids for the submitted effects associated with a trait

                size_t n_trait_effects = trait_eff_ids.size();

                size_t n_obs_trait = model.observations[trait_obs_id].size(); // number of observations for a specific trait

                n_obs.push_back(model.observations[trait_obs_id].size()); // get the number of observations for the specific trait

                std::vector<size_t> trait_levels;

                for (size_t j = 0; j < n_trait_effects; j++)
                {
                    std::vector<size_t> shape = model.all_effects[trait_eff_ids[j]].shape();

                    trait_levels.push_back(shape[0]); // on transposed effects

                    if (shape[1] != n_obs_trait)
                        throw std::string("The dimension of provided effect matrix corresponding to observations is not correct!");

                    adj_effects_order[trait_eff_ids[j]] = i_adj_effects; // map between the submitted effects ids and their recoded (consecutive) order

                    i_adj_effects++;
                }

                n_lev.push_back(trait_levels); // get the effects' levels for each trait

                // Combines a consecutive sets of incidense matrices (stack of effects matrices)
                // of different (original & unchanged) types for a specific trait
                //construct_union_z(trait_eff_ids, false); // build std::vector<std::vector<effects_storage>> z_uni;
                                                           // <- everything on disk !
                construct_union_z(trait_eff_ids, true);    // <- builds (and keeps) container directly on memory

                set_y(trait_obs_id); // moves all observations to std::vector<matrix<float>> y
            }

            set_r();
            set_g();

            // Checking if dimensions of correlated effects are matching:
            for (size_t i = 0; i < model.correlated_effects.size(); i++) // loop over provided correlation structures
            {
                matrix<int> cor_effects = model.correlated_effects[i]; // list of correlated effects IDs
                cor_effects.fread();
                size_t n_effects = cor_effects.size();                 // number of correlated effects
                
                std::vector<size_t> shape1(2); // get shape of correlation matrix
                shape1[0] = model.correlations[i].nrows();
                shape1[1] = model.correlations[i].ncols();
                //matrix<size_t> shape1 = model.correlations[i].shape();

                if (model.identity_correlations[i]) // in case of correlation matrix is identity matrix
                    shape1[0] = model.identity_dimension[i];

                for (size_t j = 0; j < n_effects; j++) // loop over correlated effects
                {
                    size_t which_effect = cor_effects(j, 0);                     // effect ID
                    std::vector<size_t> shape2 = model.all_effects[which_effect].shape(); // shape of effect matrix, it is transposed, so n_lev-by-n_obs

                    if (shape2[0] != shape1[0]) // size of n_levels should be equal to a size of correlation matrix (symmetric)
                        throw std::string("The provided correlation structure missmatch! Dimensions of a correlation matrix and correlated effects should match!");
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::process_model()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_solver::process_model()" << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::process_model()." << '\n';
            throw;
        }
    }

    size_t sparse_solver::get_levels(size_t which_trait)
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
            std::cerr << "Exception in sparse_solver::get_levels(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::get_levels(size_t)." << '\n';
            throw;
        }        
    }

    size_t sparse_solver::get_all_levels()
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
            std::cerr << "Exception in sparse_solver::get_all_levels()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::get_all_levels()." << '\n';
            throw;
        }
    }

    size_t sparse_solver::get_all_levels(size_t before_trait)
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
            std::cerr << "Exception in sparse_solver::get_all_levels(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::get_all_levels(size_t)." << '\n';
            throw;
        }        
    }

    size_t sparse_solver::num_all_levels()
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
            std::cerr << "Exception in sparse_solver::num_all_levels()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::num_all_levels()." << '\n';
            throw;
        }        
    }

    std::vector<size_t> sparse_solver::get_ordered_levels()
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
            std::cerr << "Exception in sparse_solver::get_ordered_levels()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::get_ordered_levels()." << '\n';
            throw;
        }
    }

    std::vector<std::vector<size_t>> sparse_solver::get_cov_offsets(const std::vector<size_t> &ordered_levels)
    {
        /*
        %-------- Create coordinate vectors for random effect covar. blocks ------
        %
        % Here is the strategy:
        % create a matrix 'rcov_offsets' which holds [row col] values of a very first (top left corner)
        % element of each variance/covariance block of random effects in the order
        % these blocks appear in the model (coefficient) matrix A. The matrix
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
            std::cerr << "Exception in sparse_solver::get_cov_offsets(const std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::get_cov_offsets(const std::vector<int> &)." << '\n';
            throw;
        }
    }

    void sparse_solver::construct_rhs()
    {
        // Consider trait 1 effects: x1, z1;
        //          trait 2 effects: z2;
        // Model effects stacked consecutively as [ x1 z1 z2 ]';
        // Model observations stacked as [ y1 y2 ];
        // Than RHS calculated as:
        // [ x1 ]                                            [ x1 * y1 * R11  x1 * y2 * R12 ]
        // | z1 | * [ y1 y2 ]; then summ by rows the matrix: | z1 * y1 * R11  z1 * y2 * R12 |
        // [ z2 ]                                            [ z2 * y1 * R21  z2 * y2 * R22 ]
        
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

                    matrix<float> vect;
                    z_dot_y(vect, tr_levels, i, j, which_r );

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
            std::cerr << "Exception in sparse_solver::construct_rhs()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::construct_rhs()." << '\n';
            throw;
        }
    }

    void sparse_solver::z_dot_y(matrix<float> &out_vect, size_t vect_size, size_t i_trait, size_t j_trait, size_t r_index)
    {
        // Returns a column vector which will be inserted into a specific location of the RHS vector
        //
        // Here we do matrix-by-vector multiplication:
        // get_vect_z_uni(i_trait, i) return a row of a specific effect matrix which then
        // multiplied by the observation vector v2 and then by the corresponding component
        // of an inversed residual matrix r
        
        try
        {
            out_vect.resize(vect_size, 1);

            matrix<float> v2 = y[j_trait].fget(); // observations for a specific trait

            for (size_t i = 0; i < v2.size(); i++) // multiply each v2 by relevant R(-1), corrected for missing observations
                v2(i,0) = v2(i,0) * r_map[ R_hash[i] ][r_index];

            for (size_t i = 0; i < vect_size; i++)
            {
                std::vector<float> vals;
                std::vector<size_t> keys;

                get_vect_z_uni2(i_trait, i, vals, keys);

                float result = 0;

                for( size_t j = 0; j < keys.size(); j++)
                    result = result + v2[ keys[j] ] * vals[ j ];
                
                out_vect(i, 0) = result;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::z_dot_y(size_t, size_t, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::z_dot_y(size_t, size_t, size_t, size_t)." << '\n';
            throw;
        }        
    }

    void sparse_solver::z_dot_z2(std::vector<float> &out_values, std::vector<size_t> &out_keys, size_t row, size_t vect_size, size_t i_matr, size_t j_matr, size_t r_index)
    {
        try
        {
            //out_vect.resize(1, vect_size); // (1, n_levels)

            size_t vect_dim = n_obs[i_matr];

            size_t i_row = 0;
            size_t last_col = 0;
            size_t which_eff_matr = 0; // will be changed inside z_row_to_uni_col(...)
            size_t icol[] = {0, 0}; // will be changed inside z_row_to_uni_col(...)

            float **v11 = new float *[1];
            v11[0] = new float[vect_dim];

            get_vect_z_uni2(i_matr, row, v11); // returns dense vector v11

            for (size_t i = 0; i < vect_dim; i++) // Multiply each element by correct R(-1) according to observations pattern
                v11[0][i] = v11[0][i] * r_map[ R_hash[i] ][ r_index ];

            for (size_t i = 0; i < vect_size; i++)
            {
                z_row_to_uni_col(j_matr, i, icol, which_eff_matr);

                last_col = n_obs[j_matr] - 1; // this is constant for any row, due to the same amount of observations
                i_row = icol[0];

                float t = 0.0f;

                z_uni[j_matr][which_eff_matr].row_dot_float_vect(v11, i_row, last_col, t); // returns t, the product of v11 and the row i_row

                if ( t != 0.0f )
                {
                    out_values.push_back(t);
                    out_keys.push_back(i);
                }
            }

            delete[] v11[0];
            delete[] v11;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::z_dot_z2(std::vector<float> &, std::vector<size_t> &, size_t, size_t, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::z_dot_z2(std::vector<float> &, std::vector<size_t> &, size_t, size_t, size_t, size_t)." << '\n';
            throw;
        }
    }
    void sparse_solver::set_model_matrix()
    {
        try
        {
            size_t n_all_levels = num_all_levels();

            std::vector<size_t> ordered_random_levels = get_ordered_levels();

            std::vector<std::vector<size_t>> rcov_offsets = get_cov_offsets(ordered_random_levels);

            model_matrix.resize( get_all_levels() );
            first_row_ondisk = model_matrix.size();

//auto start = std::chrono::high_resolution_clock::now();

            make_model_matrix(rcov_offsets, n_all_levels, ordered_random_levels);

//auto stop = std::chrono::high_resolution_clock::now();
//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
//std::cout <<"set_amatr() (milliseconds): "<< duration.count() << std::endl;

        }
        catch (const std::string &err)
        {
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::set_model_matrix()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::set_model_matrix()." << '\n';
            throw;
        }
    }

    /*size_t sparse_solver::get_num_of_mem_blocks()
    {
        return bin_fnames.size();
    }

    void sparse_solver::get_mem_block_range(size_t mem_blok, size_t &first, size_t &second)
    {
        first = blocks_ranges[mem_blok][0];
        second = blocks_ranges[mem_blok][1];
    }

    void sparse_solver::load_model_matrix(size_t mem_blok)
    {
        fread(bin_fnames[mem_blok], blocks_ranges[mem_blok][0], blocks_ranges[mem_blok][1]);
    }

    void sparse_solver::unload_model_matrix(size_t mem_blok)
    {
        for (size_t i = blocks_ranges[mem_blok][0]; i <= blocks_ranges[mem_blok][1]; i++)
            model_matrix[i].remove_data();
    }*/

    void sparse_solver::unload_model_matrix(size_t first_row, size_t last_row)
    {
        for (size_t i = first_row; i <= last_row; i++)
            model_matrix[i].remove_data();
    }

    void sparse_solver::fwrite(const std::string &fname, size_t first_row, size_t last_row)
    {
        std::fstream fA;
        fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        fA.open(fname, fA.binary | fA.trunc | fA.out);

        if (!fA.is_open())
            throw std::string("sparse_solver::fwrite(const std::string &, size_t, size_t): Error while opening a binary file.");

        for (size_t i = first_row; i <= last_row; i++)
            model_matrix[i].fwrite(fA);
        
        fA.close();
    }

    void sparse_solver::fwrite(std::fstream &external_fA, size_t first_row, size_t last_row)
    {
        if (!external_fA.is_open())
            throw std::string("sparse_solver::fwrite(const std::string &, size_t, size_t): Error while opening a binary file.");

        for (size_t i = first_row; i <= last_row; i++)
        {
            model_matrix[i].bin_file_read_position = external_fA.tellp(); // get the current position of the stream
            model_matrix[i].fwrite(external_fA);
        }        
    }

    void sparse_solver::fread(const std::string &fname, size_t first_row, size_t last_row)
    {
        std::fstream fA;
        fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        fA.open(fname, fA.binary | fA.in);

        if (!fA.is_open())
            throw std::string("sparse_solver::fread(const std::string &, size_t, size_t): Error while opening a binary file.");

        for (size_t i = first_row; i <= last_row; i++)
            model_matrix[i].fread(fA);
        
        fA.close();
    }

    void sparse_solver::fclear()
    {
        std::ifstream f(bin_fname.c_str());
        if (f.good())
            remove(bin_fname.c_str());
    }

    void sparse_solver::clear_model_matrix()
    {
        fclear();
        for (size_t i = 0; i < model_matrix.size(); i++)
                model_matrix[i].clear();
    }

    std::string sparse_solver::create_fname()
    {
        std::random_device rd;
        srand(rd());
        int iNum = rand() % 1000000;

        return "model_matrix_" + std::to_string(iNum);
    }

    void sparse_solver::set_memory_limit(double limit)
    {
        available_memory = limit * 0.95;
    }

    double sparse_solver::get_memory_limit()
    {
        return available_memory;
    }

    void sparse_solver::set_data_size(double dat_size)
    {
        data_size = dat_size;
    }

    void sparse_solver::set_cpu_limit(int limit)
    {
        available_cpu = limit;
    }

    void sparse_solver::get_load_per_memory_block(std::vector<std::vector<size_t>> &loads)
    {
        loads.resize(n_trait, std::vector<size_t>(2));

        row_size_upper_bound = ( (double)get_all_levels() * 4.0 ) / gb_constant; // upper bound
        double runtime_data_size = (double)available_cpu * row_size_upper_bound * (4.0); // upper bound

        double rows_per_memory = ( available_memory - 2.0 * data_size - runtime_data_size ) / row_size_upper_bound;

        if (rows_per_memory <= 0.0)
            throw std::string("Not enough memory, leading to rows_per_memory == 0. In order to fix the issue: increase a memory limit or decrease the number of available cpu.");
        
        for (size_t i_trate = 0; i_trate < n_trait; i_trate++)
        {
            size_t rows_in_trait = get_levels(i_trate);
            size_t blocks_per_trait = std::ceil( (double)rows_in_trait / rows_per_memory);

            loads[i_trate][0] = blocks_per_trait;            

            if (blocks_per_trait == 0)
            {
                loads[i_trate][0] = 1;
                loads[i_trate][1] = rows_in_trait;
            }
            else
                loads[i_trate][1] = std::floor( (double)rows_in_trait / (double)blocks_per_trait ); // rows per block            
        }    
    }

    void sparse_solver::make_model_matrix(std::vector<std::vector<size_t>> &cov_offsets, size_t num_levels, std::vector<size_t> &ordered_levels)
    {
        // Build coefficient matrix on disk
        // The used approach:
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
            bin_fname = create_fname();
            std::fstream fA;
            fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            fA.open(bin_fname, fA.binary | fA.trunc | fA.out);

            std::vector<std::vector<size_t>> loads;
            get_load_per_memory_block(loads);

            double running_data_size = (double)available_cpu * row_size_upper_bound * 4.0 + 2.0 * data_size; // upper bound
            double total_used_memory = 0.0;
            double size_of_block = row_size_upper_bound * (double)loads[0][1];
            bool calcullate_ocupied_memoty = true;

            size_t n_blocks = 0;
            size_t rows_per_block = 0;

            int all_blocks = 0;

            for (size_t i_trate = 0; i_trate < n_trait; i_trate++)
            {
                n_blocks = loads[i_trate][0]; // number of model_matrix blocks needed to fit the available memoory, sizeof(i_block) ~ sizeof(memory)
                rows_per_block = loads[i_trate][1]; // number of rows of the model_matrix that fit to a memory (block)

                size_t first_block = 0; // first processed row in the block
                size_t last_block = rows_per_block; // last processed row-1 of the block

                for (size_t b = 0; b < n_blocks; b++) // process model_matrix by blocks
                {
                    std::vector<size_t> range(2, 0); // used rows range of the current memory block
                    range[0] = get_all_levels(i_trate) + first_block; // first_processed_row of the block

                    if ( b == (n_blocks - 1) )
                        last_block = get_levels(i_trate);

#pragma omp parallel for num_threads( available_cpu )
                    for (size_t i_eff = first_block; i_eff < last_block; i_eff++) // private loop
                    {
                        size_t private_raw = get_all_levels(i_trate) + i_eff;
                        size_t all_tr_levels = get_all_levels();
                        compact_storage<float> vect;
                        get_row_cmatr2(vect, all_tr_levels, i_trate, i_eff, cov_offsets, num_levels, ordered_levels, private_raw);                        
                        model_matrix[private_raw] = vect; // the size of vect should be (0, all_tr_levels)
                    } // end of within_block loop
                    //---------------------------------------------------------------------------------------------
                    // end of memory block fitting to a memory while building a model_matrix, write it to a disck                     
                    range[1] = get_all_levels(i_trate) + last_block - 1; // last_processed_row of the block                     

                    all_blocks++;

                    if (calcullate_ocupied_memoty)
                    {
                        for (size_t i2 = range[0]; i2 <= range[1]; i2++ ) // get used memory, in GB
                            total_used_memory = total_used_memory + model_matrix[i2].size_inmem() / gb_constant;
                        
//std::cout<<" used memory: "<<total_used_memory<<" size_of_block "<<size_of_block<<" running_data_size "<<running_data_size<<" available_memory "<<available_memory<<"\n";
                        if ( ( total_used_memory + size_of_block + running_data_size ) > available_memory ) // update memory usage untill we use it all
                        {
                            calcullate_ocupied_memoty = false;
                            first_row_ondisk = range[0]; // will be the very first row writen into disck
                        }
                    }

                    if ( !calcullate_ocupied_memoty ) // update memory usage untill we use it all
                    { // if the memory completely ocupied, write the data block to the disck
            //bool t_var = true;
            //while (t_var)
            //    std::cout<<"waiting ..."<<"\r";

                        fwrite(fA, range[0], range[1]); // relocate the content of the block to a disck
                        
                        // another variant for comparison
                        //blocks_ranges.push_back(range);
                        //bin_fnames.push_back( create_fname() ); // generate and save a file name for the block
                        //fwrite(bin_fnames.back(), range[0], range[1]); // relocate the content of the block to a disck
                        //-------------------------------

                        unload_model_matrix(range[0], range[1]);
                        //unload_model_matrix(b);
//std::cout<<"writing completed, writing range: "<<range[0]<<" "<<range[1]<<"\n";
            //bool t_var = true;
            //while (t_var)
            //    std::cout<<"waiting ..."<<"\r";

                    }                    
                    //---------------------------------------------------------------------------------------------
                    // update the block range (first and last rows)
                    first_block = last_block;
                    last_block = first_block + rows_per_block;
                } // end of blocks_within_trait loop
            } // end of traits loop

            fA.close();
//std::cout<<"all blocks: "<<all_blocks<<" first_row_ondisk "<<first_row_ondisk<<" last row in matrix "<<model_matrix.size()-1<<"\n";
        }
        catch (const std::string &err)
        {
            std::cerr << err << "  In in sparse_solver::make_model_matrix(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::make_model_matrix(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::make_model_matrix(std::vector<std::vector<size_t>> &, size_t, std::vector<int> &)." << '\n';
            throw;
        }
    }

    void
    sparse_solver::get_row_cmatr2(
                        compact_storage<float> &model_matrix_row,
                        size_t rhs_size,
                        size_t i_trate,
                        size_t i_eff,
                        std::vector<std::vector<size_t>> &cov_offsets,
                        size_t num_levels,
                        std::vector<size_t> &ordered_levels,
                        size_t i_row)
    {
        try
        {
            std::vector<float> vect_a(rhs_size, 0.0f);

            model_matrix_row.resize(1, rhs_size);

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

                std::vector<float> vals;
                std::vector<size_t> keys;

                z_dot_z2(vals, keys, i_eff, tr_levels, i_trate, j_trate, which_r);

                for ( size_t i = 0; i < keys.size(); i++ )
                    vect_a[ keys[i] + (first_random_level - 1) ] = vals[i];
            }

            // adding covariance structure:
            for (size_t i1 = 0; i1 < correlations; i1++)
            {
                matrix<int> which_effects = get_corr_effects(i1); // dim = (n_eff, 0)

                matrix<size_t> shape_eff = which_effects.shape();

                for (size_t i2 = 0; i2 < shape_eff[0]; i2++) // loop over number of correlated effects
                {
                    for (size_t i3 = 0; i3 < shape_eff[0]; i3++) // correlated effects by correlated effects
                    {
                        size_t iblock_row = which_effects(i2, 0);
                        size_t iblock_col = which_effects(i3, 0);

                        std::vector<size_t> ioffset = cov_offsets[(iblock_row)*num_levels + iblock_col];
                        
                        size_t first_row = ioffset[0];
                        size_t first_col = ioffset[1];

                        size_t last_row = first_row + ordered_levels[iblock_row] - 1;

                        if (i_row >= first_row && i_row <= last_row)
                        {
                            size_t t_row = i_row - first_row;
                            size_t t_col1 = 0;
                            size_t t_col2 = ordered_levels[iblock_col] - 1;

                            bool identity = identity_correlation(i1);

                            float var = get_variance(i1, i2, i3);

                            if (identity)
                                vect_a[ first_col + t_row ] = vect_a[ first_col + t_row ] + var;
                            else
                                add_correlation(vect_a, first_col, var, i1, t_row, t_col1, t_col2);
                        }
                    }
                }
            }
            model_matrix_row.append(vect_a);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::get_row_cmatr2(size_t, size_t, size_t, std::vector<std::vector<size_t>> &, size_t, std::vector<int> &, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::get_row_cmatr2(size_t, size_t, size_t, std::vector<std::vector<size_t>> &, size_t, std::vector<int> &, size_t)." << '\n';
            throw;
        }   
    }

    void sparse_solver::memload_effects()
    {
        try
        {
            if ( z_uni.empty() )
                throw std::string("z_uni is empty!");

            if (!z_on_memory)
            {
                for (auto &e : z_uni)
                {
                    for (auto &p : e)
                        p.fread();
                }
                z_on_memory = true;
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_solver::memload_effects(): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::memload_effects()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::memload_effects()." << '\n';
            throw;
        }
    }

    void sparse_solver::diskload_effects()
    {
        try
        {
            if ( z_uni.empty() )
                throw std::string("z_uni is empty!.");

            if (z_on_memory)
            {
                for (auto &e : z_uni)
                {
                    for (auto &p : e)
                    {
                        // !!! revise, we do not need to write-out data again, only the momory needs to be cleared
                        p.fwrite();
                    }
                }
                z_on_memory = false;
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_solver::diskload_effects(): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::diskload_effects()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::diskload_effects()." << '\n';
            throw;
        }
    }

    void sparse_solver::construct_union_z(const matrix<int> &eff_ids, bool on_memory)
    {
        // Combines a consecutive sets (from left to right) of incidense
        // matrices of different (uncasted) types for a specific trait
        // (groups effects according to trait). All effects remain on disk or memory.

        // eff_ids is the indeces list of effects matrices associated to a specific trait

        try
        {
            if (on_memory) // load everything on memory
            {
                z_on_memory = true;

                std::vector<effects_storage> eff;

                effects_storage e = model.all_effects[eff_ids[0]];

                e.fread();

                eff.push_back(e); // very first effect

                e.clear();

                for (size_t i = 1; i < eff_ids.size(); i++)
                {
                    e = model.all_effects[eff_ids[i]];
                    e.fread();

                    eff.push_back(e);

                    e.clear();
                }

                z_uni.push_back(eff);
            }
            else // keep everything on disk
            {
                z_on_memory = false;

                std::vector<effects_storage> eff;

                eff.push_back(model.all_effects[eff_ids[0]]); // very first effect

                for (size_t i = 1; i < eff_ids.size(); i++)
                {
                    eff.push_back(model.all_effects[eff_ids[i]]);
                }

                z_uni.push_back(eff);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::construct_union_z(const matrix<int> &, bool)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::construct_union_z(const matrix<int> &, bool)." << '\n';
            throw;
        }
    }

    void sparse_solver::get_vect_z_uni2(const size_t &which_trait, const size_t &which_row, float **vect)
    {
        // Returns dense vector
        try
        {
            size_t which_eff_matr = 0; // will be changed inside z_row_to_uni_col(...)

            size_t irow[] = {0, n_obs[which_trait] - 1}; // this is constant for any row, due to the same amount of observations

            size_t icol[] = {0, 0}; // will be changed inside z_row_to_uni_col(...)

            z_row_to_uni_col(which_trait, which_row, icol, which_eff_matr);

            z_uni[which_trait][which_eff_matr].get_fcast(vect, icol, irow); // here is row vector extracted
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::get_vect_z_uni2(const size_t &, const size_t &, float **)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::get_vect_z_uni2(const size_t &, const size_t &, float **)." << '\n';
            throw;
        }
    }

    void sparse_solver::get_vect_z_uni2(const size_t &which_trait, const size_t &which_row, std::vector<float> &values, std::vector<size_t> &keys)
    {
        // Returns sparse vector (keys & values)
        try
        {
            size_t which_eff_matr = 0; // will be changed inside z_row_to_uni_col(...)

            size_t irow[] = {0, n_obs[which_trait] - 1}; // this is constant for any row, due to the same amount of observations

            size_t icol[] = {0, 0}; // will be changed inside z_row_to_uni_col(...)

            z_row_to_uni_col(which_trait, which_row, icol, which_eff_matr);

            z_uni[which_trait][which_eff_matr].get_fcast(values, keys, icol, irow); // here is row vector extracted
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::get_vect_z_uni2(const size_t &, const size_t &, std::vector<float> &, std::vector &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::get_vect_z_uni2(const size_t &, const size_t &, std::vector<float> &, std::vector<size_t> &)." << '\n';
            throw;
        }
    }

    void sparse_solver::z_row_to_uni_col(const size_t &which_trait, const size_t &in_z_row, size_t *out_uni_col, size_t &out_uni_matr)
    {
        // Returns the specific column and the index of an incidense (effects) matrix.

        try
        {
            std::vector<size_t> levels;

            size_t sum_levels = 0;
            levels.push_back(sum_levels);

            for (size_t i = 0; i < n_lev[which_trait].size(); i++)
            {
                sum_levels = sum_levels + n_lev[which_trait][i];
                levels.push_back(sum_levels);

                if (sum_levels - 1 >= in_z_row)
                {
                    out_uni_col[0] = out_uni_col[1] = in_z_row - levels[i];
                    out_uni_matr = i;
                    break;
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::z_row_to_uni_col(const size_t &, const size_t &, size_t *, size_t &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::z_row_to_uni_col(const size_t &, const size_t &, size_t *, size_t &)." << '\n';
            throw;
        }
    }

    void sparse_solver::set_y(const int obs_id)
    {
        try
        {
            matrix<float> yi = model.observations[obs_id].fget();

            model.observations[obs_id].fclear();
            model.observations[obs_id].clear();

            // Making missing data structure
            std::vector<bool> s_miss;
            matrix<size_t> shp = yi.shape();

            for (size_t i = 0; i < shp[0]; i++)
            {
                if (yi(i, 0) == model.missing_constant)
                {
                    s_miss.push_back(false);
                    yi(i, 0) = 0.0; // is this OK? Otherwise we'll have -999.0 in this place!
                }
                else
                    s_miss.push_back(true);
            }

            s.push_back(s_miss);

            yi.fwrite();

            y.push_back(yi);

            yi.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::set_y(const int)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::set_y(const int)." << '\n';
            throw;
        }
    }

    void sparse_solver::set_r()
    {
        try
        {
            size_t k = (n_trait * n_trait + n_trait) / 2.0; // elements in lower triangular part

            matrix<float> R = model.residuals[0].fget(); // Note! This is only for the very first submitted residual matrix.

            if (R.size() < k)
                throw std::string("The number of elements in the provided residual covariance matrix is not correspond to the number of traitsin the model!");

            std::vector<float> r_init(k, 0.0); // inverse of residual covariance for the  trivial case when all traits are unrelated
            for (size_t i = 0; i < n_trait; i++)
                for (size_t j = 0; j <= i; j++)
                    r_init[i * (i + 1) / 2 + j] = 1.0 / R(i, j); // using lower triangular part

            for (size_t i = 0; i < n_obs[0]; i++) // create hash list for each row of observations
            {
                std::vector<bool> vec_bool;
                for (size_t j = 0; j < n_trait; j++) // loop over all trait for the specific observation i
                    vec_bool.push_back(s[j][i]);

                std::hash<std::vector<bool>> hash_vector_bool;

                R_hash.push_back(hash_vector_bool(vec_bool));
            }

            std::vector<size_t> hash_keys; // unique list of hash keys
            hash_keys = R_hash;

            sort(hash_keys.begin(), hash_keys.end()); // sort and getting unique hesh keys
            hash_keys.erase(unique(hash_keys.begin(), hash_keys.end()), hash_keys.end());

            for (auto const &v : hash_keys) // initializing r_map, and modify its value's specific elements based on unique patterns
            {
                r_map[v] = r_init; // initialize by the trivial case covariance (all unrelated)

                // by looping over unique keys, get info about the specific pattern
                // find the details of the pattern by using its hash key
                size_t index_in_obs = 0;

                std::vector<size_t>::iterator it = find(R_hash.begin(), R_hash.end(), v);

                if (it != R_hash.end())
                    index_in_obs = it - R_hash.begin();
                else
                    throw std::string("Cannot find the hask key in the R_hash vector!");

                std::vector<size_t> pattern; // holds indexes needed to construct the reduced covariance matrix

                for (size_t i3 = 0; i3 < n_trait; i3++)
                    if (s[i3][index_in_obs])   // select only those indexes which have observations (true values)
                        pattern.push_back(i3); // which index should be used to construct and invert matrix

                if (pattern.size() < 2) // if the pattern is 0 everywhere or only 1 is non-missing -> is trivial case
                    continue;           // do nothing, switch to the next key

                evolm::matrix<float> iR(pattern.size(), pattern.size()); // creating reduced covariance matrix

                for (size_t i3 = 0; i3 < pattern.size(); i3++)
                    for (size_t j3 = 0; j3 <= i3; j3++)
                        iR(i3, j3) = iR(j3, i3) = R(pattern[i3], pattern[j3]); // build matrix

                iR.invert();

                for (size_t i3 = 0; i3 < pattern.size(); i3++)
                    for (size_t j3 = 0; j3 <= i3; j3++)
                        r_map[v][pattern[i3] * (pattern[i3] + 1) / 2.0 + pattern[j3]] = iR(i3, j3); // using lower triangulaar indexation
            }

            // NOTE: uncoment this in case of accepting constant Rij(-1)
            //       for all observations (not correcting for missing observations)
            // r = R;
            // r.invert();
            // End of the NOTE.

            R.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::set_r()" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_solver::set_r()" << '\n';
            std::cerr << "Reason => " << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::set_r()" << '\n';
            throw;
        }
    }

    void sparse_solver::set_g()
    {
        try
        {
            for (size_t i = 0; i < model.variances.size(); i++)
            {
                matrix<float> var;
                var = model.variances[i].fget();

                model.variances[i].fclear();
                model.variances[i].clear();

                var.invert();
                //var.fwrite(); keep variances on memory
                model.variances[i] = var;

                var.clear();
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::set_g()" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::set_g()" << '\n';
            throw;
        }
    }

    size_t sparse_solver::corr_size()
    {
        try
        {
            return model.correlations.size();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::corr_size()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::corr_size()." << '\n';
            throw;
        }
    }

    matrix<int> sparse_solver::get_corr_effects(size_t which_correlation)
    {
        try
        {
            matrix<int> which_effects = model.correlated_effects[which_correlation];

            //which_effects.fread(); loaded on memory
            matrix<size_t> shape_eff = which_effects.shape();

            for (size_t i = 0; i < shape_eff[0]; i++)
                which_effects(i, 0) = (int)adj_effects_order[which_effects(i, 0)];

            return which_effects;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::get_corr_effects(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::get_corr_effects(size_t)." << '\n';
            throw;
        }
    }

    float sparse_solver::get_variance(size_t which_correlation, size_t row, size_t col)
    {
        try
        {
            return model.variances[which_correlation](row, col);
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_solver::get_variance(size_t, size_t, size_t): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::get_variance(size_t, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::get_variance(size_t, size_t, size_t)." << '\n';
            throw;
        }
    }

    void sparse_solver::add_correlation(std::vector<float> &vect_to_add, size_t vect_first_index, float variance, size_t which_trait, size_t which_row, size_t col_1, size_t col_2)
    {
        matrix<float> cor_out;

        try
        {
            size_t icol[2] = {col_1, col_2};

            model.correlations[which_trait].add_row_to_vector(vect_to_add, vect_first_index, variance, which_row, icol);
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_solver::get_correlation(matrix<float> &, size_t, float, size_t, size_t, size_t, size_t): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::get_correlation(matrix<float> &, size_t, floatsize_t, size_t, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::get_correlation(matrix<float> &, size_t, floatsize_t, size_t, size_t, size_t)." << '\n';
            throw;
        }
    }

    matrix<float> sparse_solver::get_correlation(size_t which_trait, size_t which_row, size_t col_1, size_t col_2)
    {
        matrix<float> cor_out;

        try
        {
            size_t irow[2] = {which_row, which_row};
            size_t icol[2] = {col_1, col_2};

            matrix<float> cor;

            model.correlations[which_trait].to_dense(cor, irow, icol);

            return cor;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_solver::get_correlation(size_t, size_t, size_t, size_t): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::get_correlation(size_t, size_t, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::get_correlation(size_t, size_t, size_t, size_t)." << '\n';
            throw;
        }

        return cor_out;
    }

    void sparse_solver::memload_var()
    {
        try
        {
            if (model.variances.size() == 0)
                throw std::string("Vector of variances is empty. In sparse_solver::memload_var()");

            if (!var_onmem)
            {
                for (size_t i = 0; i < model.variances.size(); i++)
                    model.variances[i].fread();

                var_onmem = true;
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_solver::memload_var(): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::memload_var()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::memload_var()" << '\n';
            throw;
        }
    }

    void sparse_solver::diskload_var()
    {
        try
        {
            if (model.variances.size() == 0)
                throw std::string("Vector of variances is empty. In sparse_solver::memload_var()");

            if (var_onmem)
            {
                for (size_t i = 0; i < model.variances.size(); i++)
                    model.variances[i].fwrite();

                var_onmem = false;
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_solver::diskload_var(): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::diskload_var()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::diskload_var()" << '\n';
            throw;
        }
    }

    void sparse_solver::memload_cor()
    {
        try
        {
            if (model.correlations.size() == 0)
                throw std::string("Vector of correlations is empty. In sparse_solver::memload_cor()");

            if (!cor_onmem)
            {
                for (size_t i = 0; i < model.correlations.size(); i++)
                {
                    model.correlations[i].fread();
                    model.correlations[i].fread_rows_structure();
                }

                cor_onmem = true;
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_solver::memload_cor(): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::memload_cor()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::memload_cor()" << '\n';
            throw;
        }
    }

    void sparse_solver::diskload_cor()
    {
        try
        {
            if (model.correlations.size() == 0)
                throw std::string("Vector of correlations is empty. In sparse_solver::memload_cor()");

            if (cor_onmem)
            {
                for (size_t i = 0; i < model.correlations.size(); i++)
                {
                    // !!! revise here: we do not need to write out the data again, only clearing the memory is needed
                    model.correlations[i].fwrite();
                    model.correlations[i].fwrite_rows_structure();
                }

                cor_onmem = false;
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_solver::diskload_cor(): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::diskload_cor()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::diskload_cor()" << '\n';
            throw;
        }
    }

    void sparse_solver::memload_cor_effects()
    {
        try
        {
            if (model.correlated_effects.size() == 0)
                throw std::string("Vector of cor_effects is empty. In sparse_solver::memload_cor_effects()");

            for (size_t i = 0; i < model.correlated_effects.size(); i++)
                model.correlated_effects[i].fread();
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_solver::memload_cor_effects(): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::memload_cor_effects()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::memload_cor_effects()" << '\n';
            throw;
        }
    }

    void sparse_solver::diskload_cor_effects()
    {
        try
        {
            if (model.correlated_effects.size() == 0)
                throw std::string("Vector of cor_effects is empty. In sparse_solver::diskload_cor_effects()");

            for (size_t i = 0; i < model.correlated_effects.size(); i++)
                model.correlated_effects[i].fwrite();
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in sparse_solver::diskload_cor_effects(): " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::diskload_cor_effects()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::diskload_cor_effects()" << '\n';
            throw;
        }
    }

    bool sparse_solver::identity_correlation(size_t which_corr)
    {
        try
        {
            return model.identity_correlations[which_corr];
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::identity_correlation(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::identity_correlation(size_t)." << '\n';
            throw;
        }
    }

#ifdef UTEST
    matrix<float> sparse_solver::test_z_dot_y(size_t vect_size, size_t i_trait, size_t j_trait, size_t r_index)
    {
        try
        {
            matrix<float> out_vect;
            z_dot_y(out_vect, vect_size, i_trait, j_trait, r_index);
            return out_vect;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::test_z_dot_y(size_t, size_t, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::test_z_dot_y(size_t, size_t, size_t, size_t)." << '\n';
            throw;
        }
    }

    matrix<float> sparse_solver::test_rhs()
    {

        try
        {
            memload_effects();

            construct_rhs();

            return rhs;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::test_rhs()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::test_rhs()." << '\n';
            throw;
        }
    }

    size_t sparse_solver::test_num_all_levels()
    {
        try
        {
            size_t n_all_levels = num_all_levels();

            return n_all_levels;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::test_num_all_levels()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::test_num_all_levels()." << '\n';
            throw;
        }
    }

    std::vector<size_t> sparse_solver::test_ordered_levels()
    {
        try
        {
            std::vector<size_t> ordered_random_levels = get_ordered_levels();

            return ordered_random_levels;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::test_ordered_levels()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::test_ordered_levels()." << '\n';
            throw;
        }
    }

    std::vector<std::vector<size_t>> sparse_solver::test_cov_offsets()
    {
        try
        {
            std::vector<size_t> ordered_random_levels = get_ordered_levels();

            std::vector<std::vector<size_t>> rcov_offsets = get_cov_offsets(ordered_random_levels);

            return rcov_offsets;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::test_cov_offsets()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::test_cov_offsets()." << '\n';
            throw;
        }
    }

    void sparse_solver::test_vect_z_uni(const size_t &which_trait, matrix<float> &out)
    {
        try
        {
            size_t levels = 0;
            for (size_t k = 0; k < n_lev[which_trait].size(); k++)
                levels = levels + n_lev[which_trait][k];

            memload_effects();

            matrix<float> v(levels, n_obs[which_trait]);

            for (size_t j = 0; j < levels; j++)
            {
                std::vector<float> vals;
                std::vector<size_t> keys;
                get_vect_z_uni2(which_trait, j, vals, keys);

                for (size_t i = 0; i < vals.size(); i++)
                    v(j, keys[i]) = vals[i];
            }

            out = v;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::test_vect_z_uni(const size_t &, matrix<float> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::test_vect_z_uni(const size_t &, matrix<float> &)." << '\n';
            throw;
        }
    }
    std::vector<std::vector<float>> sparse_solver::test_A()
    {
        try
        {
            std::vector<std::vector<float>> A;
            memload_cor_effects();
            memload_cor();

            set_model_matrix();

            diskload_cor_effects();
            diskload_cor();

            //bin_fnames.clear();
            //blocks_ranges.clear();

            matrix<float> M(model_matrix.size(), model_matrix.size());

            for (size_t i = 0; i < first_row_ondisk; i++)
            {
                matrix<float> M0;
                model_matrix[i].to_dense(M0);
                for (size_t j = 0; j < M0.size(); j++)
                    M(i,j) = M0[j];
            }

            /*for (size_t i = 0; i < get_num_of_mem_blocks(); i++) // variant 1, on blocks with consecutive reading
            {
                size_t first_row = 0;
                size_t last_row = 0;

                load_model_matrix(i);
                get_mem_block_range(i, first_row, last_row);

                for (size_t i2 = first_row; i2 <= last_row; i2++)
                {
                    matrix<float> M0;
                    model_matrix[i2].to_dense(M0);
                    for (size_t j = 0; j < M0.size(); j++)
                        M(i2,j) = M0[j];
                }

                unload_model_matrix(i);
            }*/
            
            for (size_t i = first_row_ondisk; i < model_matrix.size(); i++) // variant 2
            {
                std::fstream fA;
                fA.open(bin_fname, fA.binary | fA.in);
                fA.seekg (model_matrix[i].bin_file_read_position, fA.beg);
                model_matrix[i].fread(fA);

                matrix<float> M0;
                model_matrix[i].to_dense(M0);
                for (size_t j = 0; j < M0.size(); j++)
                    M(i,j) = M0[j];

                fA.close();
            }

            M.to_vector2d(A);

            return A;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in sparse_solver::test_A()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in sparse_solver::test_A()." << '\n';
            throw;
        }
    }

#endif

}
