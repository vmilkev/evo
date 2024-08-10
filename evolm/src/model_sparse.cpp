#include "model_sparse.hpp"

namespace evolm
{
    model_sparse::model_sparse()
    {
    }
    //===============================================================================================================
    void model_sparse::set_missing(float val)
    {
        try
        {
            missing_constant = val;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::set_missing()" << '\n';
            std::cerr << "Reason: " << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_sparse::set_missing()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::set_missing()" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    int model_sparse::append_residual(const std::vector<float> &arr, size_t lda)
    {
        try
        {
            // lda := is leading diagonal of symmetric matrix, though arr is in a full-store format
            matrix<float> residual;

            residual.resize(lda, lda);

            for (size_t i = 0; i < lda * lda; i++)
                residual[i] = arr[i];

            residual.fwrite();

            residuals.push_back(residual);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::append_residual(const std::vector<float> &, size_t)." << '\n';
            std::cerr << "Reason: " << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_sparse::append_residual(const std::vector<float> &, size_t)." << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::append_residual(const std::vector<float> &, size_t)." << '\n';
            throw;
        }
        return 0;
    }
    //===============================================================================================================
    int model_sparse::append_residual(const std::string &fname)
    {
        try
        {
            IOInterface datstream;

            datstream.set_fname(fname);

            std::vector<std::vector<float>> arr;

            datstream.fgetdata(arr);

            size_t lda = arr.size();

            if (lda == 0)
                throw "The first dimension of the arrey is 0, expected et least 1!";

            if (arr[0].size() == 0)
                throw "The second dimension of the arrey is 0, expected et least 1!";

            if (lda != arr[0].size())
                throw "The arrey is not square, that is expected in this case!";

            matrix<float> residual;

            residual.resize(lda, lda);
            // residual.rectosym();

            for (size_t i = 0; i < lda; i++)
            {
                for (size_t j = 0; j < lda; j++)
                    residual(i, j) = arr[i][j];
            }

            residual.fwrite();
            residuals.push_back(residual);

            arr.clear();
            arr.shrink_to_fit();
        }
        catch (std::string err)
        {
            std::cerr << "Exception in model_sparse::append_residual(const std::string &)." << '\n';
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::append_residual(const std::string &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::append_residual(const std::string &)." << '\n';
            throw;
        }

        return 0;
    }
    //===============================================================================================================
    int model_sparse::append_observation(const std::vector<float> &arr, size_t lda)
    {
        try
        {
            if (arr.size() != lda)
                throw std::string("Provided vector size does not correspond to the number of elements in the observations vector!");

            matrix<float> observation;

            // lda := is size of vector

            observation.resize(lda, 1);

            for (size_t i = 0; i < lda; i++)
                observation[i] = arr[i];

            size_of_data = size_of_data + observation.size() * sizeof(float);

            observation.fwrite();

            observations.push_back(observation);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::append_observation(const std::vector<float> &, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_sparse::append_observation(const std::vector<float> &, size_t)." << '\n';
            std::cerr << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::append_observation(const std::vector<float> &, size_t)." << '\n';
            throw;
        }

        return 0;
    }
    //===============================================================================================================
    int model_sparse::append_observation(const std::vector<float> &arr, const std::vector<bool> &miss_arr, size_t lda)
    {
        try
        {
            /*size_t zeros = count(miss_arr.begin(), miss_arr.end(), 0);
            size_t non_zeros = miss_arr.size() - zeros;

            if ( (arr.size() != lda) || (non_zeros != lda)  )
                throw std::string("Provided observations  vector size does not correspond to the number of elements in the observations vector!");
            */
            matrix<float> observation;

            // lda := is size of vector

            observation.resize(lda, 1);

            for (size_t i = 0; i < lda; i++)
                observation[i] = arr[i];

            size_of_data = size_of_data + observation.size() * sizeof(float);

            observation.fwrite();
            observations.push_back(observation);

            miss_observations.push_back(miss_arr);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::append_observation(const std::vector<float> &, const std::vector<bool> &, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_sparse::append_observation(const std::vector<float> &, const std::vector<bool> &, size_t)." << '\n';
            std::cerr << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::append_observation(const std::vector<float> &, const std::vector<bool> &, size_t)." << '\n';
            throw;
        }

        return 0;
    }
    //===============================================================================================================
    int model_sparse::append_observation(const std::string &fname)
    {
        try
        {
            IOInterface datstream;

            datstream.set_fname(fname);

            std::vector<std::vector<float>> arr;

            datstream.fgetdata(arr);

            size_t lda = arr.size();

            if (lda == 0)
                throw "The first dimension of the arrey is 0, expected et least 1!";

            matrix<float> observation;

            // lda := is size of vector

            observation.resize(lda, 1);

            for (size_t i = 0; i < lda; i++)
                observation[i] = arr[i][0];

            size_of_data = size_of_data + observation.size() * sizeof(float);

            observation.fwrite();
            observations.push_back(observation);

            arr.clear();
            arr.shrink_to_fit();
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_sparse::append_observation(const std::string &)." << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::append_observation(const std::string &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::append_observation(const std::string &)." << '\n';
            throw;
        }

        return 0;
    }
    //===============================================================================================================
    template <typename T>
    void model_sparse::append_effect(compact_storage<T> &eff)
    {
        try
        {
            eff.transpose();

            eff.make_rows_list();

            size_of_data = size_of_data + eff.size_inmem();

            eff.fwrite_rows_structure();

            eff.fwrite(); // keys and values

            effects_storage e;

            e.set(eff);

            all_effects.push_back(e);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::append_effect(compact_storage<T> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_sparse::append_effect(compact_storage<T> &)." << '\n';
            std::cerr << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::append_effect(compact_storage<T> &)." << '\n';
            throw;
        }
    }
    template void model_sparse::append_effect(compact_storage<int> &eff);
    template void model_sparse::append_effect(compact_storage<float> &eff);
    template void model_sparse::append_effect(compact_storage<double> &eff);
    //===============================================================================================================
    template <typename T>
    void model_sparse::append_effect(std::vector<T> &values, size_t n_rows, size_t n_cols)
    {
        try
        {
            if (n_rows == 0 || n_cols == 0)
                throw std::string("Provided number of rows/columns is zerro!");

            if ( values.empty() )
                throw std::string("Provided all_effects vector is empty!");

            compact_storage<T> effect;

            effect.resize(n_rows, n_cols);

            effect.set_sparsity_threshold(sparsity_threshold);

            effect.append(values);

            effect.transpose();

            effect.make_rows_list();

            size_of_data = size_of_data + effect.size_inmem();

            effect.fwrite_rows_structure();

            effect.fwrite(); // keys and values

            effects_storage e;

            e.set(effect);

            all_effects.push_back(e);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::append_effect(const std::vector<T> &, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_sparse::append_effect(const std::vector<T> &, size_t, size_t)." << '\n';
            std::cerr << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::append_effect(const std::vector<T> &, size_t, size_t)." << '\n';
            throw;
        }
    }
    template void model_sparse::append_effect(std::vector<int> &values, size_t n_rows, size_t n_cols);
    template void model_sparse::append_effect(std::vector<float> &values, size_t n_rows, size_t n_cols);
    template void model_sparse::append_effect(std::vector<double> &values, size_t n_rows, size_t n_cols);
    //===============================================================================================================
    template <typename T>
    void model_sparse::append_effect(std::vector<T> &values, std::vector<size_t> &rows, std::vector<size_t> &cols, size_t n_rows, size_t n_cols)
    {
        try
        {
            if (n_rows == 0 || n_cols == 0)
                throw std::string("Provided number of rows/columns is zerro!");

            if ( (rows.size() == 0 && cols.size() != 0) || (rows.size() != 0 && cols.size() == 0))
                throw std::string("One of provided vectors of rows/cols is zerro!");

            compact_storage<T> effect;

            effect.resize(n_rows, n_cols);

            effect.set_sparsity_threshold(sparsity_threshold);

            if ( rows.empty() )
                effect.append(values);
            else
                effect.append(values, rows, cols);

            effect.transpose();

            effect.make_rows_list();

            size_of_data = size_of_data + effect.size_inmem();

            effect.fwrite_rows_structure();

            effect.fwrite(); // keys and values

            effects_storage e;

            e.set(effect);

            all_effects.push_back(e);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::append_effect(const std::vector<T> &, const std::vector<size_t> &, const std::vector<size_t> &, size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_sparse::append_effect(const std::vector<T> &, const std::vector<size_t> &, const std::vector<size_t> &, size_t, size_t)." << '\n';
            std::cerr << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::append_effect(const std::vector<T> &, const std::vector<size_t> &, const std::vector<size_t> &, size_t, size_t)." << '\n';
            throw;
        }
    }
    template void model_sparse::append_effect(std::vector<double> &values, std::vector<size_t> &rows, std::vector<size_t> &cols, size_t n_rows, size_t n_cols);
    template void model_sparse::append_effect(std::vector<float> &values, std::vector<size_t> &rows, std::vector<size_t> &cols, size_t n_rows, size_t n_cols);
    template void model_sparse::append_effect(std::vector<int> &values, std::vector<size_t> &rows, std::vector<size_t> &cols, size_t n_rows, size_t n_cols);
    //===============================================================================================================
    void model_sparse::append_effect(const std::string &fname)
    {
        try
        {
            std::fstream fA;
            size_t data_type = 3003;
            size_t B[6];

            fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            fA.open(fname, fA.binary | fA.in);

            if (!fA.is_open())
                throw std::string("Error while opening a binary file!");

            fA.read( reinterpret_cast<char *>(&B), 6 * sizeof(size_t) ); // 0. reading a storage info

            if (  B[0] != data_type )
                throw std::string("Trying to read a wrong kind of storage from file into the compact_storage format!");

            size_t var_type = B[5];

            fA.close();

            compact_storage<int> i_effect;
            compact_storage<float> f_effect;
            compact_storage<double> d_effect;

            i_effect.set_sparsity_threshold(sparsity_threshold);
            f_effect.set_sparsity_threshold(sparsity_threshold);
            d_effect.set_sparsity_threshold(sparsity_threshold);

            effects_storage e;

            switch (var_type)
            {
            case 1:
                i_effect.fread(fname);
                i_effect.transpose();
                i_effect.make_rows_list();
                size_of_data = size_of_data + i_effect.size_inmem();
                i_effect.fwrite_rows_structure();
                i_effect.fwrite(); // keys and values
                e.set(i_effect);
                all_effects.push_back(e);
                break;
            case 2:
                f_effect.fread(fname);
                f_effect.transpose();
                f_effect.make_rows_list();
                size_of_data = size_of_data + f_effect.size_inmem();
                f_effect.fwrite_rows_structure();
                f_effect.fwrite();  // keys and values
                e.set(f_effect);
                all_effects.push_back(e);
                break;
            case 3:
                d_effect.fread(fname);
                d_effect.transpose();
                d_effect.make_rows_list();
                size_of_data = size_of_data + d_effect.size_inmem();
                d_effect.fwrite_rows_structure();
                d_effect.fwrite();  // keys and values
                e.set(d_effect);
                all_effects.push_back(e);
                break;            
            default:
                throw std::string("Cannot determine the type of data in the binary file.");
            }
        }
        catch (std::string err)
        {
            std::cerr << "Exception in model_sparse::append_effect(const std::string &)." << '\n';
            std::cerr << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::append_effect(const std::string &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::append_effect(const std::string &)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void model_sparse::append_corrstruct(const std::vector<float> &var, size_t lda1, std::vector<float> &corr, size_t lda2, const std::vector<int> &which_effects)
    {
        try
        {
            matrix<float> variance;
            compact_storage<float> correlation;
            matrix<int> _effects;

            _effects.resize(which_effects.size(), 1);
            variance.resize(lda1, lda1);
            correlation.resize(lda2, lda2);

            if (lda1 != which_effects.size())
                throw std::string("The number of provided correlated all_effects does not correspond to the dimension of variance-covariance matrix!");

            for (size_t i = 0; i < which_effects.size(); i++)
                _effects[i] = which_effects[i];

            for (size_t i = 0; i < var.size(); i++)
                variance[i] = var[i];

            correlation.append(corr);

            variance.fwrite();

            variances.push_back(variance);

            correlation.make_rows_list();
            size_of_data = size_of_data + correlation.size_inmem();
            correlation.fwrite_rows_structure();

            correlation.fwrite();
            correlations.push_back(correlation);

            identity_correlations.push_back(false);
            identity_dimension.push_back(0);

            _effects.fwrite();
            correlated_effects.push_back(_effects);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::append_corrstruct(const std::vector<float> &, size_t, std::vector<float> &, size_t, const std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (std::string &err)
        {
            std::cerr << "Exception in model_sparse::append_corrstruct(const std::vector<float> &, size_t, std::vector<float> &, size_t, const std::vector<int> &)." << '\n';
            std::cerr << "Reason => " << err << '\n';
            throw err;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::append_corrstruct(const std::vector<float> &, size_t, std::vector<float> &, size_t, const std::vector<int> &)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void model_sparse::append_corrstruct(const std::vector<float> var, size_t lda1, std::string &identity, size_t lda2, const std::vector<int> which_effects)
    {
        try
        {
            std::string identity_string("I");

            if (identity_string.compare(identity) != 0)
                throw std::runtime_error("Error in the third input argument, the expected is I");

            matrix<float> variance;
            compact_storage<float> correlation;
            matrix<int> _effects;

            _effects.resize(which_effects.size(), 1);

            variance.resize(lda1, lda1);

            correlation.resize(1, 1);

            if (lda1 != which_effects.size())
                throw std::string("The number of provided correlated all_effects does not correspond to the dimension of variance-covariance matrix!");

            for (size_t i = 0; i < which_effects.size(); i++)
                _effects[i] = which_effects[i];

            for (size_t i = 0; i < var.size(); i++)
                variance[i] = var[i];

            std::vector<float> t_vect;
            t_vect.push_back(1.0f);
            correlation.append(t_vect);

            variance.fwrite();

            variances.push_back(variance);

            correlation.make_rows_list();
            size_of_data = size_of_data + correlation.size_inmem();
            correlation.fwrite_rows_structure();

            correlation.fwrite();
            correlations.push_back(correlation);

            identity_correlations.push_back(true);
            identity_dimension.push_back(lda2);

            _effects.fwrite();
            correlated_effects.push_back(_effects);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::append_corrstruct(const std::vector<float>, size_t, std::string &, size_t, const std::vector<int>)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (std::string &err)
        {
            std::cerr << "Exception in model_sparse::append_corrstruct(const std::vector<float>, size_t, std::string &, size_t, const std::vector<int>)." << '\n';
            std::cerr << "Reason => " << err << '\n';
            throw err;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::append_corrstruct(const std::vector<float>, size_t, std::string &, size_t, const std::vector<int>)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void model_sparse::append_corrstruct(const std::string &fname_var, std::string &identity, size_t lda2, const std::vector<int> which_effects)
    {
        try
        {
            std::string identity_string("I");

            if (identity_string.compare(identity) != 0)
                throw std::runtime_error("Error in the third input argument, the expected is I");

            IOInterface datstream;

            datstream.set_fname(fname_var);
            std::vector<std::vector<float>> var;
            datstream.fgetdata(var);

            size_t lda1 = var.size();

            if (lda2 == 0)
                throw "The first dimension of the arrey CORR is 0, expected et least 1!";

            if (lda1 == 0)
                throw "The first dimension of the arrey VAR is 0, expected et least 1!";

            if (var[0].size() == 0)
                throw "The second dimension of the arrey VAR is 0, expected et least 1!";

            if (lda1 != var[0].size())
                throw "The arrey VAR is not square, that is expected in this case!";

            matrix<float> variance;
            compact_storage<float> correlation;
            matrix<int> _effects;

            _effects.resize(which_effects.size(), 1);

            variance.resize(lda1, lda1);

            correlation.resize(1, 1);

            if (lda1 != which_effects.size())
                throw std::string("The number of provided correlated all_effects does not correspond to the dimension of variance-covariance matrix!");

            for (size_t i = 0; i < which_effects.size(); i++)
                _effects[i] = which_effects[i];

            for (size_t i = 0; i < lda1; i++)
            {
                for (size_t j = 0; j < lda1; j++)
                    variance(i, j) = var[i][j];
            }

            std::vector<float> t_vect;
            t_vect.push_back(1.0f);
            correlation.append(t_vect);

            variance.fwrite();

            variances.push_back(variance);

            correlation.make_rows_list();
            size_of_data = size_of_data + correlation.size_inmem();
            correlation.fwrite_rows_structure();

            correlation.fwrite();
            correlations.push_back(correlation);

            identity_correlations.push_back(true);
            identity_dimension.push_back(lda2);

            _effects.fwrite();
            correlated_effects.push_back(_effects);

            var.clear();
            var.shrink_to_fit();
        }
        catch (std::string &err)
        {
            std::cerr << "Exception in model_sparse::append_corrstruct(const std::string &, std::string &, size_t, const std::vector<int>)." << '\n';
            std::cerr << "Reason => " << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::append_corrstruct(const std::string &, std::string &, size_t, const std::vector<int>)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::append_corrstruct(const std::string &, std::string &, size_t, const std::vector<int>)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void model_sparse::append_corrstruct(const std::vector<float> &var, size_t lda1, const std::string &fname_corr, const std::vector<int> &which_effects)
    {
        try
        {
            IOInterface datstream;

            matrix<float> variance;
            compact_storage<float> correlation;
            matrix<int> _effects;

            _effects.resize(which_effects.size(), 1);
            variance.resize(lda1, lda1);

            if (lda1 != which_effects.size())
                throw std::string("The number of provided correlated all_effects does not correspond to the dimension of variance-covariance matrix!");

            for (size_t i = 0; i < which_effects.size(); i++)
                _effects[i] = which_effects[i];

            for (size_t i = 0; i < var.size(); i++)
                variance[i] = var[i];

            correlation.fread(fname_corr);

            if ( correlation.nrows() == 0 )
                throw "The first dimension of the arrey CORR is 0, expected et least 1!";

            if ( correlation.ncols() == 0 )
                throw "The second dimension of the arrey CORR is 0, expected et least 1!";

            if ( correlation.nrows() != correlation.ncols() )
                throw "The arrey CORR is not square, that is expected in this case!";

            variance.fwrite();

            variances.push_back(variance);

            correlation.make_rows_list();
            size_of_data = size_of_data + correlation.size_inmem();
            correlation.fwrite_rows_structure();

            correlation.fwrite();
            correlations.push_back(correlation);

            identity_correlations.push_back(false);
            identity_dimension.push_back(0);

            _effects.fwrite();
            correlated_effects.push_back(_effects);
        }
        catch (std::string &err)
        {
            std::cerr << "Exception in model_sparse::append_corrstruct(const std::vector<float> &, size_t, const std::string &, const std::vector<int> &)." << '\n';
            std::cerr << "Reason => " << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::append_corrstruct(const std::vector<float> &, size_t, const std::string &, const std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::append_corrstruct(const std::vector<float> &, size_t, const std::string &, const std::vector<int> &)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void model_sparse::append_corrstruct(const std::string &fname_var, const std::string &fname_corr, const std::vector<int> &which_effects)
    {
        try
        {
            IOInterface datstream;

            datstream.set_fname(fname_var);
            std::vector<std::vector<float>> var;
            datstream.fgetdata(var);

            //size_t lda2 = corr.size();
            size_t lda1 = var.size();

            if (lda1 == 0)
                throw "The first dimension of the arrey VAR is 0, expected et least 1!";

            if (var[0].size() == 0)
                throw "The second dimension of the arrey VAR is 0, expected et least 1!";

            if (lda1 != var[0].size())
                throw "The arrey VAR is not square, that is expected in this case!";

            matrix<float> variance;
            compact_storage<float> correlation;
            matrix<int> _effects;

            _effects.resize(which_effects.size(), 1);
            variance.resize(lda1, lda1);

            if (lda1 != which_effects.size())
                throw std::string("The number of provided correlated all_effects does not correspond to the dimension of variance-covariance matrix!");

            for (size_t i = 0; i < which_effects.size(); i++)
                _effects[i] = which_effects[i];

            for (size_t i = 0; i < lda1; i++)
            {
                for (size_t j = 0; j < lda1; j++)
                    variance(i, j) = var[i][j];
            }

            correlation.fread(fname_corr);

            if ( correlation.nrows() == 0 )
                throw "The first dimension of the arrey CORR is 0, expected et least 1!";

            if ( correlation.ncols() == 0 )
                throw "The second dimension of the arrey CORR is 0, expected et least 1!";

            if ( correlation.nrows() != correlation.ncols() )
                throw "The arrey CORR is not square, that is expected in this case!";

            variance.fwrite();

            variances.push_back(variance);

            correlation.make_rows_list();
            size_of_data = size_of_data + correlation.size_inmem();
            correlation.fwrite_rows_structure();

            correlation.fwrite();
            correlations.push_back(correlation);

            identity_correlations.push_back(false);
            identity_dimension.push_back(0);

            _effects.fwrite();
            correlated_effects.push_back(_effects);

            var.clear();
            var.shrink_to_fit();
        }
        catch (std::string &err)
        {
            std::cerr << "Exception in model_sparse::append_corrstruct(const std::string &, const std::string &, const std::vector<int> &)." << '\n';
            std::cerr << "Reason => " << err << '\n';
            throw err;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::append_corrstruct(const std::string &, const std::string &, const std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::append_corrstruct(const std::string &, const std::string &, const std::vector<int> &)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    int model_sparse::append_traitstruct(int obs_id, const std::vector<int> &eff_id)
    {
        try
        {
            matrix<int> eff;

            eff.resize(eff_id.size(), 1);

            for (size_t i = 0; i < eff_id.size(); i++)
                eff[i] = eff_id[i];

            observation_trait.push_back(obs_id);
            effects_trait.push_back(eff);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::append_traitstruct(int, const std::vector<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::append_traitstruct(int, const std::vector<int> &)." << '\n';
            throw;
        }

        return 0;
    }
    //===============================================================================================================
    void model_sparse::clear_residuals()
    {
        try
        {
            for (auto &e : residuals)
            {
                e.fclear();
                e.clear();
            }
            //residuals.clear();
            //residuals.shrink_to_fit();
            std::vector<matrix<float>>().swap(residuals);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::clear_residuals()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::clear_residuals()." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void model_sparse::clear_observations()
    {
        try
        {
            for (auto &e : observations)
            {
                e.fclear();
                e.clear();
            }
            //observations.clear();
            //observations.shrink_to_fit();
            std::vector<matrix<float>>().swap(observations);
            
            for (auto &e : miss_observations)
            {
                e.clear();
            }
            miss_observations.clear();
            miss_observations.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::clear_observations()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::clear_observations()." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void model_sparse::clear_effects()
    {
        try
        {
            for (auto &e : all_effects)
                e.clear();
            
            //all_effects.clear();
            //all_effects.shrink_to_fit();
            std::vector<effects_storage>().swap(all_effects);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::clear_effects()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::clear_effects()." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void model_sparse::clear_corrstruct()
    {
        try
        {
            for (auto &e : correlated_effects)
            {
                e.fclear();
                e.clear();
            }
            //correlated_effects.clear();
            //correlated_effects.shrink_to_fit();
            std::vector<matrix<int>>().swap(correlated_effects);

            for (auto &e : variances)
            {
                e.fclear();
                e.clear();
            }
            //variances.clear();
            //variances.shrink_to_fit();
            std::vector<matrix<float>>().swap(variances);

            for (auto &e : correlations)
            {
                e.fclear();
                e.clear();
            }
            //correlations.clear();
            //correlations.shrink_to_fit();
            std::vector<compact_storage<float>>().swap(correlations);

            identity_correlations.clear();
            identity_correlations.shrink_to_fit();
            identity_dimension.clear();
            identity_dimension.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::clear_corrstruct()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::clear_corrstruct()." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void model_sparse::clear_traitstruct()
    {
        try
        {
            //observation_trait.clear();
            //observation_trait.shrink_to_fit();
            //effects_trait.clear();
            //effects_trait.shrink_to_fit();
            std::vector<int>().swap(observation_trait);
            std::vector<matrix<int>>().swap(effects_trait);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::clear_traitstruct()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::clear_traitstruct()." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void model_sparse::clear()
    {
        try
        {
            clear_residuals();
            clear_observations();
            clear_effects();
            clear_corrstruct();
            clear_traitstruct();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::clear()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::clear()." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void model_sparse::set_sparsity_threshold(double threshold)
    {
        sparsity_threshold = threshold;
    }
    //===============================================================================================================
    size_t model_sparse::get_size_of_data()
    {
        return size_of_data;
    }
    //===============================================================================================================
#ifdef UTEST
    void model_sparse::print()
    {
        try
        {
            matrix<size_t> shape(2, 1);

            for (size_t i2 = 0; i2 < residuals.size(); i2++)
            {
                residuals[i2].fread();
                residuals[i2].print("residual");
                shape = residuals[i2].shape();
                residuals[i2].fwrite();
                shape.transpose();
                shape.print("Shape:");
                shape.clear();
            }
            for (size_t i2 = 0; i2 < observations.size(); i2++)
            {
                observations[i2].fread();
                observations[i2].print("observations");
                shape = observations[i2].shape();
                observations[i2].fwrite();
                shape.transpose();
                shape.print("Shape:");
                shape.clear();
            }
            for (size_t i2 = 0; i2 < all_effects.size(); i2++)
            {
                compact_storage<float> eff;
                matrix<float> M;
                std::vector<size_t> shp;

                all_effects[i2].fread();
                all_effects[i2].get(eff);
                eff.to_dense(M);
                M.print("model effect");

                shp = all_effects[i2].shape();
                all_effects[i2].fwrite();
                shape[0] = shp[0];
                shape[1] = shp[1];
                shape.transpose();
                shape.print("Shape:");
            }
            for (size_t i2 = 0; i2 < correlated_effects.size(); i2++)
            {
                correlated_effects[i2].fread();
                correlated_effects[i2].print("correlated_effects");
                shape = correlated_effects[i2].shape();
                correlated_effects[i2].fwrite();
                shape.transpose();
                shape.print("Shape:");
                shape.clear();
            }
            for (size_t i2 = 0; i2 < variances.size(); i2++)
            {
                variances[i2].fread();
                variances[i2].print("variances");
                shape = variances[i2].shape();
                variances[i2].fwrite();
                shape.transpose();
                shape.print("Shape:");
                shape.clear();
            }
            for (size_t i2 = 0; i2 < correlations.size(); i2++)
            {
                correlations[i2].fread();
                matrix<float> M;
                correlations[i2].to_dense(M);
                M.print("correlations");
                shape[0] = correlations[i2].nrows();
                shape[1] = correlations[i2].ncols();
                correlations[i2].fwrite();
                shape.transpose();
                shape.print("Shape:");
                shape.clear();
            }
            for (size_t i2 = 0; i2 < effects_trait.size(); i2++)
            {
                effects_trait[i2].print("effects_trait");
                shape = effects_trait[i2].shape();
                shape.transpose();
                shape.print("Shape:");
                shape.clear();
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::print()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::print()." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    size_t model_sparse::size_of(const std::string type)
    {
        try
        {
            if (type.compare("res") == 0)
            {
                size_t sz = 0;
                if (residuals.size() > 0)
                {
                    for (size_t i = 0; i < residuals.size(); i++)
                    {
                        sz = sz + residuals[i].size();
                    }
                }
                return sz;
            }

            if (type.compare("obs") == 0)
            {
                size_t sz = 0;
                if (observations.size() > 0)
                {
                    for (size_t i = 0; i < observations.size(); i++)
                        sz = sz + observations[i].size();
                }
                return sz;
            }

            if (type.compare("eff") == 0)
            {
                size_t sz = 0;
                if (all_effects.size() > 0)
                {
                    for (size_t i = 0; i < all_effects.size(); i++)
                    {
                        all_effects[i].fread();
                        sz = sz + all_effects[i].size();
                        all_effects[i].fwrite();
                    }
                }
                return sz;
            }

            if (type.compare("var") == 0)
            {
                size_t sz = 0;
                if (variances.size() > 0)
                {
                    for (size_t i = 0; i < variances.size(); i++)
                    {
                        sz = sz + variances[i].size();
                    }
                }
                return sz;
            }

            if (type.compare("cor") == 0)
            {
                size_t sz = 0;
                if (correlations.size() > 0)
                {
                    for (size_t i = 0; i < correlations.size(); i++)
                    {
                        correlations[i].fread();
                        sz = sz + correlations[i].size();
                        correlations[i].fwrite();
                    }
                }
                return sz;
            }

            if (type.compare("cor_eff") == 0)
            {
                size_t sz = 0;
                if (correlated_effects.size() > 0)
                {
                    for (size_t i = 0; i < correlated_effects.size(); i++)
                        sz = sz + correlated_effects[i].size();
                }
                return sz;
            }

            if (type.compare("obs_trt") == 0)
            {
                size_t sz = 0;
                if (observation_trait.size() > 0)
                {
                    for (size_t i = 0; i < observation_trait.size(); i++)
                        sz = sz + 1;
                }
                return sz;
            }

            if (type.compare("eff_trt") == 0)
            {
                size_t sz = 0;
                if (effects_trait.size() > 0)
                {
                    for (size_t i = 0; i < effects_trait.size(); i++)
                        sz = sz + effects_trait[i].size();
                }
                return sz;
            }

            return 0;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::size_of(const std::string)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::size_of(const std::string)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    std::vector<std::vector<size_t>> model_sparse::shape_of(const std::string type)
    {
        try
        {
            std::vector<size_t> sz(2, 0);
            matrix<size_t> shape(2, 1);
            std::vector<std::vector<size_t>> shapes;

            if (type.compare("res") == 0)
            {
                if (residuals.size() > 0)
                {
                    for (size_t i = 0; i < residuals.size(); i++)
                    {
                        shape = residuals[i].shape();
                        sz[0] = shape[0];
                        sz[1] = shape[1];
                        shapes.push_back(sz);
                    }
                    shape.clear();
                }
                // return shapes;
            }

            if (type.compare("obs") == 0)
            {
                if (observations.size() > 0)
                {
                    for (size_t i = 0; i < observations.size(); i++)
                    {
                        shape = observations[i].shape();
                        sz[0] = shape[0];
                        sz[1] = shape[1];
                        shapes.push_back(sz);
                    }
                    shape.clear();
                }
                // return shapes;
            }

            if (type.compare("eff") == 0)
            {
                if (all_effects.size() > 0)
                {
                    for (size_t i = 0; i < all_effects.size(); i++)
                    {
                        std::vector<size_t> shp;
                        shp = all_effects[i].shape();
                        // sz[0] = shape[0]; // on not transposed
                        // sz[1] = shape[1]; // on not transposed
                        sz[0] = shp[1]; // on transposed
                        sz[1] = shp[0]; // on transposed

                        shapes.push_back(sz);
                    }
                    shape.clear();
                }
                // return shapes;
            }

            if (type.compare("var") == 0)
            {
                if (variances.size() > 0)
                {
                    for (size_t i = 0; i < variances.size(); i++)
                    {
                        shape = variances[i].shape();
                        sz[0] = shape[0];
                        sz[1] = shape[1];
                        shapes.push_back(sz);
                    }
                    shape.clear();
                }
                // return shapes;
            }

            if (type.compare("cor") == 0)
            {
                if (correlations.size() > 0)
                {
                    for (size_t i = 0; i < correlations.size(); i++)
                    {
                        sz[0] = correlations[i].nrows();
                        sz[1] = correlations[i].ncols();
                        shapes.push_back(sz);
                    }
                }
                // return shapes;
            }

            if (type.compare("cor_eff") == 0)
            {
                if (correlated_effects.size() > 0)
                {
                    for (size_t i = 0; i < correlated_effects.size(); i++)
                    {
                        shape = correlated_effects[i].shape();
                        sz[0] = shape[0];
                        sz[1] = shape[1];
                        shapes.push_back(sz);
                    }
                    shape.clear();
                }
                // return shapes;
            }

            if (type.compare("obs_trt") == 0)
            {
                if (observation_trait.size() > 0)
                {
                    for (size_t i = 0; i < observation_trait.size(); i++)
                    {
                        // shape = observation_trait[i].shape();
                        sz[0] = 1;
                        sz[1] = 1;
                        shapes.push_back(sz);
                    }
                    shape.clear();
                }
                // return shapes;
            }

            if (type.compare("eff_trt") == 0)
            {
                if (effects_trait.size() > 0)
                {
                    for (size_t i = 0; i < effects_trait.size(); i++)
                    {
                        shape = effects_trait[i].shape();
                        sz[0] = shape[0];
                        sz[1] = shape[1];
                        shapes.push_back(sz);
                    }
                    shape.clear();
                }
                // return shapes;
            }

            return shapes;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::shape_of(const std::string)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::shape_of(const std::string)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    compact_storage<float> model_sparse::test_effects(size_t which_effect)
    {
        try
        {
            compact_storage<float> effect;

            all_effects[which_effect].fread();

            effects_storage e = all_effects[which_effect];

            all_effects[which_effect].fwrite();

            e.get(effect);

            effect.transpose(); // because we transposing when appending

            return effect;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::test_effects(size_t, T)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::test_effects(size_t, T)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    std::vector<float> model_sparse::test_observations(size_t which_observations)
    {
        try
        {
            std::vector<float> out;

            matrix<float> obs;

            obs = observations[which_observations];

            obs.fread();

            for (size_t i = 0; i < obs.size(); i++)
                out.push_back(obs[i]);

            return out;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::test_observations(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::test_observations(size_t)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    std::vector<float> model_sparse::test_residual(size_t which_residual)
    {
        try
        {
            std::vector<float> out;

            matrix<float> resid;

            resid = residuals[which_residual];

            resid.fread();

            for (size_t i = 0; i < resid.size(); i++)
                out.push_back(resid[i]);

            return out;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::test_residual(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::test_residual(size_t)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    std::vector<float> model_sparse::test_variance(size_t which_variance)
    {
        try
        {
            std::vector<float> out;

            matrix<float> var;

            var = variances[which_variance];

            var.fread();

            for (size_t i = 0; i < var.size(); i++)
                out.push_back(var[i]);

            return out;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::test_variance(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::test_variance(size_t)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    std::vector<float> model_sparse::test_correlation(size_t which_correlation)
    {
        try
        {
            std::vector<float> out;

            compact_storage<float> cor;

            cor = correlations[which_correlation];

            cor.fread();

            matrix<float> M;
            cor.to_dense(M);
            M.to_vector(out);

            return out;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_sparse::test_correlation(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_sparse::test_correlatino(size_t)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
#endif

} // end of namespace evolm
