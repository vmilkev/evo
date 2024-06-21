#include "effects_storage.hpp"

namespace evolm
{
    effects_storage::effects_storage()
    {
        is_diagonal = false;
    }
    //===============================================================================================================
    void effects_storage::set(compact_storage<int> &in_effect)
    {
        try
        {
            i_effect = in_effect;
            type = 1;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::set(compact_storage<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::set(compact_storage<int> &)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void effects_storage::set(compact_storage<float> &in_effect)
    {
        try
        {
            f_effect = in_effect;
            type = 2;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::set(compact_storage<float> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::set(compact_storage<float> &)." << '\n';
            throw;
        }
    }
    //===============================================================================================================

    void effects_storage::set(compact_storage<double> &in_effect)
    {
        try
        {
            d_effect = in_effect;
            type = 3;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::set(compact_storage<double> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::set(compact_storage<double> &)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void effects_storage::get(compact_storage<int> &out_effect)
    {
        try
        {
            switch (type)
            {
            case 1:
                out_effect = i_effect;
                break;
            case 2:
                f_effect.to_int(out_effect);
                break;
            case 3:
                d_effect.to_int(out_effect);
                break;            
            default:
                throw std::string("Undetermined or invalid data type!");
                break;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::get(compact_storage<int> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::get(compact_storage<int> &)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void effects_storage::get(compact_storage<float> &out_effect)
    {
        try
        {
            switch (type)
            {
            case 1:
                i_effect.to_float(out_effect);
                break;
            case 2:
                out_effect = f_effect;
                break;
            case 3:
                d_effect.to_float(out_effect);
                break;            
            default:
                throw std::string("Undetermined or invalid data type!");
                break;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::get(compact_storage<float> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::get(compact_storage<float> &)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void effects_storage::get(compact_storage<double> &out_effect)
    {
        try
        {           
            switch (type)
            {
            case 1:
                i_effect.to_double(out_effect);
                break;
            case 2:
                f_effect.to_double(out_effect);
                break;
            case 3:
                out_effect = d_effect;
                break;            
            default:
                throw std::string("Undetermined or invalid data type!");
                break;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::get(compact_storage<double> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::get(compact_storage<double> &)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void effects_storage::fread()
    {
        try
        {
            switch (type)
            {
            case 1:
                i_effect.fread();
                break;
            case 2:
                f_effect.fread();
                break;
            case 3:
                d_effect.fread();
                break;            
            default:
                throw std::string("Undetermined or invalid data type!");
                break;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::fread()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::fread()." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void effects_storage::fwrite()
    {
        try
        {
            switch (type)
            {
            case 1:
                i_effect.fwrite();
                break;
            case 2:
                f_effect.fwrite();
                break;
            case 3:
                d_effect.fwrite();
                break;            
            default:
                throw std::string("Undetermined or invalid data type!");
                break;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::fwrite()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::fwrite()." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void effects_storage::get_fcast(matrix<float> &out, size_t *irow, size_t *icol)
    {
        try
        {
            switch (type)
            {
            case 1:
                i_effect.cast_to_fdense(out, irow, icol);
                break;
            case 2:
                f_effect.to_dense(out, irow, icol);
                break;
            case 3:
                d_effect.cast_to_fdense(out, irow, icol);
                break;            
            default:
                throw std::string("Undetermined or invalid data type!");
                break;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::get_fcast(matrix<float> &, size_t *, size_t *)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::get_fcast(matrix<float> &, size_t *, size_t *)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void effects_storage::get_fcast(smatrix<float> &out, size_t *irow, size_t *icol)
    {
        try
        {
            switch (type)
            {
            case 1:
                i_effect.cast_to_fsparse(out, irow, icol);
                break;
            case 2:
                f_effect.to_sparse(out, irow, icol);
                break;
            case 3:
                d_effect.cast_to_fsparse(out, irow, icol);
                break;            
            default:
                throw std::string("Undetermined or invalid data type!");
                break;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::get_fcast(smatrix<float> &, size_t *, size_t *)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::get_fcast(smatrix<float> &, size_t *, size_t *)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void effects_storage::get_fcast(float **out, size_t *irow, size_t *icol)
    {
        try
        {
            switch (type)
            {
            case 1:
                i_effect.cast_to_fdense(out, irow, icol);
                break;
            case 2:
                f_effect.to_dense(out, irow, icol);
                break;
            case 3:
                d_effect.cast_to_fdense(out, irow, icol);
                break;            
            default:
                throw std::string("Undetermined or invalid data type!");
                break;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::get_fcast(float **, size_t *, size_t *)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::get_fcast(float **, size_t *, size_t *)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void effects_storage::get_fcast(std::vector<float> &vals_out, std::vector<size_t> &keys_out, size_t *irow, size_t *icol)
    {
        try
        {
            switch (type)
            {
            case 1:
                i_effect.cast_to_fsparse(vals_out, keys_out, irow, icol);
                break;
            case 2:
                f_effect.to_sparse(vals_out, keys_out, irow, icol);
                break;
            case 3:
                d_effect.cast_to_fsparse(vals_out, keys_out, irow, icol);
                break;            
            default:
                throw std::string("Undetermined or invalid data type!");
                break;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::get_fcast(std::vector<float> &, std::vector<size_t> &, size_t *, size_t *)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::get_fcast(std::vector<float> &, std::vector<size_t> &, size_t *, size_t *)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void effects_storage::row_dot_float_vect(float **in_vect, size_t irow, size_t icol, float &out_res)
    {
        try
        {
            switch (type)
            {
            case 1:
                i_effect.row_dot_float_vect(in_vect, irow, icol, out_res);
                break;
            case 2:
                f_effect.row_dot_vect(in_vect, irow, icol, out_res);
                break;
            case 3:
                d_effect.row_dot_float_vect(in_vect, irow, icol, out_res);
                break;            
            default:
                throw std::string("Undetermined or invalid data type!");
                break;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::get_fcast(float **, size_t *, size_t *)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::get_fcast(float **, size_t *, size_t *)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    std::vector<size_t> effects_storage::shape()
    {
        try
        {
            std::vector<size_t> shp(2, 0);

            switch (type)
            {
            case 1:
                shp[0] = i_effect.nrows();
                shp[1] = i_effect.ncols();
                break;
            case 2:
                shp[0] = f_effect.nrows();
                shp[1] = f_effect.ncols();
                break;
            case 3:
                shp[0] = d_effect.nrows();
                shp[1] = d_effect.ncols();
                break;
            default:
                throw std::string("Undetermined or invalid data type!");
                break;
            }

            return shp;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::shape()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::shape()." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    size_t effects_storage::size()
    {
        try
        {
            switch (type)
            {
            case 1:
                return i_effect.size();
            case 2:
                return f_effect.size();
            case 3:
                return d_effect.size();
            default:
                throw std::string("Undetermined or invalid data type!");
                break;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::size()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::size()." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void effects_storage::clear()
    {
        try
        {
            i_effect.fclear();
            i_effect.clear();
            f_effect.fclear();
            f_effect.clear();
            d_effect.fclear();
            d_effect.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::clear()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::clear()." << '\n';
            throw;
        }
    }
    //===============================================================================================================

}
