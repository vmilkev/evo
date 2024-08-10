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
                i_effect.fread_rows_structure();
                break;
            case 2:
                f_effect.fread();
                f_effect.fread_rows_structure();
                break;
            case 3:
                d_effect.fread();
                d_effect.fread_rows_structure();
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
                i_effect.fwrite_rows_structure();
                break;
            case 2:
                f_effect.fwrite();
                f_effect.fwrite_rows_structure();
                break;
            case 3:
                d_effect.fwrite();
                d_effect.fwrite_rows_structure();
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
            case 0:
                return 0;
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
            type = 0;
            is_diagonal = false;
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
    void effects_storage::extend_by(effects_storage &rhs_matr)
    {
        try
        {
            // temporal smatrix containers for LHS effects_storage
            evolm::smatrix<int> i_sm;
            evolm::smatrix<float> f_sm;
            evolm::smatrix<double> d_sm;
            // temporal smatrix containers for RHS effects_storage
            evolm::smatrix<int> i_sm_rhs;
            evolm::smatrix<float> f_sm_rhs;
            evolm::smatrix<double> d_sm_rhs;
            // temporal std::vector containers for resetting the LHS effects_storage
            std::vector<int> i_values;
            std::vector<float> f_values;
            std::vector<double> d_values;
            std::vector<size_t> keys;
            // temporal compact_storage containers for RHS effects_storage
            compact_storage<int> i_effect_rhs;    // type 1
            compact_storage<float> f_effect_rhs;  // type 2
            compact_storage<double> d_effect_rhs; // type 3

            switch (rhs_matr.type) // get the right compact_storage from the RHS effects_storage
            {
            case 1:
                rhs_matr.get(i_effect_rhs);
                break;
            case 2:
                rhs_matr.get(f_effect_rhs);
                break;
            case 3:
                rhs_matr.get(d_effect_rhs);
                break;
            default:
                throw std::string("Undetermined or invalid data type!");
                break;
            }

            if ( type == rhs_matr.type ) // if LHS type is the same as RHS type
            {
                switch (type) // do extend_by for the right type of compact_storage
                {
                case 1: // int
                    i_effect.to_sparse(i_sm); // convert to smatrix
                    i_effect_rhs.to_sparse(i_sm_rhs); // convert to smatrix
                    i_sm.extend_by(i_sm_rhs); // do extend_by                    
                    i_sm.to_vect(i_values, keys); // convert smatrix to std::vector
                    i_effect.resize(i_sm.nrows(), i_sm.ncols()); // resizing the LHS compact_storage for resetting/update
                    i_effect.append_with_keys(i_values, keys); // reset/update the LHS compact_storage
                    break;
                case 2: // float
                    f_effect.to_sparse(f_sm);
                    f_effect_rhs.to_sparse(f_sm_rhs);
                    f_sm.extend_by(f_sm_rhs);
                    f_sm.to_vect(f_values, keys);
                    f_effect.resize(f_sm.nrows(), f_sm.ncols());
                    f_effect.append_with_keys(f_values, keys);
                    break;
                case 3: // double
                    d_effect.to_sparse(d_sm);
                    d_effect_rhs.to_sparse(d_sm_rhs);
                    d_sm.extend_by(d_sm_rhs);
                    d_sm.to_vect(d_values, keys);
                    d_effect.resize(d_sm.nrows(), d_sm.ncols());
                    d_effect.append_with_keys(d_values, keys);
                    break;
                default:
                    throw std::string("Undetermined or invalid data type!");
                    break;
                }
            }
            else
            {
                // containers for converted types (LHS -> RHS or RHS -> LHS)
                compact_storage<int> i_tmp;
                compact_storage<float> f_tmp;
                compact_storage<double> d_tmp;

                if ( (type - rhs_matr.type) < 0 ) // convert LHS type to RHS type
                {
                    switch (type)
                    {
                    case 1:
                        if (rhs_matr.type == 2)
                            i_effect.to_float(f_tmp);
                        else if (rhs_matr.type == 3)
                            i_effect.to_double(d_tmp);
                        else
                            throw std::string("Unknown type of rhs_matr during compact_storage conversion!");
                        break;
                    case 2:
                        if (rhs_matr.type == 3)
                            f_effect.to_double(d_tmp);
                        else
                            throw std::string("Unknown type of rhs_matr during compact_storage conversion!");
                        break;
                    default:
                        throw std::string("Undetermined or invalid data type!");
                        break;
                    }

                    switch (type) // do extend_by for the right type of compact_storage
                    {
                    case 1:
                        if (rhs_matr.type == 2) // LHS was converted to float
                        {
                            i_effect.resize();
                            type = 2;

                            f_tmp.to_sparse(f_sm);
                            f_effect_rhs.to_sparse(f_sm_rhs);
                            f_sm.extend_by(f_sm_rhs);
                            
                            f_sm.to_vect(f_values, keys);
                            f_effect.resize(f_sm.nrows(), f_sm.ncols());
                            f_effect.append_with_keys(f_values, keys);
                        }
                        if (rhs_matr.type == 3) // LHS was converted to double
                        {
                            i_effect.resize();
                            type = 3;

                            d_tmp.to_sparse(d_sm);
                            d_effect_rhs.to_sparse(d_sm_rhs);
                            d_sm.extend_by(d_sm_rhs);
                            
                            d_sm.to_vect(d_values, keys);
                            d_effect.resize(d_sm.nrows(), d_sm.ncols());
                            d_effect.append_with_keys(d_values, keys);
                        }
                        break;
                    case 2:
                        // LHS was converted to double

                        f_effect.resize();
                        type = 3;

                        d_tmp.to_sparse(d_sm);
                        d_effect_rhs.to_sparse(d_sm_rhs);
                        d_sm.extend_by(d_sm_rhs);
                        
                        d_sm.to_vect(d_values, keys);
                        d_effect.resize(d_sm.nrows(), d_sm.ncols());
                        d_effect.append_with_keys(d_values, keys);
                        break;
                    default:
                        throw std::string("Undetermined or invalid data type!");
                        break;
                    }
                }

                if ( (type - rhs_matr.type) > 0 )
                {
                    switch (type) // convert RHS type to LHS type
                    {
                    case 3:
                        if (rhs_matr.type == 2)
                            f_effect_rhs.to_double(d_tmp);
                        else if (rhs_matr.type == 1)
                            i_effect_rhs.to_double(d_tmp);
                        else
                            throw std::string("Unknown type of rhs_matr during compact_storage conversion!");
                        break;
                    case 2:
                        if (rhs_matr.type == 1)
                            i_effect_rhs.to_float(f_tmp);
                        else
                            throw std::string("Unknown type of rhs_matr during compact_storage conversion!");
                        break;
                    default:
                        throw std::string("Undetermined or invalid data type!");
                        break;
                    }

                    switch (type) // do extend_by for the right type of compact_storage
                    {
                    case 2:
                        f_effect.to_sparse(f_sm);
                        f_tmp.to_sparse(f_sm_rhs);
                        f_sm.extend_by(f_sm_rhs);

                        f_sm.to_vect(f_values, keys);
                        f_effect.resize(f_sm.nrows(), f_sm.ncols());
                        f_effect.append_with_keys(f_values, keys);
                        break;
                    case 3:
                        d_effect.to_sparse(d_sm);
                        d_tmp.to_sparse(d_sm_rhs);
                        d_sm.extend_by(d_sm_rhs);

                        d_sm.to_vect(d_values, keys);
                        d_effect.resize(d_sm.nrows(), d_sm.ncols());
                        d_effect.append_with_keys(d_values, keys);
                        break;
                    default:
                        throw std::string("Undetermined or invalid data type!");
                        break;
                    }
                }

                i_tmp.resize();
                f_tmp.resize();
                d_tmp.resize();
            }

            i_effect_rhs.resize();
            f_effect_rhs.resize();
            d_effect_rhs.resize();

            i_values.clear(); i_values.shrink_to_fit();
            f_values.clear(); f_values.shrink_to_fit();
            d_values.clear(); d_values.shrink_to_fit();

            i_sm_rhs.clear();
            f_sm_rhs.clear();
            d_sm_rhs.clear();

            i_sm.clear();
            f_sm.clear();
            d_sm.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::extend_by(effects_storage &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::extend_by(effects_storage &)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void effects_storage::element_wise_dot(effects_storage &rhs_matr)
    {
        try
        {
            // temporal smatrix containers for LHS effects_storage
            evolm::smatrix<int> i_sm;
            evolm::smatrix<float> f_sm;
            evolm::smatrix<double> d_sm;
            // temporal smatrix containers for RHS effects_storage
            evolm::smatrix<int> i_sm_rhs;
            evolm::smatrix<float> f_sm_rhs;
            evolm::smatrix<double> d_sm_rhs;
            // temporal std::vector containers for resetting the LHS effects_storage
            std::vector<int> i_values;
            std::vector<float> f_values;
            std::vector<double> d_values;
            std::vector<size_t> keys;
            // temporal compact_storage containers for RHS effects_storage
            compact_storage<int> i_effect_rhs;    // type 1
            compact_storage<float> f_effect_rhs;  // type 2
            compact_storage<double> d_effect_rhs; // type 3

            switch (rhs_matr.type) // get the right compact_storage from the RHS effects_storage
            {
            case 1:
                rhs_matr.get(i_effect_rhs);
                break;
            case 2:
                rhs_matr.get(f_effect_rhs);
                break;
            case 3:
                rhs_matr.get(d_effect_rhs);
                break;
            default:
                throw std::string("Undetermined or invalid data type!");
                break;
            }

            if ( type == rhs_matr.type ) // if LHS type is the same as RHS type
            {
                switch (type) // do element_wise_dot for the right type of compact_storage
                {
                case 1: // int
                    i_effect.to_sparse(i_sm); // convert to smatrix
                    i_effect_rhs.to_sparse(i_sm_rhs); // convert to smatrix
                    i_sm.element_wise_dot(i_sm_rhs); // do element_wise_dot                    
                    i_sm.to_vect(i_values, keys); // convert smatrix to std::vector
                    i_effect.resize(i_sm.nrows(), i_sm.ncols()); // resizing the LHS compact_storage for resetting/update
                    i_effect.append_with_keys(i_values, keys); // reset/update the LHS compact_storage
                    break;
                case 2: // float
                    f_effect.to_sparse(f_sm);
                    f_effect_rhs.to_sparse(f_sm_rhs);
                    f_sm.element_wise_dot(f_sm_rhs);
                    f_sm.to_vect(f_values, keys);
                    f_effect.resize(f_sm.nrows(), f_sm.ncols());
                    f_effect.append_with_keys(f_values, keys);
                    break;
                case 3: // double
                    d_effect.to_sparse(d_sm);
                    d_effect_rhs.to_sparse(d_sm_rhs);
                    d_sm.element_wise_dot(d_sm_rhs);
                    d_sm.to_vect(d_values, keys);
                    d_effect.resize(d_sm.nrows(), d_sm.ncols());
                    d_effect.append_with_keys(d_values, keys);
                    break;
                default:
                    throw std::string("Undetermined or invalid data type!");
                    break;
                }
            }
            else
            {
                // containers for converted types (LHS -> RHS or RHS -> LHS)
                compact_storage<int> i_tmp;
                compact_storage<float> f_tmp;
                compact_storage<double> d_tmp;

                if ( (type - rhs_matr.type) < 0 ) // convert LHS type to RHS type
                {
                    switch (type)
                    {
                    case 1:
                        if (rhs_matr.type == 2)
                            i_effect.to_float(f_tmp);
                        else if (rhs_matr.type == 3)
                            i_effect.to_double(d_tmp);
                        else
                        {
                            std::cout<<"HERE 1, type "<<rhs_matr.type<<"\n";
                            throw std::string("Unknown type of rhs_matr during compact_storage conversion!");
                        }
                        break;
                    case 2:
                        if (rhs_matr.type == 3)
                            f_effect.to_double(d_tmp);
                        else
                        {
                            std::cout<<"HERE 2, type "<<rhs_matr.type<<"\n";
                            throw std::string("Unknown type of rhs_matr during compact_storage conversion!");
                        }
                        break;
                    default:
                        throw std::string("Undetermined or invalid data type!");
                        break;
                    }

                    switch (type) // do element_wise_dot for the right type of compact_storage
                    {
                    case 1:
                        if (rhs_matr.type == 2) // LHS was converted to float
                        {
                            i_effect.resize();
                            type = 2;

                            f_tmp.to_sparse(f_sm);
                            f_effect_rhs.to_sparse(f_sm_rhs);
                            f_sm.element_wise_dot(f_sm_rhs);
                            
                            f_sm.to_vect(f_values, keys);
                            f_effect.resize(f_sm.nrows(), f_sm.ncols());
                            f_effect.append_with_keys(f_values, keys);
                        }
                        if (rhs_matr.type == 3) // LHS was converted to double
                        {
                            i_effect.resize();
                            type = 3;

                            d_tmp.to_sparse(d_sm);
                            d_effect_rhs.to_sparse(d_sm_rhs);
                            d_sm.element_wise_dot(d_sm_rhs);
                            
                            d_sm.to_vect(d_values, keys);
                            d_effect.resize(d_sm.nrows(), d_sm.ncols());
                            d_effect.append_with_keys(d_values, keys);
                        }
                        break;
                    case 2:
                        // LHS was converted to double

                        f_effect.resize();
                        type = 3;

                        d_tmp.to_sparse(d_sm);
                        d_effect_rhs.to_sparse(d_sm_rhs);
                        d_sm.element_wise_dot(d_sm_rhs);
                        
                        d_sm.to_vect(d_values, keys);
                        d_effect.resize(d_sm.nrows(), d_sm.ncols());
                        d_effect.append_with_keys(d_values, keys);
                        break;
                    default:
                        throw std::string("Undetermined or invalid data type!");
                        break;
                    }
                }

                if ( (type - rhs_matr.type) > 0 )
                {
                    switch (type) // convert RHS type to LHS type
                    {
                    case 3:
                        //
                        if (rhs_matr.type == 2)
                            f_effect_rhs.to_double(d_tmp);
                        else if (rhs_matr.type == 1)
                            i_effect_rhs.to_double(d_tmp);
                        else
                        {std::cout<<"HERE 3, type "<<rhs_matr.type<<"\n";
                            throw std::string("Unknown type of rhs_matr during compact_storage conversion!");
                        }
                        break;
                    case 2:
                        if (rhs_matr.type == 1)
                            i_effect_rhs.to_float(f_tmp);
                        else
                        {
                            std::cout<<"HERE 4, type "<<rhs_matr.type<<" type "<<type<<"\n";
                            throw std::string("Unknown type of rhs_matr during compact_storage conversion!");
                        }
                        break;
                    default:
                        throw std::string("Undetermined or invalid data type!");
                        break;
                    }

                    switch (type) // do element_wise_dot for the right type of compact_storage
                    {
                    case 2:
                        f_effect.to_sparse(f_sm);
                        f_tmp.to_sparse(f_sm_rhs);
                        f_sm.element_wise_dot(f_sm_rhs);

                        f_sm.to_vect(f_values, keys);
                        f_effect.resize(f_sm.nrows(), f_sm.ncols());
                        f_effect.append_with_keys(f_values, keys);
                        break;
                    case 3:
                        d_effect.to_sparse(d_sm);
                        d_tmp.to_sparse(d_sm_rhs);
                        d_sm.element_wise_dot(d_sm_rhs);

                        d_sm.to_vect(d_values, keys);
                        d_effect.resize(d_sm.nrows(), d_sm.ncols());
                        d_effect.append_with_keys(d_values, keys);
                        break;
                    default:
                        throw std::string("Undetermined or invalid data type!");
                        break;
                    }
                }

                i_tmp.resize();
                f_tmp.resize();
                d_tmp.resize();
            }

            i_effect_rhs.resize();
            f_effect_rhs.resize();
            d_effect_rhs.resize();

            i_values.clear(); i_values.shrink_to_fit();
            f_values.clear(); f_values.shrink_to_fit();
            d_values.clear(); d_values.shrink_to_fit();

            i_sm_rhs.clear();
            f_sm_rhs.clear();
            d_sm_rhs.clear();

            i_sm.clear();
            f_sm.clear();
            d_sm.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::element_wise_dot(effects_storage &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::element_wise_dot(effects_storage &)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void effects_storage::print(std::string &message)
    {
        try
        {
            evolm::smatrix<int> i_m;
            evolm::smatrix<float> f_m;
            evolm::smatrix<double> d_m;

            switch (type)
            {
            case 1:
                i_effect.to_sparse(i_m);
                i_m.print(message);
                break;
            case 2:
                f_effect.to_sparse(f_m);
                f_m.print(message);
                break;
            case 3:
                d_effect.to_sparse(d_m);
                d_m.print(message);
                break;
            default:
                throw std::string("Undetermined or invalid data type!");
                break;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in effects_storage::print(std::string &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in effects_storage::print(std::string &)." << '\n';
            throw;
        }
    }
    //===============================================================================================================

}
