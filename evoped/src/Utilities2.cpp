#include "Utilities2.hpp"

namespace evoped
{
    //===============================================================================================================

    Utilities2::Utilities2()
    {
    }
    
    //===============================================================================================================

    Utilities2::~Utilities2()
    {
    }

    //===============================================================================================================

    int Utilities2::find_invect(std::vector<std::int64_t> &where, std::int64_t what)
    {
        try
        {
            std::vector<std::int64_t>::iterator it;
            it = find(where.begin(), where.end(), what);
            if (it != where.end())
                return it - where.begin();
            else
                return -1;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::find_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::find_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::find_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Utilities2::check_id(std::vector<std::int64_t> &id_list, std::vector<std::int64_t> &checked_id, std::vector<std::int64_t> &missing_id)
    {
        try
        {
            // id_list - vector of ids from where to check
            // checked_id - vector of ids to check
            // missing_id - vector of missing ids (whose not found in id vector)

            /*for (size_t i = 0; i < checked_id.size(); i++)
                if (!std::binary_search(id_list.begin(), id_list.end(), checked_id[i]))
                    missing_id.push_back(checked_id[i]);*/
            for (size_t i = 0; i < checked_id.size(); i++)
                if ( std::find(id_list.begin(), id_list.end(), checked_id[i]) == id_list.end() )
                    missing_id.push_back(checked_id[i]);

        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::check_id(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::check_id(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::check_id(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Utilities2::check_id2(std::vector<std::int64_t> &id_list, std::vector<std::int64_t> &checked_id, std::vector<std::int64_t> &is_in_id)
    {
        try
        {
            // id_list - vector of ids from where to check
            // checked_id - vector of ids to check
            // is_in_id - vector of present ids (whose are found in id vector)

            /*for (size_t i = 0; i < checked_id.size(); i++)
                if ( std::binary_search(id_list.begin(), id_list.end(), checked_id[i]) )
                    is_in_id.push_back(checked_id[i]);*/
            for (size_t i = 0; i < checked_id.size(); i++)
                if ( std::find(id_list.begin(), id_list.end(), checked_id[i]) != id_list.end() )
                    is_in_id.push_back(checked_id[i]);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::check_id2(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::check_id2(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::check_id2(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    bool Utilities2::is_value_in_vect(std::vector<std::int64_t> &where_tocheck, std::vector<std::int64_t> &what_tocheck)
    {
        // are values from what_tocheck in where_tocheck?
        try
        {
            //----------------------------------
            /*std::vector<int64_t> vec_1(where_tocheck);
            std::vector<int64_t> vec_2(what_tocheck);
            
            std::sort(vec_1.begin(), vec_1.end());
            std::sort(vec_2.begin(), vec_2.end());*/
            //----------------------------------
            bool out = true;

            std::vector<std::int64_t> missing;
            check_id(where_tocheck, what_tocheck, missing);
            //check_id(vec_1, vec_2, missing);

            if (!missing.empty())
            {
                out = false;
                missing.clear();
            }

            //vec_1.clear();
            //vec_2.clear();

            return out;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::is_value_in_vect(std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::is_value_in_vect(std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::is_value_in_vect(std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Utilities2::get_missing_in_vect(std::vector<std::int64_t> &where_tocheck, std::vector<std::int64_t> &what_tocheck, std::vector<std::int64_t> &missing)
    {
        try
        {
            std::vector<int64_t> vec_1(where_tocheck);
            std::vector<int64_t> vec_2(what_tocheck);
            
            std::sort(vec_1.begin(), vec_1.end());
            std::sort(vec_2.begin(), vec_2.end());

            check_id(vec_1, vec_2, missing);

            vec_1.clear();
            vec_2.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::get_missing_in_vect(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::get_missing_in_vect(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::get_missing_in_vect(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Utilities2::get_what_in_vect(std::vector<std::int64_t> &where_tocheck, std::vector<std::int64_t> &what_tocheck, std::vector<std::int64_t> &what_in)
    {
        try
        {
            std::vector<int64_t> vec_1(where_tocheck);
            std::vector<int64_t> vec_2(what_tocheck);
            
            std::sort(vec_1.begin(), vec_1.end());
            std::sort(vec_2.begin(), vec_2.end());

            check_id2(vec_1, vec_2, what_in);

            vec_1.clear();
            vec_2.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::get_what_in_vect(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::get_what_in_vect(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::get_what_in_vect(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    bool Utilities2::is_invect(std::vector<std::int64_t> &where, std::int64_t what)
    {
        try
        {
            std::vector<std::int64_t>::iterator it;
            it = find(where.begin(), where.end(), what);
            if (it != where.end())
                return true;
            else
                return false;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::is_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::is_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::is_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    bool Utilities2::is_unique(std::vector<std::int64_t> &x)
    {
        try
        {
            bool out = false;
            sort(x.begin(), x.end());
            if (adjacent_find(x.begin(), x.end()) == x.end())
                out = true;
            return out;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::is_unique(std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::is_unique(std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::is_unique(std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Utilities2::vect_to_binary(std::vector<std::int64_t> &vect, const std::string &fname)
    {
        try
        {
            std::fstream fA;
            fA.open(fname, fA.binary | fA.trunc | fA.out);

            if (!fA.is_open())
            {
                throw std::string("Error while opening a binary file for writing.");
            }

            size_t sz = vect.size();

            size_t B[2];
            B[0] = sz;
            B[1] = (size_t)sizeof(std::int64_t); // points to a type of data
            
            fA.write(reinterpret_cast<char *>(&B), 2 * sizeof(size_t));

            fA.write(reinterpret_cast<char *>(&vect[0]), sz * sizeof(std::int64_t));

            fA.close();

        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::vect_to_binary(std::vector<std::int64_t> &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::vect_to_binary(std::vector<std::int64_t> &, const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::vect_to_binary(std::vector<std::int64_t> &, const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Utilities2::vect_from_binary(std::vector<std::int64_t> &vect, const std::string &fname)
    {
        try
        {
            std::fstream fA;
            fA.open(fname, fA.binary | fA.in);

            if (!fA.is_open())
            {
                throw std::string("Error while opening a binary file for reading.");
            }

            size_t B[2];

            fA.read(reinterpret_cast<char *>(&B), 2 * sizeof(size_t));

            size_t sz = B[0];
            size_t var_inbytes = B[1];

            if ( var_inbytes != sizeof(std::int64_t) )
                throw std::string("The data type stored in file is not consistent with type of the vector trying to read the data.");

            if ( !vect.empty() )
                vect.clear();

            std::int64_t vect_value = 0;
            for (size_t i = 0; i < sz; i++)
            {
                fA.read(reinterpret_cast<char *>(&vect_value), sizeof(std::int64_t));
                vect.push_back(vect_value);
            }

            fA.close();

            if ( vect.size() != 0 )
                throw std::string("The amount of data received from file is zerro!.");

            if ( vect.size() != sz )
                throw std::string("The amount of data received from file is not correct!");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::vect_from_binary(std::vector<std::int64_t> &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::vect_from_binary(std::vector<std::int64_t> &, const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::vect_from_binary(std::vector<std::int64_t> &, const std::string &)" << '\n';
            throw;
        }
    }
    //===========================================================================================
    template <typename T>
    void Utilities2::fwrite_matrix(const std::string &fname, std::vector<T> &vals, std::vector<size_t> &keys, std::vector<std::int64_t> &ids, size_t n_cols)
    {
        //int integer_var;
        float float_var;
        double double_var;

        // determine type of the matrix, expectd: float or double
        //const std::type_info &type_1 = typeid(integer_var);
        const std::type_info &type_2 = typeid(float_var);
        const std::type_info &type_3 = typeid(double_var);
        const std::type_info &type_T = typeid(vals[0]);
        
        size_t var_type = 0; // type of matrix
        size_t var_type2 = 5; // int64 type of IDs

        if (type_T == type_2)
            var_type = 2;
        else if (type_T == type_3)
            var_type = 3;
        else
            throw std::string("Utilities2::fwrite_matrix(const std::string &, std::vector<T> &, std::vector<size_t> &, std::vector<T2> &, size_t): Cannot determine the type of values vector.");

        std::string name_suffix(".corbin");
        size_t find = fname.find(name_suffix);
        std::string fname2(fname);
        if( find == std::string::npos ) // there is no suffix
            fname2 = fname + name_suffix;

        std::fstream fA;
        fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        fA.open(fname2, fA.binary | fA.trunc | fA.out);

        if (!fA.is_open())
            throw std::string("Utilities2::fwrite_matrix(const std::string &, std::vector<T> &, std::vector<size_t> &, std::vector<T2> &, size_t): Error while opening a binary file.");

        size_t B[6]; // necessary data info

        // info about matrix values
        B[0] = matrix_type;
        B[1] = var_type; // expected 2 or 3
        
        // info about matrix IDs
        B[2] = var_type2; // expected 4 or 5
        B[3] = (size_t)sizeof(T);
        B[4] = ids.size();
        B[5] = n_cols;
        
        size_t vals_size = vals.size();
        size_t keys_size = keys.size();
        size_t ids_size = ids.size();

        fA.write( reinterpret_cast<const char *>(&B), 6 * sizeof(size_t) ); // 0. writing the data info

        fA.write( reinterpret_cast<const char *>( &vals_size ), sizeof(size_t) ); // 1. writing size of values
        fA.write( reinterpret_cast<const char *>( vals.data() ), vals_size * sizeof(T) ); // 2. writing all values

        fA.write( reinterpret_cast<const char *>( &keys_size ), sizeof(size_t) ); // 3. writing size of keys
        fA.write( reinterpret_cast<const char *>( keys.data() ), keys_size * sizeof(size_t) ); // 4. writing all keys

        fA.write( reinterpret_cast<const char *>( &ids_size ), sizeof(size_t) ); // 5. writing size of ids
        fA.write( reinterpret_cast<const char *>( ids.data() ), ids_size * sizeof(std::int64_t) ); // 6. writing all IDs

        fA.close();
    }
    template void Utilities2::fwrite_matrix(const std::string &fname, std::vector<float> &vals, std::vector<size_t> &keys, std::vector<std::int64_t> &ids, size_t n_cols);
    template void Utilities2::fwrite_matrix(const std::string &fname, std::vector<double> &vals, std::vector<size_t> &keys, std::vector<std::int64_t> &ids, size_t n_cols);
    //===========================================================================================
    template <typename T>
    void Utilities2::fwrite_matrix(const std::string &fname, evolm::matrix<T> &matr, std::vector<std::int64_t> &ids)
    {
        std::string name_suffix(".dmbin");
        size_t find = fname.find(name_suffix);
        std::string fname2(fname);
        if( find == std::string::npos ) // there is no suffix
            fname2 = fname + name_suffix;

        matr.fwrite(fname2, ids);
    }
    template void Utilities2::fwrite_matrix(const std::string &fname, evolm::matrix<float> &matr, std::vector<std::int64_t> &ids);
    template void Utilities2::fwrite_matrix(const std::string &fname, evolm::matrix<double> &matr, std::vector<std::int64_t> &ids);
    //===========================================================================================
    template <typename T>
    void Utilities2::fwrite_matrix(const std::string &fname, evolm::matrix<T> &matr, std::vector<std::string> &ids)
    {
        std::string name_suffix(".dmbin");
        size_t find = fname.find(name_suffix);
        std::string fname2(fname);
        if( find == std::string::npos ) // there is no suffix
            fname2 = fname + name_suffix;
        
        matr.fwrite(fname2, ids);
    }
    template void Utilities2::fwrite_matrix(const std::string &fname, evolm::matrix<float> &matr, std::vector<std::string> &ids);
    template void Utilities2::fwrite_matrix(const std::string &fname, evolm::matrix<double> &matr, std::vector<std::string> &ids);
    //===========================================================================================
    template <typename T>
    void Utilities2::fwrite_matrix(const std::string &fname, std::vector<T> &vals, std::vector<size_t> &keys, std::vector<std::string> &ids, size_t n_cols)
    {
        //int integer_var;
        float float_var;
        double double_var;

        // determine type of the matrix, expectd: float or double
        //const std::type_info &type_1 = typeid(integer_var);
        const std::type_info &type_2 = typeid(float_var);
        const std::type_info &type_3 = typeid(double_var);
        const std::type_info &type_T = typeid(vals[0]);
        
        size_t var_type = 0; // type of matrix
        size_t var_type2 = 4; // string type of IDs

        if (type_T == type_2)
            var_type = 2;
        else if (type_T == type_3)
            var_type = 3;
        else
            throw std::string("Utilities2::fwrite_matrix(const std::string &, std::vector<T> &, std::vector<size_t> &, std::vector<T2> &, size_t): Cannot determine the type of values vector.");

        std::string name_suffix(".corbin");
        size_t find = fname.find(name_suffix);
        std::string fname2(fname);
        if( find == std::string::npos ) // there is no suffix
            fname2 = fname + name_suffix;

        std::fstream fA;
        fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        fA.open(fname2, fA.binary | fA.trunc | fA.out);

        if (!fA.is_open())
            throw std::string("Utilities2::fwrite_matrix(const std::string &, std::vector<T> &, std::vector<size_t> &, std::vector<T2> &, size_t): Error while opening a binary file.");

        size_t B[6]; // necessary data info

        // info about matrix values
        B[0] = matrix_type;
        B[1] = var_type; // expected 2 or 3
        
        // info about matrix IDs
        B[2] = var_type2; // expected 4
        B[3] = (size_t)sizeof(T);
        B[4] = ids.size();
        B[5] = n_cols;
        
        size_t vals_size = vals.size();
        size_t keys_size = keys.size();
        size_t ids_size = ids.size();

        fA.write( reinterpret_cast<const char *>(&B), 6 * sizeof(size_t) ); // 0. writing the data info

        fA.write( reinterpret_cast<const char *>( &vals_size ), sizeof(size_t) ); // 1. writing size of values
        fA.write( reinterpret_cast<const char *>( vals.data() ), vals_size * sizeof(T) ); // 2. writing all values

        fA.write( reinterpret_cast<const char *>( &keys_size ), sizeof(size_t) ); // 3. writing size of keys
        fA.write( reinterpret_cast<const char *>( keys.data() ), keys_size * sizeof(size_t) ); // 4. writing all keys

        fA.write( reinterpret_cast<const char *>( &ids_size ), sizeof(size_t) ); // 5. writing size of ids
        for (size_t i = 0; i < ids_size; i++)
        {
            std::string str = ids[i];
            size_t str_size = str.size();
            fA.write( reinterpret_cast<const char *>( &str_size ), sizeof(size_t) ); // 5. writing size of ids
            fA.write( reinterpret_cast<const char *>( &str[0] ), str_size ); // 6. writing all IDs
        }

        fA.close();
    }
    template void Utilities2::fwrite_matrix(const std::string &fname, std::vector<float> &vals, std::vector<size_t> &keys, std::vector<std::string> &ids, size_t n_cols);
    template void Utilities2::fwrite_matrix(const std::string &fname, std::vector<double> &vals, std::vector<size_t> &keys, std::vector<std::string> &ids, size_t n_cols);
    //===========================================================================================
    template <typename T>
    void Utilities2::fread_matrix(const std::string &fname, std::vector<T> &vals, std::vector<size_t> &keys, std::vector<std::int64_t> &ids)
    {
        std::fstream fA;
        fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        std::string name_suffix(".corbin");
        size_t find = fname.find(name_suffix);
        std::string fname2(fname);
        if( find == std::string::npos ) // there is no suffix
            fname2 = fname + name_suffix;
        
        fA.open(fname2, fA.binary | fA.in);

        if (!fA.is_open())
            throw std::string("Utilities2::fread_matrix(const std::string &, std::vector<T> &, std::vector<size_t> &, std::vector<T2> &): Error while opening a binary file.");

        size_t B[6];

        fA.read( reinterpret_cast<char *>(&B), 6 * sizeof(size_t) ); // 0. reading a storage info

        // info about matrix values
        size_t var_inbytes = B[3];
        size_t var_type = B[1]; // expected 2 or 3
        
        // info about matrix IDs
        size_t var_type2 = B[2]; // expected 4 or 5

        size_t vals_size;
        size_t keys_size;
        size_t ids_size;

        fA.read( reinterpret_cast<char *>( &vals_size ), sizeof(size_t) ); // 1. reading size of values
        vals.resize(vals_size);
        fA.read( reinterpret_cast<char *>( vals.data() ), vals_size * var_inbytes ); // 2. reading all values

        fA.read( reinterpret_cast<char *>( &keys_size ), sizeof(size_t) ); // 3. reading size of keys
        keys.resize(keys_size);
        fA.read( reinterpret_cast<char *>( keys.data() ), keys_size * sizeof(size_t) ); // 4. reading all keys

        fA.read( reinterpret_cast<char *>( &ids_size ), sizeof(size_t) ); // 3. reading size of keys
        ids.resize(ids_size);
        fA.read( reinterpret_cast<char *>( ids.data() ), ids_size * sizeof(std::int64_t) ); // 4. reading all keys

        fA.close();
    }
    template void Utilities2::fread_matrix(const std::string &fname, std::vector<float> &vals, std::vector<size_t> &keys, std::vector<std::int64_t> &ids);
    template void Utilities2::fread_matrix(const std::string &fname, std::vector<double> &vals, std::vector<size_t> &keys, std::vector<std::int64_t> &ids);
    //===========================================================================================
    template <typename T>
    void Utilities2::fread_matrix(const std::string &fname, evolm::matrix<T> &matr, std::vector<std::int64_t> &ids)
    {
        std::string name_suffix(".dmbin");
        size_t find = fname.find(name_suffix);
        std::string fname2(fname);
        if( find == std::string::npos ) // there is no suffix
            fname2 = fname + name_suffix;
        
        matr.fread(fname2, ids);
    }
    template void Utilities2::fread_matrix(const std::string &fname, evolm::matrix<float> &matr, std::vector<std::int64_t> &ids);
    template void Utilities2::fread_matrix(const std::string &fname, evolm::matrix<double> &matr, std::vector<std::int64_t> &ids);
    //===========================================================================================
    template <typename T>
    void Utilities2::fread_matrix(const std::string &fname, evolm::matrix<T> &matr, std::vector<std::string> &ids)
    {
        std::string name_suffix(".dmbin");
        size_t find = fname.find(name_suffix);
        std::string fname2(fname);
        if( find == std::string::npos ) // there is no suffix
            fname2 = fname + name_suffix;
        
        matr.fread(fname2, ids);
    }
    template void Utilities2::fread_matrix(const std::string &fname, evolm::matrix<float> &matr, std::vector<std::string> &ids);
    template void Utilities2::fread_matrix(const std::string &fname, evolm::matrix<double> &matr, std::vector<std::string> &ids);
    //===========================================================================================
    template <typename T>
    void Utilities2::fread_matrix(const std::string &fname, std::vector<T> &vals, std::vector<size_t> &keys, std::vector<std::string> &ids)
    {
        std::fstream fA;
        fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        std::string name_suffix(".corbin");
        size_t find = fname.find(name_suffix);
        std::string fname2(fname);
        if( find == std::string::npos ) // there is no suffix
            fname2 = fname + name_suffix;

        fA.open(fname2, fA.binary | fA.in);

        if (!fA.is_open())
            throw std::string("Utilities2::fread_matrix(const std::string &, std::vector<T> &, std::vector<size_t> &, std::vector<T2> &): Error while opening a binary file.");

        size_t B[6];

        fA.read( reinterpret_cast<char *>(&B), 6 * sizeof(size_t) ); // 0. reading a storage info

        // info about matrix values
        size_t var_inbytes = B[3];
        size_t var_type = B[1]; // expected 2 or 3
        
        // info about matrix IDs
        size_t var_type2 = B[2]; // expected 4 or 5

        size_t vals_size;
        size_t keys_size;
        size_t ids_size;

        fA.read( reinterpret_cast<char *>( &vals_size ), sizeof(size_t) ); // 1. reading size of values
        vals.resize(vals_size);
        fA.read( reinterpret_cast<char *>( vals.data() ), vals_size * var_inbytes ); // 2. reading all values

        fA.read( reinterpret_cast<char *>( &keys_size ), sizeof(size_t) ); // 3. reading size of keys
        keys.resize(keys_size);
        fA.read( reinterpret_cast<char *>( keys.data() ), keys_size * sizeof(size_t) ); // 4. reading all keys

        fA.read( reinterpret_cast<char *>( &ids_size ), sizeof(size_t) ); // 3. reading size of keys
        ids.resize(ids_size);
        for (size_t i = 0; i < ids_size; i++)
        {
            size_t str_size;
            std::string str;
            fA.read( reinterpret_cast<char *>( &str_size ), sizeof(size_t) ); // 5. writing size of ids
            str.resize(str_size);
            fA.read( reinterpret_cast<char *>( &str[0] ), str_size ); // 6. reading specific id
            ids[i] = str;
        }

        fA.close();
    }
    template void Utilities2::fread_matrix(const std::string &fname, std::vector<float> &vals, std::vector<size_t> &keys, std::vector<std::string> &ids);
    template void Utilities2::fread_matrix(const std::string &fname, std::vector<double> &vals, std::vector<size_t> &keys, std::vector<std::string> &ids);
    //===========================================================================================
    void Utilities2::fread_matrix_info(const std::string &fname, size_t &info)
    {
        std::fstream fA;
        fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        std::string name_suffix(".corbin");
        size_t find = fname.find(name_suffix);
        std::string fname2(fname);
        if( find == std::string::npos ) // there is no suffix
            fname2 = fname + name_suffix;

        fA.open(fname2, fA.binary | fA.in);

        if (!fA.is_open())
            throw std::string("Utilities2::fread_matrix_info(const std::string &, size_t &): Error while opening a binary file.");

        fA.read( reinterpret_cast<char *>(info), 6 * sizeof(size_t) ); // reading a storage info

        fA.close();
    }
    //===============================================================================================================

    size_t Utilities2::fget_matrix_kind(const std::string &fname)
    {
        try
        {
            size_t mat_kind = 0;

            std::fstream fA;
            fA.open(fname, fA.binary | fA.in);

            if (!fA.is_open())
                throw std::string("Error while opening a binary file for reading.");

            fA.read(reinterpret_cast<char *>(&mat_kind), sizeof(size_t));
            
            fA.close();

            return mat_kind;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::fget_matrix_kind(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::fget_matrix_kind(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::fget_matrix_kind(const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Utilities2::get_RecodedIdMap(std::map<std::int64_t, std::int64_t> &id_map,
                                std::vector<std::int64_t> &idVect)
    {
        try
        {
            size_t code_id = 1;
            for (auto const &elem : idVect)
            {
                id_map[elem] = code_id;
                code_id++;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::get_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::get_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::get_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Utilities2::find_RecodedIdMap(std::map<std::int64_t, std::int64_t> &id_map,
                                 std::vector<std::int64_t> &whereIidList,
                                 std::vector<std::int64_t> &whatIdList)
    {
        try
        {
            for (size_t i = 0; i < whatIdList.size(); i++)
                id_map[ whatIdList[i] ] = find_invect( whereIidList, whatIdList[i] ) + 1; // +1 because we start indexing from 1
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::find_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::find_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::find_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    template <typename T>
    void Utilities2::dense_to_sparse(evolm::matrix<T> &from, evolm::smatrix<T> &to)
    {
        try
        {
            evolm::matrix<size_t> shp;
            shp = from.shape();

            T zerro_value = (T)0;

            for (size_t i = 0; i < shp[0]; i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    if (from(i, j) != zerro_value)
                        to(i, j) = from(i, j);
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::dense_to_sparse(evolm::matrix<T> &, evolm::smatrix<T> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::dense_to_sparse(evolm::matrix<T> &, evolm::smatrix<T> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::dense_to_sparse(evolm::matrix<T> &, evolm::smatrix<T> &)" << '\n';
            throw;
        }
    }

    template void Utilities2::dense_to_sparse(evolm::matrix<float> &from, evolm::smatrix<float> &to);
    template void Utilities2::dense_to_sparse(evolm::matrix<double> &from, evolm::smatrix<double> &to);

    //===============================================================================================================

    template <typename T>
    void Utilities2::sparse_to_dense(evolm::smatrix<T> &from, evolm::matrix<T> &to)
    {
        try
        {
            std::vector<size_t> key_list;            
            from.get_keyslist( key_list );

            for (size_t i = 0; i < key_list.size(); i++)
                to[ key_list[i] ] = from[ key_list[i] ];

            key_list.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::sparse_to_dense(evolm::smatrix<T> &, evolm::matrix<T> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::sparse_to_dense(evolm::smatrix<T> &, evolm::matrix<T> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::sparse_to_dense(evolm::smatrix<T> &, evolm::matrix<T> &)" << '\n';
            throw;
        }
    }

    template void Utilities2::sparse_to_dense(evolm::smatrix<float> &from, evolm::matrix<float> &to);
    template void Utilities2::sparse_to_dense(evolm::smatrix<double> &from, evolm::matrix<double> &to);

    //===============================================================================================================

    template <typename T>
    void Utilities2::solve_ls(evolm::matrix<T> &L, T *b, T *x)
    {
        // Solving the linear system L*Lt*x = b by back substitution
        try
        {
            evolm::matrix<size_t> shp_L;
            shp_L = L.shape();

            // Forward iteration
            T *y = new T [ shp_L[1] ];
            y[0] = b[0] / L(0,0);
            for (size_t i = 1; i < shp_L[1]; i++)
            {
                T l_sum = (T)0.0;
                for (size_t j = 0; j < i; j++)
                    l_sum = l_sum + L(i,j) * y[j];
                y[i] = ( b[i] - l_sum ) / L(i,i);
            }

            // Backward iteration
            x[ shp_L[1]-1 ] = y[ shp_L[1]-1 ] / L(shp_L[1]-1,shp_L[1]-1);
            for (std::int64_t i = shp_L[1]-2; i >= 0; i--)
            {
                T l_sum = (T)0.0;
                for (std::int64_t j = shp_L[1]-1; j > i; j--)
                {
                    l_sum = l_sum + L(i,j) * x[j];
                }
                x[i] = ( y[i] - l_sum ) / L(i,i);
            }

            delete [] y;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::solve_ls(evolm::matrix<T> &, T *, T *)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::solve_ls(evolm::matrix<T> &, T *, T *)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::solve_ls(evolm::matrix<T> &, T *, T *)" << '\n';
            throw;
        }
    }

    template void Utilities2::solve_ls(evolm::matrix<float> &L, float *b, float *x);
    template void Utilities2::solve_ls(evolm::matrix<double> &L, double *b, double *x);

    //===============================================================================================================

    template <typename T>
    bool Utilities2::is_float(T val)
    {
        bool isFlt = false;

        try
        {
            float flt_var;
            double dbl_var;
            const std::type_info &ti1 = typeid(flt_var);
            const std::type_info &ti2 = typeid(dbl_var);
            const std::type_info &ti3 = typeid(val);

            if (ti3 == ti1)
                isFlt = true;

            if ( (ti3 != ti1) && (ti3 != ti2) )
                throw std::string("Cannot detect the type of variable.");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilities2::is_float(T)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Utilities2::is_float(T)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilities2::is_float(T)" << '\n';
            throw;
        }

        return isFlt;
    }

    template bool Utilities2::is_float(float val);
    template bool Utilities2::is_float(double val);

    //===============================================================================================================

} // end of namespace evoped