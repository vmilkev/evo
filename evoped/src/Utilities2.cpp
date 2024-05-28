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

            for (size_t i = 0; i < checked_id.size(); i++)
                if (!std::binary_search(id_list.begin(), id_list.end(), checked_id[i]))
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

            for (size_t i = 0; i < checked_id.size(); i++)
                if ( std::binary_search(id_list.begin(), id_list.end(), checked_id[i]) )
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
            std::vector<int64_t> vec_1(where_tocheck);
            std::vector<int64_t> vec_2(what_tocheck);
            
            std::sort(vec_1.begin(), vec_1.end());
            std::sort(vec_2.begin(), vec_2.end());
            //----------------------------------
            bool out = true;

            std::vector<std::int64_t> missing;
            //check_id(where_tocheck, what_tocheck, missing);
            check_id(vec_1, vec_2, missing);

            if (!missing.empty())
            {
                /*std::cout<<"missing size: "<<missing.size()<<"\n";
                size_t sz = std::min((size_t)100, missing.size());
                for (size_t i = 0; i<sz; i++)
                    std::cout<<"missing: "<<missing[i]<<"\n";*/

                out = false;
                missing.clear();
            }

            vec_1.clear();
            vec_2.clear();

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

} // end of namespace evoped