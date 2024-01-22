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

    bool Utilities2::is_value_in_vect(std::vector<std::int64_t> &where_tocheck, std::vector<std::int64_t> &what_tocheck)
    {
        // are values from what_tocheck in where_tocheck?
        try
        {
            bool out = true;

            std::vector<std::int64_t> missing;
            check_id(where_tocheck, what_tocheck, missing);

            if (!missing.empty())
            {
                out = false;
                missing.clear();
            }

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

} // end of namespace evoped