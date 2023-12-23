#include "Hmat.hpp"

namespace evoped
{
    //===============================================================================================================

    Hmat::Hmat()
    {
    }

    //===============================================================================================================

    Hmat::~Hmat()
    {
    }
    //===============================================================================================================

    void Hmat::check_id(std::vector<std::int64_t> &id_list, std::vector<std::int64_t> &checked_id, std::vector<std::int64_t> &missing_id)
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
            std::cerr << "Exception in Hmat::check_id(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::check_id(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::check_id(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    bool Hmat::is_value_in_vect(std::vector<std::int64_t> &where_tocheck, std::vector<std::int64_t> &what_tocheck)
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
            std::cerr << "Exception in Hmat::is_value_in_vect(std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::is_value_in_vect(std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::is_value_in_vect(std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

}