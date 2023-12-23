#include "Gmat.hpp"

namespace evoped
{
    //===============================================================================================================

    Gmat::Gmat()
    {
    }

    //===============================================================================================================

    Gmat::~Gmat()
    {
    }

    //===============================================================================================================
    template <typename T>
    void Gmat::read_gmatrix(const std::string &gmat_file, std::vector<std::int64_t> &g_row, std::vector<std::int64_t> &g_col, std::vector<T> &g_val)
    {
        // reads G matrix from gmat.dat file into the linked list container

        std::ifstream ped;

        try
        {
            std::string where("Error during reading G matrix");
            std::string line;

            char *end;
            const char *p;
            std::vector<T> tmp_list;

            ped.exceptions(std::ifstream::failbit | std::ifstream::badbit);

            size_t diagonals = 0; /* for debugging */

            ped.open(gmat_file, std::fstream::in);

            if (ped)
            {
                while (getline(ped, line))
                {
                    p = line.c_str();
                    for (double f = std::strtod(p, &end); p != end; f = std::strtod(p, &end))
                    {
                        p = end;
                        if (errno == ERANGE)
                        {
                            throw std::string("Range error during reading G matrix");
                            errno = 0;
                        }
                        tmp_list.push_back(static_cast<T>(f));
                    }
                    g_row.push_back(static_cast<std::int64_t>(tmp_list[0]));
                    g_col.push_back(static_cast<std::int64_t>(tmp_list[1]));
                    g_val.push_back(tmp_list[2]);
                    gmatID.push_back(int(tmp_list[0]));
                    gmatID.push_back(int(tmp_list[1]));

                    if (static_cast<std::int64_t>(tmp_list[0]) == static_cast<std::int64_t>(tmp_list[1]))
                        diagonals++;

                    tmp_list.erase(tmp_list.begin(), tmp_list.end());

                    if (ped.eof())
                    {
                        ped.close();

                        if (!is_unique(gmatID))
                            gmatID.erase(unique(gmatID.begin(), gmatID.end()), gmatID.end()); // here the vector should be sorted and unique

                        if (diagonals != gmatID.size())
                            throw std::string("There are missing diagonals in G-matrix file.");
                        else
                            return;
                    }
                }
                ped.close();

                if (!is_unique(gmatID))
                    gmatID.erase(unique(gmatID.begin(), gmatID.end()), gmatID.end()); // here the vector should be sorted and unique
            }
            else
                throw std::string("Cannot open G-matrix file!");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::read_gmatrix(const std::string &, std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<T> &)" << '\n';
            std::cerr << e.what() << '\n';
            if (ped.eof())
                std::cerr << "Check for empty line(s) at the endd of the file!"
                          << "\n";
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::read_gmatrix(const std::string &, std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<T> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::read_gmatrix(const std::string &, std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<T> &)" << '\n';
            throw;
        }
    }

    template void Gmat::read_gmatrix(const std::string &fname, std::vector<std::int64_t> &g_row, std::vector<std::int64_t> &g_col, std::vector<float> &g_val);
    template void Gmat::read_gmatrix(const std::string &fname, std::vector<std::int64_t> &g_row, std::vector<std::int64_t> &g_col, std::vector<double> &g_val);

    //===============================================================================================================

    size_t Gmat::find_invect(std::vector<std::int64_t> &where, std::int64_t what)
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
            std::cerr << "Exception in Gmat::find_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::find_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::find_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::find_RecodedIdMap(std::map<std::int64_t, std::int64_t> &id_map,
                                 std::vector<std::int64_t> &whereIidList,
                                 std::vector<std::int64_t> &whatIdList)
    {
        try
        {
            for (size_t i = 0; i < whatIdList.size(); i++)
                id_map[whatIdList[i]] = find_invect(whereIidList, whatIdList[i]) + 1;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::find_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::find_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::find_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::get_RecodedIdMap(std::map<std::int64_t, std::int64_t> &id_map,
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
            std::cerr << "Exception in Gmat::get_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::get_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::get_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    bool Gmat::is_unique(std::vector<std::int64_t> &x)
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
            std::cerr << "Exception in Gmat::is_unique(std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::is_unique(std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::is_unique(std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================
    //===============================================================================================================
    //===============================================================================================================
}