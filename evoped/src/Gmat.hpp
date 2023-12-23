#ifndef Gmat_hpp__
#define Gmat_hpp__

#include <vector>
#include <map>
#include <iostream>
#include <fstream>

namespace evoped
{
    class Gmat
    {
    public:
        Gmat();
        ~Gmat();

        template <typename T>
        void read_gmatrix(const std::string &gmat_file, std::vector<std::int64_t> &g_row, std::vector<std::int64_t> &g_col, std::vector<T> &g_val);

    private:
        
        std::vector<std::int64_t> gmatID;      /* container for the list of G matrix IDs */

        void find_RecodedIdMap(std::map<std::int64_t, std::int64_t> &id_map,
                               std::vector<std::int64_t> &whereIidList,
                               std::vector<std::int64_t> &whatIdList); // not sure if this is needed

        void get_RecodedIdMap(std::map<std::int64_t, std::int64_t> &id_map,
                              std::vector<std::int64_t> &idVect); // not sure if this is needed
        
        size_t find_invect(std::vector<std::int64_t> &where, std::int64_t what); // was find_invect2

         bool is_unique(std::vector<std::int64_t> &x);

    };

} // end of namespace evoped

#endif // Gmat_hpp__