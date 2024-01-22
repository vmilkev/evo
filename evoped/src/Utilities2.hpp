#ifndef Utilities2_hpp__
#define Utilities2_hpp__

#include <vector>
#include <map>
#include <iostream>

namespace evoped
{
    class Utilities2
    {
    private:
        
    public:
        Utilities2();
        ~Utilities2();

        void check_id(std::vector<std::int64_t> &id_list,
                      std::vector<std::int64_t> &checked_id,
                      std::vector<std::int64_t> &missing_id);
        bool is_value_in_vect(std::vector<std::int64_t> &where_tocheck,
                              std::vector<std::int64_t> &what_tocheck);
        int find_invect(std::vector<std::int64_t> &where,
                           std::int64_t what);
        bool is_invect(std::vector<std::int64_t> &where,
                       std::int64_t what);               // was find_invect
        bool is_unique(std::vector<std::int64_t> &x);
        void get_RecodedIdMap(std::map<std::int64_t,std::int64_t> &id_map,
                              std::vector<std::int64_t> &idVect);
        void find_RecodedIdMap(std::map<std::int64_t, std::int64_t> &id_map,
                               std::vector<std::int64_t> &whereIidList,
                               std::vector<std::int64_t> &whatIdList);
    };

} // end of namespace evoped

#endif // Utilities2_hpp__