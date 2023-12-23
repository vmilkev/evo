#ifndef Hmat_hpp__
#define Hmat_hpp__

#include <vector>
#include <map>
#include <iostream>

namespace evoped
{
    class Hmat
    {
    public:
        Hmat();
        ~Hmat();

    private:
        void check_id(std::vector<std::int64_t> &id_list,
                      std::vector<std::int64_t> &checked_id,
                      std::vector<std::int64_t> &missing_id);
        /* methods to check if IDs are in pedigree and in G matrix */        
        bool is_value_in_vect(std::vector<std::int64_t> &where_tocheck, std::vector<std::int64_t> &what_tocheck);       // Not sure if this is needed!
       
    };

} // end of namespace evoped

#endif // Hmat_hpp__