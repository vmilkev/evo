#ifndef group_hpp__
#define group_hpp__

#include "Population.hpp"
#include <vector>
#include <functional>

namespace evogen
{
    class Group
    {
    public:
        Group();
        ~Group();

        size_t size();
        size_t size_at(size_t at);

        void clear();
        void remove();
        void move(Population &pop);
        void add(Population &pop); // add entire population to the group
        void add(Population &pop, size_t which_one);
        void add(Group &grp);

        friend class Trait;

    private:
        std::vector< std::reference_wrapper<Population> > pop_list;
        std::vector< std::vector<size_t> > individuals_list;

        Population & get_population( size_t which_pop );
        std::vector<size_t> get_individuals( size_t which_pop );
        
    protected:
    };

} // end of namespace evogen

#endif // group_hpp__