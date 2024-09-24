#ifndef group_hpp__
#define group_hpp__

#include "Population.hpp"
#include "Trait.hpp"
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
        void add(Population &pop);          // add entire population to the group
        void add(Population &pop,
                 size_t which_one);
        void add(Group &grp);
        void mate();                        // default: sexual reproduction, 1 - offspring with 1.0 success rate
        void mate(bool sexual_reproduction, // true - sexual, false - asexsual;
                  int max_offspring,        // max number of offspring from the same mating event
                  float success_rate);     // probability of getting max_offspring value (binomial distribution)
        void mate(bool sexual_reproduction,
                  float max_offspring,
                  float success_rate);
        void regroup_newborn(Group &grp);

        void aging(int delta_t);
        void genotype();
        void kill();

#ifdef PYBIND
        void make_observation(Trait &trt,
                              pybind11::array_t<float> py_env);

        void make_observation(Trait &trt,
                              pybind11::array_t<float> py_env,
                              const std::string &out_trvalues);

        void make_observation(Trait &trt,
                              pybind11::array_t<float> py_env,
                              const std::string &out_trvalues,
                              const std::string &out_genotypes);

        void make_observation(Trait &trt,
                              pybind11::array_t<float> py_env,
                              pybind11::array_t<float> py_trvalues);

        void make_observation(Trait &trt,
                              pybind11::array_t<float> py_env,
                              pybind11::array_t<float> py_trvalues,
                              pybind11::array_t<int> py_genotypes);
#else
        void make_observation(Trait &trt,
                              std::vector<float> &env);

        void make_observation(Trait &trt,
                              std::vector<float> &env,
                              const std::string &out_trvalues);

        void make_observation(Trait &trt,
                              std::vector<float> &env,
                              const std::string &out_trvalues,
                              const std::string &out_genotypes);

        void make_observation(Trait &trt,
                              std::vector<float> &env,
                              std::vector<std::vector<float>> &out_trvalues);

        void make_observation(Trait &trt,
                              std::vector<float> &env,
                              std::vector<std::vector<float>> &out_trvalues,
                              std::vector<std::vector<short>> &out_genotypes);
#endif

    private:
        std::vector<std::reference_wrapper<Population>> pop_list;
        std::vector<std::vector<size_t>> individuals_list;

        std::vector<std::reference_wrapper<Population>> pop_list_newborn;
        std::vector<std::vector<size_t>> individuals_list_newborn;

        Population &get_population(size_t which_pop);
        std::vector<size_t> get_individuals(size_t which_pop);
        void add_newborn(Population &pop, size_t which_one);

    protected:
    };

} // end of namespace evogen

#endif // group_hpp__