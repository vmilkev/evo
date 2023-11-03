#ifndef utilites_hpp__
#define utilites_hpp__

#include <iostream>
#include <random>
#include <chrono>

namespace evogen
{
    class Utilites
    {
    public:
        Utilites();

        unsigned long long rdtsc();
        unsigned long get_randi(unsigned long range);

        template <typename T>
        T get_randi(T low_bound,
                    T upp_bound); /* returns a uniformly distributed random integer within the specified interval */

        template <typename T>
        std::vector<T> get_uni_rand(size_t n_values,
                                    T low_bound,
                                    T upp_bound,
                                    bool reset); /* returns a vector of uniformly distributed random integer within the specified interval */

        std::vector<double> get_gamma_rand(size_t n_values,
                                           double alpha,
                                           double beta,
                                           bool reset); /* returns a vector of gamma distributed random value */

    private:

    protected:
    };

} // end of namespace evogen

#endif // utilites_hpp__