#include "Utilites.hpp"

namespace evogen
{
    //===============================================================================================================

    Utilites::Utilites()
    {
        try
        {
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilites::Utilites()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilites::Utilites()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    unsigned long long Utilites::rdtsc()
    {
        unsigned long long number = 0;

        try
        {
            /* Seed for random number generator. */

            unsigned int lo, hi;
            __asm__ __volatile__("rdtsc"
                                 : "=a"(lo), "=d"(hi));

            number = ((unsigned long long)hi << 32) | lo;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilites::rdtsc()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilites::rdtsc()." << '\n';
            throw;
        }

        return number;
    }

    //===============================================================================================================

    unsigned long Utilites::get_randi(unsigned long range)
    {
        unsigned long rnum = 0;

        try
        {
            srand(static_cast<unsigned int>(rdtsc()));
            rnum = std::rand() % range;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilites::get_randi( unsigned long )." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilites::get_randi( unsigned long )." << '\n';
            throw;
        }

        return rnum;
    }

    //===============================================================================================================

    template <typename T>
    T Utilites::get_randi(T low_bound, T upp_bound)
    {
        T rnum = 0;

        try
        {
            std::mt19937 generator(static_cast<unsigned int>(rdtsc()));

            std::uniform_int_distribution<T> distribution(low_bound, upp_bound);

            // distribution.reset();

            rnum = distribution(generator);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilites::get_randi( T, T )." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilites::get_randi( T, T )." << '\n';
            throw;
        }

        return rnum;
    }

    template short Utilites::get_randi(short low_bound, short upp_bound);
    template int Utilites::get_randi(int low_bound, int upp_bound);
    template ushort Utilites::get_randi(ushort low_bound, ushort upp_bound);
    template uint Utilites::get_randi(uint low_bound, uint upp_bound);
    template long Utilites::get_randi(long low_bound, long upp_bound);
    template ulong Utilites::get_randi(ulong low_bound, ulong upp_bound);
    template unsigned long long Utilites::get_randi(unsigned long long low_bound, unsigned long long upp_bound);

    //===============================================================================================================

    template <typename T>
    std::vector<T> Utilites::get_uni_rand(size_t n_values, T low_bound, T upp_bound, bool reset)
    {
        std::vector<T> rnum_vect;

        try
        {
            std::mt19937 generator(static_cast<unsigned int>(rdtsc()));

            std::uniform_int_distribution<T> distribution(low_bound, upp_bound);

            for (size_t i = 0; i < n_values; i++)
            {
                if (reset)
                    distribution.reset();

                T number = distribution(generator);
                rnum_vect.push_back(number);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilites::get_randi( size_t, T, T, bool )." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilites::get_randi( size_t, T, T, bool )." << '\n';
            throw;
        }

        return rnum_vect;
    }

    template std::vector<short> Utilites::get_uni_rand(size_t n_values, short low_bound, short upp_bound, bool reset);
    template std::vector<int> Utilites::get_uni_rand(size_t n_values, int low_bound, int upp_bound, bool reset);
    template std::vector<ushort> Utilites::get_uni_rand(size_t n_values, ushort low_bound, ushort upp_bound, bool reset);
    template std::vector<uint> Utilites::get_uni_rand(size_t n_values, uint low_bound, uint upp_bound, bool reset);
    template std::vector<long> Utilites::get_uni_rand(size_t n_values, long low_bound, long upp_bound, bool reset);
    template std::vector<ulong> Utilites::get_uni_rand(size_t n_values, ulong low_bound, ulong upp_bound, bool reset);
    template std::vector<unsigned long long> Utilites::get_uni_rand(size_t n_values, unsigned long long low_bound, unsigned long long upp_bound, bool reset);

    //===============================================================================================================

    std::vector<double> Utilites::get_gamma_rand(size_t n_values, double alpha, double beta, bool reset)
    {
        /*
           n_values - number of random numbers generated;
           alpha - the shape parameter, alpha > 0;
           beta - scale parameter, beta > 0;
           reset - resets the distribution, if true - all the numbers will be from different distributions )
                   independent, if false - all numbers are from the same distribution.
        */

        std::vector<double> rnum_vect;

        try
        {
            std::mt19937 generator(static_cast<unsigned int>(rdtsc()));

            std::gamma_distribution<double> distribution(alpha, beta);

            for (size_t i = 0; i < n_values; i++)
            {
                if (reset)
                    distribution.reset();

                double number = distribution(generator);
                rnum_vect.push_back(number);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilites::get_gamma_rand(size_t, double, double, bool)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilites::get_gamma_rand(size_t, double, double, bool)." << '\n';
            throw;
        }

        return rnum_vect;
    }

    //===============================================================================================================

}