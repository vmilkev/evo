#include "Utilites.hpp"

namespace evogen
{
    //===============================================================================================================

    Utilites::Utilites()
    {
        try
        {
            // std::random_device rd;
            // std::mt19937 generator( rd() );
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

    unsigned long Utilites::get_randi(unsigned long range)
    {
        unsigned long rnum = 0;

        try
        {
            std::random_device rd;
            srand(rd());

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
            std::random_device rd;
            std::mt19937 generator(rd());

            std::uniform_int_distribution<T> distribution(low_bound, upp_bound);

            distribution.reset();

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
    template unsigned short Utilites::get_randi(unsigned short low_bound, unsigned short upp_bound);
    template unsigned char Utilites::get_randi(unsigned char low_bound, unsigned char upp_bound);
    template unsigned int Utilites::get_randi(unsigned int low_bound, unsigned int upp_bound);
    template long Utilites::get_randi(long low_bound, long upp_bound);
    template unsigned long Utilites::get_randi(unsigned long low_bound, unsigned long upp_bound);
    template unsigned long long Utilites::get_randi(unsigned long long low_bound, unsigned long long upp_bound);

    //===============================================================================================================

    template <typename T>
    std::vector<T> Utilites::get_uni_rand(size_t n_values, T low_bound, T upp_bound, bool reset)
    {
        std::vector<T> rnum_vect;

        try
        {
            std::random_device rd;
            std::mt19937 generator(rd());

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
            std::cerr << "Exception in Utilites::get_uni_rand( size_t, T, T, bool )." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilites::get_uni_rand( size_t, T, T, bool )." << '\n';
            throw;
        }

        return rnum_vect;
    }

    template std::vector<short> Utilites::get_uni_rand(size_t n_values, short low_bound, short upp_bound, bool reset);
    template std::vector<int> Utilites::get_uni_rand(size_t n_values, int low_bound, int upp_bound, bool reset);
    template std::vector<ushort> Utilites::get_uni_rand(size_t n_values, ushort low_bound, ushort upp_bound, bool reset);
    template std::vector<uint> Utilites::get_uni_rand(size_t n_values, uint low_bound, uint upp_bound, bool reset);
    template std::vector<long> Utilites::get_uni_rand(size_t n_values, long low_bound, long upp_bound, bool reset);
    template std::vector<unsigned long> Utilites::get_uni_rand(size_t n_values, unsigned long low_bound, unsigned long upp_bound, bool reset);
    template std::vector<unsigned long long> Utilites::get_uni_rand(size_t n_values, unsigned long long low_bound, unsigned long long upp_bound, bool reset);

    //===============================================================================================================

    template <typename T>
    std::vector<T> Utilites::get_runi_rand(size_t n_values, T low_bound, T upp_bound, bool reset)
    {
        std::vector<T> rnum_vect;

        try
        {
            std::random_device rd;
            std::mt19937 generator(rd());

            std::uniform_real_distribution<T> distribution(low_bound, upp_bound);

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
            std::cerr << "Exception in Utilites::get_runi_rand( size_t, T, T, bool )." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilites::get_runi_rand( size_t, T, T, bool )." << '\n';
            throw;
        }

        return rnum_vect;
    }

    template std::vector<double> Utilites::get_runi_rand(size_t n_values, double low_bound, double upp_bound, bool reset);
    template std::vector<float> Utilites::get_runi_rand(size_t n_values, float low_bound, float upp_bound, bool reset);

    //===============================================================================================================

    template <typename T>
    std::vector<T> Utilites::get_norm_rand(size_t n_values, T dist_mean, T dist_std, bool reset)
    {
        std::vector<T> rnum_vect;

        try
        {
            std::random_device rd;
            std::mt19937 generator(rd());

            std::normal_distribution<T> distribution(dist_mean, dist_std);

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
            std::cerr << "Exception in Utilites::get_norm_rand( size_t, T, T, bool )." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilites::get_norm_rand( size_t, T, T, bool )." << '\n';
            throw;
        }

        return rnum_vect;
    }

    template std::vector<double> Utilites::get_norm_rand(size_t n_values, double dist_mean, double dist_std, bool reset);
    template std::vector<float> Utilites::get_norm_rand(size_t n_values, float dist_mean, float dist_std, bool reset);

    //===============================================================================================================

    template <typename T>
    std::vector<T> Utilites::get_gamma_rand(size_t n_values, T alpha, T beta, bool reset)
    {
        /*
           n_values - number of random numbers generated;
           alpha - the shape parameter, alpha > 0;
           beta - scale parameter, beta > 0;
           reset - resets the distribution, if true - all the numbers will be from different distributions )
                   independent, if false - all numbers are from the same distribution.
        */

        std::vector<T> rnum_vect;

        try
        {
            std::random_device rd;
            std::mt19937 generator(rd());

            std::gamma_distribution<T> distribution(alpha, beta);

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
            std::cerr << "Exception in Utilites::get_gamma_rand(size_t, T, T, bool)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilites::get_gamma_rand(size_t, T, T, bool)." << '\n';
            throw;
        }

        return rnum_vect;
    }

    template std::vector<double> Utilites::get_gamma_rand(size_t n_values, double alpha, double beta, bool reset);
    template std::vector<float> Utilites::get_gamma_rand(size_t n_values, float alpha, float beta, bool reset);

    //===============================================================================================================

    std::vector<int> Utilites::get_bin_rand(size_t n_values, int n, double p, bool reset)
    {
        /*
           n_values - number of random numbers generated;
           n - number of trials, n > 0;
           p - probability of success, 0 <= p <= 1;
           reset - resets the distribution, if true - all the numbers will be from different distributions )
                   independent, if false - all numbers are from the same distribution.
        */

        std::vector<int> rnum_vect;

        try
        {
            std::random_device rd;
            std::mt19937 generator(rd());

            std::binomial_distribution<int> distribution(n, p);

            for (size_t i = 0; i < n_values; i++)
            {
                if (reset)
                    distribution.reset();

                int number = distribution(generator);
                rnum_vect.push_back(number);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilites::get_bin_rand(size_t, int, double, bool)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilites::get_bin_rand(size_t, int, double, bool)." << '\n';
            throw;
        }

        return rnum_vect;
    }

    //===============================================================================================================

    int Utilites::bin_rand(int n, double p)
    {
        /*
            n - number of trials, n > 0;
            p - probability of success, 0 <= p <= 1;
            reset - resets the distribution, if true - all the numbers will be from different distributions )
                    independent, if false - all numbers are from the same distribution.
        */

        int rnum_vect = 0;

        try
        {
            std::random_device rd;
            std::mt19937 generator(rd());

            std::binomial_distribution<int> distribution(n, p);

            //if (reset)
            //    distribution.reset();

            rnum_vect = distribution(generator);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Utilites::get_bin_rand(size_t, int, double, bool)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Utilites::get_bin_rand(size_t, int, double, bool)." << '\n';
            throw;
        }

        return rnum_vect;
    }
    
    //===============================================================================================================

    void Utilites::fremove(std::string file_name)
    {
        try
        {
            std::filesystem::remove(file_name);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::fremove(std::string)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::fremove(std::string)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::fremove(std::string)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

}