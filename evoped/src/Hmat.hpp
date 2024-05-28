#ifndef Hmat_hpp__
#define Hmat_hpp__

#include "Utilities2.hpp"

namespace evoped
{
    template <typename T>
    class Hmat
    {
    public:
        Hmat();
        ~Hmat();

        void make_matrix(const std::string &a_matr, // supposed to be python interface, needs rewriting!
                         const std::string &a_ids,
                         const std::string &a_red_matr,
                         const std::string &a_red_ids,
                         const std::string &g_matr,
                         const std::string &g_ids);

        void make_matrix(evolm::matrix<T>& a_matr,
                         std::vector<std::int64_t>& a_ids,
                         evolm::matrix<T>& a_red_matr,
                         std::vector<std::int64_t>& a_red_ids,
                         evolm::matrix<T>& g_matr,
                         std::vector<std::int64_t>& g_ids);

        void make_matrix(evolm::smatrix<T>& a_matr,
                         std::vector<std::int64_t>& a_ids,
                         evolm::matrix<T>& a_red_matr,
                         std::vector<std::int64_t>& a_red_ids,
                         evolm::matrix<T>& g_matr,
                         std::vector<std::int64_t>& g_ids);

        void make_matrix(evolm::matrix<T>& a_matr,
                         std::vector<std::int64_t>& a_ids,
                         evolm::smatrix<T>& a_red_matr,
                         std::vector<std::int64_t>& a_red_ids,
                         evolm::matrix<T>& g_matr,
                         std::vector<std::int64_t>& g_ids);

        void make_matrix(evolm::smatrix<T>& a_matr,
                         std::vector<std::int64_t>& a_ids,
                         evolm::smatrix<T>& a_red_matr,
                         std::vector<std::int64_t>& a_red_ids,
                         evolm::matrix<T>& g_matr,
                         std::vector<std::int64_t>& g_ids);

        void get_matrix(evolm::matrix<T>& arr,
                        std::vector<std::int64_t>& ids);

        void get_matrix(evolm::smatrix<T>& arr,
                        std::vector<std::int64_t>& ids);

        void save_matrix(const std::string & arr,
                        const std::string & ids);

        void clear();

    private:

        struct EmptyDataStatus
        {
            bool H;
            bool H_s;
        };
        EmptyDataStatus IsEmpty;

        evolm::matrix<T> H;
        evolm::smatrix<T> H_s;
        std::vector<std::int64_t> hmat_id;
    };

} // end of namespace evoped

#endif // Hmat_hpp__