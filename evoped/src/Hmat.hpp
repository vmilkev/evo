#ifndef Hmat_hpp__
#define Hmat_hpp__

#include "cs_matrix.hpp"
#include "Utilities2.hpp"

namespace evoped
{
    class Hmat
    {
    public:
        Hmat();
        ~Hmat();

        void make_matrix(std::vector<double>& a_matr,
                         std::vector<std::int64_t>& a_ids,
                         std::vector<double>& a_red_matr,
                         std::vector<std::int64_t>& a_red_ids,
                         std::vector<double>& g_matr,
                         std::vector<std::int64_t>& g_ids);
        void make_matrix(evolm::matrix<double>& a_matr,
                         std::vector<std::int64_t>& a_ids,
                         evolm::matrix<double>& a_red_matr,
                         std::vector<std::int64_t>& a_red_ids,
                         evolm::matrix<double>& g_matr,
                         std::vector<std::int64_t>& g_ids);
        void get_matrix(evolm::matrix<double>& arr,
                        std::vector<std::int64_t>& ids,
                        bool keep_ondisk);
        void get_matrix(std::vector<double>& arr,
                        std::vector<std::int64_t>& ids,
                        bool keep_ondisk);
        void clear();

    private:

        evolm::matrix<double> H;
        std::vector<std::int64_t> hmat_id;
    };

} // end of namespace evoped

#endif // Hmat_hpp__