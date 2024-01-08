#ifndef Hmat_hpp__
#define Hmat_hpp__

#include <vector>
#include <map>
#include <iostream>

#include "cs_matrix.hpp"

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

        void get_submatrix(std::vector<double>& full_matr,
                           std::vector<std::int64_t>& matr_ids,
                           std::vector<std::int64_t>& selected_ids,
                           bool use_invverse);
        void get_submatrix(evolm::matrix<double>& full_matr,
                           std::vector<std::int64_t>& matr_ids,
                           std::vector<std::int64_t>& selected_ids,
                           bool use_invverse);
        void bin_write();
        void bin_read();
        void get_matrix(evolm::matrix<double>& arr);
        void get_matrix(std::vector<double>& arr);
        void get_ids(std::vector<std::int64_t>& ids);
        void clear();

    private:

        evolm::matrix<double> H;
        std::vector<std::int64_t> hmat_id;

        void check_id(std::vector<std::int64_t> &id_list,
                      std::vector<std::int64_t> &checked_id,
                      std::vector<std::int64_t> &missing_id);
        bool is_value_in_vect(std::vector<std::int64_t> &where_tocheck,
                              std::vector<std::int64_t> &what_tocheck);
        int find_invect(std::vector<std::int64_t> &where,
                           std::int64_t what);       
    };

} // end of namespace evoped

#endif // Hmat_hpp__