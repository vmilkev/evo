#ifndef Utilities2_hpp__
#define Utilities2_hpp__

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "dense_matrix.hpp"
#include "sparse_matrix.hpp"

namespace evoped
{
    class Utilities2
    {
    private:
        const size_t matrix_type = 3003;
    public:
        Utilities2();
        ~Utilities2();

        void check_id(std::vector<std::int64_t> &id_list,
                      std::vector<std::int64_t> &checked_id,
                      std::vector<std::int64_t> &missing_id);
        void check_id2(std::vector<std::int64_t> &id_list,
                      std::vector<std::int64_t> &checked_id,
                      std::vector<std::int64_t> &is_in_id);
        bool is_value_in_vect(std::vector<std::int64_t> &where_tocheck,
                              std::vector<std::int64_t> &what_tocheck);
        void get_missing_in_vect(std::vector<std::int64_t> &where_tocheck,
                                 std::vector<std::int64_t> &what_tocheck,
                                 std::vector<std::int64_t> &missing);
        void get_what_in_vect(std::vector<std::int64_t> &where_tocheck,
                                 std::vector<std::int64_t> &what_tocheck,
                                 std::vector<std::int64_t> &what_in);
        int find_invect(std::vector<std::int64_t> &where,
                        std::int64_t what);
        bool is_invect(std::vector<std::int64_t> &where,
                       std::int64_t what); // was find_invect
        bool is_unique(std::vector<std::int64_t> &x);
        void get_RecodedIdMap(std::map<std::int64_t, std::int64_t> &id_map,
                              std::vector<std::int64_t> &idVect);
        void find_RecodedIdMap(std::map<std::int64_t, std::int64_t> &id_map,
                               std::vector<std::int64_t> &whereIidList,
                               std::vector<std::int64_t> &whatIdList);
        void vect_to_binary(std::vector<std::int64_t> &vect, const std::string &fname);
        void vect_from_binary(std::vector<std::int64_t> &vect, const std::string &fname);

        template <typename T>
        void dense_to_sparse(evolm::matrix<T> &from, evolm::smatrix<T> &to);
        template <typename T>
        void sparse_to_dense(evolm::smatrix<T> &from, evolm::matrix<T> &to);

        size_t fget_matrix_kind(const std::string &fname);

        template <typename T>
        void fwrite_matrix(const std::string &fname, std::vector<T> &vals, std::vector<size_t> &keys, std::vector<std::int64_t> &ids, size_t n_cols = 0);
        template <typename T>
        void fwrite_matrix(const std::string &fname, std::vector<T> &vals, std::vector<size_t> &keys, std::vector<std::string> &ids, size_t n_cols = 0);
        template <typename T>
        void fwrite_matrix(const std::string &fname, evolm::matrix<T> &matr, std::vector<std::int64_t> &ids);
        template <typename T>
        void fwrite_matrix(const std::string &fname, evolm::matrix<T> &matr, std::vector<std::string> &ids);

        template <typename T>
        void fread_matrix(const std::string &fname, std::vector<T> &vals, std::vector<size_t> &keys, std::vector<std::int64_t> &ids);
        template <typename T>
        void fread_matrix(const std::string &fname, std::vector<T> &vals, std::vector<size_t> &keys, std::vector<std::string> &ids);
        template <typename T>
        void fread_matrix(const std::string &fname, evolm::matrix<T> &matr, std::vector<std::int64_t> &ids);
        template <typename T>
        void fread_matrix(const std::string &fname, evolm::matrix<T> &matr, std::vector<std::string> &ids);
        void fread_matrix_info(const std::string &fname, size_t &info);

        template <typename T>
        void solve_ls(evolm::matrix<T> &L, T *b, T *x);

        template <typename T>
        bool is_float(T val);

    };

} // end of namespace evoped

#endif // Utilities2_hpp__