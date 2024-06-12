/**
 * @file Gmat.hpp
 * @author Viktor Milkevych
 * @brief 
 * @version 0.1
 * @date 2024-05-14
 * 
 * @copyright Copyright (c) 2024
 * 
 */
#ifndef Gmat_hpp__
#define Gmat_hpp__

#include <fstream>
#include <unordered_map>

#include "dense_matrix.hpp"
#include "Utilities2.hpp"
#include <plinkio/plinkio.h>

#define workload 100000

namespace evoped
{
    template <typename T>
    class Gmat
    {
    public:
        Gmat();
        ~Gmat();

        void read_matrix(const std::string &gmat_file);
        
        void make_matrix(const std::string &fname);
        void make_matrix(const std::string &fname, const std::string &fname_ids);
        
        void scale_genotypes(const std::string &fname);
        void scale_genotypes(const std::string &fname, const std::string &fname_ids);        
        
        void scale_diag(T scale_coef);        
        void scale_matrix(evolm::matrix<T>& scale_matr, T scaling_weight);
        
        void invert_matrix();
        void invert_matrix(bool full_store);        
        void invert_matrix(std::vector<std::int64_t>& core_id); // sparse inverse (APY)        
        
        void get_matrix(evolm::matrix<T> &arr, std::vector<std::int64_t> &ids);
        void get_matrix(evolm::matrix<T> &arr);
        void get_ids(std::vector<std::int64_t> &ids);

        void save_matrix(const std::string &arr, const std::string &ids);
        void save_matrix(const std::string &arr);
        void save_ids(const std::string &ids);
        
        void clear();

        // Python specific interfaces
        void scale_matrix(std::vector<T>& scale_matr, T scaling_weight);
        void get_matrix(std::vector<T>& arr, std::vector<std::int64_t>& ids);

#ifdef UTEST
        void get_alpha_beta(T &alpha,
                            T &beta,
                            T &a_diag,
                            T &a_ofd,
                            T &g_diag,
                            T &g_ofd);
#endif

    private:
        evolm::matrix<T> G;          // G or inv. of G matrix, permanent object container        
        evolm::matrix<T> Z;          // Z (snp) matrix, temporal
        std::vector<std::int64_t> gmatID; // container for the list of G matrix IDs, initiated while reading pre-built G-matrix from file
        T freq;
        std::map<std::int64_t, std::string> snp_map; // initial markers data, temporal
        std::map<size_t, std::int64_t> anim_id_map;  // key: consecutive index, value: animal ID
        std::unordered_map<size_t, char *> samples_id_map; // key: assigned sample id, value: original id (as in .fam file)

        void read_snp(const std::string &snp_file);
        void read_snp(const std::string &snp_file, const std::string& ids_file);
        void parse_string(std::string &snp_str, std::vector<int> &markers);
        void make_zmatrix();
        void make_zmatrix( evolm::matrix<int> &M );
        void make_matrix();
        void read_matrix(const std::string &gmat_file, std::vector<std::int64_t> &g_row, std::vector<std::int64_t> &g_col, std::vector<T> &g_val);
        
        bool is_plink_file(const std::string &fname);
        void get_m_matrix(const std::string &fname, evolm::matrix<int> &M);


#ifdef UTEST
        T scaling_a;
        T scaling_b;
        T scaling_a_diag;
        T scaling_a_ofd;
        T scaling_g_diag;
        T scaling_g_ofd;
#endif

    };

} // end of namespace evoped

#endif // Gmat_hpp__