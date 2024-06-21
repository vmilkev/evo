#ifndef sparse_solver_hpp__
#define sparse_solver_hpp__

#include "model_sparse.hpp"

namespace evolm
{
    class sparse_solver
    {
    public:
        
        void append_model(const model_sparse &m);
        void remove_model();
        void set_memory_limit(double limit);
        
        sparse_solver();
        ~sparse_solver();

        //virtual void solve() = 0;
        //virtual ~sparse_solver();

        size_t get_num_of_mem_blocks();
        void load_model_matrix(size_t mem_blok);
        void get_mem_block_range(size_t mem_blok, size_t &first, size_t &second);
        void unload_model_matrix(size_t mem_blok);

#ifdef UTEST
        void test_vect_z_uni(const size_t &which_trait, matrix<float> &out);
        matrix<float> test_z_dot_y(size_t vect_size, size_t i_trait, size_t j_trait, size_t r_index);
        matrix<float> test_rhs();
        size_t test_num_all_levels();
        std::vector<size_t> test_ordered_levels();
        std::vector<std::vector<size_t>> test_cov_offsets();
        std::vector<std::vector<float>> test_A(); // model matrix
#endif

    private:

        void process_model();
        void construct_union_z(const matrix<int> &eff_ids, bool on_memory);
        void set_y(const int obs_id);
        void set_g();
        void set_r();
        void z_row_to_uni_col(const size_t &which_trait,
                              const size_t &in_z_row,
                              size_t *out_uni_col,
                              size_t &out_uni_matr);
        
        void get_load_per_memory_block(std::vector<std::vector<size_t>> &loads);
        std::string create_fname();
        void fwrite(const std::string &fname, size_t first_row, size_t last_row);
        void fread(const std::string &fname, size_t first_row, size_t last_row);

        double available_memory = 100; // default available memory in GB
        const double gb_constant = 1073741824; // num of bytes in 1 GB

        std::fstream fA;
        std::vector<std::string> bin_fnames;
         std::vector<std::vector<size_t>> blocks_ranges;

    protected:

        model_sparse model;                             // instance of the model to be solved
        std::vector<size_t> n_obs;               // number of observations for each trait
        size_t n_trait = 0;                          // number of traits
        std::vector<std::vector<size_t>> n_lev;  // num of effects' levels for each trait
        std::vector<matrix<float>> y;            // Observations for each trait.
        std::vector<std::vector<effects_storage>> z_uni; // Combines a consecutive sets of incidense matrices.

        std::vector<std::vector<bool>> s;           // boolean container indicating which data on which trait is missing (false), which is not (true)
        std::vector<size_t> R_hash;                 // the hash keys corresponding each observation pattern, the size of n_obs; the size and values are the same for all traits
        std::map<size_t, std::vector<float>> r_map; // hash values are keys, and specific (according to the observation pattern) covar matrices are values of the map
                                                    // the usage: r_map[ R_hash[i] ][j], where i is observation, j is indexing j = r*(r+1)/2 +c pointing to an element of R(-1)

        //---------------------------------------------
        matrix<float> rhs;

        size_t get_levels(size_t which_trait);
        size_t get_all_levels();
        size_t get_all_levels(size_t before_trait);
        size_t num_all_levels();
        std::vector<size_t> get_ordered_levels();
        std::vector<std::vector<size_t>> get_cov_offsets(const std::vector<size_t> &ordered_levels);
        void construct_rhs();
        void z_dot_y(matrix<float> &out_vect, size_t vect_size, size_t i_trait, size_t j_trait, size_t r_index);

        void z_dot_z(matrix<float> &out_vect, size_t row, size_t vect_size, size_t i_matr, size_t j_matr, size_t r_index);
        void z_dot_z2(std::vector<float> &out_values, std::vector<size_t> &out_keys, size_t row, size_t vect_size, size_t i_matr, size_t j_matr, size_t r_index);
        //---------------------------------------------
        std::vector<compact_storage<float>> model_matrix;

        matrix<float> amatr; // temporal !

        void make_model_matrix(std::vector<std::vector<size_t>> &cov_offsets, size_t num_levels, std::vector<size_t> &ordered_levels);
        void set_model_matrix();

        void get_row_cmatr2(
                            compact_storage<float> &model_matrix_row,
                            size_t rhs_size,
                            size_t i_trate,
                            size_t i_eff,
                            std::vector<std::vector<size_t>> &cov_offsets,
                            size_t num_levels,
                            std::vector<size_t> &ordered_levels,
                            size_t i_row);

        matrix<float> get_row_cmatr(
                            size_t rhs_size,
                            size_t i_trate,
                            size_t i_eff,
                            std::vector<std::vector<size_t>> &cov_offsets,
                            size_t num_levels,
                            std::vector<size_t> &ordered_levels,
                            size_t i_row);

        //---------------------------------------------

        bool z_on_memory = false;
        bool var_onmem = false;
        bool cor_onmem = false;

        std::map<size_t, size_t> adj_effects_order; // map between the submitted effects ids and their recoded (consecutive) order
        
        void get_vect_z_uni2(const size_t &which_trait, const size_t &which_row, std::vector<float> &values, std::vector<size_t> &keys);
        void get_vect_z_uni2(const size_t &which_trait, const size_t &which_row, float **vect);
                
        size_t corr_size();
        matrix<int> get_corr_effects(size_t which_correlation);
        float get_variance(size_t which_correlation, size_t row, size_t col);
        matrix<float> get_correlation(size_t which_trait, size_t which_row, size_t col_1, size_t col_2);
        void add_correlation(std::vector<float> &vect_to_add, size_t vect_first_index, float variance, size_t which_trait, size_t which_row, size_t col_1, size_t col_2);
        bool identity_correlation(size_t which_corr);

        void memload_effects();
        void diskload_effects();
        void memload_var();
        void diskload_var();
        void memload_cor();
        void diskload_cor();
        void memload_cor_effects();
        void diskload_cor_effects();
    };

} // end of namespace evolm

#endif // solver_hpp__
