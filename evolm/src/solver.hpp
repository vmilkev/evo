#ifndef solver_hpp__
#define solver_hpp__

#include <vector>
#include <map>
#include "model.hpp"

namespace evolm
{

    class Solver
    {

    public:
        virtual void solve() = 0;

        void append_model(const Model &m);

        void remove_model();

        Solver()
        {
            n_trait = 0;
            z_on_memory = false;
            is_model_added = false;
            var_onmem = false;
            cor_onmem = false;
        };

        virtual ~Solver();

#ifdef UTEST
        void print();
        void test_vect_z_uni(const size_t &which_trait, matrix<float> &out);
#endif

    private:
        void process_model();
        void set_complete_z(const matrix<int> &eff_ids); // build full incedense matrix, z_dat
        void construct_union_z(const matrix<int> &eff_ids);
        void construct_union_z(const matrix<int> &eff_ids, bool on_memory);
        void set_y(const int obs_id);
        void set_g();
        void set_r();
        void z_row_to_uni_col(const size_t &which_trait,
                              const size_t &in_z_row,
                              size_t *out_uni_col,
                              size_t &out_uni_matr);

    protected:
        Model model;                             // instance of the model to be solved
        std::vector<size_t> n_obs;               // number of observations for each trait
        size_t n_trait;                          // number of traits
        std::vector<std::vector<size_t>> n_lev;  // num of effects' levels for each trait
        std::vector<matrix<float>> z_dat;        // combined incidence matrix for each trait.
        std::vector<matrix<float>> y;            // Observations for each trait.
        std::vector<std::vector<Effects>> z_uni; // Combines a consecutive sets of incidense matrices.

        std::vector<std::vector<bool>> s;             // boolean container indicating which data on which trait is missing (false), which is not (true)
        std::vector<size_t> R_hash;                   // the hash keys corresponding each observation pattern, the size of n_obs; the size and values are the same for all traits
        std::map< size_t, std::vector<float> > r_map; // hash values are keys, and specific (according to the observation pattern) covar matrices are values of the map
                                                      // the usage: r_map[ R_hash[i] ][j], where i is observation, j is indexing j = r*(r+1)/2 +c pointing to an element of R(-1)
        
        //matrix<float> r;
        //evolm::matrix<float> T;  // lower tringular matrix, the result of T = chol(r)
        //evolm::matrix<float> iT; // inverse of T
        //matrix<float> iT_tr; // transposed iT

        bool z_on_memory;
        bool is_model_added;
        bool var_onmem;
        bool cor_onmem;

        std::map<size_t, size_t> adj_effects_order; // map between the submitted effects ids and their recoded (consecutive) order

        matrix<float> get_vect_z_uni(const size_t &which_trait, const size_t &which_row);
        void get_vect_z_uni2(const size_t &which_trait, const size_t &which_row, std::vector<std::vector<float>> &vect);
        void get_vect_z_uni2(const size_t &which_trait, const size_t &which_row, float **vect);
        size_t corr_size();
        matrix<int> get_corr_effects(size_t which_correlation);
        matrix<float> get_variance(size_t which_correlation, size_t row, size_t col);
        matrix<float> get_correlation(size_t which_trait, size_t which_row, size_t col_1, size_t col_2);
        bool identity_correlation(size_t which_corr);

        void memload_effects();
        void diskload_effects();
        void memload_var();
        void diskload_var();
        void memload_cor();
        void diskload_cor();
    };

} // end of namespace evolm

#endif // solver_hpp__
