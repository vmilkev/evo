#ifndef model_sparse_hpp__
#define model_sparse_hpp__

#include <string>
#include <vector>
#include <typeinfo>

#include "effects_storage.hpp"
#include "iointerface.hpp"

namespace evolm
{
        class model_sparse
        {
        public:
                friend class sparse_solver;
                // ----------------------------------------------
                model_sparse();
                // ----------------------------------------------
                void set_missing(float val);
                // ----------------------------------------------
                int append_residual(const std::vector<float> &arr, size_t lda);
                int append_residual(const std::string &fname);
                // ----------------------------------------------
                int append_observation(const std::vector<float> &arr, size_t lda);
                int append_observation(const std::vector<float> &arr, const std::vector<bool> &miss_arr, size_t lda);
                int append_observation(const std::string &fname);
                // ----------------------------------------------
                template <typename T>
                void append_effect(compact_storage<T> &eff);
                template <typename T>
                void append_effect(std::vector<T> &values, size_t n_rows, size_t n_cols);
                template <typename T>
                void append_effect(std::vector<T> &values, std::vector<size_t> &rows, std::vector<size_t> &cols, size_t n_rows, size_t n_cols);
                void append_effect(const std::string &fname);
                // ----------------------------------------------
                int append_traitstruct(int obs_id, const std::vector<int> &eff_id);
                // ----------------------------------------------
                void append_corrstruct(const std::vector<float> &var, size_t lda1, std::vector<float> &corr, size_t lda2, const std::vector<int> &which_effects);
                void append_corrstruct(const std::vector<float> &var, size_t lda1, const std::string &fname_corr, const std::vector<int> &which_effects);
                void append_corrstruct(const std::string &fname_var, const std::string &fname_corr, const std::vector<int> &which_effects);
                
                void append_corrstruct(const std::vector<float> var, size_t lda1, std::string &identity, size_t lda2, const std::vector<int> which_effects);
                void append_corrstruct(const std::string &fname_var, std::string &identity, size_t lda2, const std::vector<int> which_effects);
                // ----------------------------------------------
                void clear_residuals();
                void clear_observations();
                void clear_effects();
                void clear_corrstruct();
                void clear_traitstruct();
                void clear(); // clears all
                // ----------------------------------------------
                void set_sparsity_threshold(double threshold);
                size_t get_size_of_data();
#ifdef UTEST
                size_t size_of(const std::string type);
                std::vector<std::vector<size_t>> shape_of(const std::string type);
                void print();
                compact_storage<float> test_effects(size_t which_effect);
                std::vector<float> test_observations(size_t which_observations);
                std::vector<float> test_residual(size_t which_residual);
                std::vector<float> test_variance(size_t which_variance);
                std::vector<float> test_correlation(size_t which_correlation);
#endif
                // ----------------------------------------------
        private:
                double sparsity_threshold = 0.5;
                size_t size_of_data = 0.0;

                // Data
                std::vector<matrix<float>> residuals;             // assumed to full matrix
                std::vector<matrix<float>> observations;          // assumed to be (nx1) vectors
                std::vector<std::vector<bool>> miss_observations; // assumed to be (nx1) vectors
                std::vector<effects_storage> all_effects;

                // Correlation structure
                std::vector<matrix<int>> correlated_effects; // assumed to be (nx1) vectors. In matlab: which_correlated_effects
                std::vector<matrix<float>> variances;        // assumed to be full matrix. In matlab: G
                std::vector<compact_storage<float>> correlations;
                std::vector<bool> identity_correlations;     // indicates which provided correlations should be treated as identity matrices
                std::vector<size_t> identity_dimension;

                // Trait structure
                std::vector<int> observation_trait;     // assumed to be (nx1) vectors
                std::vector<matrix<int>> effects_trait; // assumed to be (nx1) vectors

                float missing_constant = -999.0; // variable indicating which number recogniised as missing in observations
        };

} // end of namespace evolm

#endif // model_sparse_hpp__
