#ifndef trait_hpp__
#define trait_hpp__

#include "Genome.hpp"
#include "Population.hpp"
#include "cs_matrix.hpp"
#include "Utilites.hpp"
#include <iostream>
#include <algorithm>
#include <math.h>

namespace evogen
{
    class Trait
    {
    public:
        Trait();
        ~Trait();

        void set_trait(Population &pop,                          // reference to the base population
                       std::vector<double> &trmean,              // expected values of trait in base population, size=(n_trait,1)
                       std::vector<double> &qtl_prop_chrom,      // proportion of snps used as qtls for each chromosome, size=(n_chr,1)
                       std::vector<std::vector<double>> &corr_g, // correlations of genetic effects, size=(n_trait, n_trait)
                       std::vector<double> &varr_g,              // variances of genetic effects, size=(n_trait, 1)
                       std::vector<std::vector<double>> &corr_e, // correlations of environmental (residual) effects, size=(n_trait, n_trait)
                       std::vector<double> &varr_e,              // variances of environmental effects, size=(n_trait, 1)
                       std::vector<double> &envr,                // Control of an environment (accounting its changes), size=(n_trait,1)
                       size_t dist_model,                        // type of distribution used to sample dominance
                       std::vector<double> &dist_par);           // dominance distribution parameters, size=(2,1)

        void reset_trait(Population &pop,
                         std::vector<double> &trmean,
                         std::vector<double> &qtl_prop_chrom,
                         std::vector<std::vector<double>> &corr_g,
                         std::vector<double> &varr_g,
                         std::vector<std::vector<double>> &corr_e,
                         std::vector<double> &varr_e,
                         std::vector<double> &envr,
                         size_t dist_model,
                         std::vector<double> &dist_par);

        void get_observations(Population &in_pop,                             // reference to the base population for which making observations
                              std::vector<double> &out_trvalues,              // container with observed trait values for each individual in in_pop
                              std::vector<std::vector<bool>> &out_genotypes); // container with observed genotypes for each individual in in_pop

        void get_observations(Population &in_pop,
                              const std::string &out_trvalues,   // output file name where observed trait values for each individual in in_pop will be writen
                              const std::string &out_genotypes); // output file name where observed genotypes for each individual in in_pop will be writen

        void get_observations(Population &in_pop,
                              std::vector<double> &out_trvalues);

        void get_observations(Population &in_pop,
                              std::vector<double> &env,
                              const std::string &out_t);

    private:
        evolm::matrix<double> a;    // additive effects, size=(qtls,n_trait)
        evolm::matrix<double> e;    // environmental effects, size=(qtls,n_trait)
        evolm::matrix<double> k;    // dominance effects, size=(qtls,1)
        evolm::matrix<double> qtls; // trait responsible qtls locations, size=(n_qtls,1)
        evolm::matrix<double> t_mean; // correction mean for the base population
        evolm::matrix<double> ta;   // genotype-determined trait
        evolm::matrix<double> te;   // environment-determined trait

        void clear();
        void sample_dom(size_t which_dist, std::vector<double> &dist_param);
        void sample_genes(std::vector<size_t> &n_qtl,
                          std::vector<std::vector<unsigned long>> &stable,
                          std::vector<double> &qtl_prop);
        void sample_effects(size_t n_trate);
        void calculate_trait(Population &in_pop, std::vector<double> &envr, size_t n_trate);

        double ploidy_effect(size_t n_ploidy, double degree);

        std::vector<int> get_locus_state(Population &in_pop,
                                       size_t individ,
                                       size_t qtl_position);

        evolm::matrix<double> var_diag( evolm::matrix<double> &arr );
        evolm::matrix<double> get_scaler( std::vector<double> &in_var, evolm::matrix<double> &in_arr );
        void realloc_traits(Population &in_pop,size_t n_trait);
        void calculate_correction_mean(std::vector<double> &in_mean);

    protected:
    };

} // end of namespace evogen

#endif // trait_hpp__