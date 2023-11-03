/*
    Abstract class
*/

#ifndef genome_hpp__
#define genome_hpp__

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "Utilites.hpp"

namespace evogen
{
    class Genome
    {
    public:
        Genome();
        ~Genome();

        void set_genome(std::vector<std::vector<bool>> &snp,
                        std::vector<std::vector<unsigned long>> &gstructure);   /* uses prepared SNP variants either from file or from reproduction gamete */
        
        void set_genome(std::vector<std::vector<unsigned long>> &gstructure,
                        float ref_allele_probability,
                        short nploidy);                                         /* simulate SNP variants */
        
        void get_genome( std::vector<std::vector<bool>> &haplotypes,
                         std::vector<std::vector<unsigned long>> &gstructure ); /* returns genome information: variants and structure */
        
        short get_ploidy();

        void get_reproduction_gamete( size_t cross_per_chr, size_t mut_per_genome, std::vector<bool> &out_gamete, short &out_sex_chr_id ); /* NOTE! we allow possibility for polyploid gametes */
        
        // For tests
        void show_genome( int max_variants ); // this is for debugging
        void show_recombination();

    private:
        std::vector<std::vector<bool>> markers;   /* SNP (markers) variants for entire genome, the structure descriped in 'structure' */
        std::vector<std::vector<unsigned long>> structure; /* length of each chromosome (bp) and markers density (constant, bp) */
        std::vector<std::vector<unsigned long>> snp_table; /* for each chromosome indicates the first and the last snp in consecutive order */

        short asign_snp_variant( float ref_allele_probability );

        // part of reproduction process
        void mutation( std::vector<bool> &in_gamete, size_t events_genome );
        void recombination( std::vector<bool> &out_gamete, short &out_sex_chr_id, size_t n_crosses );
        std::vector<bool> crossover( size_t in_chr, size_t out_strand, size_t n_crossovers );
        std::vector<bool> get_chromatide( size_t which_chr, size_t which_strand );
        void def_snp_table();

    protected:
    };
} // end of namespace evo

#endif // genome_hpp__