#ifndef genome_hpp__
#define genome_hpp__

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
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

        void get_reproduction_gamete(std::vector<std::vector<bool>> &out_gamete,
                                     size_t cross_per_chr,
                                     float mut_freq);                   /* NOTE! we allow possibility for polyploid gametes */

        std::vector<short> get_genome();
        
        // For tests
        void show_genome( int max_variants );
        void show_recombination();

        std::vector<short> get_genome_at(size_t locus);

        std::vector<std::vector<unsigned long>> get_snp_table();
        std::vector<std::vector<unsigned long>> get_genome_structure();

        void clear();

        void operator=(const Genome& g)
        {
            markers = g.markers;
            structure = g.structure;
            snp_table = g.snp_table;
            hotspots = g.hotspots;
            std_hotspots = g.std_hotspots;
        }
        
    private:
    
        std::vector<std::vector<bool>> markers;   /* SNP (markers) variants for entire genome, the structure descriped in 'structure' */
        std::vector<std::vector<unsigned long>> structure; /* length of each chromosome (bp) and markers density (constant, bp) */
        std::vector<std::vector<unsigned long>> snp_table; /* for each chromosome indicates the first and the last snp in consecutive order */
        std::vector<std::vector<size_t>> hotspots; // locations of recombination hotspots for each chromosome
        std::vector<float> std_hotspots; // allowed standart deviation of recombination events around specific hotspot

        const float freq_hotspot = 0.2; // frequency of recombination hotspots; n_hotspot = (freq_hotspot/mark_density) * n_snp_chr, number of equally distributed hotspots along each chromosome
        const float step_devider = 5.0;

        short asign_snp_variant( float ref_allele_probability );

        // part of reproduction process
        void mutation( std::vector<std::vector<bool>> &in_gamete, float mut_freq );
        void recombination( std::vector<std::vector<bool>> &out_gamete, size_t n_crosses );
        std::vector<bool> crossover( size_t in_chr, size_t out_strand, size_t n_crossovers, size_t which_strands_set );
        std::vector<bool> get_chromatide( size_t which_chr, size_t which_strand, size_t which_strands_set );
        
        void crossover_model_gamma(size_t in_chr, size_t out_strand, size_t n_crossovers, size_t which_strands_set, std::vector<bool> &out_chromatide);
        void crossover_model_uniform(size_t in_chr, size_t out_strand, size_t n_crossovers, size_t which_strands_set, std::vector<bool> &out_chromatide);
        
        void def_snp_table();
        void def_hotspot_table();

    protected:
    };
} // end of namespace evo

#endif // genome_hpp__