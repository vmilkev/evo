#ifndef Population_hpp__
#define Population_hpp__

#include "Animal.hpp"
#include "Iointerface.hpp"
#include <vector>

namespace evogen
{
    class Population
    {
    public:
        Population();
        ~Population();

        void set_population(size_t nindividuals,
                            const std::string &genotype_structure,
                            float ref_allele_probability);         /* set a 'wild' Population with simulated (or empty?) genome */

        void set_population(const std::string &haplotypes_fname,
                            const std::string &genotype_structure,
                            bool with_pedigree);                   /* sets Population with predefined genotypes which provided from file */
        size_t get_size();

        size_t get_ploidy();
        //size_t get_nmarkers();

        std::vector<short> get_genome_at(size_t which_genome, size_t locus);

        std::vector<std::vector<unsigned long>> get_genome_table();

        void show_animals( size_t max_animals, size_t max_snps ); // for testing

    private:
        std::vector<evogen::Animal> individuals;

    protected:
    };
} // end of namespace evo

#endif // Population_hpp__