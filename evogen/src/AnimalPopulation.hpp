#ifndef AnimalPopulation_hpp__
#define AnimalPopulation_hpp__

#include "Animal.hpp"
#include "Iointerface.hpp"
#include <vector>

namespace evogen
{
    class AnimalPopulation
    {
    public:
        AnimalPopulation();
        ~AnimalPopulation();

        void set_population(size_t nindividuals,
                            const std::string &genotype_structure,
                            float ref_allele_probability);         /* set a 'wild' AnimalPopulation with simulated (or empty?) genome */

        void set_population(const std::string &haplotypes_fname,
                            const std::string &genotype_structure,
                            bool with_pedigree);                   /* sets AnimalPopulation with predefined genotypes which provided from file */
        size_t get_size();

        void show_animals( int max_animals, int max_snps ); // for testing

    private:
        std::vector<evogen::Animal> individuals;

    protected:
    };
} // end of namespace evo

#endif // AnimalPopulation_hpp__