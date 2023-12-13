#ifndef population_hpp__
#define population_hpp__

#include "Animal.hpp"
#include "Iointerface.hpp"
#include <vector>
#include <algorithm>
#include <memory>
#include <numeric>

namespace evogen
{
    class Population
    {
    public:
        Population();
        ~Population();

        void set_population(size_t nindividuals,
                            const std::string &genotype_structure,
                            float ref_allele_probability,
                            int n_ploidy); /* set a 'wild' Population with simulated (or empty?) genome */

        void set_population(const std::string &haplotypes_fname,
                            const std::string &genotype_structure,
                            bool with_pedigree); /* sets Population with predefined genotypes which provided from file */
        
        size_t size();     // number of active individuals based on the active_individuals list
        size_t capacity(); // number of all (active + disabled) individuals, based on the individuals list
        void clear();
        void reshape();    // reduce capaciity to fit size (hence, memory)

        void aging(int delta_t);

        void id_at(size_t at, unsigned long id);
        unsigned long id_at(size_t at);

        void sire_at(size_t at, unsigned long id);
        unsigned long sire_at(size_t at);

        void dame_at(size_t at, unsigned long id);
        unsigned long dame_at(size_t at);

        void age_at(size_t at, int age);
        int age_at(size_t at);

        void alive_at(size_t at, bool alive);
        bool alive_at(size_t at);
        
        void isgenotyped_at(size_t at, bool genotyped);
        bool isgenotyped_at(size_t at);

        void sex_at(size_t at, int sex);
        short sex_at(size_t at);

        void phenotype_at(size_t at, std::vector<double> &phen); // do not make this public ?
        std::vector<double> phenotype_at(size_t at);

        void breedingvalue_at(size_t at, std::vector<double> &bv); // do not make this public ?
        std::vector<double> breedingvalue_at(size_t at);

        // ------------- Not publik in Python ------------------------
        std::vector<short> get_genome_at(size_t which_genome, size_t locus);
        void get_all_genotypes(const std::string &file_out);
        void get_all_genotypes(std::vector<std::vector<short>> &vect_out);

        // ------------- Required in Trait class ---------------------
        std::vector<std::vector<unsigned long>> get_genome_table();
        std::vector<std::vector<unsigned long>> get_genome_structure();
        size_t get_ploidy();
        // ------------------------------------------------------------

        friend class Group;

        // ------------------ TESTING --------------------------------
        void show_animals(size_t max_animals, size_t max_snps);

        void show_pop()
        {
            for (size_t i = 0; i < individuals.size(); i++)
            {
                std::cout<<"from individuals: i = "<<i<<", id = "<<individuals[i].get_id()<<", active = "<<individuals[i].get_active()<<"\n";
            }
            std::cout<<"\n";
            for (size_t i = 0; i < active_individuals.size(); i++)
            {
                std::cout<<"from active_individuals: i = "<<i<<", value = "<<active_individuals[i]<<"\n";
            }
        }

    private:
        std::vector<Animal> individuals;
        std::vector<size_t> active_individuals;

        void remove_at(size_t which_list_position);
        void add(Animal &in_a); // adds existent entiety from another population

    protected:
    };

    inline bool operator==(const Population &rhs, const Population &lhs)
    {
        return std::addressof(rhs) == std::addressof(lhs);
    }

    inline bool operator!=(const Population &rhs, const Population &lhs)
    {
        return std::addressof(rhs) != std::addressof(lhs);
    }

} // end of namespace evo

#endif // population_hpp__