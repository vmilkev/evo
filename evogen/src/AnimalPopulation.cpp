#include "AnimalPopulation.hpp"

namespace evogen
{
    //===============================================================================================================

    AnimalPopulation::AnimalPopulation()
    {
        //
    }

    //===============================================================================================================

    AnimalPopulation::~AnimalPopulation()
    {
        //
    }

    //===============================================================================================================

    void AnimalPopulation::set_population(size_t nindividuals, const std::string &genotype_structure, float ref_allele_probability)
    {
        try
        {
            if (ref_allele_probability > 1.0)
            {
                std::string err("The reference allele probability should be less then 1.0!");
                throw err;
            }

            std::vector<std::vector<unsigned long>> genome_structure; /* length of each chromosome (bp) and markers density (constant, bp) */

            evogen::IOInterface data;
            data.set_fname(genotype_structure);
            data.fgetdata(genome_structure);

            for (size_t i = 0; i < nindividuals; i++)
            {
                Animal a;
                a.genome.set_genome(genome_structure, ref_allele_probability, 2);
                individuals.push_back(a);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in AnimalPopulation::set_population(size_t, const std::string &, float)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &err)
        {
            std::cerr << "Exception in AnimalPopulation::set_population(size_t, const std::string &, float)." << '\n';
            std::cerr << err << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in AnimalPopulation::set_population(size_t, const std::string &, float)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    size_t AnimalPopulation::get_size()
    {
        size_t out = 0;

        try
        {
            out = individuals.size();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in AnimalPopulation::get_size()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in AnimalPopulation::get_size()." << '\n';
            throw;
        }

        return out;
    }

    //===============================================================================================================

    void AnimalPopulation::set_population(const std::string &haplotypes_fname, const std::string &genotype_structure, bool with_pedigree)
    {
        /*
        Create population using haplotypes file. This variant reads unnamed haplotypes, hence,
        assign random IDs to individuals with unknown parents (unrelated population).
        */
        try
        {
            /* because it is the animal population we have diploid genotype,
            hence, read piars of haplotype for each individual */

            // 1. Read the data files
            evogen::IOInterface data;

            data.set_fname(haplotypes_fname);

            std::vector<std::vector<unsigned long>> pedigree;
            std::vector<std::vector<bool>> haplotypes;                /* haplotypes of all individuals, rows = 2 * no.individuals */
            std::vector<std::vector<unsigned long>> genome_structure; /* length of each chromosome (bp) and markers density (constant, bp) */

            if (with_pedigree)
                data.fgetdata(pedigree, haplotypes);
            else
                data.fgetdata(haplotypes);

            data.set_fname(genotype_structure);
            data.fgetdata(genome_structure);

            // 2. Creating animals with a specific genotype

            size_t n_haplotypes = haplotypes.size();
            size_t n_markers = haplotypes[0].size();

            std::cout << "haplotypes, markers: " << n_haplotypes << ". " << n_markers << "\n";

            if (n_haplotypes % 2 != 0)
            {
                std::string err("The number of haplotypes in the data file is not even! Exit.");
                throw err;
            }

            for (int i = 0; i < n_haplotypes;)
            {
                std::vector<std::vector<bool>> markers(2, std::vector<bool>(n_markers));

                markers[0] = haplotypes[i];
                markers[1] = haplotypes[i + 1];

                Animal a;
                a.genome.set_genome(markers, genome_structure);

                if (with_pedigree)
                {
                    unsigned long id = pedigree[i][0];
                    short sex = pedigree[i][1];
                    unsigned long sire = pedigree[i][2];
                    unsigned long dame = pedigree[i + 1][2];

                    a.set_id(id);
                    a.set_sex(sex);
                    a.set_sire(sire);
                    a.set_dame(dame);
                }

                individuals.push_back(a);

                i = i + 2;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in AnimalPopulation::set_population(const std::string &, const std::string &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in AnimalPopulation::set_population(const std::string &, const std::string &)." << '\n';
            std::cerr << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in AnimalPopulation::set_population(const std::string &, const std::string &)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void AnimalPopulation::show_animals( int max_animals, int max_snps )
    {
        try
        {
            std::cout << "Animals showup:"
                      << "\n";
            size_t icount = 0;
            for (auto &e : individuals)
            {
                std::cout << "genome:"
                          << "\n";
                e.genome.show_genome(max_snps);
                e.genome.show_recombination();

                std::cout << "id: " << e.get_id() << "\n";
                std::cout << "sire: " << e.get_sire() << "\n";
                std::cout << "dame: " << e.get_dame() << "\n";
                std::cout << "sex: " << e.get_sex() << "\n";
                std::cout << "age: " << e.get_age() << "\n";
                std::cout << "alive: " << e.get_alive() << "\n";

                if (icount >= max_animals)
                    break;

                icount++;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in AnimalPopulation::show_animals(int, int)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in AnimalPopulation::show_animals(int, int)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    //===============================================================================================================

    //===============================================================================================================

    //===============================================================================================================

    //===============================================================================================================

    //===============================================================================================================

    //===============================================================================================================

    //===============================================================================================================

    //===============================================================================================================

    //===============================================================================================================

}