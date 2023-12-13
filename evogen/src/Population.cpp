#include "Population.hpp"

namespace evogen
{
    //===============================================================================================================

    Population::Population()
    {
        //
    }

    //===============================================================================================================

    Population::~Population()
    {
        //
    }

    //===============================================================================================================

    void Population::aging(int delta_t)
    {
        try
        {
            for (size_t i = 0; i < size(); i++)
            {
                int new_age = age_at(i) + delta_t;
                age_at(i,new_age);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::aging(int)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &err)
        {
            std::cerr << "Exception in Population::aging(int)" << '\n';
            std::cerr << err << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::aging(int)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::set_population(size_t nindividuals, const std::string &genotype_structure, float ref_allele_probability, int n_ploidy)
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
                a.genome.set_genome(genome_structure, ref_allele_probability, n_ploidy);
                individuals.push_back(a);
                active_individuals.push_back(i);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::set_population(size_t, const std::string &, float, int)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &err)
        {
            std::cerr << "Exception in Population::set_population(size_t, const std::string &, float, int)" << '\n';
            std::cerr << err << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::set_population(size_t, const std::string &, float, int)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    size_t Population::size()
    {
        size_t out = 0;

        try
        {
            out = active_individuals.size();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::size()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &err)
        {
            std::cerr << "Exception in Population::size()" << '\n';
            std::cerr << err << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::size()." << '\n';
            throw;
        }

        return out;
    }

    //===============================================================================================================

    size_t Population::capacity()
    {
        size_t out = 0;

        try
        {
            out = individuals.size();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::capacity()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &err)
        {
            std::cerr << "Exception in Population::capacity()" << '\n';
            std::cerr << err << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::capacity()." << '\n';
            throw;
        }

        return out;
    }

    //===============================================================================================================

    void Population::remove_at(size_t which_list_position)
    {
        try
        {
            // removes by value in the active_individuals list
            /*std::vector<size_t>::iterator position = std::find(active_individuals.begin(), active_individuals.end(), which_animal_position);
            if (position != active_individuals.end()) // if the element was actually found
            {
                active_individuals.erase(position);
                individuals[ active_individuals[which_animal_position] ].clear();
            }*/

            // NOTE: remember, we are operating always through active_individuals list;
            // individuals determined by value of specific position in this list;
            // the indexes obviously consecutive but the values associated with these indexes
            // (which are pointing to the individuals list) because of removal on request.
            // When we retriave the position in the list we use this to index the individuals list
            // which is always consecutive because we are not removing anything from it, but
            // rather clear memory of a specific object (class Animal) in that list.

            // removes by position in the active_individuals list
            //std::cout<<"removing id_at: "<<individuals[active_individuals[which_list_position]].get_id()<<", position: "<<which_list_position<<", value: "<<active_individuals[which_list_position]<<"\n";
            individuals[active_individuals[which_list_position]].clear();
            active_individuals.erase(active_individuals.begin() + which_list_position); 
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::remove_at(size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &err)
        {
            std::cerr << "Exception in Population::remove_at(size_t)" << '\n';
            std::cerr << err << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::remove_at(size_t)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::clear()
    {
        try
        {
            for (size_t i = 0; i < active_individuals.size(); i++)
            {
                individuals[active_individuals[i]].clear();
            }
            active_individuals.clear();
            active_individuals.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::clear()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &err)
        {
            std::cerr << "Exception in Population::clear()" << '\n';
            std::cerr << err << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::clear()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::set_population(const std::string &haplotypes_fname, const std::string &genotype_structure, bool with_pedigree)
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

            size_t i_indiv = 0;
            for (size_t i = 0; i < n_haplotypes;)
            {
                std::vector<std::vector<bool>> markers(2, std::vector<bool>(n_markers));

                markers[0] = haplotypes[i];
                markers[1] = haplotypes[i + 1];

                Animal a;
                a.genome.set_genome(markers, genome_structure);

                if (with_pedigree)
                {
                    unsigned long id = pedigree[i][0];
                    int sex = (int)pedigree[i][1];
                    unsigned long sire = pedigree[i][2];
                    unsigned long dame = pedigree[i + 1][2];

                    a.set_id(id);
                    a.set_sex(sex);
                    a.set_sire(sire);
                    a.set_dame(dame);
                }

                individuals.push_back(a);
                active_individuals.push_back(i_indiv);

                i = i + 2;
                i_indiv++;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::set_population(const std::string &, const std::string &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::set_population(const std::string &, const std::string &)." << '\n';
            std::cerr << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::set_population(const std::string &, const std::string &)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<std::vector<unsigned long>> Population::get_genome_table()
    {
        try
        {
            if (active_individuals.size() >= 1)
                return individuals[active_individuals[0]].genome.get_snp_table();
            else
                throw std::string("Cannot return snp_table for the empty population");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::get_genome_table()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::get_genome_table()." << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::get_genome_table()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<std::vector<unsigned long>> Population::get_genome_structure()
    {
        try
        {
            if (active_individuals.size() >= 1)
                return individuals[active_individuals[0]].genome.get_genome_structure();
            else
                throw std::string("Cannot return genome structure for the empty population");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::get_genome_structure()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::get_genome_structure()." << '\n';
            std::cerr << "Reason: " << e << '\n';
            //std::cerr<<"id: "<<individuals[active_individuals[0]].get_id()<<", active: "<<individuals[active_individuals[0]].get_active()<<"\n";
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::get_genome_structure()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    size_t Population::get_ploidy()
    {
        try
        {
            if (active_individuals.size() >= 1)
                return (size_t)individuals[active_individuals[0]].genome.get_ploidy();
            else
                throw std::string("Cannot return ploidy for the empty population");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::get_ploidy()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::get_ploidy()." << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::get_ploidy()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::get_all_genotypes(const std::string &file_out)
    {
        try
        {
            if (active_individuals.size() >= 1)
            {
                std::vector<short> vect_out = individuals[active_individuals[0]].genome.get_genome();
                size_t n_snp = vect_out.size();

                vect_out.clear();
                vect_out.shrink_to_fit();

                FILE *pFile;
                pFile = fopen(file_out.c_str(), "a");
                if (pFile != NULL)
                {
                    for (size_t i = 0; i < active_individuals.size(); i++)
                    {
                        vect_out = individuals[active_individuals[i]].genome.get_genome();

                        for (size_t j = 0; j < n_snp; j++)
                        {
                            fprintf(pFile, "%3d", vect_out[j]);
                        }
                        fprintf(pFile, "\n");
                    }
                    fclose(pFile);
                }
                else
                    throw std::string("Cannot open file for writing out genotypes!");
            }
            else
                throw std::string("Cannot return the genotypes for the empty population");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::get_all_genotypes(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::get_all_genotypes(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::get_all_genotypes(const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::get_all_genotypes(std::vector<std::vector<short>> &vect_out)
    {
        try
        {
            if (active_individuals.size() >= 1)
            {
                for (size_t i = 0; i < active_individuals.size(); i++)
                    vect_out.push_back(individuals[active_individuals[i]].genome.get_genome());
            }
            else
                throw std::string("Cannot return the genotypes for the empty population");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::get_all_genotypes(std::vector<std::vector<short>> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::get_all_genotypes(std::vector<std::vector<short>> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::get_all_genotypes(std::vector<std::vector<short>> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<short> Population::get_genome_at(size_t which_genome, size_t locus)
    {
        try
        {
            if (active_individuals.size() < which_genome + 1)
                throw std::string("Cannot return the markers value for the requested individual, which is not exists!");
            else
                return individuals[active_individuals[which_genome]].genome.get_genome_at(locus);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::get_genome_at(size_t, size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::get_genome_at(size_t, size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::get_genome_at(size_t, size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::show_animals(size_t max_animals, size_t max_snps)
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
            std::cerr << "Exception in Population::show_animals(int, int)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::show_animals(int, int)." << '\n';
            throw;
        }
    }
    //===============================================================================================================

    unsigned long Population::id_at(size_t at)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            return individuals[active_individuals[at]].get_id();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::id_at(size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::id_at(size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::id_at(size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::id_at(size_t at, unsigned long id)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            individuals[active_individuals[at]].set_id(id);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::id_at(size_t, unsigned long)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::id_at(size_t, unsigned long)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::id_at(size_t, unsigned long)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    unsigned long Population::sire_at(size_t at)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            return individuals[active_individuals[at]].get_sire();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::sire_at(size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::sire_at(size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::sire_at(size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::sire_at(size_t at, unsigned long sire)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            individuals[active_individuals[at]].set_sire(sire);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::sire_at(size_t, unsigned long)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::sire_at(size_t, unsigned long)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::sire_at(size_t, unsigned long)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    unsigned long Population::dame_at(size_t at)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            return individuals[active_individuals[at]].get_dame();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::dame_at(size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::dame_at(size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::dame_at(size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::dame_at(size_t at, unsigned long dame)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            individuals[active_individuals[at]].set_dame(dame);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::dame_at(size_t, unsigned long)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::dame_at(size_t, unsigned long)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::dame_at(size_t, unsigned long)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    int Population::age_at(size_t at)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            return individuals[active_individuals[at]].get_age();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::age_at(size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::age_at(size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::age_at(size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::age_at(size_t at, int age)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            individuals[active_individuals[at]].set_age(age);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::age_at(size_t, int)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::age_at(size_t, int)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::age_at(size_t, int)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    short Population::sex_at(size_t at)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            return individuals[active_individuals[at]].get_sex();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::sex_at(size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::sex_at(size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::sex_at(size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::sex_at(size_t at, int sex)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            individuals[active_individuals[at]].set_sex(sex);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::sex_at(size_t, int)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::sex_at(size_t, int)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::sex_at(size_t, int)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    bool Population::alive_at(size_t at)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            return individuals[active_individuals[at]].get_alive();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::alive_at(size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::alive_at(size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::alive_at(size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::alive_at(size_t at, bool alive)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            individuals[active_individuals[at]].set_alive(alive);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::alive_at(size_t, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::alive_at(size_t, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::alive_at(size_t, bool)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    bool Population::isgenotyped_at(size_t at)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            return individuals[active_individuals[at]].get_isgenotyped();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::isgenotyped_at(size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::isgenotyped_at(size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::isgenotyped_at(size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::isgenotyped_at(size_t at, bool genotyped)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            individuals[active_individuals[at]].set_isgenotyped(genotyped);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::isgenotyped_at(size_t, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::isgenotyped_at(size_t, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::isgenotyped_at(size_t, bool)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<double> Population::phenotype_at(size_t at)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            return individuals[active_individuals[at]].get_phenotype();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::phenotype_at(size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::phenotype_at(size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::phenotype_at(size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::phenotype_at(size_t at, std::vector<double> &phen)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            return individuals[active_individuals[at]].set_phenotype(phen);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::phenotype_at(size_t, std::vector<double> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::phenotype_at(size_t, std::vector<double> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::phenotype_at(size_t, std::vector<double> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<double> Population::breedingvalue_at(size_t at)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            return individuals[active_individuals[at]].get_breeding_value();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::breedingvalue_at(size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::breedingvalue_at(size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::breedingvalue_at(size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::breedingvalue_at(size_t at, std::vector<double> &bv)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            return individuals[active_individuals[at]].set_breeding_value(bv);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::breedingvalue_at(size_t, std::vector<double> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::breedingvalue_at(size_t, std::vector<double> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::breedingvalue_at(size_t, std::vector<double> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::add(Animal &in_a)
    {
        try
        {
            size_t next_index = individuals.size();
            Animal a;
            a = in_a;
            individuals.push_back(a);
            active_individuals.push_back(next_index);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::add(Animal &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::add(Animal &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::add(Animal &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::reshape()
    {
        try
        {
            size_t i = 0;
            while ( i < individuals.size() )
            {
                if ( !individuals[i].get_active() )
                    individuals.erase( individuals.begin() + i );
                else
                    i = i + 1;
            }
            
            active_individuals.clear();
            active_individuals.shrink_to_fit();
            
            for (size_t j = 0; j < individuals.size(); j++)
                active_individuals.push_back(j);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::reshape()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &err)
        {
            std::cerr << "Exception in Population::reshape()" << '\n';
            std::cerr << err << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::reshape()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    //===============================================================================================================

    //===============================================================================================================

    //===============================================================================================================

}