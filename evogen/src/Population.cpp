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
        clear();
    }

    //===============================================================================================================

    void Population::aging(int delta_t)
    {
        //    Increase an age of every individual in a population:
        //    new_age = current_age + delta_t.
        //    delta_t - the amount of time added to the current age of a particular individual;
        //              the time unit is arbitrary.
        
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
            if (ref_allele_probability > 1.0f)
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
                individuals[active_individuals[i]].clear();

            individuals.shrink_to_fit();

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
        // Create population using haplotypes file. This variant reads unnamed haplotypes, hence,
        // assign random IDs to individuals with unknown parents (unrelated population).
        
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

            //std::cout << "haplotypes, markers: " << n_haplotypes << ". " << n_markers << "\n";

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
            if (active_individuals.size() == 0)
                throw std::string("Cannot return the genotypes for the empty population");

            for (size_t i = 0; i < active_individuals.size(); i++)
                vect_out.push_back(individuals[active_individuals[i]].genome.get_genome());
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

    void Population::get_all_haplotypes(std::vector<std::vector<bool>> &vect_out, std::vector<std::vector<unsigned long>> &out_snp_table)
    {
        try
        {
            if (active_individuals.size() == 0)
                throw std::string("Cannot return the halpotypes for the empty population");

            std::vector<std::vector<bool>> i_haplotypes;
            std::vector<std::vector<unsigned long>> i_gstructure;

            for (size_t i = 0; i < active_individuals.size(); i++)
            {
                individuals[active_individuals[i]].genome.get_genome(i_haplotypes, i_gstructure); // get all haplotypes for individual i
                
                for ( size_t j = 0; j < i_haplotypes.size(); j++) // write haplotypes to the population haplotypes pool
                    vect_out.push_back( i_haplotypes[j] );
                                
                i_haplotypes.clear();
                i_haplotypes.shrink_to_fit();
                i_gstructure.clear();
                i_gstructure.shrink_to_fit();
            }

            out_snp_table = individuals[active_individuals[0]].genome.get_snp_table();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::get_all_haplotypes(std::vector<std::vector<bool>> &, std::vector<std::vector<unsigned long>> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::get_all_haplotypes(std::vector<std::vector<bool>> &, std::vector<std::vector<unsigned long>> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::get_all_haplotypes(std::vector<std::vector<bool>> &, std::vector<std::vector<unsigned long>> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::get_ld(const std::string &out_file)
    {
        try
        {
            std::vector<std::vector<bool>> haplotypes;
            std::vector<std::vector<unsigned long>> snp_table;

            get_all_haplotypes(haplotypes, snp_table);

            // do some file management ------------------
            std::string prefix(out_file + "_chr_");
            std::string suffix(".ld");
            std::vector<std::string> fnames(snp_table.size());
            for(size_t i = 0; i < snp_table.size(); i++)
                fnames[i] = prefix + std::to_string(i) + suffix;
            // ------------------------------------------

            size_t n_variants = haplotypes[0].size();
            size_t n_haplotypes = haplotypes.size();

            for (size_t chr_i = 0; chr_i < snp_table.size(); chr_i++) // calculate and write-out LD chromosome-by-chromosome
            {
                size_t snp_0 = snp_table[chr_i][0]; // first snp in chromosome
                size_t snp_1 = snp_table[chr_i][1]; // last snp in shromosome

                size_t chr_variants = snp_1 - snp_0 + 1; // number of variants
                std::vector<double> freq(chr_variants, 0.0); // frequency of individual variants
                std::vector<std::vector<double>> freq_paired(chr_variants, std::vector<double> (chr_variants, 0.0));
                std::vector<std::vector<double>> corr(chr_variants, std::vector<double> (chr_variants, 0.0));

#pragma omp parallel for
                for (size_t snp_i = snp_0; snp_i <= snp_1; snp_i++) // for every chromosomal variant
                {
                    for (size_t hap_i = 0; hap_i < n_haplotypes; hap_i++) // sum over all haplotypes
                        freq[snp_i-snp_0] = freq[snp_i-snp_0] + (double)haplotypes[hap_i][snp_i];
                    
                    freq[snp_i-snp_0] = freq[snp_i-snp_0] / (double)n_haplotypes; // calculate individual frequecy

                    for (size_t snp_j = snp_i; snp_j <= snp_1; snp_j++) // calculate paired frequencies
                    {
                        for (size_t hap_i = 0; hap_i < n_haplotypes; hap_i++) // sum over all haplotypes
                            if ( haplotypes[hap_i][snp_i] == haplotypes[hap_i][snp_j] )
                                freq_paired[snp_i-snp_0][snp_j-snp_0] = freq_paired[snp_i-snp_0][snp_j-snp_0] + (double)haplotypes[hap_i][snp_j];
                        
                        freq_paired[snp_i-snp_0][snp_j-snp_0] = freq_paired[snp_i-snp_0][snp_j-snp_0] / (double)n_haplotypes; // calculate individual frequecy
                    }
                }

#pragma omp parallel for
                for (size_t snp_i = snp_0; snp_i <= snp_1; snp_i++)
                {
                    for (size_t snp_j = snp_i; snp_j <= snp_1; snp_j++)
                    {
                        double top_val = (freq_paired[snp_i-snp_0][snp_j-snp_0] - freq[snp_i-snp_0] * freq[snp_j-snp_0]) * (freq_paired[snp_i-snp_0][snp_j-snp_0] - freq[snp_i-snp_0] * freq[snp_j-snp_0]);
                        double bot_val = freq[snp_i-snp_0] * (1 - freq[snp_i-snp_0]) * freq[snp_j-snp_0] * (1 - freq[snp_j-snp_0]);

                        //corr[snp_i-snp_0][snp_j-snp_0] = freq_paired[snp_i-snp_0][snp_j-snp_0] - freq[snp_i-snp_0] * freq[snp_j-snp_0]; // unscaled LD
                        
                        if (bot_val == 0.0f)
                            corr[snp_i-snp_0][snp_j-snp_0] = 0.0f;
                        else
                            corr[snp_i-snp_0][snp_j-snp_0] = top_val / bot_val;
                    }
                }
                
                std::ofstream fout(fnames[chr_i]);
                fout << std::setprecision(10);
                for (size_t i = 0; i < chr_variants; i++)
                {
                    for (size_t j = 0; j < chr_variants; j++)
                        fout<<corr[i][j]<<" ";
                    fout<<'\n';
                }
            }

            // write haplotypes
            std::ofstream fout("ld_related_haplotypes.hpt");
            for (size_t i = 0; i < haplotypes.size(); i++)
            {
                for (size_t j = 0; j < haplotypes[i].size(); j++)
                    fout << haplotypes[i][j] << " ";
                fout << '\n';
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::get_ld(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::get_ld(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::get_ld(const std::string &)" << '\n';
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

#ifdef UTEST
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
#endif

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
            if ( at >= active_individuals.size() )
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

#ifdef PYBIND

    pybind11::array_t<float> Population::phenotype_at(size_t at)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            //return individuals[active_individuals[at]].get_phenotype();

            std::vector<float> vals = individuals[active_individuals[at]].get_phenotype();

            size_t N = vals.size();

            pybind11::array_t<float, pybind11::array::c_style> arr({N, (size_t)1});

            auto arr_obj = arr.mutable_unchecked();

            for (size_t i = 0; i < N; i++)
            {
                arr_obj(i,0) = vals[i];
            }

            return arr;
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

    void Population::phenotype_at(size_t at, pybind11::array_t<float> phen)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            std::vector<float> phen_vect;

            pybind11::buffer_info buf1 = phen.request();

            if (buf1.ndim != 1)
                throw std::runtime_error("Number of dimensions must be one");

            float *ptr1 = static_cast<float *>(buf1.ptr);

            for (pybind11::ssize_t i = 0; i < buf1.shape[0]; i++)
                phen_vect.push_back(ptr1[i]);

            individuals[active_individuals[at]].set_phenotype(phen_vect);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::phenotype_at(size_t, std::vector<float> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::phenotype_at(size_t, std::vector<float> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::phenotype_at(size_t, std::vector<float> &)" << '\n';
            throw;
        }
    }

#endif

    //===============================================================================================================

    std::vector<float> Population::phenotype_at_cpp(size_t at)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            return individuals[active_individuals[at]].get_phenotype();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::phenotype_at_cpp(size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::phenotype_at_cpp(size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::phenotype_at_cpp(size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::phenotype_at_cpp(size_t at, std::vector<float> &phen)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            individuals[active_individuals[at]].set_phenotype(phen);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::phenotype_at_cpp(size_t, std::vector<float> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::phenotype_at_cpp(size_t, std::vector<float> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::phenotype_at_cpp(size_t, std::vector<float> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

#ifdef PYBIND

    pybind11::array_t<float> Population::breedingvalue_at(size_t at)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            //return pybind11::cast( individuals[active_individuals[at]].get_breeding_value() );
            
            std::vector<float> vals = individuals[active_individuals[at]].get_breeding_value();

            size_t N = vals.size();

            pybind11::array_t<float, pybind11::array::c_style> arr({N, (size_t)1});

            auto arr_obj = arr.mutable_unchecked();

            for (size_t i = 0; i < N; i++)
            {
                arr_obj(i, 0) = vals[i];
            }

            return arr;

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

    void Population::breedingvalue_at(size_t at, pybind11::array_t<float> bv)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            std::vector<float> bv_vect;

            pybind11::buffer_info buf1 = bv.request();

            if (buf1.ndim != 1)
                throw std::runtime_error("Number of dimensions must be one");

            float *ptr1 = static_cast<float *>(buf1.ptr);

            for (pybind11::ssize_t i = 0; i < buf1.shape[0]; i++)
                bv_vect.push_back(ptr1[i]);

            individuals[active_individuals[at]].set_breeding_value(bv_vect);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::breedingvalue_at(size_t, std::vector<float> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::breedingvalue_at(size_t, std::vector<float> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::breedingvalue_at(size_t, std::vector<float> &)" << '\n';
            throw;
        }
    }

#endif

    //===============================================================================================================

    std::vector<float> Population::breedingvalue_at_cpp(size_t at)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            return individuals[active_individuals[at]].get_breeding_value();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::breedingvalue_at_cpp(size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::breedingvalue_at_cpp(size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::breedingvalue_at_cpp(size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::breedingvalue_at_cpp(size_t at, std::vector<float> &bv)
    {
        try
        {
            if ( at >= active_individuals.size() || at < 0 )
                throw std::string("Illegal value of passed parameter!");

            individuals[active_individuals[at]].set_breeding_value(bv);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::breedingvalue_at_cpp(size_t, std::vector<float> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::breedingvalue_at_cpp(size_t, std::vector<float> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::breedingvalue_at_cpp(size_t, std::vector<float> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::add(Animal &in_a)
    {
        try
        {
            size_t next_index = individuals.size();
            individuals.push_back(in_a);
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

    void Population::resize(size_t n_elements)
    {
        try
        {
            size_t sz_1 = individuals.size();
            size_t sz_2 = active_individuals.size();
            
            individuals.resize( sz_1 + n_elements );            
            active_individuals.resize( sz_2 + n_elements );
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::resize(size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::resize(size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::resize(size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::add_at(Animal &in_a, size_t position, size_t position2)
    {
        try
        {
            if ( position >= individuals.size() )
                throw std::string("position >= individuals.size()");

            if ( position2 >= active_individuals.size() )
                throw std::string("position2 >= active_individuals.size()");

            individuals[position] = in_a;
            active_individuals[position2] = position;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::add_at(Animal &, size_t, size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::add_at(Animal &, size_t, size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::add_at(Animal &, size_t, size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Population::reshape()
    {
        try
        {
            //auto start = std::chrono::high_resolution_clock::now();

            std::vector<Animal> t_individuals( active_individuals.size() );

#pragma omp parallel for
            for (size_t i = 0; i < active_individuals.size(); i++)
                t_individuals[i] = individuals[ active_individuals[i] ];

            individuals.clear();
            individuals.shrink_to_fit();

            individuals = t_individuals;

            t_individuals.clear();
            t_individuals.shrink_to_fit();
            
            active_individuals.clear();
            active_individuals.shrink_to_fit();
            
            for (size_t j = 0; j < individuals.size(); j++)
                active_individuals.push_back(j);

            //auto stop = std::chrono::high_resolution_clock::now();
            //auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
            //std::cout <<"==> reshaping duration, sec "<< duration.count() << std::endl;
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