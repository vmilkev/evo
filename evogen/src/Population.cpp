#include "Population.hpp"

namespace evogen
{
    //===============================================================================================================
    Population::Population()
    {
        origin_id = assign_origin_id();
    }
    //===============================================================================================================
    Population::~Population()
    {
        clear();
    }
    //===============================================================================================================
    popid_t Population::assign_origin_id()
    {
        popid_t id = 0;
        try
        {
            Utilites u;
            popid_t lower_bound = std::numeric_limits<popid_t>::min();
            popid_t upper_bound = std::numeric_limits<popid_t>::max();
            id = u.get_randi( lower_bound, upper_bound );
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::assign_origin_id()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::assign_origin_id()." << '\n';
            throw;
        }
        return id;
    }
    //===============================================================================================================
    size_t Population::get_popid()
    {
        return (size_t)origin_id;
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
            for (size_t i = 0; i < active_individuals.size(); i++)
            {
                if ( active_individuals[i] == -1 )
                    continue;
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
            if ( (ref_allele_probability > 1.0f) || (ref_allele_probability < 0.0f) )
            {
                std::string err("The reference allele probability should be in the range [0.0, 1.0].");
                throw err;
            }

            std::vector<std::vector<unsigned long>> genome_structure; /* length of each chromosome (bp) and markers density (constant, bp) */
            std::vector<double> gen_distance; /* genetic distance of each marker for all chromosomes (cM) */
            std::vector<float> gen_dist;
            std::vector<double> chromosome; /* chromosome id of each marker, will be used to make the genome_structure data table */

            evogen::IOInterface data;
            data.set_fname(genotype_structure);

            // determine if we work with MAP or STRUCTURE files
            if ( data.is_var_in_header(gdistance_name) ) // is a map file
            {
                data.fgetvar(gdistance_name, gen_distance); // reading as double
                std::copy(gen_distance.begin(), gen_distance.end(), std::back_inserter(gen_dist)); // converting to float

                data.fgetvar(chr_name, chromosome);

                // convert cromosome vector to genome_structure, where first col is num of snp in a cromosome and the second col is 1 (snps spaced by dist 1)
                // 1. get unique chromosomes
                std::vector<double> chr(chromosome);
                std::sort(chr.begin(), chr.end());
                auto it = std::unique(chr.begin(), chr.end());  
                chr.resize(std::distance(chr.begin(),it));

                genome_structure.resize( chr.size(), std::vector<unsigned long>(2,0) );

                // 2. make map between chromosome id and consecutive chromosome in genome_structure
                std::map<double, unsigned long> chr_ids;
                for (size_t i = 0; i < chr.size(); i++)
                    chr_ids[chr[i]] = i;

                // 3. create genome_structure table
                for (size_t i = 0; i < chromosome.size(); i++) // filling the first column of the table
                    genome_structure[ chr_ids[chromosome[i]] ][0] = genome_structure[ chr_ids[chromosome[i]] ][0] + 1;

                for (size_t i = 0; i < genome_structure.size(); i++) // filling the second column of the table
                    genome_structure[i][1] = 1;

                chr.clear();
                chr.shrink_to_fit();
                chromosome.clear();
                chromosome.shrink_to_fit();
                chr_ids.clear();
                gen_distance.clear();
                gen_distance.shrink_to_fit();
            }
            else // is a structure file
                data.fgetdata(genome_structure);

            for (size_t i = 0; i < nindividuals; i++)
            {
                Animal a;
                a.genome.set_genome(genome_structure, ref_allele_probability, (short)n_ploidy, gen_dist, origin_id); // also initializes gene ancestry
                individuals.push_back(a);
                active_individuals.push_back((poplen_t)i);
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
    void Population::set_population(const std::string &haplotypes_fname, const std::string &genotype_structure, bool with_pedigree)
    {
        // Create population using haplotypes file. This variant reads unnamed haplotypes, hence,
        // assign random IDs to individuals with unknown parents (unrelated population).
        
        try
        {
            std::vector<std::vector<unsigned long>> pedigree;
            std::vector<std::vector<bool>> haplotypes;                /* haplotypes of all individuals, rows = 2 * no.individuals */
            std::vector<std::vector<unsigned long>> genome_structure; /* length of each chromosome (bp) and markers density (constant, bp) */
            std::vector<double> gen_distance; /* genetic distance of each marker for all chromosomes (cM) */
            std::vector<float> gen_dist;
            std::vector<double> chromosome; /* chromosome id of each marker, will be used to make the genome_structure data table */

            /* because it is the animal population we have diploid genotype,
            hence, read piars of haplotype for each individual */

            evogen::IOInterface data;
            
            // Read the data files

            data.set_fname(haplotypes_fname);
            if (with_pedigree)
                data.fgetdata(pedigree, haplotypes);
            else
                data.fgetdata(haplotypes);

            data.set_fname(genotype_structure);

            // determine if we work with MAP or STRUCTURE files
            if ( data.is_var_in_header(gdistance_name) ) // is a map file
            {
                data.fgetvar(gdistance_name, gen_distance); // reading as double

                std::copy(gen_distance.begin(), gen_distance.end(), std::back_inserter(gen_dist)); // converting to float

                data.fgetvar(chr_name, chromosome);

                // convert cromosome vector to genome_structure, where first col is num of snp in a cromosome and the second col is 1 (snps spaced by dist 1)
                // 1. get unique chromosomes
                std::vector<double> chr(chromosome);
                std::sort(chr.begin(), chr.end());
                auto it = std::unique(chr.begin(), chr.end());  
                chr.resize(std::distance(chr.begin(),it));

                genome_structure.resize( chr.size(), std::vector<unsigned long>(2,0) );

                // 2. make map between chromosome id and consecutive chromosome in genome_structure
                std::map<double, unsigned long> chr_ids;
                for (size_t i = 0; i < chr.size(); i++)
                    chr_ids[chr[i]] = i;

                // 3. create genome_structure table
                for (size_t i = 0; i < chromosome.size(); i++) // filling the first column of the table
                    genome_structure[ chr_ids[chromosome[i]] ][0] = genome_structure[ chr_ids[chromosome[i]] ][0] + 1;

                for (size_t i = 0; i < genome_structure.size(); i++) // filling the second column of the table
                    genome_structure[i][1] = 1;

                chr.clear();
                chr.shrink_to_fit();
                chromosome.clear();
                chromosome.shrink_to_fit();
                chr_ids.clear();
                gen_distance.clear();
                gen_distance.shrink_to_fit();
            }
            else // is a structure file
                data.fgetdata(genome_structure);

            // 2. Creating animals with a specific genotype

            size_t n_haplotypes = haplotypes.size();
            size_t n_markers = haplotypes[0].size();

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
                a.genome.set_genome(markers, genome_structure, gen_dist, origin_id); // also initializes gene ancestry

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
                active_individuals.push_back((poplen_t)i_indiv);

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
    size_t Population::size()
    {
        size_t len = 0;
        for (auto const &e: active_individuals)
        {
            if (e != -1) len++;
        }
        return len;
        //return active_individuals.size();
    }
    //===============================================================================================================
    size_t Population::capacity()
    {
        return individuals.size();
    }
    //===============================================================================================================
    void Population::remove_at(size_t which_list_position)
    {
        try
        {
            // NOTE: we are operating always through active_individuals list;
            // individuals determined by value of specific position in this list;
            // the indexes obviously consecutive but the values associated with these indexes
            // (which are pointing to the individuals list) because of removal on request.
            // When we retriave the position in the list we use this to index the individuals list
            // which is always consecutive because we are not removing anything from it, but
            // rather clear memory of a specific object (class Animal) in that list.

            // removes by position in the active_individuals list
            if ( active_individuals[which_list_position] == -1 )
                throw std::string("Cannot remove individual due to the pointing position is out of range. This can happen if the individual was already removed/relocated by operation on different groups.");

            individuals[ (size_t)active_individuals[ which_list_position ] ].clear();
            active_individuals[which_list_position] = -1;
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
        for (size_t i = 0; i < active_individuals.size(); i++)
        {
            if ( active_individuals[i] == -1 )
                continue;
            individuals[(size_t)active_individuals[i]].clear();
        }
        individuals.shrink_to_fit();
        active_individuals.clear();
        active_individuals.shrink_to_fit();
    }
    //===============================================================================================================
    std::vector<std::vector<unsigned long>> Population::get_genome_table()
    {
        try
        {
            size_t count = 0;
            for (auto const &e: active_individuals) // looking for the first non-negative value
            {
                if ( e != -1 ) break;
                count++;
            }

            if ( count != active_individuals.size() )
                return individuals[(size_t)active_individuals[count]].genome.get_snp_table();
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
            size_t count = 0;
            for (auto const &e: active_individuals) // looking for the first non-negative value
            {
                if ( e != -1 ) break;
                count++;
            }

            if ( count != active_individuals.size() )
                return individuals[(size_t)active_individuals[count]].genome.get_genome_structure();
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
            size_t count = 0;
            for (auto const &e: active_individuals) // looking for the first non-negative value
            {
                if ( e != -1 ) break;
                count++;
            }

            if ( count != active_individuals.size() )
                return (size_t)individuals[(size_t)active_individuals[count]].genome.get_ploidy();
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
            size_t count = 0;
            for (auto const &e: active_individuals) // looking for the first non-negative value
            {
                if ( e != -1 ) break;
                count++;
            }

            if ( count != active_individuals.size() )
            {
                std::vector<short> vect_out = individuals[(size_t)active_individuals[count]].genome.get_genome();
                size_t n_snp = vect_out.size();

                vect_out.clear();
                vect_out.shrink_to_fit();

                FILE *pFile;
                pFile = fopen(file_out.c_str(), "a");
                if (pFile != NULL)
                {
                    for (size_t i = 0; i < active_individuals.size(); i++)
                    {
                        if (active_individuals[i] == -1)
                            continue;
                        vect_out = individuals[(size_t)active_individuals[i]].genome.get_genome();
                        unsigned long id = individuals[(size_t)active_individuals[i]].get_id();
                        
                        fprintf(pFile, "%lu %c", id, ' ');

                        for (size_t j = 0; j < n_snp; j++)
                            fprintf(pFile, "%3d", vect_out[j]);
                        
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
    void Population::get_ancestry(const std::string &file_out)
    {
        try
        {
            size_t count = 0;
            for (auto const &e: active_individuals) // looking for the first non-negative value
            {
                if ( e != -1 ) break;
                count++;
            }
            if ( count == active_individuals.size() )
                throw std::string("Cannot return ancestry genotypes for the empty population");

            std::ofstream fout(file_out);
            if (fout.is_open())
            {
                for (size_t i = 0; i < active_individuals.size(); i++)
                {
                    if (active_individuals[i] == -1)
                        continue;

                    std::vector<std::vector<ancestry_segment>> ancestry;
                    individuals[(size_t)active_individuals[i]].genome.get_ancestry(ancestry);
                    unsigned long id = individuals[(size_t)active_individuals[i]].get_id();

                    for (size_t i2 = 0; i2 < ancestry.size(); i2++) // expected/required ancestry.size() = n_ploidy 
                    {
                        fout << id;
                        for (size_t j = 0; j < ancestry[i2].size(); j++)
                            fout << " " << "[" << std::get<0>(ancestry[i2][j]) << ", " << std::get<1>(ancestry[i2][j]) << "]";
                        fout << '\n';
                    }
                }
                fout.close();
            }
            else
                std::cerr << "Error while opening file " + file_out << "\n";    
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::get_ancestry(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::get_ancestry(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::get_ancestry(const std::string &)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Population::get_genotypes(const std::string &file_out)
    {
        try
        {
            get_all_genotypes(file_out);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::get_genotypes(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::get_genotypes(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::get_genotypes(const std::string &)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Population::get_haplotypes(const std::string &file_out)
    {
        try
        {
            std::vector<std::vector<bool>> haplotypes;
            std::vector<std::vector<unsigned long>> snp_table;
            std::vector<float> gen_distance;

            get_all_haplotypes(haplotypes, snp_table, gen_distance);

            snp_table.clear();
            gen_distance.clear();

            std::vector<unsigned long> id_list;
            for (size_t i = 0; i < active_individuals.size(); i++)
            {
                if (active_individuals[i] == -1)
                    continue;
                id_list.push_back( individuals[(size_t)active_individuals[i]].get_id() );
            }

            size_t n_ploidy = get_ploidy();

            if ( (size_t)( id_list.size() * n_ploidy ) != haplotypes.size() )
                throw std::string("(size_t)( id_list.size() * get_ploidy() ) != haplotypes.size()");

            std::ofstream fout(file_out); // writing haplotypes to file
            if (fout.is_open())
            {
                int id_count = -1;
                for (size_t i = 0; i < haplotypes.size(); i++)
                {
                    if ( ! (i % n_ploidy) )
                        id_count++;
                    fout << id_list[(size_t)id_count] << " ";
                    for (size_t j = 0; j < haplotypes[i].size(); j++)
                        fout << haplotypes[i][j] << " ";
                    fout << '\n';
                }
                fout.close();
            }
            else
                std::cerr << "Error while opening file " + file_out << "\n";
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::get_haplotypes(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::get_haplotypes(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::get_haplotypes(const std::string &)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Population::get_all_genotypes(std::vector<std::vector<short>> &vect_out)
    {
        try
        {
            size_t count = 0;
            for (auto const &e: active_individuals) // looking for the first non-negative value
            {
                if ( e != -1 ) break;
                count++;
            }
            if ( count == active_individuals.size() )
                throw std::string("Cannot return the genotypes for the empty population");

            for (size_t i = 0; i < active_individuals.size(); i++)
            {
                if (active_individuals[i] == -1)
                    continue;
                vect_out.push_back(individuals[(size_t)active_individuals[i]].genome.get_genome());
            }
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
    void Population::get_all_haplotypes(std::vector<std::vector<bool>> &vect_out, std::vector<std::vector<unsigned long>> &out_snp_table, std::vector<float> &out_gen_distance)
    {
        try
        {
            size_t count = 0;
            for (auto const &e: active_individuals) // looking for the first non-negative value
            {
                if ( e != -1 ) break;
                count++;
            }
            if ( count == active_individuals.size() )
                throw std::string("Cannot return the halpotypes for the empty population");

            std::vector<std::vector<bool>> i_haplotypes;
            std::vector<std::vector<unsigned long>> i_gstructure;

            for (size_t i = 0; i < active_individuals.size(); i++)
            {
                if (active_individuals[i] == -1)
                    continue;
                individuals[(size_t)active_individuals[i]].genome.get_genome(i_haplotypes, i_gstructure); // get all haplotypes for individual i
                
                for ( size_t j = 0; j < i_haplotypes.size(); j++) // write haplotypes to the population haplotypes pool
                    vect_out.push_back( i_haplotypes[j] );
                                
                i_haplotypes.clear();
                i_haplotypes.shrink_to_fit();
                i_gstructure.clear();
                i_gstructure.shrink_to_fit();
            }

            out_snp_table = individuals[(size_t)active_individuals[count]].genome.get_snp_table();
            out_gen_distance = individuals[(size_t)active_individuals[count]].genome.get_gen_distance();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::get_all_haplotypes(std::vector<std::vector<bool>> &, std::vector<std::vector<unsigned long>> &, std::vector<float> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::get_all_haplotypes(std::vector<std::vector<bool>> &, std::vector<std::vector<unsigned long>> &, std::vector<float> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::get_all_haplotypes(std::vector<std::vector<bool>> &, std::vector<std::vector<unsigned long>> &, std::vector<float> &)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Population::get_ld(const std::string &out_file, bool full_info, int which_chr, unsigned int snp_step)
    {
        try
        {
            std::vector<std::vector<bool>> haplotypes;
            std::vector<std::vector<unsigned long>> snp_table;
            std::vector<float> gen_distance;
            std::vector<std::vector<unsigned long>> structure; // need to know distance between snps

            get_all_haplotypes(haplotypes, snp_table, gen_distance);

            if ( gen_distance.empty() ) // only if we are not using genetic distance
                structure = get_genome_structure();

            // do some file management ------------------
            std::string hapt_file(out_file + "_related_haplotypes.hpt");
            std::string prefix(out_file + "_chr_");
            std::string suffix(".matrld");
            std::string suffix2(".fq");
            std::string suffix3(".ld");
            std::vector<std::string> fnames(snp_table.size()); // LD matrix results for specific chromosome
            std::vector<std::string> fnames2(snp_table.size()); // Frequency results for specific chromosome
            std::vector<std::string> fnames3(snp_table.size()); // averaged LD results for specific chromosome
            for(size_t i = 0; i < snp_table.size(); i++)
            {
                fnames[i] = prefix + std::to_string(i) + suffix;
                fnames2[i] = prefix + std::to_string(i) + suffix2;
                fnames3[i] = prefix + std::to_string(i) + suffix3;
            }
            // ------------------------------------------

            size_t n_variants = haplotypes[0].size();
            size_t n_haplotypes = haplotypes.size();

            // ------- make bitset copy of vector<bool> ------------
            const size_t bts_const = 1000; // num of haplotyppes at one loci
            std::vector< std::vector<std::bitset<bts_const>> > haplotypes3;
            std::bitset<bts_const> bts;
            size_t l = 0;
            for (size_t snp_i = 0; snp_i < n_variants; snp_i++) // loop ove loci
            {
                std::vector<std::bitset<bts_const>> loci_i;
                for (size_t hap_i = 0; hap_i < n_haplotypes; hap_i++) // loop over all haplotypes
                {
                    l++;
                    if ( haplotypes[hap_i][snp_i] == 1 )
                        bts.set(l-1);
                    if ( l == bts_const || hap_i == n_haplotypes-1 )
                    {
                        loci_i.push_back(bts);
                        l = 0;
                        bts.reset();
                    }
                }
                haplotypes3.push_back(loci_i);
            }
            size_t n_bitsets = haplotypes3[0].size();

            haplotypes.clear();
            haplotypes.shrink_to_fit();
            // ------------------------------------------

            size_t chr_first = 0;
            size_t chr_last = snp_table.size()-1;

            if ( which_chr != -1 ) // if not for all chromosomes
            {
                if ( which_chr < (int)chr_first || which_chr > (int)chr_last )
                    throw std::string("Provided chromosome id is not correct!");
                chr_first = which_chr;
                chr_last = which_chr;
            }

            for (size_t chr_i = chr_first; chr_i <= chr_last; chr_i++) // calculate and write-out LD chromosome-by-chromosome
            {
                size_t snp_0 = snp_table[chr_i][0]; // first snp in chromosome
                size_t snp_1 = snp_table[chr_i][1]; // last snp in shromosome

                if ( snp_step < 1 || snp_step > (snp_1 - snp_0 + 1) )
                    throw std::string("Provided snp step value is not correct!");
                
                std::vector<size_t> variants_list;
                for (size_t snp_i = snp_0; snp_i <= snp_1;)
                {
                    variants_list.push_back(snp_i);
                    snp_i = snp_i + snp_step;
                }

                size_t chr_variants = variants_list.size(); // number of variants
                //size_t chr_variants = snp_1 - snp_0 + 1; // number of variants
                std::vector<float> freq(chr_variants, 0.0f); // frequency of individual variants
                std::vector<std::vector<float>> freq_paired(chr_variants, std::vector<float> (chr_variants, 0.0f));
                std::vector<std::vector<float>> corr(chr_variants, std::vector<float> (chr_variants, 0.0f));

#pragma omp parallel for
                for (size_t snp_i = 0; snp_i < chr_variants; snp_i++) // for every chromosomal variant
                {
                    for (size_t i_bts = 0; i_bts < n_bitsets; i_bts++)
                       freq[snp_i] = freq[snp_i] + (float)haplotypes3[ variants_list[snp_i] ][i_bts].count();

                    freq[snp_i] = freq[snp_i] / (float)n_haplotypes; // calculate individual frequecy

                    for (size_t snp_j = snp_i; snp_j < chr_variants; snp_j++) // calculate paired frequencies
                    {
                        for (size_t i_bts = 0; i_bts < n_bitsets; i_bts++)
                        {
                            std::bitset<bts_const> b = haplotypes3[ variants_list[snp_i] ][i_bts] & haplotypes3[ variants_list[snp_j] ][i_bts];
                            freq_paired[snp_i][snp_j] = freq_paired[snp_i][snp_j] + (float)b.count();
                        }

                        freq_paired[snp_i][snp_j] = freq_paired[snp_i][snp_j] / (float)n_haplotypes; // calculate individual frequecy
                    }
                }

                if ( full_info )
                {
                    std::ofstream fout(fnames2[chr_i]); // writing frequencies to file
                    if (fout.is_open())
                    {
                        for(auto const &v: freq)
                            fout<<v<<"\n";
                        fout.close();
                    }
                    else
                        std::cerr << "Error while opening file " + fnames2[chr_i] << "\n";
                }

#pragma omp parallel for
                for (size_t snp_i = 0; snp_i < chr_variants; snp_i++)
                {
                    float top_val = 0.0f;
                    float bot_val = 0.0f;
                    for (size_t snp_j = snp_i; snp_j < chr_variants; snp_j++)
                    {
                        top_val = (freq_paired[snp_i][snp_j] - freq[snp_i] * freq[snp_j]) * (freq_paired[snp_i][snp_j] - freq[snp_i] * freq[snp_j]);
                        bot_val = freq[snp_i] * (1 - freq[snp_i]) * freq[snp_j] * (1 - freq[snp_j]);
                        
                        if (bot_val != 0.0f)
                            corr[snp_i][snp_j] = top_val / bot_val;
                    }
                }

                if ( full_info )
                {
                    std::ofstream fout(fnames[chr_i]); // writing LD matrix to file
                    if (fout.is_open())
                    {
                        fout << std::setprecision(10);
                        for (size_t i = 0; i < chr_variants; i++)
                        {
                            for (size_t j = 0; j < chr_variants; j++)
                                fout<<corr[i][j]<<" ";
                            fout<<'\n';
                        }
                        fout.close();
                    }
                    else
                        std::cerr << "Error while opening file " + fnames[chr_i] << "\n";
                }

                // calculating averaged LD
                std::vector<float> non_zeros(chr_variants, 1.0f);
                for (size_t i = 0; i < chr_variants; i++)
                {
                    for (size_t j = i; j < chr_variants; j++)
                    {
                        if (corr[i][j] != 0.0f)
                        {
                            corr[0][j-i] = corr[0][j-i] + corr[i][j];
                            non_zeros[j-i]++;
                        }
                    }
                }
                for (size_t i = 0; i < chr_variants; i++)
                    corr[0][i] = corr[0][i]/non_zeros[i];

                std::ofstream fout(fnames3[chr_i]); // writing averaged LD to file
                if (fout.is_open())
                {
                    if ( !gen_distance.empty() )
                    {
                        fout<<"distance_(cM)"<<" "<<"correlation"<<"\n";
                        for(size_t i = 0; i < chr_variants; i++)
                            fout<<gen_distance[ variants_list[i] ]<<" "<<corr[0][i]<<"\n";
                    }
                    else
                    {
                        float distance = (float)structure[chr_i][1];
                        fout<<"distance_(bp)"<<" "<<"correlation"<<"\n";
                        for(size_t i = 0; i < chr_variants; i++)
                            fout<< (float)( variants_list[i] - snp_0 ) * distance <<" "<<corr[0][i]<<"\n";
                    }
                    
                    fout.close();
                }
                else
                    std::cerr << "Error while opening file " + fnames3[chr_i] << "\n";
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::get_ld(const std::string &, bool, int, unsigned int)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::get_ld(const std::string &, bool, int, unsigned int)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::get_ld(const std::string &, bool, int, unsigned int)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    std::vector<short> Population::get_genome_at(size_t which_genome, size_t locus)
    {
        try
        {
            if (size() < which_genome + 1)
                throw std::string("Cannot return the markers value for the requested individual, which is not exists!");
            
            if (active_individuals[which_genome] == -1)
                throw std::string("Cannot return the markers value for the requested individual, which is not alive!");
            
            return individuals[(size_t)active_individuals[which_genome]].genome.get_genome_at(locus);
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
    size_t Population::get_valid_pos(poplen_t in_pos)
    {
        // in_pos is position in the active_individuals if all -1 are ignored;
        // when calling pop.size() we get a number of non-negative values (active individuals),
        // this method allow access to i-th (in_pos) non-negative number 
        try
        {
            if ( in_pos >= (poplen_t)active_individuals.size() || in_pos < 0 )
                throw std::string("Illegal value of passed parameter!");
            
            poplen_t count = -1;
            size_t out_pos = 0;
            for (auto const &e: active_individuals) // looking for the first non-negative value
            {
                if ( e != -1 ) count++;
                if ( count == in_pos ) break;
                out_pos++;
            }
            if ( count == -1 )
                throw std::string("Cannot return a valid position for the empty population");
    
            return out_pos;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Population::get_valid_pos(size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Population::get_valid_pos(size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Population::get_valid_pos(size_t)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    bool Population::is_valid_pos(poplen_t in_pos)
    {
        if ( in_pos >= (poplen_t)active_individuals.size() || in_pos < 0 )
            return false;
        if ( active_individuals[in_pos] == -1 )
            return false;
        return true;
    }
    //===============================================================================================================
    unsigned long Population::id_at(size_t at)
    {
        try
        {
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");
                
            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            return individuals[(size_t)active_individuals[at]].get_id();
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            individuals[(size_t)active_individuals[at]].set_id(id);
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            return individuals[(size_t)active_individuals[at]].get_sire();
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            individuals[(size_t)active_individuals[at]].set_sire(sire);
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            return individuals[(size_t)active_individuals[at]].get_dame();
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            individuals[(size_t)active_individuals[at]].set_dame(dame);
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            return individuals[(size_t)active_individuals[at]].get_age();
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            individuals[(size_t)active_individuals[at]].set_age(age);
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            return (short)individuals[(size_t)active_individuals[at]].get_sex();
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            individuals[(size_t)active_individuals[at]].set_sex(sex);
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

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            return individuals[(size_t)active_individuals[at]].get_alive();
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            individuals[(size_t)active_individuals[at]].set_alive(alive);
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            return individuals[(size_t)active_individuals[at]].get_isgenotyped();
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            individuals[(size_t)active_individuals[at]].set_isgenotyped(genotyped);
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            std::vector<float> vals = individuals[(size_t)active_individuals[at]].get_phenotype();

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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            std::vector<float> phen_vect;

            pybind11::buffer_info buf1 = phen.request();

            if (buf1.ndim != 1)
                throw std::runtime_error("Number of dimensions must be one");

            float *ptr1 = static_cast<float *>(buf1.ptr);

            for (pybind11::ssize_t i = 0; i < buf1.shape[0]; i++)
                phen_vect.push_back(ptr1[i]);

            individuals[(size_t)active_individuals[at]].set_phenotype(phen_vect);
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            return individuals[(size_t)active_individuals[at]].get_phenotype();
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            individuals[(size_t)active_individuals[at]].set_phenotype(phen);
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            std::vector<float> vals = individuals[(size_t)active_individuals[at]].get_breeding_value();

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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            std::vector<float> bv_vect;

            pybind11::buffer_info buf1 = bv.request();

            if (buf1.ndim != 1)
                throw std::runtime_error("Number of dimensions must be one");

            float *ptr1 = static_cast<float *>(buf1.ptr);

            for (pybind11::ssize_t i = 0; i < buf1.shape[0]; i++)
                bv_vect.push_back(ptr1[i]);

            individuals[(size_t)active_individuals[at]].set_breeding_value(bv_vect);
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            return individuals[(size_t)active_individuals[at]].get_breeding_value();
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
            if ( at >= active_individuals.size() )
                throw std::string("Illegal value of passed parameter!");

            if (active_individuals[at] == -1)
                throw std::string("Cannot access the requested individual, it is not alive!");

            individuals[(size_t)active_individuals[at]].set_breeding_value(bv);
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
            poplen_t next_index = (poplen_t)individuals.size();
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
            size_t sz_1 = individuals.size(); // == active_individuals.size()
            individuals.resize( sz_1 + n_elements );
            active_individuals.resize( sz_1 + n_elements );
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
    void Population::add_at(Animal &in_a, size_t position)
    {
        try
        {
            if ( position >= individuals.size() )
                throw std::string("position >= individuals.size()");

            if ( position >= active_individuals.size() )
                throw std::string("position >= active_individuals.size()");

            individuals[position] = in_a;
            active_individuals[position] = (poplen_t)position;
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
            std::vector<size_t> valid_index;
            for (size_t i = 0; i < active_individuals.size(); i++)
            {
                if (active_individuals[i] == -1)
                    continue;
                valid_index.push_back(i);
            }

            std::vector<Animal> t_individuals( valid_index.size() );

#pragma omp parallel for
            for (size_t i = 0; i < valid_index.size(); i++)
                t_individuals[ i ] = individuals[ (size_t)active_individuals[ valid_index[i] ] ];

            individuals.clear();
            individuals.shrink_to_fit();

            individuals = t_individuals;

            t_individuals.clear();
            t_individuals.shrink_to_fit();
            
            active_individuals.clear();
            active_individuals.shrink_to_fit();

            for (size_t j = 0; j < individuals.size(); j++)
                active_individuals.push_back((poplen_t)j);
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