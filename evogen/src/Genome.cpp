#include "Genome.hpp"

namespace evogen
{
    //===============================================================================================================

    Genome::Genome()
    {
    }

    //===============================================================================================================

    Genome::~Genome() {}

    //===============================================================================================================

    void Genome::clear()
    {
        try
        {
            markers.clear();
            markers.shrink_to_fit();
            structure.clear();
            structure.shrink_to_fit();
            snp_table.clear();
            snp_table.shrink_to_fit();
            hotspots.clear();
            hotspots.shrink_to_fit();
            std_hotspots.clear();
            std_hotspots.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::clear()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Genome::clear()." << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::clear()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::def_snp_table()
    {
        try
        {
            size_t n_chr = structure.size();

            unsigned long previous = 0;

            for (size_t i = 0; i < n_chr; i++)
            {
                std::vector<unsigned long> range;

                range.push_back(previous);
                range.push_back(previous + structure[i][0] / structure[i][1] - 1);
                previous = previous + structure[i][0] / structure[i][1];

                snp_table.push_back(range);
            }

            def_hotspot_table();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::def_snp_table()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Genome::def_snp_table()." << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::def_snp_table()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::def_hotspot_table()
    {
        try
        {
            for (size_t i = 0; i < structure.size(); i++) // for each chromosome
            {
                float resolution = (float)structure[i][1];
                float n_snp_chr = snp_table[i][1] - snp_table[i][0] + 1; // length of chromosome in terms of number of markers
                size_t n_hotspot = std::ceil(n_snp_chr * freq_hotspot / resolution);

                if (n_hotspot < 2)
                    throw std::string("The chromosome resolution/density id too high leading to zerro number of recombination hotspots in the chromosome!");

                std::vector<size_t> locations;
                size_t step = n_snp_chr / n_hotspot;
                size_t pos = snp_table[i][0] + step;

                std_hotspots.push_back( (float)step / step_devider );

                size_t centromere_pos = n_hotspot / 4;
//std::cout<<"resolution "<<resolution<<" n_snp_chr "<<n_snp_chr<<" n_hotspot "<<n_hotspot<<" step "<<step<<" centromere_pos "<<centromere_pos<<'\n';                
                for (size_t j = 0; j < n_hotspot - 1; j++) // assign hotspots locations for the i-th chromosome
                {
                    if ( j == centromere_pos ) // there is no hotspot at centromere
                    {
                        locations.push_back( 0 );
                        pos = pos + step;
                    }
                    
                    locations.push_back( pos );
                    pos = pos + step;
                }

                hotspots.push_back(locations);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::def_hotspot_table()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Genome::def_hotspot_table()." << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::def_hotspot_table()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<std::vector<unsigned long>> Genome::get_snp_table()
    {
        try
        {
            if ( !markers.empty() )
                return snp_table;
            else
                throw std::string("The genome for this individual is empty!");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::get_snp_table()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Genome::get_snp_table()." << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_snp_table()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<std::vector<unsigned long>> Genome::get_genome_structure()
    {
        try
        {
            if ( !markers.empty() )
                return structure;
            else
                throw std::string("The genome for this individual is empty!");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::get_genome_structure()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Genome::get_genome_structure()." << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_genome_structure()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<short> Genome::get_genome()
    {
        // Reads (multiple) haplotypes and returns single genotype
        try
        {
            if ( markers.empty() )
                throw std::string("The genome for this individual is empty!");
            
            std::vector<short> out;
            size_t strands = markers.size();
            for ( size_t i = 0; i < markers[0].size(); i++ )
            {
                short strand_j = 0;
                for (size_t j = 0; j < strands; j++){
                    strand_j += (short)markers[j][i];
                }
                out.push_back( strand_j );
            }
            return out;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::get_genome()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Genome::get_genome()." << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_genome()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::set_genome(std::vector<std::vector<bool>> &snp, std::vector<std::vector<unsigned long>> &gstructure)
    {
        // uses prepared SNP variants either from file or from reproduction gamete

        try
        {
            markers = snp;
            structure = gstructure;
            int ploidy = snp.size();

            if (ploidy % 2 != 0)
                throw std::string("The number of ploidy is odd!");

            def_snp_table();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::set_genome(std::vector<std::vector<short>> &, std::vector<std::vector<float>> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Genome::set_genome(std::vector<std::vector<short>> &, std::vector<std::vector<float>> &)." << '\n';
            std::cerr << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::set_genome(std::vector<std::vector<short>> &, std::vector<std::vector<float>> &)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::set_genome(std::vector<std::vector<unsigned long>> &gstructure, float ref_allele_probability, short nploidy)
    {
        // Simulation of genome (snp variants)

        try
        {
            structure = gstructure;
            int ploidy = (int)nploidy;

            if (ploidy % 2 != 0)
                throw std::string("The number of provided ploidy is odd!");

            if ( ref_allele_probability > 1.0f || ref_allele_probability < 0.0f )
                throw std::string("ref_allele_probability > 1.0f || ref_allele_probability < 0.0f");

            def_snp_table();

            size_t snp_variants = 0;
            
            for (size_t i = 0; i < structure.size(); i++) // calculate snp variants in the genome
                snp_variants = snp_variants + std::floor(structure[i][0] / structure[i][1]);

            if ( ref_allele_probability == 0.0f || ref_allele_probability == 1.0f) // if we want each individual have a first haplotype with all 1, and the second with all 0
            {
                for (size_t i = 0; i < (size_t)nploidy; i++)
                {                    
                    if ( i % 2 )
                        markers.emplace_back(snp_variants, true);
                    else
                        markers.emplace_back(snp_variants, false);
                }
            }
            else // normal case
            {
                for (size_t i = 0; i < (size_t)nploidy; i++)
                {
                    std::vector<bool> variants;
                    for (size_t j = 0; j < snp_variants; j++)
                        variants.push_back(asign_snp_variant(ref_allele_probability));
                    markers.push_back(variants);
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::set_genome(std::vector<std::vector<unsigned long>> &, float, short)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Genome::set_genome(std::vector<std::vector<unsigned long>> &, float, short)." << '\n';
            std::cerr << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::set_genome(std::vector<std::vector<unsigned long>> &, float, short)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    short Genome::asign_snp_variant(float ref_allele_probability)
    {
        
        //    ref_allele_probability - is the probability of appearance of a reference allele in a loci
        //    snp_variant is 0 => num. of ref. alleles is 1
        //    snp_variant is 1 => num. of ref. alleles is 0
        
        short snp_variant = 0;

        try
        {
            Utilites u;

            int id = u.get_randi(1, 100);

            if (id <= ref_allele_probability * 100)
                snp_variant = 0; // in this loci we observe a reference allele, the snp_variant is 0
            else
                snp_variant = 1; // in this loci we observe an alternative allele, the snp_variant is 1
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::asign_snp_variant(float)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::asign_snp_variant(float)." << '\n';
            throw;
        }

        return snp_variant;
    }

    //===============================================================================================================

    short Genome::get_ploidy()
    {
        short ploidy = 0;

        try
        {
            if ( markers.empty() )
                throw std::string("The genome for this individual is empty!");

            ploidy = markers.size();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::get_ploidy()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_ploidy()." << '\n';
            throw;
        }

        return ploidy;
    }

    //===============================================================================================================

    std::vector<short> Genome::get_genome_at(size_t locus)
    {
        try
        {
            if ( markers.empty() )
                throw std::string("The genome for this individual is empty!");

            if ( markers.size() < 2 )
                throw std::string("The genome is empty, cannot retreat the markers values!");

            size_t n_ploidy = (size_t)get_ploidy();
            
            std::vector<short> v;

            for (size_t i = 0; i < n_ploidy; i++)
                v.push_back( (short)markers[i][locus] );
            
            return v;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::get_genome_at(size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Genome::get_genome_at(size_t)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_genome_at(size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::show_genome(int max_variants)
    {
        try
        {
            if ( markers.empty() )
                throw std::string("The genome for this individual is empty!");

            size_t n_chr = structure.size();

            std::cout << "n. chromosomes = " << n_chr << "\n";

            for (const auto &e : markers)
            {
                size_t max_markers = std::min(max_variants, (int)e.size());

                for (size_t i = 0; i < max_markers; i++)
                    std::cout << e[i] << " ";
                std::cout << "\n";
            }

            for (const auto &e : structure)
            {
                for (size_t i = 0; i < e.size(); i++)
                    std::cout << e[i] << " ";
                std::cout << "\n";
            }

            std::cout << "snp table:"
                      << "\n";

            for (const auto &e : snp_table)
            {
                for (size_t i = 0; i < e.size(); i++)
                    std::cout << e[i] << " ";
                std::cout << "\n";
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::get_genome()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_genome()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::get_genome(std::vector<std::vector<bool>> &haplotypes, std::vector<std::vector<unsigned long>> &gstructure)
    {
        try
        {
            if ( markers.empty() )
                throw std::string("The genome for this individual is empty!");

            haplotypes = markers;
            gstructure = structure;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::get_genome(std::vector<std::vector<bool>> &, std::vector<std::vector<unsigned long>> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_genome(std::vector<std::vector<bool>> &, std::vector<std::vector<unsigned long>> &)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::get_reproduction_gamete(std::vector<std::vector<bool>> &out_gamete, size_t cross_per_chr, float mut_freq)
    {
        // out_sex_chr_id: for each produced gamete indicates where does sex chromosome (the last one) comes from in the gamete: 0 - paternal, 1 - maternal.

        try
        {
            if ( markers.empty() )
                throw std::string("The genome for this individual is empty!");

            recombination(out_gamete, cross_per_chr);
            mutation(out_gamete, mut_freq);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::get_reproduction_gamete( std::vector<std::vector<bool>> &, size_t, float )." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_reproduction_gamete( std::vector<std::vector<bool>> &, size_t, float )." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::mutation(std::vector<std::vector<bool>> &in_gamete, float mut_freq)
    {
        // mut_freq - probability of mutation of a single snp during meyosis
        try
        {
            if ( mut_freq >= 1.0f )
                throw std::string("Mutation frequency should be in the range [0.0, 1.0].");

            for (size_t i_set = 0; i_set < in_gamete.size(); i_set++) // for each specific gamete
            {
                // Sample mutation points
                Utilites u;

                std::vector<int> events = u.get_bin_rand(1, in_gamete[i_set].size(), (double)mut_freq, false); // based on the mut_freq probability, calculate the number of mutations (modified snps) for entire gamete
                size_t events_genome = events[0]; // number of mutations
                
//std::cout<<"mutation events: "<<events_genome<<'\n';
                std::vector<size_t> locations = u.get_uni_rand(events_genome, (size_t)0, in_gamete[i_set].size() - 1, false); // sample of mutation locations (specific snps subject to modification)

                for (size_t i = 0; i < events_genome; i++)
                {
                    size_t point = locations[i];
//std::cout << "Mutation point: " << point << "; current val: " << in_gamete[i_set][point] << "; obtained val: ";
                    if ( !in_gamete[i_set][point] ) // if 0
                        in_gamete[i_set][point] = true; // change to 1
                    else
                        in_gamete[i_set][point] = false; // change to 0
//std::cout << in_gamete[i_set][point] << "\n";
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::mutation( std::vector<std::vector<bool>> &, float )" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::mutation( std::vector<std::vector<bool>> &, float )" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::recombination(std::vector<std::vector<bool>> &out_gamete, size_t n_crosses)
    {
        try
        {
            size_t i_set = 0;
            for (size_t n_sets = 0; n_sets < markers.size()/2; n_sets++) // loop over ploidy sets of DNA pairs
            {
                Utilites u;

                size_t n_chr = structure.size();

                //......................................................

                /* Table A.
                Random table of 2 cols to select between paternal and maternal chromosomes
                which aligned (distributed) between two cells after crossing over. Two produced cells
                will consist of mixture of paternal and maternal chromosomes, among which
                there are one - original, and one - crossed.
                */

                std::vector<std::vector<short>> is_paternal;

                for (size_t i = 0; i < n_chr; i++)
                {
                    std::vector<short> i_chr;

                    int rnum = u.get_randi(1, 100);

                    short paternal = 0;
                    short maternal = 1;

                    if (rnum >= 50)
                    {
                        paternal = 1;
                        maternal = 0;
                    }
                    i_chr.push_back(paternal);
                    i_chr.push_back(maternal);

                    is_paternal.push_back(i_chr);
                }

                //......................................................

                /* Table B.
                Random table of 4 cols to select between original and crossed chromosoms
                distributed among 4 gamete cells. Two cells will consist original (paternal or maternal)
                chromosomes, and two other cells will consist crossed chromosomes.
                */

                std::vector<std::vector<short>> is_crossed;

                for (size_t i = 0; i < n_chr; i++)
                {
                    // Utilites u;

                    std::vector<short> i_chr;

                    for (size_t j = 0; j < 2; j++)
                    {
                        int rnum = u.get_randi(1, 100);

                        short crossed_0 = 0;
                        short crossed_1 = 1;

                        if (rnum >= 50)
                        {
                            crossed_0 = 1;
                            crossed_1 = 0;
                        }
                        i_chr.push_back(crossed_0);
                        i_chr.push_back(crossed_1);
                    }

                    is_crossed.push_back(i_chr);
                }
                //......................................................
                // Sample which cols from the previous two tables will be used:
                //......................................................

                /* Sample the col of the table B:
                among 4 final gamete cells define (sample) which cell is selected for reproduction;
                the variable which_cell will be used as index in the table is_crossed.
                */

                short which_cell = 0;
                // Utilites u;
                int sample_cell = u.get_randi(1, 100);

                if (sample_cell <= 25)
                    which_cell = 1;
                if (sample_cell > 25 && sample_cell <= 50)
                    which_cell = 2;
                if (sample_cell > 50 && sample_cell <= 75)
                    which_cell = 3;
                if (sample_cell > 75 && sample_cell <= 100)
                    which_cell = 0;

                //......................................................

                /* Sample the col of the table A:
                among 2 pre-gamete cells define (sample) which cell line
                should be used (is selected) for reproduction;
                the var which_pat will be used as index in the table is_paternal.
                */

                short which_pat = 0;
                int sample_pat = u.get_randi(1, 100);

                if (sample_pat <= 50)
                    which_pat = 1;
                if (sample_pat > 50)
                    which_pat = 0;

                //......................................................
                // Here we know where does sex chromosome will come: from sire or dame

                //out_sex_chr_id = is_paternal[n_chr - 1][which_pat];

                //......................................................

                std::vector<bool> gamete;

                for (size_t i = 0; i < n_chr; i++)
                {
                    std::vector<bool> chromatide;
                    short apply_cross = is_crossed[i][which_cell];  // defines the sequence of original/crossed chromosomes
                    short which_strand = is_paternal[i][which_pat]; // defines the sequence of chromosome origin: paternal/maternal

                    // the values of which_cell & which_pat are constant for a specific gamete, and expected to be different for every new generated gamete.

                    //std::cout << "chr.no: " << i << "; apply_cross & which strand: " << apply_cross << ", " << which_strand << "; which_pat & which_cell: " << which_pat << ", " << which_cell << "\n";

                    if (apply_cross)
                        chromatide = crossover(i, which_strand, n_crosses, i_set);
                    else
                        chromatide = get_chromatide(i, which_strand, i_set);

                    for (size_t i = 0; i < chromatide.size(); i++)
                        gamete.push_back(chromatide[i]);
                }

                out_gamete.push_back(gamete);

                i_set = i_set + 2; // select the next set of DNA strands, for ploidy > 2

            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::recombination( std::vector<std::vector<bool>> &, short &, size_t )." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::recombination( std::vector<std::vector<bool>> &, short &, size_t )." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<bool> Genome::crossover(size_t in_chr, size_t out_strand, size_t n_crossovers, size_t which_strands_set)
    {
        std::vector<bool> chromatide;

        try
        {
            //crossover_model_gamma(in_chr, out_strand, n_crossovers, which_strands_set, chromatide);
            crossover_model_uniform(in_chr, out_strand, n_crossovers, which_strands_set, chromatide);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::crossover(size_t, size_t, size_t, size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::crossover(size_t, size_t, size_t, size_t)" << '\n';
            throw;
        }

        return chromatide;
    }

    //===============================================================================================================

    void Genome::crossover_model_gamma(size_t in_chr, size_t out_strand, size_t n_crossovers, size_t which_strands_set, std::vector<bool> &out_chromatide)
    {
        try
        {
            Utilites u;

            size_t first_snp = snp_table[in_chr][0];
            size_t last_snp = snp_table[in_chr][1];

            std::vector<long> cross_locations;

            long pos = (long)u.get_randi( (size_t)0, hotspots[in_chr].size()-1 ); // sample position in hotspots table for a very first crossover event

            if ( (pos != (size_t)0) && (hotspots[in_chr][pos] == (size_t)0) ) // if we catch the centromere, change the position
                pos--;

            std::vector<float> location = u.get_norm_rand(1, (float)hotspots[in_chr][pos], std_hotspots[in_chr], false); // the actual crossover position sampled from normal distr using pos as a mean

            if ( (size_t)location[0] > first_snp && (size_t)location[0] < last_snp )
                cross_locations.push_back( (size_t)location[0] ); // very first crossover location
            else
                cross_locations.push_back( hotspots[in_chr][pos] );

            size_t chr_length = last_snp - first_snp + 1; // in number of snps (markers)
            float beta = 1.0 / ((float)n_crossovers / (float)chr_length );
            float alpha = 2.0;

            size_t step = std_hotspots[in_chr] * step_devider;

            // Left move along chromosome
            size_t index = 0;
            for (size_t i = 0; i < 2 * n_crossovers; i++)
            {
                std::vector<float> left_cross_point = u.get_gamma_rand(1, alpha, beta, false);

                if ( ( cross_locations[index] - (long)std::round(left_cross_point[0]) ) <= (long)first_snp )
                    break;
                size_t position = ( cross_locations[index] - (long)std::round(left_cross_point[0]) - first_snp ) / step; // position in the hotspots table
                if ( hotspots[in_chr][position] == (size_t)0 ) // this is centromere, where there is no crossover
                    continue;
                std::vector<float> locat = u.get_norm_rand(1, (float)hotspots[in_chr][position], std_hotspots[in_chr], false); // sample actual recombination position
                if ( (size_t)locat[0] > first_snp && (size_t)locat[0] < last_snp )
                {
                    cross_locations.push_back( (size_t)locat[0] ); // crossover location
                    index++;
                }
                else
                    break;
            }

            std::vector<float> right_cross_point = u.get_gamma_rand(1, alpha, beta, false);
            size_t position0 = ( cross_locations[0] + (long)std::round(right_cross_point[0]) - first_snp ) / step; // position in the hotspots table
            
            if ( (position0 != (size_t)0) && (hotspots[in_chr][position0] == (size_t)0) ) // if we catch the centromere, resample
                position0--;

            if ( position0 > 0 && position0 < hotspots[in_chr].size() )
            {
                std::vector<float> locat = u.get_norm_rand(1, (float)hotspots[in_chr][position0], std_hotspots[in_chr], false); // sample actual recombination position
                
                if ( (size_t)locat[0] > first_snp && (size_t)locat[0] < last_snp )
                {
                    cross_locations.push_back( (size_t)locat[0] ); // crossover location
                    size_t processed_crosses = cross_locations.size() - 1;

                    index = 0;
                    for (size_t i = 1; i < 2 * n_crossovers - 1; i++)
                    {
                        right_cross_point.clear(); right_cross_point.shrink_to_fit();

                        right_cross_point = u.get_gamma_rand(1, alpha, beta, false);

                        if ( ( cross_locations[processed_crosses + index] + (long)std::round(right_cross_point[0]) ) >= (long)last_snp )
                            break;
                        
                        size_t position = ( cross_locations[processed_crosses + index] + (long)std::round(right_cross_point[0]) - first_snp ) / step; // position in the hotspots table
                        if ( hotspots[in_chr][position] == (size_t)0 )
                            continue;

                        if ( (position > 0) && (position < hotspots[in_chr].size()) )
                        {                        
                            //std::vector<float> locat = u.get_norm_rand(1, (float)hotspots[in_chr][position], std_hotspots[in_chr], false); // sample actual recombination position
                            
                            if ( (size_t)locat[0] > first_snp && (size_t)locat[0] < last_snp )
                            {
                                cross_locations.push_back( (size_t)locat[0] ); // crossover location
                                index++;
                            }
                            else
                                break;
                        }
                        break;
                    }
                }
            }

            std::sort(cross_locations.begin(), cross_locations.end());

            /*std::cout << "hotspots on the chromosome " << in_chr << ": ";
            for (size_t i = 0; i < hotspots[in_chr].size(); i++)
                std::cout<<hotspots[in_chr][i]<<" ";
            std::cout << "\n";
            //std::cout << "\n";
            std::cout << "cross locations in chromosome " << in_chr << ": ";
            for (size_t i = 0; i < cross_locations.size(); i++)
            {
                std::cout << cross_locations[i] << " ";
            }
            std::cout << "\n";*/
            

            // .........................................................
            // Staring exchnging of sections of chromosome
            // .........................................................

            // temporal containers
            std::vector<bool> chromatide_pat;
            std::vector<bool> chromatide_mat;

            for (size_t i = first_snp; i <= last_snp; i++)
            {
                chromatide_pat.push_back(markers[0 + which_strands_set][i]);
                chromatide_mat.push_back(markers[1 + which_strands_set][i]);
            }

            for (size_t i = 0; i < cross_locations.size(); i++)
            {
                std::vector<bool> chromatide_tmp(1, 0);

                for (size_t j = 0; j < cross_locations[i] - first_snp; j++)
                {
                    chromatide_tmp[0] = chromatide_pat[j];

                    chromatide_pat[j] = chromatide_mat[j];
                    chromatide_mat[j] = chromatide_tmp[0];
                }
            }
            // .........................................................

            // Return gamete's chromosom

            if (out_strand == 0)
                out_chromatide = chromatide_pat;
            else
                out_chromatide = chromatide_mat;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::crossover_model_gamma(size_t, size_t, size_t, size_t, std::vector<bool> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::crossover_model_gamma(size_t, size_t, size_t, size_t, std::vector<bool> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::crossover_model_uniform(size_t in_chr, size_t out_strand, size_t n_crossovers, size_t which_strands_set, std::vector<bool> &out_chromatide)
    {
        try
        {
            Utilites u;

            size_t first_snp = snp_table[in_chr][0];
            size_t last_snp = snp_table[in_chr][1];

            std::vector<long> cross_locations;

            std::vector<int> sampled_crossovers = u.get_bin_rand(1, (int)n_crossovers, 0.8, false); // sample the actual number of crossovers for this chromosome

            long pos = 0;
            std::vector<float> location;

            for (size_t i = 0; i < (size_t)sampled_crossovers[0]; i++)
            {
                pos = (long)u.get_randi( (size_t)0, hotspots[in_chr].size()-1 ); // sample position in hotspots table for a crossover event
                if ( (pos != (size_t)0) && (hotspots[in_chr][pos] == (size_t)0) ) // if we catch the centromere, change the position
                    pos--;
                location = u.get_norm_rand(1, (float)hotspots[in_chr][pos], std_hotspots[in_chr], false); // the actual crossover position within the hotspot sampled from normal distr using pos as a mean
                
                if ( (size_t)location[0] > first_snp && (size_t)location[0] < last_snp )
                    cross_locations.push_back( (size_t)location[0] ); // very first crossover location
                else
                    cross_locations.push_back( hotspots[in_chr][pos] );
            }

            std::sort(cross_locations.begin(), cross_locations.end());

            /*std::cout << "hotspots on the chromosome " << in_chr << ": ";
            for (size_t i = 0; i < hotspots[in_chr].size(); i++)
                std::cout<<hotspots[in_chr][i]<<" ";
            std::cout << "\n";
            //std::cout << "\n";
            std::cout << "cross locations in chromosome " << in_chr << ": ";
            for (size_t i = 0; i < cross_locations.size(); i++)
            {
                std::cout << cross_locations[i] << " ";
            }
            std::cout << "\n";*/            

            // .........................................................
            // Staring exchnging of sections of chromosome
            // .........................................................

            // temporal containers
            std::vector<bool> chromatide_pat;
            std::vector<bool> chromatide_mat;

            for (size_t i = first_snp; i <= last_snp; i++)
            {
                chromatide_pat.push_back(markers[0 + which_strands_set][i]);
                chromatide_mat.push_back(markers[1 + which_strands_set][i]);
            }

            for (size_t i = 0; i < cross_locations.size(); i++)
            {
                std::vector<bool> chromatide_tmp(1, 0);

                for (size_t j = 0; j < cross_locations[i] - first_snp; j++)
                {
                    chromatide_tmp[0] = chromatide_pat[j];

                    chromatide_pat[j] = chromatide_mat[j];
                    chromatide_mat[j] = chromatide_tmp[0];
                }
            }
            // .........................................................

            // Return gamete's chromosom

            if (out_strand == 0)
                out_chromatide = chromatide_pat;
            else
                out_chromatide = chromatide_mat;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::crossover_model_uniform(size_t, size_t, size_t, size_t, std::vector<bool> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::crossover_model_uniform(size_t, size_t, size_t, size_t, std::vector<bool> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<bool> Genome::get_chromatide(size_t which_chr, size_t which_strand, size_t which_strands_set)
    {
        std::vector<bool> chromatide;

        try
        {
            size_t first_snp = snp_table[which_chr][0];
            size_t last_snp = snp_table[which_chr][1];
            for (size_t i = first_snp; i <= last_snp; i++)
            {
                chromatide.push_back(markers[which_strand + which_strands_set][i]);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::get_chromatide(size_t, size_t, size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_chromatide(size_t, size_t, size_t)" << '\n';
            throw;
        }

        return chromatide;
    }

    //===============================================================================================================

    void Genome::show_recombination()
    {
        try
        {
            if ( markers.empty() )
                throw std::string("The genome for this individual is empty!");

            std::vector<std::vector<bool>> chromatide;

            get_reproduction_gamete(chromatide, 3, 2);

            std::cout << "\n";

            std::cout << "showing gamete with "<<chromatide.size()<<" strands:"<< "\n";

            for (size_t j = 0; j < chromatide.size(); j++)
            {
                size_t snps = std::min((size_t)100, chromatide[j].size());

                for (size_t i = 0; i < snps; i++)
                {
                    std::cout << chromatide[j][i] << " ";
                }
                std::cout << "\n";
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::show_recombination()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::show_recombination()." << '\n';
            throw;
        }
    }

    //===============================================================================================================


    //===============================================================================================================

}
