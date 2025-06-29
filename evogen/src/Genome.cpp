#include "Genome.hpp"

namespace evogen
{
    //===============================================================================================================

    Genome::Genome(){}

    //===============================================================================================================

    Genome::~Genome(){}

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
            gen_distance.clear();
            gen_distance.shrink_to_fit();
            gen_ancestry.clear();
            gen_ancestry.shrink_to_fit();
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
            std::cerr << "Reason: " << e << '\n';
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

            if ( gen_distance.empty() )
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
            std::cerr << "Reason: " << e << '\n';
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
                // float freq_hotspot = 0.05f;
                // float step_devider = 5.0f;

                float n_snp_chr = snp_table[i][1] - snp_table[i][0] + 1; // length of chromosome in terms of number of markers
                // size_t n_hotspot = 10;//std::ceil(n_snp_chr * freq_hotspot);

                if (n_snp_chr < n_hotspot)
                    throw std::string("The chromosome resolution/density is too coarse leading to zerro number of recombination hotspots in the chromosome!");

                std::vector<size_t> locations;
                size_t step = n_snp_chr / n_hotspot;
                size_t pos = snp_table[i][0] + step;

                // std_hotspots.push_back( (float)0.0f );
                std_hotspots.push_back((float)step / 10);
                // std_hotspots.push_back( (float)step );

                size_t centromere_pos_start = snp_table[i][1] - (size_t)std::round(n_snp_chr / 2.5);
                size_t centromere_pos_end = snp_table[i][1] - (size_t)std::round(n_snp_chr / 2.7);

                while (pos < snp_table[i][1])
                {
                    if ((pos >= centromere_pos_start) && (pos <= centromere_pos_end)) // there is no hotspot at centromere
                    {
                        pos = pos + step;
                        continue;
                    }
                    locations.push_back(pos);
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
            std::cerr << "Reason: " << e << '\n';
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
            if (!markers.empty())
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
            std::cerr << "Reason: " << e << '\n';
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
            if (!markers.empty())
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
            std::cerr << "Reason: " << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_genome_structure()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<float> Genome::get_gen_distance()
    {
        try
        {
            return gen_distance;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::get_gen_distance()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Genome::get_gen_distance()." << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_gen_distance()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<short> Genome::get_genome()
    {
        // Reads (multiple) haplotypes and returns single genotype
        try
        {
            if (markers.empty())
                throw std::string("The genome for this individual is empty!");

            std::vector<short> out;
            size_t strands = markers.size();
            for (size_t i = 0; i < markers[0].size(); i++)
            {
                short strand_j = 0;
                for (size_t j = 0; j < strands; j++)
                {
                    strand_j += (short)markers[j][i];
                }
                out.push_back(strand_j);
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
            std::cerr << "Reason: " << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_genome()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::set_genome(std::vector<std::vector<bool>> &snp, std::vector<std::vector<unsigned long>> &gstructure, std::vector<float> &gen_dist, popid_t pop_origin_id)
    {
        // uses prepared SNP variants either from file or from reproduction gamete; initializes gen_ancestry

        try
        {
            markers = snp;
            structure = gstructure;
            int ploidy = snp.size();

            if (ploidy % 2 != 0)
                throw std::string("The number of ploidy is odd!");

            if (!gen_dist.empty())
                std::copy(gen_dist.begin(), gen_dist.end(), std::back_inserter(gen_distance));

            def_snp_table();

            for (size_t i = 0; i < (size_t)ploidy; i++) // add very first tuples pointing that the entire genome originates from the pop_origin_id 
            {
                std::vector<ancestry_segment> ancestry;
                ancestry.push_back( std::make_tuple(0,pop_origin_id) );
                gen_ancestry.push_back(ancestry);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::set_genome(std::vector<std::vector<short>> &, std::vector<std::vector<float>> &, std::vector<float> &, popid_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Genome::set_genome(std::vector<std::vector<short>> &, std::vector<std::vector<float>> &, std::vector<float> &, popid_t)." << '\n';
            std::cerr << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::set_genome(std::vector<std::vector<short>> &, std::vector<std::vector<float>> &, std::vector<float> &, popid_t)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::set_genome(std::vector<std::vector<bool>> &snp, std::vector<std::vector<ancestry_segment>> &ancestry, std::vector<std::vector<unsigned long>> &gstructure, std::vector<float> &gen_dist)
    {
        // uses prepared SNP variants either from file or from reproduction gamete; does not initialize gen_ancestry

        try
        {
            markers = snp;
            gen_ancestry = ancestry;
            structure = gstructure;
            int ploidy = snp.size();

            if (ploidy % 2 != 0)
                throw std::string("The number of ploidy is odd!");

            if (!gen_dist.empty())
                std::copy(gen_dist.begin(), gen_dist.end(), std::back_inserter(gen_distance));

            def_snp_table();

            if ( gen_ancestry[0].empty() || gen_ancestry[1].empty() )
                throw std::string("The gen_ancestry vector is empty. Expected to be at least initialized, hence have the size = 1.");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::set_genome(std::vector<std::vector<short>> &, std::vector<std::vector<ancestry_segment>> &, std::vector<std::vector<float>> &, std::vector<float> &, popid_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Genome::set_genome(std::vector<std::vector<short>> &, std::vector<std::vector<ancestry_segment>> &, std::vector<std::vector<float>> &, std::vector<float> &, popid_t)." << '\n';
            std::cerr << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::set_genome(std::vector<std::vector<short>> &, std::vector<std::vector<ancestry_segment>> &, std::vector<std::vector<float>> &, std::vector<float> &, popid_t)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::set_genome(std::vector<std::vector<unsigned long>> &gstructure, float ref_allele_probability, short nploidy, std::vector<float> &gen_dist, popid_t pop_origin_id)
    {
        // Simulation of genome (snp variants)

        try
        {
            structure = gstructure;
            int ploidy = (int)nploidy;

            if (ploidy % 2 != 0)
                throw std::string("The number of provided ploidy is odd!");

            if (ref_allele_probability > 1.0f || ref_allele_probability < 0.0f)
                throw std::string("ref_allele_probability > 1.0f || ref_allele_probability < 0.0f");

            if (!gen_dist.empty())
                std::copy(gen_dist.begin(), gen_dist.end(), std::back_inserter(gen_distance));

            def_snp_table();

            for (size_t i = 0; i < (size_t)nploidy; i++) // add very first tuples pointing that the entire genome originates from the pop_origin_id 
            {
                std::vector<ancestry_segment> ancestry;
                ancestry.push_back( std::make_tuple(0,pop_origin_id) ); // only one segment because the entire genome has the same origin
                gen_ancestry.push_back(ancestry);
            }

            size_t snp_variants = 0;

            for (size_t i = 0; i < structure.size(); i++) // calculate snp variants in the genome
                snp_variants = snp_variants + std::floor(structure[i][0] / structure[i][1]);

            if (ref_allele_probability == 0.0f || ref_allele_probability == 1.0f) // if we want each individual have a first haplotype with all 1, and the second with all 0
            {
                for (size_t i = 0; i < (size_t)nploidy; i++)
                {
                    if (i % 2)
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
            std::cerr << "Exception in Genome::set_genome(std::vector<std::vector<unsigned long>> &, float, short, std::vector<float> &, popid_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Genome::set_genome(std::vector<std::vector<unsigned long>> &, float, short, std::vector<float> &, popid_t)." << '\n';
            std::cerr << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::set_genome(std::vector<std::vector<unsigned long>> &, float, short, std::vector<float> &, popid_t)." << '\n';
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
            if (markers.empty())
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
            if (markers.empty())
                throw std::string("The genome for this individual is empty!");

            if (markers.size() < 2)
                throw std::string("The genome is empty, cannot retreat the markers values!");

            size_t n_ploidy = (size_t)get_ploidy();

            std::vector<short> v;

            for (size_t i = 0; i < n_ploidy; i++)
                v.push_back((short)markers[i][locus]);

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
            std::cerr << "Reason: " << e << '\n';
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
            if (markers.empty())
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
            if (markers.empty())
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

    void Genome::get_ancestry(std::vector<std::vector<ancestry_segment>> &ancestry_haplotypes)
    {
        try
        {
            if (gen_ancestry.empty())
                throw std::string("The ancestry_haplotypes for this individual is empty!");

                ancestry_haplotypes = gen_ancestry;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::get_ancestry(std::vector<std::vector<ancestry_segment>> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_ancestry(std::vector<std::vector<ancestry_segment>> &)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::get_reproduction_gamete(std::vector<std::vector<bool>> &out_gamete, std::vector< std::vector<ancestry_segment>> &out_ancestry, size_t cross_per_chr, double mut_freq)
    {
        try
        {
            if (markers.empty())
                throw std::string("The genome for this individual is empty!");

            recombination(out_gamete, out_ancestry, cross_per_chr);
            mutation(out_gamete, mut_freq);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::get_reproduction_gamete( std::vector<std::vector<bool>> &, std::vector< std::vector<ancestry_segment>> &, size_t, double )." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_reproduction_gamete( std::vector<std::vector<bool>> &, std::vector< std::vector<ancestry_segment>> &, size_t, double )." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::mutation(std::vector<std::vector<bool>> &in_gamete, double mut_freq)
    {
        // mut_freq - probability of mutation of a single snp during meyosis
        try
        {
            if (mut_freq >= 1.0)
                throw std::string("Mutation frequency should be in the range [0.0, 1.0].");

            for (size_t i_set = 0; i_set < in_gamete.size(); i_set++) // for each specific gamete
            {
                // Sample mutation points
                Utilites u;

                std::vector<int> events = u.get_bin_rand(1, in_gamete[i_set].size(), mut_freq, false); // based on the mut_freq probability, calculate the number of mutations (modified snps) for entire gamete
                size_t events_genome = events[0];                                                      // number of mutations

                std::vector<size_t> locations = u.get_uni_rand(events_genome, (size_t)0, in_gamete[i_set].size() - 1, false); // sample of mutation locations (specific snps subject to modification)

                for (size_t i = 0; i < events_genome; i++)
                {
                    size_t point = locations[i];

                    if (!in_gamete[i_set][point])       // if 0
                        in_gamete[i_set][point] = true; // change to 1
                    else
                        in_gamete[i_set][point] = false; // change to 0
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::mutation( std::vector<std::vector<bool>> &, double )" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::mutation( std::vector<std::vector<bool>> &, double )" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::recombination(std::vector<std::vector<bool>> &out_gamete, std::vector< std::vector<ancestry_segment>> &out_ancestry, size_t n_crosses)
    {
        try
        {
            Utilites u;
            size_t i_set = 0;

            for (size_t n_sets = 0; n_sets < markers.size() / 2; n_sets++) // loop over ploidy sets of DNA pairs
            {
                size_t n_chr = structure.size();

                size_t snp_variants = 0;
                for (size_t i = 0; i < structure.size(); i++) // calculate snp variants in the genome
                    snp_variants = snp_variants + std::floor(structure[i][0] / structure[i][1]);
        
                std::vector<bool> gamete(snp_variants,false);
                std::vector<ancestry_segment> ancestry_gamete;

                for (size_t i = 0; i < n_chr; i++)
                {
                    std::vector<bool> chromatide;
                    std::vector<ancestry_segment> ancestry_chromatide;

                    short which_case = 0;
                    int sample_cell = u.get_randi(1, 100);
                    if (sample_cell <= 25)
                        which_case = 1;
                    if (sample_cell > 25 && sample_cell <= 50)
                        which_case = 2;
                    if (sample_cell > 50 && sample_cell <= 75)
                        which_case = 3;
                    if (sample_cell > 75 && sample_cell <= 100)
                        which_case = 0;

                    // Agreement:
                    // the first set of markers/variants is maternal
                    // the second set of markers/variants is paternal

                    switch (which_case)
                    {
                    case 0: // maternal with no recombination
                        get_chromatide(i, 0, i_set, chromatide, ancestry_chromatide);
                        break;
                    case 1: // maternal with recombination
                        crossover(i, 0, n_crosses, i_set, chromatide, ancestry_chromatide);
                        break;
                    case 2: // paternal with no recombination
                        get_chromatide(i, 1, i_set, chromatide, ancestry_chromatide);
                        break;
                    case 3: // paternal with recombination
                        crossover(i, 1, n_crosses, i_set, chromatide, ancestry_chromatide);
                        break;
                    default:
                        break;
                    }

                    for (size_t i2 = 0; i2 < chromatide.size(); i2++)
                        gamete[ snp_table[i][0] + i2 ] = chromatide[i2];

                    for (size_t i2 = 0; i2 < ancestry_chromatide.size(); i2++)
                        ancestry_gamete.push_back( ancestry_chromatide[i2] );
                }

                remove_redundant_segments(ancestry_gamete);
            
                out_gamete.push_back(gamete);
                out_ancestry.push_back(ancestry_gamete);

                i_set = i_set + 2; // select the next set of DNA strands, for ploidy > 2
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::recombination( std::vector<std::vector<bool>> &, std::vector< std::vector<ancestry_segment>> &, short &, size_t )." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::recombination( std::vector<std::vector<bool>> &, std::vector< std::vector<ancestry_segment>> &, short &, size_t )." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::crossover(size_t in_chr, size_t out_strand, size_t n_crossovers, size_t which_strands_set, std::vector<bool> &chromatide, std::vector<ancestry_segment> &ancestry_chromatide)
    {
        try
        {
            if ( !gen_distance.empty() )
                crossover_model_genmap(in_chr, out_strand, n_crossovers, which_strands_set, chromatide, ancestry_chromatide);
            else
                crossover_model_uniform(in_chr, out_strand, n_crossovers, which_strands_set, chromatide, ancestry_chromatide);
                // crossover_model_gamma(in_chr, out_strand, n_crossovers, which_strands_set, chromatide, ancestry_chromatide);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::crossover(size_t, size_t, size_t, size_t, std::vector<bool> &, std::vector<ancestry_segment> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::crossover(size_t, size_t, size_t, size_t, std::vector<bool> &, std::vector<ancestry_segment> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::crossover_model_gamma(size_t in_chr, size_t out_strand, size_t n_crossovers, size_t which_strands_set, std::vector<bool> &chromatide, std::vector<ancestry_segment> &ancestry_chromatide)
    {
        try
        {
            Utilites u;

            size_t first_snp = snp_table[in_chr][0];
            size_t last_snp = snp_table[in_chr][1];

            std::vector<long> cross_locations;

            long pos = (long)u.get_randi((size_t)0, hotspots[in_chr].size() - 1); // sample position in hotspots table for a very first crossover event

            if ((pos != (size_t)0) && (hotspots[in_chr][pos] == (size_t)0)) // if we catch the centromere, change the position
                pos--;

            std::vector<float> location = u.get_norm_rand(1, (float)hotspots[in_chr][pos], std_hotspots[in_chr], false); // the actual crossover position sampled from normal distr using pos as a mean

            if ((size_t)location[0] > first_snp && (size_t)location[0] < last_snp)
                cross_locations.push_back((size_t)location[0]); // very first crossover location
            else
                cross_locations.push_back(hotspots[in_chr][pos]);

            size_t chr_length = last_snp - first_snp + 1; // in number of snps (markers)
            float beta = 1.0 / ((float)n_crossovers / (float)chr_length);
            float alpha = 2.0;

            // size_t step = std_hotspots[in_chr] * step_devider;
            size_t step = std_hotspots[in_chr] * 1.0f;

            // Left move along chromosome
            size_t index = 0;
            for (size_t i = 0; i < 2 * n_crossovers; i++)
            {
                std::vector<float> left_cross_point = u.get_gamma_rand(1, alpha, beta, false);

                if ((cross_locations[index] - (long)std::round(left_cross_point[0])) <= (long)first_snp)
                    break;
                size_t position = (cross_locations[index] - (long)std::round(left_cross_point[0]) - first_snp) / step; // position in the hotspots table
                if (hotspots[in_chr][position] == (size_t)0)                                                           // this is centromere, where there is no crossover
                    continue;
                std::vector<float> locat = u.get_norm_rand(1, (float)hotspots[in_chr][position], std_hotspots[in_chr], false); // sample actual recombination position
                if ((size_t)locat[0] > first_snp && (size_t)locat[0] < last_snp)
                {
                    cross_locations.push_back((size_t)locat[0]); // crossover location
                    index++;
                }
                else
                    break;
            }

            std::vector<float> right_cross_point = u.get_gamma_rand(1, alpha, beta, false);
            size_t position0 = (cross_locations[0] + (long)std::round(right_cross_point[0]) - first_snp) / step; // position in the hotspots table

            if ((position0 != (size_t)0) && (hotspots[in_chr][position0] == (size_t)0)) // if we catch the centromere, resample
                position0--;

            if (position0 > 0 && position0 < hotspots[in_chr].size())
            {
                std::vector<float> locat = u.get_norm_rand(1, (float)hotspots[in_chr][position0], std_hotspots[in_chr], false); // sample actual recombination position

                if ((size_t)locat[0] > first_snp && (size_t)locat[0] < last_snp)
                {
                    cross_locations.push_back((size_t)locat[0]); // crossover location
                    size_t processed_crosses = cross_locations.size() - 1;

                    index = 0;
                    for (size_t i = 1; i < 2 * n_crossovers - 1; i++)
                    {
                        right_cross_point.clear();
                        right_cross_point.shrink_to_fit();

                        right_cross_point = u.get_gamma_rand(1, alpha, beta, false);

                        if ((cross_locations[processed_crosses + index] + (long)std::round(right_cross_point[0])) >= (long)last_snp)
                            break;

                        size_t position = (cross_locations[processed_crosses + index] + (long)std::round(right_cross_point[0]) - first_snp) / step; // position in the hotspots table
                        if (hotspots[in_chr][position] == (size_t)0)
                            continue;

                        if ((position > 0) && (position < hotspots[in_chr].size()))
                        {
                            // std::vector<float> locat = u.get_norm_rand(1, (float)hotspots[in_chr][position], std_hotspots[in_chr], false); // sample actual recombination position

                            if ((size_t)locat[0] > first_snp && (size_t)locat[0] < last_snp)
                            {
                                cross_locations.push_back((size_t)locat[0]); // crossover location
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

            // .........................................................
            // Staring exchnging of sections of chromosome
            // .........................................................

            // temporal containers
            std::vector<bool> chromatide_pat;
            std::vector<bool> chromatide_mat;

            for (size_t i = first_snp; i <= last_snp; i++)
            {
                chromatide_mat.push_back(markers[0 + which_strands_set][i]);
                chromatide_pat.push_back(markers[1 + which_strands_set][i]);
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
                chromatide = chromatide_mat;
            else
                chromatide = chromatide_pat;

            // ------------------------------------------
            // temporal containers
            std::vector<ancestry_segment> segments_pat;
            std::vector<ancestry_segment> segments_mat;
            
            // copy data to the temporal containers
            get_ancestry_segment(first_snp, last_snp, 0, which_strands_set, segments_mat); // get maternal ancestry list
            get_ancestry_segment(first_snp, last_snp, 1, which_strands_set, segments_pat); // get paternal ancestry list
            for (size_t i = 0; i < cross_locations.size(); i++) // add new recombination locations as segments to the lists
            {
                if ( !in_segment(segments_mat, cross_locations[i]) )
                {
                    append_segment(segments_mat, cross_locations[i]);
                    std::sort(segments_mat.begin(), segments_mat.end());
                }
                if ( !in_segment(segments_pat, cross_locations[i]) )
                {
                    append_segment(segments_pat, cross_locations[i]);
                    std::sort(segments_pat.begin(), segments_pat.end());
                }
            }
            for (size_t i = 0; i < cross_locations.size(); i++) // perform segments exchange (similar to recombinational exchange)
            {
                // temporal storages
                std::vector<ancestry_segment> segments_tmp1;
                std::vector<ancestry_segment> segments_tmp2;
    
                int mat_recomb_pos = get_segment_pos(segments_mat, cross_locations[i]);
                int pat_recomb_pos = get_segment_pos(segments_pat, cross_locations[i]);
                
                // --- this section is temporal -----
                if ( mat_recomb_pos == -1 || pat_recomb_pos == -1 )
                    std::string ("There is the negative index, this should not happen!");
                // ----------------------------------
                
                // mat segments exchange
                for (size_t l = 0; l < (size_t)mat_recomb_pos; l++) // copy all mat items before recombination point to tmp storage
                    segments_tmp1.push_back(segments_mat[l]);
                for (size_t l = (size_t)mat_recomb_pos; l < segments_mat.size(); l++) // all mat items after recombination point to tmp storage
                    segments_tmp2.push_back(segments_mat[l]);
                segments_mat.clear();
                for (size_t l = 0; l < (size_t)pat_recomb_pos; l++) // copy pat segments before recombination point to mat list
                    segments_mat.push_back(segments_pat[l]);
                for (size_t l = 0; l < segments_tmp2.size(); l++) // copy all mat segments after recombination point back to mat list
                    segments_mat.push_back(segments_tmp2[l]);
                segments_tmp2.clear();
                
                // pat segments exchange
                for (size_t l = (size_t)pat_recomb_pos; l < segments_pat.size(); l++) // copy all pat items after recombination point to tmp storage
                    segments_tmp2.push_back(segments_pat[l]);
                segments_pat.clear();
                for (size_t l = 0; l < segments_tmp1.size(); l++) // copy all mat segments before recombination points from tmp storage to pat list
                    segments_pat.push_back(segments_tmp1[l]);
                for (size_t l = 0; l < segments_tmp2.size(); l++) // copy all pat segments after recombination points from tmp storage to pat list
                    segments_pat.push_back(segments_tmp2[l]);
            }

            if (out_strand == 0)
                ancestry_chromatide = segments_mat;
            else
                ancestry_chromatide = segments_pat;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::crossover_model_gamma(size_t, size_t, size_t, size_t, std::vector<bool> &std::vector<ancestry_segment> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::crossover_model_gamma(size_t, size_t, size_t, size_t, std::vector<bool> &std::vector<ancestry_segment> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::crossover_model_uniform(size_t in_chr, size_t out_strand, size_t n_crossovers, size_t which_strands_set, std::vector<bool> &chromatide, std::vector<ancestry_segment> &ancestry_chromatide)
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
                pos = (long)u.get_randi((size_t)0, hotspots[in_chr].size() - 1);                          // sample position in hotspots table for a crossover event
                location = u.get_norm_rand(1, (float)hotspots[in_chr][pos], std_hotspots[in_chr], false); // the actual crossover position within the hotspot sampled from normal distr using pos as a mean

                if ((size_t)std::round(location[0]) > first_snp && (size_t)std::round(location[0]) < last_snp)
                    cross_locations.push_back((size_t)std::round(location[0])); // very first crossover location
                else
                    cross_locations.push_back(hotspots[in_chr][pos]);
            }

            std::sort(cross_locations.begin(), cross_locations.end());

            // .........................................................
            // Staring exchnging of sections of chromosome
            // .........................................................

            // temporal containers
            std::vector<bool> chromatide_pat;
            std::vector<bool> chromatide_mat;

            for (size_t i = first_snp; i <= last_snp; i++)
            {
                chromatide_mat.push_back(markers[0 + which_strands_set][i]);
                chromatide_pat.push_back(markers[1 + which_strands_set][i]);
            }
            
            for (size_t i = 0; i < cross_locations.size(); i++)
            {
                std::vector<bool> chromatide_tmp(1, false);

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
                chromatide = chromatide_mat;
            else
                chromatide = chromatide_pat;
            
            // ------------------------------------------
            // temporal containers
            std::vector<ancestry_segment> segments_pat;
            std::vector<ancestry_segment> segments_mat;
            
            // copy data to the temporal containers
            get_ancestry_segment(first_snp, last_snp, 0, which_strands_set, segments_mat); // get maternal ancestry list
            get_ancestry_segment(first_snp, last_snp, 1, which_strands_set, segments_pat); // get paternal ancestry list
            for (size_t i = 0; i < cross_locations.size(); i++) // add new recombination locations as segments to the lists
            {
                if ( !in_segment(segments_mat, cross_locations[i]) )
                {
                    append_segment(segments_mat, cross_locations[i]);
                    std::sort(segments_mat.begin(), segments_mat.end());
                }
                if ( !in_segment(segments_pat, cross_locations[i]) )
                {
                    append_segment(segments_pat, cross_locations[i]);
                    std::sort(segments_pat.begin(), segments_pat.end());
                }
            }
            for (size_t i = 0; i < cross_locations.size(); i++) // perform segments exchange (similar to recombinational exchange)
            {
                // temporal storages
                std::vector<ancestry_segment> segments_tmp1;
                std::vector<ancestry_segment> segments_tmp2;
    
                int mat_recomb_pos = get_segment_pos(segments_mat, cross_locations[i]);
                int pat_recomb_pos = get_segment_pos(segments_pat, cross_locations[i]);
                
                // --- this section is temporal -----
                if ( mat_recomb_pos == -1 || pat_recomb_pos == -1 )
                    std::string ("There is the negative index, this should not happen!");
                // ----------------------------------
                
                // mat segments exchange
                for (size_t l = 0; l < (size_t)mat_recomb_pos; l++) // copy all mat items before recombination point to tmp storage
                    segments_tmp1.push_back(segments_mat[l]);
                for (size_t l = (size_t)mat_recomb_pos; l < segments_mat.size(); l++) // all mat items after recombination point to tmp storage
                    segments_tmp2.push_back(segments_mat[l]);
                segments_mat.clear();
                for (size_t l = 0; l < (size_t)pat_recomb_pos; l++) // copy pat segments before recombination point to mat list
                    segments_mat.push_back(segments_pat[l]);
                for (size_t l = 0; l < segments_tmp2.size(); l++) // copy all mat segments after recombination point back to mat list
                    segments_mat.push_back(segments_tmp2[l]);
                segments_tmp2.clear();
                
                // pat segments exchange
                for (size_t l = (size_t)pat_recomb_pos; l < segments_pat.size(); l++) // copy all pat items after recombination point to tmp storage
                    segments_tmp2.push_back(segments_pat[l]);
                segments_pat.clear();
                for (size_t l = 0; l < segments_tmp1.size(); l++) // copy all mat segments before recombination points from tmp storage to pat list
                    segments_pat.push_back(segments_tmp1[l]);
                for (size_t l = 0; l < segments_tmp2.size(); l++) // copy all pat segments after recombination points from tmp storage to pat list
                    segments_pat.push_back(segments_tmp2[l]);
            }

            if (out_strand == 0)
                ancestry_chromatide = segments_mat;
            else
                ancestry_chromatide = segments_pat;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::crossover_model_uniform(size_t, size_t, size_t, size_t, std::vector<bool> &, std::vector<ancestry_segment> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::crossover_model_uniform(size_t, size_t, size_t, size_t, std::vector<bool> &, std::vector<ancestry_segment> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::crossover_model_genmap(size_t in_chr, size_t out_strand, size_t n_crossovers, size_t which_strands_set, std::vector<bool> &chromatide, std::vector<ancestry_segment> &ancestry_chromatide)
    {
        try
        {
            Utilites u;

            size_t first_snp = snp_table[in_chr][0];
            size_t last_snp = snp_table[in_chr][1];

            std::vector<genlen_t> cross_locations;

            genlen_t centromere_pos = first_snp + (float)(last_snp - first_snp + 1) * 0.6f;
            std::vector<genlen_t> recomb_list;
            recomb_list.push_back(centromere_pos); // centromere is at pos 0 here will be excluded from recombination
            
            std::vector<int> sampled_crossovers = u.get_bin_rand(1, (int)n_crossovers, 0.8, false); // sample the actual number of crossovers for this chromosome
            size_t n_samples = (size_t)sampled_crossovers[0];
            while ( n_samples > 0 )
            {
                std::vector<float> disstances;
                genlen_t sampled_pos = (genlen_t)u.get_randi( first_snp, last_snp );

                for (size_t j = 0; j < recomb_list.size(); j++) // all distances btw existing recombination pos and the new potential one
                    disstances.emplace_back( std::abs(gen_distance[sampled_pos] - gen_distance[recomb_list[j]]) );

                std::vector<float>::iterator min_result = std::min_element(disstances.begin(), disstances.end()); // get minimal distance
                double p = map_function_haldane( (double)(*min_result) ); // get probability of recombination
                if ( u.bin_rand(1, p) == 1 ) // in the case of recombination
                {
                    cross_locations.push_back( sampled_pos );
                    recomb_list.push_back( sampled_pos );
                    n_samples--;
                }
            }

            std::sort(cross_locations.begin(), cross_locations.end());

            // .........................................................
            // Staring exchnging of sections of chromosome
            // .........................................................

            // temporal containers
            std::vector<bool> chromatide_pat;
            std::vector<bool> chromatide_mat;

            // auto start = std::chrono::high_resolution_clock::now();

            for (size_t i = first_snp; i <= last_snp; i++)
            {
                chromatide_mat.push_back(markers[0 + which_strands_set][i]);
                chromatide_pat.push_back(markers[1 + which_strands_set][i]);
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
                chromatide = chromatide_mat;
            else
                chromatide = chromatide_pat;

            // ------------------------------------------
            // temporal containers
            std::vector<ancestry_segment> segments_pat;
            std::vector<ancestry_segment> segments_mat;
            
            // copy data to the temporal containers
            get_ancestry_segment(first_snp, last_snp, 0, which_strands_set, segments_mat); // get maternal ancestry list
            get_ancestry_segment(first_snp, last_snp, 1, which_strands_set, segments_pat); // get paternal ancestry list
            for (size_t i = 0; i < cross_locations.size(); i++) // add new recombination locations as segments to the lists
            {
                if ( !in_segment(segments_mat, cross_locations[i]) )
                {
                    append_segment(segments_mat, cross_locations[i]);
                    std::sort(segments_mat.begin(), segments_mat.end());
                }
                if ( !in_segment(segments_pat, cross_locations[i]) )
                {
                    append_segment(segments_pat, cross_locations[i]);
                    std::sort(segments_pat.begin(), segments_pat.end());
                }
            }
            for (size_t i = 0; i < cross_locations.size(); i++) // perform segments exchange (similar to recombinational exchange)
            {
                // temporal storages
                std::vector<ancestry_segment> segments_tmp1;
                std::vector<ancestry_segment> segments_tmp2;
    
                int mat_recomb_pos = get_segment_pos(segments_mat, cross_locations[i]);
                int pat_recomb_pos = get_segment_pos(segments_pat, cross_locations[i]);
                
                // --- this section is temporal -----
                if ( mat_recomb_pos == -1 || pat_recomb_pos == -1 )
                    std::string ("There is the negative index, this should not happen!");
                // ----------------------------------
                
                // mat segments exchange
                for (size_t l = 0; l < (size_t)mat_recomb_pos; l++) // copy all mat items before recombination point to tmp storage
                    segments_tmp1.push_back(segments_mat[l]);
                for (size_t l = (size_t)mat_recomb_pos; l < segments_mat.size(); l++) // all mat items after recombination point to tmp storage
                    segments_tmp2.push_back(segments_mat[l]);
                segments_mat.clear();
                for (size_t l = 0; l < (size_t)pat_recomb_pos; l++) // copy pat segments before recombination point to mat list
                    segments_mat.push_back(segments_pat[l]);
                for (size_t l = 0; l < segments_tmp2.size(); l++) // copy all mat segments after recombination point back to mat list
                    segments_mat.push_back(segments_tmp2[l]);
                segments_tmp2.clear();
                
                // pat segments exchange
                for (size_t l = (size_t)pat_recomb_pos; l < segments_pat.size(); l++) // copy all pat items after recombination point to tmp storage
                    segments_tmp2.push_back(segments_pat[l]);
                segments_pat.clear();
                for (size_t l = 0; l < segments_tmp1.size(); l++) // copy all mat segments before recombination points from tmp storage to pat list
                    segments_pat.push_back(segments_tmp1[l]);
                for (size_t l = 0; l < segments_tmp2.size(); l++) // copy all pat segments after recombination points from tmp storage to pat list
                    segments_pat.push_back(segments_tmp2[l]);
            }

            if (out_strand == 0)
                ancestry_chromatide = segments_mat;
            else
                ancestry_chromatide = segments_pat;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::crossover_model_genmap(size_t, size_t, size_t, size_t, std::vector<bool> &, std::vector<ancestry_segment> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::crossover_model_genmap(size_t, size_t, size_t, size_t, std::vector<bool> &, std::vector<ancestry_segment> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    double Genome::map_function_haldane(double distance)
    {
        double recomb_probablility = 0.0;

        try
        {
            recomb_probablility = 1.0 - exp( -2.0 * distance / 100.0 );
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::map_function_haldane(double)" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::map_function_haldane(double)" << '\n';
            throw;
        }

        return recomb_probablility;
    }

    //===============================================================================================================

    double Genome::map_function_kosambi(double distance)
    {
        double recomb_probablility = 0.0;

        try
        {
            recomb_probablility = ( exp( 4.0 * distance / 100.0 ) - 1.0 ) / ( exp( 4.0 * distance / 100.0 ) + 1.0 );
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::map_function_kosambi(double)" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::map_function_kosambi(double)" << '\n';
            throw;
        }

        return recomb_probablility;
    }
    
    //===============================================================================================================

    void Genome::get_chromatide(size_t which_chr, size_t which_strand, size_t which_strands_set, std::vector<bool> &chromatide, std::vector<ancestry_segment> &ancestry_chromatide)
    {
        try
        {
            size_t first_snp = snp_table[which_chr][0];
            size_t last_snp = snp_table[which_chr][1];
            
            for (size_t i = first_snp; i <= last_snp; i++)
                chromatide.push_back(markers[which_strand + which_strands_set][i]);

            get_ancestry_segment(first_snp, last_snp, which_strand, which_strands_set, ancestry_chromatide);
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
    }

    //===============================================================================================================

    void Genome::get_ancestry_segment(size_t first_snp, size_t last_snp, size_t which_strand, size_t which_strands_set, std::vector<ancestry_segment> &ancestry_chromatide)
    {
        try
        {
            // because the vector gen_ancestry supposed to haave much less elements than markers vector,
            // there is high chance at there will be no records in the search interval [first_snp, last_snp],
            // in this case in order to not return empty vector ancestry_chromatide we fill it with the closest
            // related data; than, whan this segment will be combined with other segments (upstream)
            // a redundant data will be removed

            size_t first_index = 0;
            size_t last_index = 0;

            auto it_first = std::find_if(gen_ancestry[which_strand + which_strands_set].begin(), gen_ancestry[which_strand + which_strands_set].end(),
                                   [&](const ancestry_segment& e) {return std::get<0>(e) >= first_snp && std::get<0>(e) <= last_snp;} );
            if (it_first != gen_ancestry[which_strand + which_strands_set].end())
                first_index = std::distance(std::begin(gen_ancestry[which_strand + which_strands_set]), it_first);
            else
            {
                // cannot find anything; find which pop id should be at the range [first_snp, last_snp] using previous segments; looking in revercing order
                auto it_first2 = std::find_if(gen_ancestry[which_strand + which_strands_set].rbegin(), gen_ancestry[which_strand + which_strands_set].rend(),
                                              [&](const ancestry_segment& e) {return std::get<0>(e) >= 0 && std::get<0>(e) <= last_snp;} );
                if (it_first2 != gen_ancestry[which_strand + which_strands_set].rend())
                    first_index = std::distance(it_first2, std::rend(gen_ancestry[which_strand + which_strands_set])) - 1;
                    ancestry_segment sgm = gen_ancestry[which_strand + which_strands_set][first_index];
                    std::get<0>(sgm) = first_snp;
                    ancestry_chromatide.push_back(sgm);
                    return;
            }

            auto it_last = std::find_if(gen_ancestry[which_strand + which_strands_set].rbegin(), gen_ancestry[which_strand + which_strands_set].rend(),
                                   [&](const ancestry_segment& e) {return std::get<0>(e) >= first_snp && std::get<0>(e) <= last_snp;} );
            if (it_last != gen_ancestry[which_strand + which_strands_set].rend())
                last_index = std::distance(it_last, std::rend(gen_ancestry[which_strand + which_strands_set])) - 1;

            for (size_t i = first_index; i <= last_index; i++)
                ancestry_chromatide.push_back(gen_ancestry[which_strand + which_strands_set][i]);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::get_ancestry_segment(size_t, size_t, size_t, size_t, std::vector<ancestry_segment> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_ancestry_segment(size_t, size_t, size_t, size_t, std::vector<ancestry_segment> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    bool Genome::in_segment(std::vector<ancestry_segment> &which_segment, genlen_t which_item)
    {
        bool exists = false;

        try
        {
            auto it = std::find_if( which_segment.begin(), which_segment.end(), [&](const ancestry_segment& e) { return std::get<0>(e) == which_item; } );
            if (it != which_segment.end())
                exists = true;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::in_segment(std::vector<ancestry_segment> &, genlen_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::in_segment(std::vector<ancestry_segment> &, genlen_t)" << '\n';
            throw;
        }
        return exists;
    }

    //===============================================================================================================

    void Genome::remove_redundant_segments(std::vector<ancestry_segment> &which_segment)
    {
        try
        {
            // assume we do not need to sort vector
            //std::sort(which_segment.begin(), which_segment.end());

            std::vector<size_t> remove_list;
            for (size_t i = 0; i < which_segment.size()-1; i++)
            {
                popid_t id = std::get<1>(which_segment[i]);
                if ( std::get<1>(which_segment[i+1]) == id )
                    remove_list.push_back(i+1);
            }

            if ( !remove_list.empty() ) // there is something to remove
            {
                std::vector<ancestry_segment> tmp; // temporal container
                for (size_t i = 0; i < which_segment.size(); i++) // copy to tmp those elements which should not be removed
                {
                    if ( std::find(remove_list.begin(), remove_list.end(), i) == remove_list.end() ) // index is not in the remove list
                        tmp.push_back(which_segment[i]);
                }
                which_segment.clear();
                which_segment = tmp;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::remove_redundant_segments(std::vector<ancestry_segment> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::remove_redundant_segments(std::vector<ancestry_segment> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::append_segment(std::vector<ancestry_segment> &which_segment, genlen_t which_item)
    {
        try
        {
            auto it = std::find_if( which_segment.begin(), which_segment.end(), [&](const ancestry_segment& e) { return std::get<0>(e) > which_item; } );
            
            int index = -1;
            
            if (it != which_segment.end())
                index = (int)std::distance(std::begin(which_segment), it) - 1; // we need index of previous element
            
            if (index == -1) // either the item is has max value (it == which_segment.end()) or the found value is in the position 0, which gives -1
                index = 0;
            
            ancestry_segment new_segment(which_item, std::get<1>(which_segment[(size_t)index]));
            which_segment.push_back( new_segment );
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::in_segment(std::vector<ancestry_segment> &, genlen_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::in_segment(std::vector<ancestry_segment> &, genlen_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    int Genome::get_segment_pos(std::vector<ancestry_segment> &which_segment, genlen_t which_item)
    {
        int pos = -1;
        try
        {
            auto it = std::find_if( which_segment.begin(), which_segment.end(), [&](const ancestry_segment& e) { return std::get<0>(e) == which_item; } );
            if (it != which_segment.end())
                pos = (int)std::distance(std::begin(which_segment), it);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::in_segment(std::vector<ancestry_segment> &, genlen_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::in_segment(std::vector<ancestry_segment> &, genlen_t)" << '\n';
            throw;
        }
        return pos;
    }

    //===============================================================================================================

    void Genome::show_recombination()
    {
        try
        {
            if (markers.empty())
                throw std::string("The genome for this individual is empty!");

            std::vector<std::vector<bool>> chromatide;
            std::vector<std::vector<ancestry_segment>> ancestry_chromatide;

            get_reproduction_gamete(chromatide, ancestry_chromatide, 3, 2);

            std::cout << "\n";

            std::cout << "showing gamete with " << chromatide.size() << " strands:" << "\n";

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
