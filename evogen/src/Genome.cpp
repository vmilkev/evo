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

    void Genome::def_snp_table()
    {
        try
        {
            size_t n_chr = structure.size();

            ulong previous = 0;

            for (size_t i = 0; i < n_chr; i++)
            {
                std::vector<unsigned long> range;

                range.push_back(previous);
                range.push_back(previous + structure[i][0] / structure[i][1] - 1);
                previous = previous + structure[i][0] / structure[i][1];

                snp_table.push_back(range);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::def_snp_table()." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::def_snp_table()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::set_genome(std::vector<std::vector<bool>> &snp, std::vector<std::vector<unsigned long>> &gstructure)
    {
        /* uses prepared SNP variants either from file or from reproduction gamete */

        try
        {
            markers = snp;
            structure = gstructure;
            int ploidy = snp.size();

            if (ploidy % 2 != 0)
            {
                throw std::string("The number of ploidy is odd!");
            }

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
        /* Simulation of genome (snp variants) */

        try
        {
            structure = gstructure;
            int ploidy = (int)nploidy;

            if (ploidy % 2 != 0)
            {
                throw std::string("The number of provided ploidy is odd!");
            }

            def_snp_table();

            size_t snp_variants = 0;

            // calculate snp variants in the genome
            for (size_t i = 0; i < structure.size(); i++)
            {
                snp_variants = snp_variants + std::floor(structure[i][0] / structure[i][1]);
            }

            for (size_t i = 0; i < nploidy; i++)
            {
                std::vector<bool> variants;
                for (size_t j = 0; j < snp_variants; j++)
                {
                    variants.push_back(asign_snp_variant(ref_allele_probability));
                }
                markers.push_back(variants);
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
        short snp_variant = 0;

        try
        {
            Utilites u;

            int id = u.get_randi(1, 100);

            if (id <= ref_allele_probability * 100)
                snp_variant = 1;
            else
                snp_variant = 0;
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

    void Genome::show_genome(int max_variants)
    {
        try
        {
            size_t n_chr = structure.size();

            std::cout << "n. chromosomes = " << n_chr << "\n";

            for (const auto &e : markers)
            {
                size_t max_markers = std::min(max_variants, (int)e.size());

                for (auto i = 0; i < max_markers; i++)
                    std::cout << e[i] << " ";
                std::cout << "\n";
            }

            for (const auto &e : structure)
            {
                for (auto i = 0; i < e.size(); i++)
                    std::cout << e[i] << " ";
                std::cout << "\n";
            }

            std::cout << "snp table:"
                      << "\n";

            for (const auto &e : snp_table)
            {
                for (auto i = 0; i < e.size(); i++)
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

    void Genome::get_reproduction_gamete(size_t cross_per_chr, size_t mut_per_genome, std::vector<bool> &out_gamete, short &out_sex_chr_id)
    {
        /* out_sex_chr_id: for each produced gamete indicates where does sex chromosome (the last one) comes from in the gamete: 0 - paternal, 1 - maternal. */

        try
        {
            recombination(out_gamete, out_sex_chr_id, cross_per_chr);
            mutation(out_gamete, mut_per_genome);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::get_reproduction_gamete( size_t, size_t, std::vector<bool> &, short & )." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_reproduction_gamete( size_t, size_t, std::vector<bool> &, short & )." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::mutation(std::vector<bool> &in_gamete, size_t events_genome)
    {
        try
        {
            // Sample mutation points
            Utilites u;
            std::vector<size_t> locations = u.get_uni_rand(events_genome, (size_t)0, in_gamete.size() - 1, false);

            for (size_t i = 0; i < events_genome; i++)
            {
                size_t point = locations[i];
                std::cout << "Mutation point: " << point << "; current val: " << in_gamete[point] << "; obtained val: ";
                if (in_gamete[point] == 0)
                    in_gamete[point] = 1;
                else
                    in_gamete[point] = 0;
                std::cout << in_gamete[point] << "\n";
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::mutation( std::vector<bool> &, size_t )." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::mutation( std::vector<bool> &, size_t )." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Genome::recombination(std::vector<bool> &out_gamete, short &out_sex_chr_id, size_t n_crosses)
    {
        try
        {
            Utilites u;

            size_t n_chr = structure.size();

            //......................................................

            /* Random table of 2 cols to select between paternal and maternal chromosomes
               which alogned (distributed) between two cells after crossing over. Two produced cells
               will consist of mixture of paternal and maternal chromosomes, among which
               there are one - original, and one - crossed.
            */

            std::vector<std::vector<short>> is_paternal;

            for (size_t i = 0; i < n_chr; i++)
            {
                // Utilites u;

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

            /* Random table of 4 cols to select between original and crossed chromosoms
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

            /* among 4 final gamete cells define (sample) which cell is selected for reproduction;
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

            /* among 2 pre-gamete cells define (sample) which cell line
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
            /* Here we know where does sex chromosome will come: from sire or dame */

            out_sex_chr_id = is_paternal[n_chr - 1][which_pat];

            //......................................................

            for (size_t i = 0; i < n_chr; i++)
            {
                std::vector<bool> chromatide;
                short apply_cross = is_crossed[i][which_cell];  // defines the sequence of original/crossed chromosomes
                short which_strand = is_paternal[i][which_pat]; // defines the sequence of chromosome origin: paternal/maternal

                /* the values of which_cell & which_pat are constant for a specific gamete, and expected to be different for every new generated gamete. */

                std::cout << "chr.no: " << i << "; apply_cross & which strand: " << apply_cross << ", " << which_strand << "; which_pat & which_cell: " << which_pat << ", " << which_cell << "\n";

                if (apply_cross)
                    chromatide = crossover(i, which_strand, n_crosses);
                else
                    chromatide = get_chromatide(i, which_strand);

                for (size_t i = 0; i < chromatide.size(); i++)
                    out_gamete.push_back(chromatide[i]);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::recombination( std::vector<bool> &, short &, size_t )." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::recombination( std::vector<bool> &out_gamete, short &, size_t )." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<bool> Genome::crossover(size_t in_chr, size_t out_strand, size_t n_crossovers)
    {
        std::vector<bool> chromatide;

        try
        {
            size_t first_snp = snp_table[in_chr][0];
            size_t last_snp = snp_table[in_chr][1];

            std::vector<long> cross_locations;

            Utilites u;

            size_t chr_length = last_snp - first_snp + 1;

            cross_locations.push_back((long)u.get_randi(first_snp, last_snp)); /* sample the very first crossover location */

            double beta = 1.0 / ((double)n_crossovers / (double)chr_length);
            double alpha = 2.0;

            // std::cout<<"first & last snp: "<<first_snp<<", "<<last_snp<<"; beta: "<<beta <<"; chr_length: "<<chr_length<<"; n_crosses: "<<n_crossovers<<"\n";

            std::vector<double> left_cross_points = u.get_gamma_rand(2 * n_crossovers, alpha, beta, false);
            std::vector<double> right_cross_points = u.get_gamma_rand(2 * n_crossovers, alpha, beta, false);

            // std::cout << "chromosome: " << in_chr << "; first location: " << cross_locations[0] << "\n";

            // Left move along chromosome
            for (size_t i = 0; i < 2 * n_crossovers; i++)
            {
                long next_location = cross_locations[i] - (long)std::round(left_cross_points[i]);

                if (next_location > (long)first_snp)
                {
                    cross_locations.push_back(next_location);
                    // std::cout << "left previous point: " << cross_locations[i] << ", next point: " << next_location << ", shift: " << std::round(left_cross_points[i]) << "\n";
                }
                else
                    break;
            }

            // right move along chromosome
            long next_location = cross_locations[0] + (long)std::round(right_cross_points[0]);
            size_t processed_crosses = cross_locations.size();

            if (next_location < (long)last_snp)
            {
                cross_locations.push_back(next_location);
                processed_crosses = cross_locations.size();
                // std::cout << "right previous point: " << cross_locations[0] << ", next point: " << next_location << ", shift: " << std::round(right_cross_points[0]) << "\n";

                for (size_t i = 0; i < 2 * n_crossovers; i++)
                {
                    long next_location = cross_locations[i + processed_crosses - 1] + (long)std::round(right_cross_points[i + 1]);

                    if (next_location < (long)last_snp)
                    {
                        cross_locations.push_back(next_location);
                        // std::cout << "right previous point: " << cross_locations[i + processed_crosses - 1] << ", next point: " << next_location << ", shift: " << std::round(right_cross_points[i+1]) << "\n";
                    }
                    else
                        break;
                }
            }

            std::sort(cross_locations.begin(), cross_locations.end());

            std::cout << "\n";
            std::cout << "cross locations in chromosome " << in_chr << ": ";
            for (size_t i = 0; i < cross_locations.size(); i++)
            {
                std::cout << cross_locations[i] << " ";
            }
            std::cout << "\n";

            // .........................................................
            // Staring exchnging of sections of chromosome
            // .........................................................

            // temporal containers
            std::vector<bool> chromatide_pat;
            std::vector<bool> chromatide_mat;

            for (size_t i = first_snp; i <= last_snp; i++)
            {
                chromatide_pat.push_back(markers[0][i]);
                chromatide_mat.push_back(markers[1][i]);
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
                chromatide = chromatide_pat;
            else
                chromatide = chromatide_mat;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::crossover(size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::crossover(size_t, size_t)." << '\n';
            throw;
        }

        return chromatide;
    }

    //===============================================================================================================

    std::vector<bool> Genome::get_chromatide(size_t which_chr, size_t which_strand)
    {
        std::vector<bool> chromatide;

        try
        {
            size_t first_snp = snp_table[which_chr][0];
            size_t last_snp = snp_table[which_chr][1];
            for (size_t i = first_snp; i <= last_snp; i++)
            {
                chromatide.push_back(markers[which_strand][i]);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Genome::get_chromatide(size_t, size_t)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Genome::get_chromatide(size_t, size_t)." << '\n';
            throw;
        }

        return chromatide;
    }

    //===============================================================================================================

    void Genome::show_recombination()
    {

        try
        {
            std::vector<bool> chromatide;
            short which_sex_chr = 0;

            get_reproduction_gamete(3, 2, chromatide, which_sex_chr);

            std::cout << "\n";

            std::cout << "showing gamete:"
                      << "\n";

            size_t snps = std::min((size_t)100, chromatide.size());

            for (size_t i = 0; i < snps; i++)
            {
                std::cout << chromatide[i] << " ";
            }
            std::cout << "\n";
            std::cout << "sex chromosome from: " << which_sex_chr
                      << "\n";
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

}
