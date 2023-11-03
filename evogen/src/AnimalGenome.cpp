#include "AnimalGenome.hpp"

namespace evogen
{
    AnimalGenome::AnimalGenome()
    {
    }
    AnimalGenome::~AnimalGenome()
    {
    }

/*     void AnimalGenome::set_genome(std::vector<std::vector<short>> &snp, std::vector<std::vector<float>> &gstructure)
    {
        try
        {
            markers = snp;
            structure = gstructure;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in AnimalGenome::set_genome(std::vector<std::vector<short>> &, std::vector<std::vector<float>> &)." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (...)
        {
            std::cerr << "Exception in AnimalGenome::set_genome(std::vector<std::vector<short>> &, std::vector<std::vector<float>> &)." << '\n';
            throw;
        }
    }

    void AnimalGenome::get_genome()
    {
        try
        {
            for (const auto &e : markers)
            {
                for (auto i = 0; i < e.size(); i++)
                    std::cout << e[i] << " ";
                std::cout << "\n";
            }

            for (const auto &e : structure)
            {
                for (auto i = 0; i < e.size(); i++)
                    std::cout << e[i] << " ";
                std::cout << "\n";
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in AnimalGenome::get_genome()." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (...)
        {
            std::cerr << "Exception in AnimalGenome::get_genome()." << '\n';
            throw;
        }
    }
 */
}