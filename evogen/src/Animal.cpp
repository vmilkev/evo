#include "Animal.hpp"

namespace evogen
{
    //===============================================================================================================
    Animal::Animal()
    {
        try
        {
            // genotype2 = new Genome();
            properties.id = asign_id();
            properties.sire = 0;
            properties.dame = 0;
            properties.age = 0;
            properties.alive = true;
            properties.isgenotyped = false;
            properties.active = true;
            properties.sex = asign_sex();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::Animal()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::Animal()" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::Animal()." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    Animal::Animal(Animal &a_sire, Animal &b_dame, double mutation_freq, size_t num_crossovers)
    {
        try
        {
            size_t n_crossovers = num_crossovers;
            double mutat_frequency = mutation_freq;

            std::vector<std::vector<bool>> gamete_a;
            std::vector<std::vector<ancestry_segment>> ancestry_a;
            a_sire.genome.get_reproduction_gamete(gamete_a, ancestry_a, n_crossovers, mutat_frequency);

            std::vector<std::vector<bool>> gamete_b;
            std::vector<std::vector<ancestry_segment>> ancestry_b;
            b_dame.genome.get_reproduction_gamete(gamete_b, ancestry_b, n_crossovers, mutat_frequency);

            std::vector<std::vector<unsigned long>> gstructure = a_sire.genome.get_genome_structure();

            if (gstructure != b_dame.genome.get_genome_structure())
                    throw std::string("The genome structures of parents are not the same!");

            if (gamete_a.size() != gamete_b.size())
                    throw std::string("The ploidy of parents are not the same!");

            if (ancestry_a.size() != ancestry_b.size())
                    throw std::string("The ploidy of ancestry vectors are not the same!");
//std::cout<<"ancestry_a.size() = "<<ancestry_a.size()<<", gamete_a.size() = "<<gamete_a.size()<<'\n';
            if (ancestry_a.size() != gamete_a.size())
                    throw std::string("The ploidy of markers and ancestry vectors are not the same!");
                    
            if (gamete_b.size() != ancestry_b.size())
                    throw std::string("The ploidy of markers and ancestry vectors are not the same!");

            std::vector<float> g_dist = a_sire.genome.get_gen_distance();

            std::vector<std::vector<bool>> gamete;
            std::vector<std::vector<ancestry_segment>> ancestry;

            for (size_t i = 0; i < gamete_b.size(); i++) // in case of poly-ploidy
            {
                // adding mat and pat gametes by pairs
                gamete.push_back( gamete_b[i] ); // first, maternal gamete
                gamete.push_back( gamete_a[i] ); // second, paternal gamete
                ancestry.push_back( ancestry_b[i] ); // first, maternal gamete
                ancestry.push_back( ancestry_a[i] ); // second, paternal gamete
                // ... than the next pair of gametes in case of ploidy > 2
            }
            
            genome.set_genome(gamete, ancestry, gstructure, g_dist); // here we do not initialize gene ancestry, the gen_ancestry should already be non-empty

            properties.id = asign_id();
            properties.sire = a_sire.get_id();
            properties.dame = b_dame.get_id();
            properties.age = 0;
            properties.alive = true;
            properties.isgenotyped = false;
            properties.active = true;
            properties.sex = asign_sex();

            //std::cout<<"NEW Animal "<<properties.id<<" sire "<<a_sire.properties.id<<", dame "<<b_dame.properties.id<<'\n';
            
            //std::vector<std::vector<bool>> haplotypes;
            //std::vector<std::vector<unsigned long>> gstructure2;

            /*a_sire.genome.get_genome(haplotypes, gstructure2);

            std::cout<<"sire "<<properties.sire<<" genome: "<<'\n';
            for (size_t i = 0; i < haplotypes.size(); i++)
            {
                for (size_t j = 0; j < haplotypes[i].size(); j++)
                    std::cout<< haplotypes[i][j] << " ";
                std::cout << '\n';
            }

            haplotypes.clear(); haplotypes.shrink_to_fit(); gstructure2.clear(); gstructure2.shrink_to_fit();

            b_dame.genome.get_genome(haplotypes, gstructure2);

            std::cout<<"dame "<<properties.dame<<" genome: "<<'\n';
            for (size_t i = 0; i < haplotypes.size(); i++)
            {
                for (size_t j = 0; j < haplotypes[i].size(); j++)
                    std::cout<< haplotypes[i][j] << " ";
                std::cout << '\n';
            }

            haplotypes.clear(); haplotypes.shrink_to_fit(); gstructure2.clear(); gstructure2.shrink_to_fit();*/

            //std::cout<<"offspring "<<properties.id<<" gamete (size & dame): "<<'\n';

            /*for (size_t i = 0; i < gamete_a.size(); i++)
            {
                for (size_t j = 0; j < gamete_a[i].size(); j++)
                    std::cout<< gamete_a[i][j] << " ";
                std::cout << '\n';
            }            
            //std::cout<<"dame gamete: "<<'\n';
            for (size_t i = 0; i < gamete_b.size(); i++)
            {
                for (size_t j = 0; j < gamete_b[i].size(); j++)
                    std::cout<< gamete_b[i][j] << " ";
                std::cout << '\n';
            }*/
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::Animal(Animal &, Animal &, double, size_t, popid_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::Animal(Animal &, Animal &, double, size_t, popid_t)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::Animal(Animal &, Animal &, double, size_t, popid_t)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    Animal::~Animal()
    {
        clear();
    }
    //===============================================================================================================
    void Animal::set_id(unsigned long id)
    {
        properties.id = id;
    }
    //===============================================================================================================
    void Animal::set_sire(unsigned long sire)
    {
        properties.sire = sire;
    }
    //===============================================================================================================
    void Animal::set_dame(unsigned long dame)
    {
        properties.dame = dame;
    }
    //===============================================================================================================
    void Animal::set_age(int age)
    {
        properties.age = age;
    }
    //===============================================================================================================
    void Animal::set_alive(bool alive)
    {
        properties.alive = alive;
    }
    //===============================================================================================================
    void Animal::set_isgenotyped(bool genotyped)
    {
        properties.isgenotyped = genotyped;
    }
    //===============================================================================================================
    void Animal::set_active(bool active)
    {
        properties.active = active;
    }
    //===============================================================================================================
    void Animal::set_sex(int sex)
    {
        properties.sex = sex;
    }
    //===============================================================================================================
    void Animal::set_phenotype(std::vector<float> &phen)
    {
        try
        {
            // clear the old data
            properties.phenotype.clear();
            properties.phenotype.shrink_to_fit();
            // write the new records
            properties.phenotype = phen;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::set_phenotype(std::vector<float> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::set_phenotype(std::vector<float> &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::set_phenotype(std::vector<float> &)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Animal::set_breeding_value(std::vector<float> &bv)
    {
        try
        {
            // clear the old data
            properties.breeding_value.clear();
            properties.breeding_value.shrink_to_fit();
            // write the new records
            properties.breeding_value = bv;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::set_breeding_value(std::vector<float> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::set_breeding_value(std::vector<float> &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::set_breeding_value(std::vector<float> &)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    unsigned long Animal::get_id()
    {
        return properties.id;
    }
    //===============================================================================================================
    unsigned long Animal::get_sire()
    {
        return properties.sire;
    }
    //===============================================================================================================
    unsigned long Animal::get_dame()
    {
        return properties.dame;
    }
    //===============================================================================================================
    int Animal::get_age()
    {
        return properties.age;
    }
    //===============================================================================================================
    bool Animal::get_alive()
    {
        return properties.alive;
    }
    //===============================================================================================================
    bool Animal::get_isgenotyped()
    {
        return properties.isgenotyped;
    }
    //===============================================================================================================
    bool Animal::get_active()
    {
        return properties.active;
    }
    //===============================================================================================================
    int Animal::get_sex()
    {
        return properties.sex;
    }
    //===============================================================================================================
    std::vector<float> Animal::get_phenotype()
    {
        return properties.phenotype;
    }
    //===============================================================================================================
    std::vector<float> Animal::get_breeding_value()
    {
        return properties.breeding_value;
    }
    //===============================================================================================================
    void Animal::clear()
    {
        try
        {
            properties.active = false;
            genome.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::clear()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::clear()" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::clear()." << '\n';
            throw;
        }

    }
    //===============================================================================================================
    unsigned long Animal::asign_id()
    {
        unsigned long id = 0;
        try
        {
            Utilites u;
            id = u.get_randi( 1, 1000000000 );
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::asign_id()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::asign_id()." << '\n';
            throw;
        }
        return id;
    }
    //===============================================================================================================
    int Animal::asign_sex()
    {
        int sex = 0;
        try
        {
            Utilites u;
            sex = u.bin_rand(1,0.5);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::asign_sex()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::asign_sex()." << '\n';
            throw;
        }
        return sex;
    }
    //===============================================================================================================
}