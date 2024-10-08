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

    Animal::Animal(Animal &a_sire, Animal &b_dame)
    {
        try
        {
            size_t n_crossovers = 4;
            float mutat_frequency = 0.005;

            std::vector<std::vector<bool>> gamete_a;
            a_sire.genome.get_reproduction_gamete(gamete_a, n_crossovers, mutat_frequency);

            std::vector<std::vector<bool>> gamete_b;
            b_dame.genome.get_reproduction_gamete(gamete_b, n_crossovers, mutat_frequency);

            std::vector<std::vector<unsigned long>> gstructure = a_sire.genome.get_genome_structure();

            if (gstructure != b_dame.genome.get_genome_structure())
                    throw std::string("The genome structures of parents are not the same!");

            if (gamete_a.size() != gamete_b.size())
                    throw std::string("The ploidy of parents are not the same!");

            //std::vector<std::vector<bool>> gamete( gamete_a );
            std::vector<std::vector<bool>> gamete;

            for (size_t i = 0; i < gamete_b.size(); i++)
            {
                gamete.push_back( gamete_b[i] );
                gamete.push_back( gamete_a[i] );
            }
                        
            genome.set_genome(gamete, gstructure);

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
            std::cerr << "Exception in Animal::Animal(Animal &, Animal &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::Animal(Animal &, Animal &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::Animal(Animal &, Animal &)" << '\n';
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
        try
        {
            properties.id = id;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::set_id(unsigned long)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::set_id(unsigned long)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::set_id(unsigned long)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Animal::set_sire(unsigned long sire)
    {
        try
        {
            properties.sire = sire;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::set_sire(unsigned long)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::set_sire(unsigned long)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::set_sire(unsigned long)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Animal::set_dame(unsigned long dame)
    {
        try
        {
            properties.dame = dame;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::set_sire(unsigned long)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::set_dame(unsigned long)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::set_sire(unsigned long)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Animal::set_age(int age)
    {
        try
        {
            properties.age = age;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::set_age(int)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::set_age(int)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::set_age(int)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Animal::set_alive(bool alive)
    {
        try
        {
            properties.alive = alive;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::set_alive(bool)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::set_alive(bool)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::set_alive(bool)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Animal::set_isgenotyped(bool genotyped)
    {
        try
        {
            properties.isgenotyped = genotyped;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::set_isgenotyped(bool)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::set_isgenotyped(bool)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::set_isgenotyped(bool)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Animal::set_active(bool active)
    {
        try
        {
            properties.active = active;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::set_active(bool)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::set_active(bool)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::set_active(bool)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Animal::set_sex(int sex)
    {
        try
        {
            properties.sex = sex;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::set_sex(int)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::set_sex(int)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::set_sex(int)." << '\n';
            throw;
        }
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
        unsigned long id = 0;

        try
        {
            id = properties.id;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::get_id()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::get_id()" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::get_id()." << '\n';
            throw;
        }

        return id;
    }

    //===============================================================================================================

    unsigned long Animal::get_sire()
    {
        unsigned long sire = 0;

        try
        {
            sire = properties.sire;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::get_sire()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::get_sire()" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::get_sire()." << '\n';
            throw;
        }

        return sire;
    }

    //===============================================================================================================

    unsigned long Animal::get_dame()
    {
        unsigned long dame = 0;

        try
        {
            dame = properties.dame;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::get_dame()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::get_dame()" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::get_dame()." << '\n';
            throw;
        }

        return dame;
    }

    //===============================================================================================================

    int Animal::get_age()
    {
        int age = 0;

        try
        {
            age = properties.age;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::get_age()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::get_age()" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::get_age()." << '\n';
            throw;
        }

        return age;
    }

    //===============================================================================================================

    bool Animal::get_alive()
    {
        bool alive = false;

        try
        {
            alive = properties.alive;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::get_alive()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::get_alive()" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::get_alive()." << '\n';
            throw;
        }

        return alive;
    }

    //===============================================================================================================

    bool Animal::get_isgenotyped()
    {
        try
        {
            return properties.isgenotyped;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::get_isgenotyped()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::get_isgenotyped()" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::get_isgenotyped()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    bool Animal::get_active()
    {
        bool active = false;

        try
        {
            active = properties.active;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::get_active()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::get_active()" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::get_active()." << '\n';
            throw;
        }

        return active;
    }

    //===============================================================================================================

    int Animal::get_sex()
    {
        int sex = 1;

        try
        {
            sex = properties.sex;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::get_sex()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::get_sex()" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::get_sex()." << '\n';
            throw;
        }

        return sex;
    }

    //===============================================================================================================

    std::vector<float> Animal::get_phenotype()
    {
        std::vector<float> out;

        try
        {
            out = properties.phenotype;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::get_phenotype()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::get_phenotype()" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::get_phenotype()" << '\n';
            throw;
        }

        return out;
    }

    //===============================================================================================================

    std::vector<float> Animal::get_breeding_value()
    {
        std::vector<float> out;

        try
        {
            out = properties.breeding_value;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::get_breeding_value()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Animal::get_breeding_value()" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::get_breeding_value()" << '\n';
            throw;
        }

        return out;
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

            unsigned long id = u.get_randi(1, 100);

            if (id <= 50)
                sex = 1;
            else
                sex = 0;
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