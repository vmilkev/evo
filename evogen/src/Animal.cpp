#include "Animal.hpp"

namespace evogen
{
    //===============================================================================================================

    Animal::Animal()
    {
        try
        {
            // genotype2 = new AnimalGenome();

            properties.id = asign_id();
            properties.sire = 0;
            properties.dame = 0;
            properties.age = 0;
            properties.alive = true;
            properties.sex = asign_sex();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::Animal()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::Animal()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    Animal::~Animal()
    {
    }

    //===============================================================================================================


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
        catch (...)
        {
            std::cerr << "Exception in Animal::set_id(unsigned long)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Animal::set_sire(unsigned long id)
    {
        try
        {
            properties.sire = id;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::set_sire(unsigned long)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::set_sire(unsigned long)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Animal::set_dame(unsigned long id)
    {
        try
        {
            properties.dame = id;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::set_sire(unsigned long)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::set_sire(unsigned long)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Animal::set_age(unsigned long age)
    {
        try
        {
            properties.dame = age;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::set_age(unsigned long)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::set_age(unsigned long)." << '\n';
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
        catch (...)
        {
            std::cerr << "Exception in Animal::set_alive(bool)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Animal::set_sex(int sex)
    {
        try
        {
            properties.sex = (short)sex;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Animal::set_sex(int)." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Animal::set_sex(int)." << '\n';
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
        catch (...)
        {
            std::cerr << "Exception in Animal::get_dame()." << '\n';
            throw;
        }

        return dame;
    }

    //===============================================================================================================

    unsigned long Animal::get_age()
    {
        unsigned long age = 0;

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
        catch (...)
        {
            std::cerr << "Exception in Animal::get_alive()." << '\n';
            throw;
        }

        return alive;
    }

    //===============================================================================================================

    short Animal::get_sex()
    {
        short sex = 1;

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
        catch (...)
        {
            std::cerr << "Exception in Animal::get_sex()." << '\n';
            throw;
        }

        return sex;
    }

    //===============================================================================================================

}