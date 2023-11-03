#include "Organism.hpp"

namespace evogen
{
    //===============================================================================================================

    Organism::Organism()
    {
        try
        {
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Organism::Organism()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Organism::Organism()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    unsigned long Organism::asign_id()
    {
        unsigned long id = 0;

        try
        {
            Utilites u;

            id = u.get_randi(1000000);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Organism::asign_id()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Organism::asign_id()." << '\n';
            throw;
        }

        return id;
    }

    //===============================================================================================================

    short Organism::asign_sex()
    {
        short sex = 0;

        try
        {
            Utilites u;

            unsigned long id = u.get_randi(1000000);

            if (id <= 500000)
                sex = 1;
            else
                sex = 0;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Organism::asign_sex()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Organism::asign_sex()." << '\n';
            throw;
        }

        return sex;
    }

    //===============================================================================================================

    //===============================================================================================================

    //===============================================================================================================

}