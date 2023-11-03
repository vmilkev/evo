/*
    Abstract class
*/

#ifndef organism_hpp__
#define organism_hpp__

#include "Utilites.hpp"
#include <iostream>

namespace evogen
{
    class Organism
    {
    public:
        Organism();
        unsigned long asign_id();
        short asign_sex();

    private:
    protected:
    };

} // end of namespace evogen

#endif // organism_hpp__