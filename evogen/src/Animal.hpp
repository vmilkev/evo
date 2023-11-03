#ifndef animal_hpp__
#define animal_hpp__

#include "Organism.hpp"
#include "AnimalGenome.hpp"

namespace evogen
{
    class Animal : public Organism
    {
    public:
        Animal();
        ~Animal();

        AnimalGenome genome;
        // AnimalGenome *genotype2 = nullptr;

        void set_id(unsigned long id);
        void set_sire(unsigned long id);
        void set_dame(unsigned long id);
        void set_age(unsigned long age);
        void set_alive(bool alive);
        void set_sex(int sex);

        unsigned long get_id();
        unsigned long get_sire();
        unsigned long get_dame();
        unsigned long get_age();
        bool get_alive();
        short get_sex();

    private:
        struct Property
        {
            unsigned long id;
            unsigned long sire;
            unsigned long dame;
            unsigned long age;
            bool alive;
            short sex;
        };

        Animal::Property properties;

    protected:
    };

} // end of namespace evo

#endif // animal_hpp__