#ifndef animal_hpp__
#define animal_hpp__

#include "Genome.hpp"

namespace evogen
{
    class Animal
    {
    public:
        Animal();
        Animal(Animal &a, Animal &b);
        ~Animal();

        Genome genome;
        // Genome *genotype2 = nullptr;

        void set_id(unsigned long id);
        void set_sire(unsigned long id);
        void set_dame(unsigned long id);
        void set_age(int age);
        void set_alive(bool alive);
        void set_isgenotyped(bool genotyped);
        void set_sex(int sex);
        void set_phenotype(std::vector<float> &phen);
        void set_breeding_value(std::vector<float> &bv);

        unsigned long get_id();
        unsigned long get_sire();
        unsigned long get_dame();
        int get_age();
        bool get_alive();
        bool get_isgenotyped();
        bool get_active();
        int get_sex();
        std::vector<float> get_phenotype();
        std::vector<float> get_breeding_value();

        void clear();

        void operator=(const Animal& a)
        {
            genome = a.genome;
            properties.id = a.properties.id;
            properties.sire = a.properties.sire;
            properties.dame = a.properties.dame;
            properties.age = a.properties.age;
            properties.alive = a.properties.alive;
            properties.active = a.properties.active;
            properties.isgenotyped = a.properties.isgenotyped;
            properties.sex = a.properties.sex;
            properties.phenotype = a.properties.phenotype;
            properties.breeding_value = a.properties.breeding_value;
        }

    private:
        struct Property
        {
            unsigned long id;
            unsigned long sire;
            unsigned long dame;
            int age;
            bool alive;
            bool active;
            bool isgenotyped;
            int sex;

            std::vector<float> phenotype;
            std::vector<float> breeding_value;
        };
        Animal::Property properties;

        void set_active(bool active); // maybe delete this because clear() is defined!
        unsigned long asign_id();
        int asign_sex();

    protected:
    };

} // end of namespace evo

#endif // animal_hpp__