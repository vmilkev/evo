#include <iostream>
#include <fstream>
#include <vector>
#include "Animal.hpp"
#include "AnimalPopulation.hpp"
#include "Utilites.hpp"

int main(int argc, char *argv[])
{
	try
	{
		/*
		evogen::Animal *cow = new evogen::Animal;

		cow->get_id();

		evogen::AnimalPopulation a, b; // simulated

		a.set_population(5, "tests/data/struct_haplotypes_pop1.dat", 0.4);
		b.set_population(5, "tests/data/struct_haplotypes_pop1.dat", 0.6);

		std::cout << "showing wild population A:"
				  << "\n";
		a.show_animals(5, 50);

		std::cout << "showing wild population B:"
				  << "\n";
		b.show_animals(5, 50);

		std::cout << "creating population from files:"
				  << "\n";
*/
		evogen::AnimalPopulation c, d; // using data from files

		c.set_population("tests/data/haplotypes_pop1.dat", "tests/data/struct_haplotypes_pop1.dat", true);				// true: haplotypes with pedigree and sex
		
		d.set_population("tests/data/haplotypes_pop_large2.dat", "tests/data/struct_haplotypes_pop_large2.dat", false); // false: only haplotypes

		std::cout << "showing population C:"
				  << "\n";

		c.show_animals(10, 60);

		std::cout << "showing population D:"
				  << "\n";

		d.show_animals(10, 50);

		std::cout << "Random numbers" << std::endl;
		evogen::Utilites u;

		std::ofstream myfile_u("tests/data/uniform_random.txt");
		std::ofstream myfile_g("tests/data/gamma_random.txt");

		std::vector<double> gamma_numbers = u.get_gamma_rand(1000, 2.0, 2.0, false);

		for (auto i = 0; i < 1000; i++)
		{
			//std::cout << "numbers in range 100: " << u.get_randi(100) << "; in range 1000: " << u.get_randi(1000) << "\n";
			myfile_u << u.get_randi(100) << "\n";
			myfile_g << gamma_numbers[i] << "\n";
		}
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	return 0;
}
