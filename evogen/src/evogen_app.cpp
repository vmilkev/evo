#include <iostream>
#include <fstream>
#include <vector>
#include "Animal.hpp"
#include "Population.hpp"
#include "Utilites.hpp"
#include "Trait.hpp"
#include "cs_matrix.hpp"

int main(int argc, char *argv[])
{
	try
	{
		/*
		evogen::Animal *cow = new evogen::Animal;

		cow->get_id();*/

		evogen::Population a, b; // simulated

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

		evogen::Population c, d; // using data from files
		std::cout << "setting c."
				  << "\n";
		c.set_population("tests/data/haplotypes_pop1.dat", "tests/data/struct_haplotypes_pop1.dat", true); // true: haplotypes with pedigree and sex
		std::cout << "setting d."
				  << "\n";
		d.set_population("tests/data/haplotypes_pop2.dat", "tests/data/struct_haplotypes_pop2.dat", false); // false: only haplotypes

		std::cout << "showing population C:"
				  << "\n";

		c.show_animals(10, 60);

		std::cout << "showing population D:"
				  << "\n";

		d.show_animals(10, 50);

		std::cout << "Random numbers" << std::endl;
		evogen::Utilites u;

		std::ofstream myfile_u("tests/data/uniform_random.txt");
		std::ofstream myfile_u2("tests/data/uniform_random2.txt");
		std::ofstream myfile_g("tests/data/gamma_random.txt");

		// std::vector<int> uni_numbers = u2.get_uni_rand(1000, 1, 100, false);
		std::vector<double> gamma_numbers = u.get_gamma_rand(1000, 2.0, 2.0, false);

		for (auto i = 0; i < 1000; i++)
		{
			// std::cout << "numbers in range 100: " << u.get_randi(100) << "; in range 1000: " << u.get_randi(1000) << "\n";
			myfile_u << u.get_randi(100) << "\n";
			myfile_u2 << u.get_randi(1, 100) << "\n";
			// std::vector<int> uni_numbers = u2.get_uni_rand(1, 1, 100, false);
			// myfile_u2 << uni_numbers[0] << "\n";
			myfile_g << gamma_numbers[i] << "\n";
		}

		std::cout << "\n";
		std::cout << "Testing trait:"
				  << "\n";
		std::cout << "\n";

		evogen::Population pop; // simulated

		pop.set_population(500, "tests/data/struct_haplotypes_pop3.dat", 0.7); // (1) population

		std::cout << "pop is ready."<<"\n";

		std::vector<double> tr_mean{40.0, 5.0, 0.5};							   // (2) trait means
		std::vector<double> qtl_prop{0.65, 0.65, 0.65, 0.65}; // (3) proportion of snps selected as qtls
		//std::vector<double> qtl_prop{1.0, 1.0, 1.0, 1.0};
		std::vector<std::vector<double>> cor_g; // (4) genomic correlations
		std::vector<std::vector<double>> cor_e; // (5) rsidual correlations
		std::vector<double> g1{1.0, 0.5, 0.7};
		std::vector<double> g2{0.5, 1.0, 0.2};
		std::vector<double> g3{0.7, 0.2, 1.0};
		std::vector<double> e1{1.0, 0.3, 0.5};
		std::vector<double> e2{0.3, 1.0, 0.4};
		std::vector<double> e3{0.5, 0.4, 1.0};
		cor_g.push_back(g1);
		cor_g.push_back(g2);
		cor_g.push_back(g3);
		cor_e.push_back(e1);
		cor_e.push_back(e2);
		cor_e.push_back(e3);

		std::vector<double> var_g{100.0, 10.0, 0.1}; // (6) genomic variances
		std::vector<double> var_e{200.0, 20.0, 0.3}; // (7) residual variances

		std::vector<double> env{0.0, 0.0, 0.0}; // (8) enviironment

		// (i) For uniform distribution, model 1
		std::vector<double> k_range_U{-1.0, 1.0}; // range of k parameter
		// (ii) For normal distribution, model 2
		std::vector<double> k_range_N{0.0, 0.5}; // mean & std
		// (iii) For gamma distribution, mdodel 3
		std::vector<double> k_range_G{0.05}; // expected value of k in loci; is '1/b' param. in the distribution

		size_t which_model = 1;

		evogen::Trait T;
		T.set_trait(pop, tr_mean, qtl_prop, cor_g, var_g, cor_e, var_e, env, which_model, k_range_U);

		T.get_observations(pop, env, "trait.dat");

	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	return 0;
}
