#include <iostream>
#include <fstream>
#include <vector>

#include "Animal.hpp"
#include "Population.hpp"
#include "Group.hpp"
#include "Utilites.hpp"
#include "Trait.hpp"
#include "dense_matrix.hpp"

int main(int argc, char *argv[])
{
	try
	{
		/*
		evogen::Animal *cow = new evogen::Animal;

		cow->get_id();*/

		evogen::Population a, b; // simulated

		a.set_population(5, "tests/data/struct_haplotypes_pop1.dat", 0.4, 4);
		b.set_population(5, "tests/data/struct_haplotypes_pop1.dat", 0.6, 6);

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

		/*std::cout << "Random numbers" << std::endl;
		evogen::Utilites u;

		std::ofstream myfile_u("tests/data/uniform_random.txt");
		std::ofstream myfile_u2("tests/data/uniform_random2.txt");
		std::ofstream myfile_g("tests/data/gamma_random.txt");

		// std::vector<int> uni_numbers = u2.get_uni_rand(1000, 1, 100, false);
		std::vector<float> gamma_numbers = u.get_gamma_rand(1000, 2.0, 2.0, false);

		for (auto i = 0; i < 1000; i++)
		{
			// std::cout << "numbers in range 100: " << u.get_randi(100) << "; in range 1000: " << u.get_randi(1000) << "\n";
			myfile_u << u.get_randi(100) << "\n";
			myfile_u2 << u.get_randi(1, 100) << "\n";
			// std::vector<int> uni_numbers = u2.get_uni_rand(1, 1, 100, false);
			// myfile_u2 << uni_numbers[0] << "\n";
			myfile_g << gamma_numbers[i] << "\n";
		}*/

		// ---------------------------------------------------------------
		std::cout << "\n";
		std::cout << "Testing trait:"
				  << "\n";
		std::cout << "\n";
		// ---------------------------------------------------------------

		evogen::Population pop; // simulated

		pop.set_population(20, "tests/data/struct_haplotypes_pop3.dat", 0.7, 4); // (1) population

		std::cout << "pop is ready."<<"\n";

		std::vector<float> tr_mean{40.0, 5.0, 0.5};							   // (2) trait means
		std::vector<float> qtl_prop{0.65, 0.65, 0.65, 0.65}; // (3) proportion of snps selected as qtls
		//std::vector<float> qtl_prop{1.0, 1.0, 1.0, 1.0};
		std::vector<std::vector<float>> cor_g; // (4) genomic correlations
		std::vector<std::vector<float>> cor_e; // (5) rsidual correlations
		std::vector<float> g1{1.0, 0.5, 0.7};
		std::vector<float> g2{0.5, 1.0, 0.2};
		std::vector<float> g3{0.7, 0.2, 1.0};
		std::vector<float> e1{1.0, 0.3, 0.5};
		std::vector<float> e2{0.3, 1.0, 0.4};
		std::vector<float> e3{0.5, 0.4, 1.0};
		cor_g.push_back(g1);
		cor_g.push_back(g2);
		cor_g.push_back(g3);
		cor_e.push_back(e1);
		cor_e.push_back(e2);
		cor_e.push_back(e3);

		std::vector<float> var_g{100.0, 10.0, 0.1}; // (6) genomic variances
		std::vector<float> var_e{200.0, 20.0, 0.3}; // (7) residual variances

		std::vector<float> env{0.0, 0.0, 0.0}; // (8) enviironment

		// (i) For uniform distribution, model 1
		std::vector<float> k_range_U{-1.0, 1.0}; // range of k parameter
		// (ii) For normal distribution, model 2
		std::vector<float> k_range_N{0.0, 0.5}; // mean & std
		// (iii) For gamma distribution, mdodel 3
		std::vector<float> k_range_G{0.05}; // expected value of k in loci; is '1/b' param. in the distribution

		size_t which_model = 1;

		std::cout << "Setting trait:"<<"\n";
		evogen::Trait T(pop, tr_mean, qtl_prop, cor_g, var_g, cor_e, var_e, env, which_model, k_range_U);

		evogen::Trait T2;
		T2.set_trait(pop, tr_mean, qtl_prop, cor_g, var_g, cor_e, var_e, env, which_model, k_range_U);

		std::cout << "Make observations on pop:"<<"\n";
		T.get_observations(pop, env, "T_trait_pop.dat");
		T2.get_observations(pop, env, "T2_trait.dat", "T2_genotypes.dat");

		std::cout << "Creating groups:"<<"\n";
		evogen::Group G, G2;
		std::cout<<"G size: "<<G.size()<<"\n";

		std::cout << "Adding pop to G:"<<"\n";
		G.add(pop);
		std::cout<<"G size: "<<G.size()<<", size at: "<<G.size_at(0)<<"\n";
		
		std::cout << "Make observations on G:"<<"\n";
		
		//T.get_observations(G, env, "trait_G_from_pop.dat");
		G.make_observation(T,env);
		
		std::cout << "Show pop:"<<"\n";
		pop.show_pop();

		std::cout<<"pop size: "<<pop.size()<<", capacity: "<<pop.capacity()<<"\n";

		std::cout << "Select in pop ids < 500000000 to G2:"<<"\n";
		std::cout<<"group G2 size: "<<G2.size_at(0)<<"... before selection."<<"\n";
		size_t count = 0;
		for (size_t i = 0; i < pop.size(); i++)
		{
			if (pop.id_at(i) < 500000000)
			{
				count++;
				std::cout<<"selecting ids: "<<pop.id_at(i)<<"\n";
				G2.add(pop, i);
			}
		}

		std::cout<<"Number of selected ids: "<<count<<"\n";
		std::cout<<"group G2 size: "<<G2.size_at(0)<<"\n";
		
		std::cout << "Testing G2.remove()"<<"\n";
		//G2.remove();
		//pop.show_pop();

		//pop.clear();		
		//G.clear();
		std::cout<<"group G size: "<<G.size_at(0)<<"\n";
		
		for (size_t i =0; i < pop.size(); i++)
			std::cout<<"all pop ids: "<<pop.id_at(i)<<", active: "<<pop.alive_at(i)<<"\n";

		std::cout<<"Create empty population pop2:"<<"\n";
		evogen::Population pop2;

		std::cout << "Show pop2:"<<"\n";
		pop2.show_pop();

		std::cout<<"pop2 size: "<<pop2.size()<<", capacity: "<<pop2.capacity()<<"\n";
		
		//T.get_observations(G2, env, "trait_G2.dat");
		//T.get_observations(G2, env);
		std::cout<<"Move the content of G2 into pop2, moved ids should dissapiar from pop and G2 should be empty.:"<<"\n";
		G2.move(pop2);
		
		std::cout<<"pop2 size: "<<pop2.size()<<", capacity: "<<pop2.capacity()<<"\n";
		std::cout<<"pop size: "<<pop.size()<<", capacity: "<<pop.capacity()<<"\n";
		std::cout<<"group G2 size: "<<G2.size()<<", size_at: "<<G2.size_at(0)<<"\n";

		std::cout << "Show pop:"<<"\n";
		pop.show_pop();

		std::cout<<"Tracking 1: The conteent of pop:"<<"\n";
		for (size_t i =0; i < pop.size(); i++)
			std::cout<<"pop ids: position "<<i<<", "<<pop.id_at(i)<<", active: "<<pop.alive_at(i)<<"\n";
        
        std::cout<<"pop size: "<<pop.size()<<", capacity: "<<pop.capacity()<<"\n";
        
        std::cout<<"disabling 0, 1, and 2 ..."<<'\n';
        
        pop.alive_at(0, false);
        pop.alive_at(1, false);
        pop.alive_at(2, false);

		std::cout<<"Tracking 2: The conteent of pop:"<<"\n";
		for (size_t i =0; i < pop.size(); i++)
			std::cout<<"pop ids: position "<<i<<", "<<"pop ids: "<<pop.id_at(i)<<", active: "<<pop.alive_at(i)<<"\n";

        std::cout<<"pop size: "<<pop.size()<<", capacity: "<<pop.capacity()<<"\n";

        std::cout << "Tracking 3: reshaping"<<"\n";

        pop.reshape();

		std::cout<<"Tracking 3: After reshaping:"<<"\n";
		for (size_t i =0; i < pop.size(); i++)
			std::cout<<"pop ids: position "<<i<<", "<<"pop ids: "<<pop.id_at(i)<<", active: "<<pop.alive_at(i)<<"\n";

        std::cout<<"pop size: "<<pop.size()<<", capacity: "<<pop.capacity()<<"\n";

		std::cout << "Show pop2:"<<"\n";
		pop2.show_pop();

		std::cout<<"The conteent of pop2:"<<"\n";
		for (size_t i =0; i < pop2.size(); i++)
			std::cout<<"pop2 ids: "<<pop2.id_at(i)<<", active: "<<pop2.alive_at(i)<<"\n";

		std::cout<<"Reshaping the pop:"<<"\n";
		pop.reshape();

		std::cout<<"The conteent of pop after reshaping:"<<"\n";
		std::cout << "Show pop:"<<"\n";
		pop.show_pop();

		std::cout << "Show pop2:"<<"\n";
		pop2.show_pop();

		std::cout<<"pop size: "<<pop.size()<<", capacity: "<<pop.capacity()<<"\n";
		std::cout<<"pop2 size: "<<pop2.size()<<", capacity: "<<pop2.capacity()<<"\n";
		
		std::cout << "After reshaping. Calculating trait in pop2:"<<"\n";
		T.get_observations(pop2, env, "trait_pop2.dat");
		
		std::cout << "After reshaping. Calculating trait in pop:"<<"\n";
		T.get_observations(pop, env, "trait_pop_reduced.dat");
		
		T.get_observations(pop, env);
		
		G.clear();
		G2.clear();

		pop.aging(3);
		pop2.aging(5);

		G.add(pop);
		G.add(pop2);

		G.mate();
		//G.mate(false, 2, 0.8);

		pop.show_pop();
		pop2.show_pop();

		std::cout<<"\n";
		std::cout<<"Relocating new-borns:"<<"\n";

		evogen::Group G_newborn,G3;
		evogen::Population pop_newborn, pop3;
		
		G.regroup_newborn(G_newborn);

		std::cout<<"new-borns group, size:"<<G_newborn.size()<<"\n";
		for (size_t i = 0; i < G_newborn.size(); i++)
		{
			std::cout<<"new-borns group, size at:"<<i<<", "<<G_newborn.size_at(i)<<"\n";
		}

		G_newborn.aging(4);
		G_newborn.genotype();
		G_newborn.kill();

		G_newborn.move(pop_newborn);

		std::cout<<"\n";
		std::cout<<"Showing populations:"<<"\n";
		std::cout<<"Showing pop:"<<"\n";
		pop.show_pop();
		std::cout<<"Showing pop2:"<<"\n";
		pop2.show_pop();
		std::cout<<"Showing newborn_pop:"<<"\n";
		pop_newborn.show_pop();

		G3.add(pop);
		G3.add(pop2);
		G3.add(pop_newborn);
		G3.move(pop3);

		std::cout<<"\n";
		std::cout<<"Showing populations:"<<"\n";
		std::cout<<"Showing pop:"<<"\n";
		pop.show_pop();
		std::cout<<"Showing pop2:"<<"\n";
		pop2.show_pop();
		std::cout<<"Showing newborn_pop:"<<"\n";
		pop_newborn.show_pop();
		std::cout<<"Showing pop3:"<<"\n";
		pop3.show_pop();

		T.clear();/**/
	}
	catch (const std::exception &e)
	{
		std::cerr << e.what() << '\n';
	}

	return 0;
}
