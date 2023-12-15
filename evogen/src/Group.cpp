#include "Group.hpp"

namespace evogen
{
    //===============================================================================================================

    Group::Group()
    {
        //
    }

    //===============================================================================================================

    Group::~Group()
    {
        //
    }

    //===============================================================================================================

    Population & Group::get_population( size_t which_pop )
    {
        try
        {
            return pop_list[which_pop].get();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in &Group::get_population( size_t )" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in &Group::get_population( size_t )" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in &Group::get_population( size_t )" << '\n';
            throw;
        }

    }

    //===============================================================================================================

    std::vector<size_t> Group::get_individuals( size_t which_pop )
    {
        try
        {
            return individuals_list[which_pop];
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::get_individuals( size_t )" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::get_individuals( size_t )" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::get_individuals( size_t )" << '\n';
            throw;
        }

    }

    //===============================================================================================================

    void Group::clear()
    {
        try
        {
            pop_list.clear();
            pop_list.shrink_to_fit();

            individuals_list.clear();
            individuals_list.shrink_to_fit();

            pop_list_newborn.clear();
            pop_list_newborn.shrink_to_fit();

            individuals_list_newborn.clear();
            individuals_list_newborn.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::clear()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::clear()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::clear()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    size_t Group::size()
    {
        size_t out = 0;

        try
        {
            out = pop_list.size();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::size()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::size()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::size()" << '\n';
            throw;
        }

        return out;
    }

    //===============================================================================================================

    size_t Group::size_at(size_t at)
    {
        size_t out = 0;

        try
        {
            if ( at < pop_list.size() )
                return individuals_list[at].size();
            else
                return 0;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::size_at(size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::size_at(size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::size_at(size_t)" << '\n';
            throw;
        }

        return out;
    }

    //===============================================================================================================

    void Group::remove()
    {
        try
        {
            if ( individuals_list.size() > 0 && pop_list.size() > 0 )
            {
                for (size_t i = 0; i < pop_list.size(); i++)
                {
                    if ( pop_list[i].get().size() != 0 )
                    {
                        if ( individuals_list[i].size() == 0 )
                            throw std::string("Cannot remove individuals from the empty list in the Group!");
                        
                        for (size_t j = 0; j < individuals_list[i].size(); j++)
                        {
                            // remember: individuals_list[i][j] points to the position (index)
                            // in the active_individuals list
                            //std::cout<<"removed id: "<<pop_list[i].get().id_at( individuals_list[i][j] - j )<<"\n";
                            pop_list[i].get().remove_at( individuals_list[i][j] - j );
                        }
                    }
                }
                clear();
            }
            else
                throw std::string("Cannot remove individuals from the empty Group!");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::remove()." << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::remove()." << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::remove()." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::add(Population &pop)
    {
        try
        {
            std::vector<std::reference_wrapper<Population>>::iterator position = std::find( pop_list.begin(), pop_list.end(), std::ref(pop) );
            
            if (position == pop_list.end()) // if the population is not in the host group
            {
                pop_list.push_back(pop);
                std::vector<size_t> list;
                for (size_t i = 0; i < pop.active_individuals.size(); i++ )
                {
                    list.push_back(i);
                }                
                individuals_list.push_back(list);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::add(Population &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::add(Population &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::add(Population &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::add(Group &grp)
    {
        try
        {
            for (size_t i = 0; i < grp.size(); i++)
            {
                std::vector<std::reference_wrapper<Population>>::iterator position = std::find( pop_list.begin(), pop_list.end(), std::ref(grp.pop_list[i]) );
                if (position != pop_list.end()) // if the population was actually found
                {
                    // the population is already in the list,
                    // therefore just adding the id to correct position
                    size_t pos = position - pop_list.begin();

                    // avoid duplication of ids, do not add
                    // the ones which are already in the host group;                        
                    // (1) sort to allow duplicates appiars consecutively,
                    // (2) unique, to remove all duplicates,
                    // (3) adjust the container size.
                    std::vector<size_t> list = grp.individuals_list[i];
                    std::sort(list.begin(), list.end());
                    std::vector<size_t>::iterator last = std::unique(list.begin(), list.end());
                    list.erase(last, list.end());

                    for (size_t j = 0; j < list.size(); j++)
                        individuals_list[pos].push_back( list[j] );
                }
                else
                {
                    // we do not need to check for the unique ids here
                    // because we assume all ids coming from the same population
                    // are unique by default
                    pop_list.push_back(grp.pop_list[i]);
                    std::vector<size_t> list;
                    for (size_t j = 0; j < grp.size_at(i); j++)
                    {
                        list.push_back(grp.individuals_list[i][j]);
                    }
                    individuals_list.push_back(list);
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::add(Group &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::add(Group &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::add(Group &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::add(Population &pop, size_t which_one)
    {
        try
        {
            if ( which_one < 0 || which_one > pop.active_individuals.size()-1 )
                throw std::string("Adding ID is illegal!");

            std::vector<std::reference_wrapper<Population>>::iterator position = std::find( pop_list.begin(), pop_list.end(), std::ref(pop) );
            
            if (position != pop_list.end()) // if the population is already in the group
            {
                // the population is already in the list,
                // therefore just adding the id to correct position
                size_t pos = position - pop_list.begin();
                
                // check for possible duplicates
                std::vector<size_t>::iterator list_position = std::find( individuals_list[pos].begin(), individuals_list[pos].end(), which_one );
                if (list_position == individuals_list[pos].end()) // if the element was not found
                    individuals_list[pos].push_back( which_one );
            }
            else
            {
                // we do not need to check for the unique ids here
                // because we assume all ids coming from the same population
                // are unique by default
                pop_list.push_back(pop);
                std::vector<size_t> list;
                list.push_back(which_one);
                individuals_list.push_back(list);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::add(Group &, size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::add(Group &, size_t)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::add(Group &, size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::add_newborn(Population &pop, size_t which_one)
    {
        try
        {
            if ( which_one < 0 || which_one > pop.active_individuals.size()-1 )
                throw std::string("Adding ID is illegal!");

            std::vector<std::reference_wrapper<Population>>::iterator position = std::find( pop_list_newborn.begin(), pop_list_newborn.end(), std::ref(pop) );
            
            if (position != pop_list_newborn.end()) // if the population is already in the group
            {
                // the population is already in the list,
                // therefore just adding the id to correct position
                size_t pos = position - pop_list_newborn.begin();
                
                // check for possible duplicates
                std::vector<size_t>::iterator list_position = std::find( individuals_list_newborn[pos].begin(), individuals_list_newborn[pos].end(), which_one );
                if (list_position == individuals_list_newborn[pos].end()) // if the element was not found
                    individuals_list_newborn[pos].push_back( which_one );
            }
            else
            {
                // we do not need to check for the unique ids here
                // because we assume all ids coming from the same population
                // are unique by default
                pop_list_newborn.push_back(pop);
                std::vector<size_t> list;
                list.push_back(which_one);
                individuals_list_newborn.push_back(list);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::add_newborn(Group &, size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::add_newborn(Group &, size_t)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::add_newborn(Group &, size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::move(Population &pop)
    {
        try
        {
            for (size_t i = 0; i < pop_list.size(); i++)
            {
                if (pop != pop_list[i])
                {
                    for (size_t j = 0; j < individuals_list[i].size(); j++)
                    {
                        // remember: individuals_list[i][j] points to the position (index)
                        // in the active_individuals list inside population
                        pop.add( pop_list[i].get().individuals[ pop_list[i].get().active_individuals[ individuals_list[i][j] - j ] ] );
                        pop_list[i].get().remove_at( individuals_list[i][j] - j );
                    }
                }
            }
            clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::move(Population &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::move(Population &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::move(Population &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::mate(bool sexual_reproduction, int max_offspring, double success_rate)
    {
        try
        {
            if (sexual_reproduction)
            {
                // association lists
                std::vector<size_t> female_list;
                std::vector<size_t> male_list;
                std::vector<size_t> female_pops;
                std::vector<size_t> male_pops;

                // initialize lists
                for (size_t i = 0; i < pop_list.size(); i++)
                {
                    for (size_t j = 0; j < individuals_list[i].size(); j++)
                    {
                        if ( pop_list[i].get().sex_at( individuals_list[i][j] ) == 1 )
                        {
                            male_list.push_back(j);
                            male_pops.push_back(i);
                        }
                        else
                        {
                            female_list.push_back(j);
                            female_pops.push_back(i);
                        }
                    }
                }

                Utilites u;

                size_t n_male = male_list.size()-1;

                // consecutively selecting females to mate them onece with available mails (randomly sampled)
                for (size_t i = 0; i < female_list.size(); i++)
                {
                    std::vector<int> n_offsprings = u.get_bin_rand(1,max_offspring,success_rate, false); // sample number of offsprings for this female
                    std::vector<unsigned long> which_mail = u.get_uni_rand(1,(unsigned long)0,(unsigned long)n_male,false); // sample male to mate

                    size_t i_pop = female_pops[i];
                    size_t i_female = female_list[i];
                    size_t i_anim = pop_list[i_pop].get().active_individuals[i_female];
                    //Animal a_dame =  pop_list[i_pop].get().individuals[i_anim]; // Female individual

                    // select the sampled male
                    size_t i_pop2 = male_pops[ which_mail[0] ];
                    size_t i_male = male_list[ which_mail[0] ];
                    size_t i_anim2 = pop_list[i_pop2].get().active_individuals[i_male];
                    //Animal b_sire =  pop_list[i_pop2].get().individuals[i_anim2]; // Male individual

                    for (size_t j = 0; j < (size_t)n_offsprings[0]; j++)
                    {
                        Animal c( pop_list[i_pop2].get().individuals[i_anim2], pop_list[i_pop].get().individuals[i_anim] ); // constructing new-born individual
                        pop_list[i_pop].get().add(c); // add to female's population
                        add(pop_list[i_pop],pop_list[i_pop].get().size()-1); // add to the current group
                        add_newborn(pop_list[i_pop],pop_list[i_pop].get().size()-1); // add to the current group, new-born lists
                    }
                }
            }
            else
            {
                std::vector<size_t> all_list;
                std::vector<size_t> all_pops;

                // initialize lists
                for (size_t i = 0; i < pop_list.size(); i++)
                {
                    for (size_t j = 0; j < individuals_list[i].size(); j++)
                    {
                        all_list.push_back(j);
                        all_pops.push_back(i);
                    }
                }

                Utilites u;

                size_t n_individuals = all_list.size()-1;

                for (size_t i = 0; i < all_list.size()/2; i++)
                {
                    std::vector<int> n_offsprings = u.get_bin_rand(1,max_offspring,success_rate, false); // sample number of offsprings for the specific mating
                    std::vector<unsigned long> which_indiv = u.get_uni_rand(2,(unsigned long)0,(unsigned long)n_individuals,false); // sample two individuals to mate

                    size_t i_pop = all_pops[ which_indiv[0] ];
                    size_t i_female = all_list[ which_indiv[0] ];
                    size_t i_anim = pop_list[i_pop].get().active_individuals[i_female];
                    //Animal a_dame =  pop_list[i_pop].get().individuals[i_anim]; // Individual 1

                    // select the sampled male
                    size_t i_pop2 = all_pops[ which_indiv[1] ];
                    size_t i_male = all_list[ which_indiv[1] ];
                    size_t i_anim2 = pop_list[i_pop2].get().active_individuals[i_male];
                    //Animal b_sire =  pop_list[i_pop2].get().individuals[i_anim2]; // Male individual

                    for (size_t j = 0; j < (size_t)n_offsprings[0]; j++)
                    {
                        Animal c( pop_list[i_pop2].get().individuals[i_anim2], pop_list[i_pop].get().individuals[i_anim] ); // constructing new-born individual
                        pop_list[i_pop].get().add(c); // add to female's population
                        add(pop_list[i_pop],pop_list[i_pop].get().size()-1); // add to the current group
                        add_newborn(pop_list[i_pop],pop_list[i_pop].get().size()-1); // add to the current group, new-born lists
                    }
                }

            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::mate(bool, int, double)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::mate(bool, int, double)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::mate(bool, int, double)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::mate()
    {
        try
        {
            mate(true, 2, 0.8);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::mate()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::mate()" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::mate()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::regroup_newborn( Group &grp )
    {
        try
        {
            Group tmp;
            tmp.pop_list = pop_list_newborn;
            tmp.individuals_list = individuals_list_newborn;

            grp.add(tmp);

            tmp.clear();

            // remove new-borns from the group
            for (size_t i = 0; i < pop_list_newborn.size(); i++)
            {
                std::vector<std::reference_wrapper<Population>>::iterator position = std::find( pop_list.begin(), pop_list.end(), std::ref(pop_list_newborn[i]) );
                
                if (position != pop_list.end()) // if the population is already in the group
                {
                    // the population is already in the list,
                    // therefore just remove id from correct position
                    size_t pos = position - pop_list.begin();
                    
                    for (size_t j = 0; j < individuals_list_newborn[i].size(); j++)
                    {
                        // check for ids position to remove
                        size_t id = individuals_list_newborn[i][j];
                        std::vector<size_t>::iterator list_position = std::find( individuals_list[pos].begin(), individuals_list[pos].end(), id );
                        size_t id_pos = list_position - individuals_list[pos].begin();
                        // remove from the group
                        if (list_position != individuals_list[pos].end()) // if the element was found
                            individuals_list[pos].erase(individuals_list[pos].begin() + id_pos);
                    }
                }
            }
            // ---------------------------------------

            pop_list_newborn.clear();
            pop_list_newborn.shrink_to_fit();

            individuals_list_newborn.clear();
            individuals_list_newborn.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::regroup_newborn(Group &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::regroup_newborn(Group &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::regroup_newborn(Group &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::aging(int delta_t)
    {
        try
        {
            for (size_t i = 0; i < pop_list.size(); i++)
            {                    
                for (size_t j = 0; j < individuals_list[i].size(); j++)
                {
                    int new_age = pop_list[i].get().age_at( individuals_list[i][j] ) + delta_t;
                    pop_list[i].get().age_at( individuals_list[i][j], new_age );
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::aging(int)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::aging(int)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::aging(int)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::genotype()
    {
        try
        {
            for (size_t i = 0; i < pop_list.size(); i++)
            {                    
                for (size_t j = 0; j < individuals_list[i].size(); j++)
                {
                    pop_list[i].get().isgenotyped_at( individuals_list[i][j], true );
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::genotype()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::genotype()" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::genotype()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::kill()
    {
        try
        {
            for (size_t i = 0; i < pop_list.size(); i++)
            {                    
                for (size_t j = 0; j < individuals_list[i].size(); j++)
                {
                    pop_list[i].get().alive_at( individuals_list[i][j], false );
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::kill()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::kill()" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::kill()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::make_observation(Trait &trt, std::vector<double> &env)
    {
        try
        {
            if (trt.cleared)
                throw std::string("The trait is configured, hence the call of set_trait() is required!");

            evolm::matrix<size_t> shape(2, 1);
            shape = trt.a.shape();

            size_t gr_size = size();

            if (gr_size < 1)
                throw std::string("Cannot provide observation on empty group!");

            for (size_t g = 0; g < gr_size; g++)
            {
                Population in_pop = get_population(g);
                std::vector<size_t> individuals_list = get_individuals(g);

                size_t n_individuals = individuals_list.size();

                if (n_individuals < 1)
                    throw std::string("Cannot provide observation on empty group => there are no selected individuals in the group!");

                size_t n_traits = shape[1];

                if (env.size() != n_traits)
                    throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

                if (trt.base_genome_structure != in_pop.get_genome_structure())
                    throw std::string("The genome structure of the base population is not the same as in the population being observed!");

                trt.realloc_traits(n_individuals, n_traits);

                trt.calculate_trait(in_pop, individuals_list, env, n_traits);

                evolm::matrix<double> t(n_individuals, n_traits);

                for (size_t i = 0; i < n_traits; i++)
                {
                    for (size_t j = 0; j < n_individuals; j++)
                    {
                        t(j, i) = trt.ta(j, i) + trt.te(j, i) + trt.t_mean(i, 0);
                    }
                }

                // register the observed phenotypes for corresponding individual
                for (size_t i = 0; i < n_individuals; i++)
                {
                    std::vector<double> obs;
                    for (size_t j = 0; j < n_traits; j++)
                    {
                        obs.push_back( t(i,j) );
                    }
                    in_pop.phenotype_at(individuals_list[i],obs);
                }

                trt.ta.clear();
                trt.te.clear();
                t.clear();
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<double> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<double> &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<double> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::make_observation(Trait &trt, std::vector<double> &env, const std::string &out_trvalues)
    {
        try
        {
            if (trt.cleared)
                throw std::string("The trait is configured, hence the call of set_trait() is required!");

            Utilites u;
            u.fremove(out_trvalues);

            evolm::matrix<size_t> shape(2, 1);
            shape = trt.a.shape();

            size_t gr_size = size();

            if (gr_size < 1)
                throw std::string("Cannot provide observation on empty group!");

            for (size_t g = 0; g < gr_size; g++)
            {
                Population in_pop = get_population(g);
                std::vector<size_t> individuals_list = get_individuals(g);

                size_t n_individuals = individuals_list.size();

                if (n_individuals < 1)
                    throw std::string("Cannot provide observation on empty group => there are no selected individuals in the group!");

                size_t n_traits = shape[1];

                if (env.size() != n_traits)
                    throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

                if (trt.base_genome_structure != in_pop.get_genome_structure())
                    throw std::string("The genome structure of the base population is not the same as in the population being observed!");

                trt.realloc_traits(n_individuals, n_traits);

                trt.calculate_trait(in_pop, individuals_list, env, n_traits);

                evolm::matrix<double> t(n_individuals, n_traits);

                for (size_t i = 0; i < n_traits; i++)
                {
                    for (size_t j = 0; j < n_individuals; j++)
                    {
                        t(j, i) = trt.ta(j, i) + trt.te(j, i) + trt.t_mean(i, 0);
                    }
                }

                // register the observed phenotypes for corresponding individual
                for (size_t i = 0; i < n_individuals; i++)
                {
                    std::vector<double> obs;
                    for (size_t j = 0; j < n_traits; j++)
                    {
                        obs.push_back( t(i,j) );
                    }
                    in_pop.phenotype_at(individuals_list[i],obs);
                }

                t.printf(out_trvalues, true); // in append mode

                trt.ta.clear();
                trt.te.clear();
                t.clear();
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<double> &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<double> &, const std::string &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<double> &, const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::make_observation(Trait &trt, std::vector<double> &env, const std::string &out_trvalues, const std::string &out_genotypes)
    {
        try
        {
            if (trt.cleared)
                throw std::string("The trait is configured, hence the call of set_trait() is required!");

            Utilites u;
            u.fremove(out_trvalues);
            u.fremove(out_genotypes);

            evolm::matrix<size_t> shape(2, 1);
            shape = trt.a.shape();

            size_t gr_size = size();

            if (gr_size < 1)
                throw std::string("Cannot provide observation on empty group!");

            for (size_t g = 0; g < gr_size; g++)
            {
                Population in_pop = get_population(g);
                std::vector<size_t> individuals_list = get_individuals(g);

                size_t n_individuals = individuals_list.size();

                if (n_individuals < 1)
                    throw std::string("Cannot provide observation on empty group => there are no selected individuals in the group!");

                size_t n_traits = shape[1];

                if (env.size() != n_traits)
                    throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

                if (trt.base_genome_structure != in_pop.get_genome_structure())
                    throw std::string("The genome structure of the base population is not the same as in the population being observed!");

                trt.realloc_traits(n_individuals, n_traits);

                trt.calculate_trait(in_pop, individuals_list, env, n_traits);

                evolm::matrix<double> t(n_individuals, n_traits);

                for (size_t i = 0; i < n_traits; i++)
                {
                    for (size_t j = 0; j < n_individuals; j++)
                    {
                        t(j, i) = trt.ta(j, i) + trt.te(j, i) + trt.t_mean(i, 0);
                    }
                }

                // register the observed phenotypes for corresponding individual
                for (size_t i = 0; i < n_individuals; i++)
                {
                    std::vector<double> obs;
                    for (size_t j = 0; j < n_traits; j++)
                    {
                        obs.push_back( t(i,j) );
                    }
                    in_pop.phenotype_at(individuals_list[i],obs);
                }

                t.printf(out_trvalues, true); // in append mode

                trt.ta.clear();
                trt.te.clear();
                t.clear();

                in_pop.get_all_genotypes(out_genotypes); // in append mode
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<double> &, const std::string &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<double> &, const std::string &, const std::string &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<double> &, const std::string &, const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::make_observation(Trait &trt, std::vector<double> &env, std::vector<std::vector<double>> &out_trvalues)
    {
        try
        {
            if (trt.cleared)
                throw std::string("The trait is configured, hence the call of set_trait() is required!");

            evolm::matrix<size_t> shape(2, 1);
            shape = trt.a.shape();

            size_t gr_size = size();

            if (gr_size < 1)
                throw std::string("Cannot provide observation on empty group!");

            for (size_t g = 0; g < gr_size; g++)
            {
                Population in_pop = get_population(g);
                std::vector<size_t> individuals_list = get_individuals(g);

                size_t n_individuals = individuals_list.size();

                if (n_individuals < 1)
                    throw std::string("Cannot provide observation on empty group => there are no selected individuals in the group!");

                size_t n_traits = shape[1];

                if (env.size() != n_traits)
                    throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

                if (trt.base_genome_structure != in_pop.get_genome_structure())
                    throw std::string("The genome structure of the base population is not the same as in the population being observed!");

                trt.realloc_traits(n_individuals, n_traits);

                trt.calculate_trait(in_pop, individuals_list, env, n_traits);

                evolm::matrix<double> t(n_individuals, n_traits);

                for (size_t i = 0; i < n_traits; i++)
                {
                    for (size_t j = 0; j < n_individuals; j++)
                    {
                        t(j, i) = trt.ta(j, i) + trt.te(j, i) + trt.t_mean(i, 0);
                    }
                }

                // register the observed phenotypes for corresponding individual
                for (size_t i = 0; i < n_individuals; i++)
                {
                    std::vector<double> obs;
                    for (size_t j = 0; j < n_traits; j++)
                    {
                        obs.push_back( t(i,j) );
                    }
                    in_pop.phenotype_at(individuals_list[i],obs);
                }

                t.to_vector2d(out_trvalues);

                trt.ta.clear();
                trt.te.clear();
                t.clear();
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<double> &, std::vector<std::vector<double>> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<double> &, std::vector<std::vector<double>> &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<double> &, std::vector<std::vector<double>> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::make_observation(Trait &trt, std::vector<double> &env, std::vector<std::vector<double>> &out_trvalues, std::vector<std::vector<short>> &out_genotypes)
    {
        try
        {
            if (trt.cleared)
                throw std::string("The trait is configured, hence the call of set_trait() is required!");

            evolm::matrix<size_t> shape(2, 1);
            shape = trt.a.shape();

            size_t gr_size = size();

            if (gr_size < 1)
                throw std::string("Cannot provide observation on empty group!");

            for (size_t g = 0; g < gr_size; g++)
            {
                Population in_pop = get_population(g);
                std::vector<size_t> individuals_list = get_individuals(g);

                size_t n_individuals = individuals_list.size();

                if (n_individuals < 1)
                    throw std::string("Cannot provide observation on empty group => there are no selected individuals in the group!");

                size_t n_traits = shape[1];

                if (env.size() != n_traits)
                    throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

                if (trt.base_genome_structure != in_pop.get_genome_structure())
                    throw std::string("The genome structure of the base population is not the same as in the population being observed!");

                trt.realloc_traits(n_individuals, n_traits);

                trt.calculate_trait(in_pop, individuals_list, env, n_traits);

                evolm::matrix<double> t(n_individuals, n_traits);

                for (size_t i = 0; i < n_traits; i++)
                {
                    for (size_t j = 0; j < n_individuals; j++)
                    {
                        t(j, i) = trt.ta(j, i) + trt.te(j, i) + trt.t_mean(i, 0);
                    }
                }

                // register the observed phenotypes for corresponding individual
                for (size_t i = 0; i < n_individuals; i++)
                {
                    std::vector<double> obs;
                    for (size_t j = 0; j < n_traits; j++)
                    {
                        obs.push_back( t(i,j) );
                    }
                    in_pop.phenotype_at(individuals_list[i],obs);
                }

                t.to_vector2d(out_trvalues);

                trt.ta.clear();
                trt.te.clear();
                t.clear();

                in_pop.get_all_genotypes(out_genotypes);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<double> &, std::vector<std::vector<double>> &, std::vector<std::vector<short>> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<double> &, std::vector<std::vector<double>> &, std::vector<std::vector<short>> &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<double> &, std::vector<std::vector<double>> &, std::vector<std::vector<short>> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================


} // end of namespace evogen