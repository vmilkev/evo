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
        clear();
    }
    //===============================================================================================================
    Population & Group::get_population( size_t which_pop )
    {
        return pop_list[which_pop].get();
    }
    //===============================================================================================================
    std::vector<size_t> Group::get_individuals( size_t which_pop )
    {
        return individuals_list[which_pop];
    }
    //===============================================================================================================
    void Group::clear()
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
    //===============================================================================================================
    size_t Group::size()
    {
        // Return the number of different populations consisting the group.
        return pop_list.size();
    }
    //===============================================================================================================
    size_t Group::size_at(size_t at)
    {
        //Returns the number of individuals in the group belonnging to a population at the position at in the group.
        if ( at < pop_list.size() )
            return individuals_list[at].size();
        else
            return 0;
    }
    //===============================================================================================================
    void Group::remove()
    {
        // For all individuals in the group: Remove individuals from their original populations and from the group.
        // This is also clears the group.
        try
        {
            if ( individuals_list.empty() && pop_list.empty() )
            {
                std::cout<< "WARNING: Cannot remove individuals from the empty Group!" <<'\n';
                return;
            }

            for (size_t i = 0; i < pop_list.size(); i++)
            {
                if ( pop_list[i].get().size() == 0 )
                    continue;

                if ( individuals_list[i].empty() )
                    throw std::string("Cannot remove individuals from the empty list in the Group!");
                
                for (size_t j = 0; j < individuals_list[i].size(); j++)
                    pop_list[i].get().remove_at( individuals_list[i][j] ); // individuals_list[i][j] points to the position (index) in the active_individuals list
            }
            clear();
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
        // Adds all individuals from the population pop to the group.
        try
        {
            if (pop.size() == 0)
            {
                std::cout << "WARNING: Adding the empty population to the group!" <<'\n';
                return;
            }
            
            std::vector<std::reference_wrapper<Population>>::iterator position = std::find( pop_list.begin(), pop_list.end(), std::ref(pop) );
            
            if (position != pop_list.end()) // the population is in the host group
            {
                std::cout<<"WARNING: Adding the population which is already in the group. This may lead to the segmentation fault error later on if the adding population was modified since the last group's usage!" << '\n';
                std::cout<<"         To update the group with a modified population, the method clean() should be called on the group before the method can be re-applied." << '\n';
                std::cout<<"         Otherwise, add each specific individual from the population to the group directly." << '\n';
                return;
            }

            pop_list.push_back(pop);
            std::vector<size_t> list;
            for (size_t i = 0; i < pop.active_individuals.size(); i++ )
            {
                if (pop.active_individuals[i] == -1)
                    continue;
                list.push_back(i);
            }
            individuals_list.push_back(list);
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
        // Adds all individuals from the group grp to the calling group.
        try
        {
            for (size_t i = 0; i < grp.size(); i++)
            {
                std::vector<std::reference_wrapper<Population>>::iterator position = std::find( pop_list.begin(), pop_list.end(), std::ref(grp.pop_list[i]) );
                if (position != pop_list.end()) // if the population was actually found
                {
                    // the population is already in the list, therefore just adding the id to correct position
                    size_t pos = position - pop_list.begin();

                    // avoid duplication of ids, do not add the ones which are already in the host group;
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
                    // we do not need to check for the unique ids here because we assume
                    // all ids coming from the same population are unique by default
                    pop_list.push_back(grp.pop_list[i]);
                    std::vector<size_t> list;
                    for (size_t j = 0; j < grp.size_at(i); j++)
                        list.push_back(grp.individuals_list[i][j]);
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
        // Adds the specific individual (determined as the position which_one) from the population pop to the group.
        try
        {
            if ( which_one > pop.active_individuals.size()-1 )
                throw std::string("Trying to add illegal ID: out of range!");

            if (pop.active_individuals[which_one] == -1)
                throw std::string("Trying to add illegal ID: not alive!");

            std::vector<std::reference_wrapper<Population>>::iterator position = std::find( pop_list.begin(), pop_list.end(), std::ref(pop) );
            
            if (position != pop_list.end()) // if the population is already in the group
            {
                // the population is already in the list, therefore just adding the id to correct position
                size_t pos = position - pop_list.begin();
                // check for possible duplicates
                std::vector<size_t>::iterator list_position = std::find( individuals_list[pos].begin(), individuals_list[pos].end(), which_one );
                if (list_position == individuals_list[pos].end()) // if the element was not found
                    individuals_list[pos].push_back( which_one );
            }
            else
            {
                // we do not need to check for the unique ids here because we assume
                // all ids coming from the same population are unique by default
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
            if ( which_one > pop.active_individuals.size()-1 )
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
                // we do not need to check for the unique ids here because we assume
                // all ids coming from the same population are unique by default
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
        // Relocate all individuals in the group to the population pop.
        // The relocated individuals will be removed from their original populations.
        try
        {
            for (size_t i = 0; i < pop_list.size(); i++)
            {
                if (pop == pop_list[i])
                    continue;

                for (size_t j = 0; j < individuals_list[i].size(); j++)
                {
                    poplen_t pos_in_activelist = (poplen_t)individuals_list[i][j]; // position in the individuals_list[i][j], points to the position (index) in the active_individuals list inside population
                    poplen_t pos_animal = pop_list[i].get().active_individuals[pos_in_activelist]; // position in individuals vector (container in population)

                    if ( pos_animal == -1 )
                        throw std::string("Trying to access already removed individual!");
                    
                    pop.add( pop_list[i].get().individuals[ (size_t)pos_animal ] );
                    pop_list[i].get().remove_at( (size_t)pos_in_activelist );
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
    void Group::mate(bool sexual_reproduction, int max_offspring, float success_rate, double mutation_frequency, size_t num_crossovers)
    {
        try
        {
            // auto start = std::chrono::high_resolution_clock::now();

            Utilites u;

            // association lists
            std::vector<size_t> male_list;
            std::vector<size_t> male_pops;

            std::vector<std::vector<size_t>> fem_by_pop; // !!!
            std::vector<int> fem_pop; // !!!

            // initialize lists
            for (size_t i = 0; i < pop_list.size(); i++) // loop over all populations in the group
            {
                int fem_ipop = -1;
                std::vector<size_t> fem_list;

                for (size_t j = 0; j < individuals_list[i].size(); j++) // loop over all individuals from the specific population
                {
                    short sex = 1;

                    if ( !sexual_reproduction ) // mix group of males and females to let them mate regardless of sex (though, this does not allow a self-mate)
                    {
                        std::vector<int> which_sex = u.get_uni_rand( 1, 1, 100, false );
                        if ( which_sex[0] < 50 )
                            sex = 0;
                    }

                    if ( pop_list[i].get().sex_at( individuals_list[i][j] ) == sex ) // if male
                    {
                        male_list.push_back( individuals_list[i][j] );
                        male_pops.push_back(i);
                    }
                    else // if female
                    {
                        fem_list.push_back( individuals_list[i][j] );
                        fem_ipop = (int)i;
                    }
                }

                if ( !fem_list.empty()  )
                {
                    fem_by_pop.push_back(fem_list);
                    fem_pop.push_back(fem_ipop);
                }
            }

            // auto stop = std::chrono::high_resolution_clock::now();
            // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            // std::cout <<"==> maating part 1 duration, millisec "<< duration.count() << std::endl;


            if ( male_list.empty() )
            {
                std::cout << "WARNING: The list of males available for mating is empty." <<'\n';
                return;
            }

            if ( fem_by_pop.empty() )
            {
                std::cout << "WARNING: The list of females available for mating is empty." <<'\n';
                return;
            }

            size_t n_male = male_list.size()-1;

            // start = std::chrono::high_resolution_clock::now();

            for (size_t i = 0; i < fem_by_pop.size(); i++) // loop over different populations present in the group
            {
                size_t i_pop = (size_t)fem_pop[i];

                std::vector<int> n_offsprings = u.get_bin_rand( fem_by_pop[i].size(), max_offspring, success_rate, false ); // sample number of offsprings for each female in pop i
                std::vector<unsigned long> which_male = u.get_uni_rand( fem_by_pop[i].size(),(unsigned long)0,n_male,false ); // sample males to mate with each specific female

                int all_expected_offs = std::accumulate(n_offsprings.begin(), n_offsprings.end(), 0); // calculate number of elements the population i_pop will be etxended

                // before resizing we need to know the starting indeces for population containers
                std::vector<size_t> in_pos2(fem_by_pop[i].size()); // writing position2 for each female
                in_pos2[0] = pop_list[ i_pop ].get().capacity(); // at active_individuals == at individuals

                for (size_t j = 1; j < fem_by_pop[i].size(); j++)
                    in_pos2[j] = in_pos2[j-1] + n_offsprings[j];

                pop_list[ i_pop ].get().resize( (size_t)all_expected_offs); // resize population

#pragma omp parallel for
                for (size_t l = 0; l < fem_by_pop[i].size(); l++) // loop over females from the specific population
                {
                    size_t i_female = fem_by_pop[i][l];
                    size_t i_anim = pop_list[i_pop].get().active_individuals[i_female];

                    // select the sampled male
                    size_t i_pop2 = male_pops[ which_male[l] ];
                    size_t i_male = male_list[ which_male[l] ];
                    size_t i_anim2 = pop_list[i_pop2].get().active_individuals[i_male];
                    for (size_t j = 0; j < (size_t)n_offsprings[l]; j++)
                    {
                        Animal c( pop_list[i_pop2].get().individuals[i_anim2], pop_list[i_pop].get().individuals[i_anim], mutation_frequency, num_crossovers ); // constructing new-born individual
                        pop_list[i_pop].get().add_at(c, in_pos2[l]+j); // add to female's population
                    }
                }

                for (size_t l = 0; l < fem_by_pop[i].size(); l++)
                {
                    for (size_t j = 0; j < (size_t)n_offsprings[l]; j++)
                    {
                        add( pop_list[i_pop], in_pos2[l]+j ); // add to the current group
                        add_newborn( pop_list[i_pop], in_pos2[l]+j ); // add to the current group, new-born lists
                    }
                }
            }
            // stop = std::chrono::high_resolution_clock::now();
            // duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            // std::cout <<"==> maating part 2 duration, millisec "<< duration.count() << std::endl;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::mate(bool, int, float, double, size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::mate(bool, int, float, double, size_t)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::mate(bool, int, float, double, size_t)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Group::mate()
    {
        try
        {
            mate(true, 2, 0.8f, 0.001, 2);
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
    void Group::mate(bool sexual_reproduction, float max_offspring, float success_rate, double mutation_frequency, size_t num_crossovers)
    {
        // this is just to allow the float type for max_offspring variable (if someone will pass it in Python)
        try
        {
            mate(sexual_reproduction, (int)max_offspring, success_rate, mutation_frequency, num_crossovers);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::mate(bool, float, float, double, size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::mate(bool, float, float, double, size_t)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::mate(bool, float, float, double, size_t)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Group::regroup_newborn( Group &grp )
    {
        // Move all new-born individuals from the calling group to another group grp.
        // The individuals will be cleared from the calling group
        // but will retain the connections to their original populations.
        try
        {
            Group tmp;
            tmp.pop_list = pop_list_newborn;
            tmp.individuals_list = individuals_list_newborn;

            grp.add(tmp);

            tmp.clear();

            for (size_t i = 0; i < pop_list_newborn.size(); i++) // remove new-borns from the group's main container
            {
                std::vector<std::reference_wrapper<Population>>::iterator position = std::find( pop_list.begin(), pop_list.end(), std::ref(pop_list_newborn[i]) );
                
                if (position == pop_list.end()) // the population is not in the group
                    continue;

                size_t pos = position - pop_list.begin(); // population's position
                
                for (size_t j = 0; j < individuals_list_newborn[i].size(); j++)
                {
                    size_t id = individuals_list_newborn[i][j];
                    std::vector<size_t>::iterator list_position = std::find( individuals_list[pos].begin(), individuals_list[pos].end(), id );
                    size_t id_pos = list_position - individuals_list[pos].begin(); // id's position to remove
                    if (list_position != individuals_list[pos].end()) // if the element was found
                        individuals_list[pos].erase(individuals_list[pos].begin() + id_pos); // remove from the group
                }
            }

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
    void Group::select_basic( Group &selected, Group &notselected, const std::string category, const std::string mode, int num, float val )
    {
        try
        {
            int _category = 0;
            int _mode = 0;

            if ( category == "phen") _category = 1;
            if ( category == "bv") _category = 2;
            if ( category == "age") _category = 3;
            if ( category == "sex") _category = 4;
            if ( category == "alive") _category = 5;
            if ( category == "active") _category = 6;
            if ( category == "genot") _category = 7;
            if ( category == "id") _category = 8;
            
            if ( mode == "rand") _mode = 1;
            if ( mode == "best") _mode = 2;
            if ( mode == "worst") _mode = 3;

            switch (_category)
            {
            case 4: case 5: case 6: case 7:
                if ( (val != 0.0f) && (val != 1.0f) )
                    throw std::string("The expeted last input parameter's value for the category " + category + " should be either 0 or 1.");
                break;
            default:
                break;
            }

            if ( num <= 0 )
                throw std::string("The parameter indicating number of selecting candidates should be positive value greater than 0.");

            for (size_t i = 0; i < pop_list.size(); i++) // loop over added (in advance) populations in the group
            {
                std::vector<selection_candidate> candidates; // container for selection candidates

                for (size_t j = 0; j < individuals_list[i].size(); j++) // looop over individuals of population i in the group
                {
                    selection_candidate cand;
                    size_t pos_in_activelist = individuals_list[i][j]; // position in the individuals_list[i][j]
                    size_t pos_animal = pop_list[i].get().active_individuals[pos_in_activelist]; // position in individuals vector (container in population)
                    
                    std::vector<float> vect;
                    float t_var1;
                    float t_var2;

                    bool empty_candidate = true;
                    
                    switch (_category)
                    {
                    case 1:
                        vect = pop_list[i].get().individuals[pos_animal].get_phenotype(); // vector of phenotypes
                        if ( val == -1.0f ) // calculate phenotype mean
                        {
                            t_var1 = std::accumulate(vect.begin(), vect.end(), 0.0f);
                            t_var2 = t_var1 / (float)vect.size();
                        }
                        else // use values of specific phenotype
                        {
                            if ( val >= (float)vect.size() )
                                throw std::string("The expeted last input parameter's value for the category " + category + " should be in the range [0, num_of_traits-1].");
                            t_var2 = vect[ (size_t)val ];
                        }
                        std::get<0>(cand) = t_var2;
                        std::get<1>(cand) = pos_in_activelist;
                        empty_candidate = false;
                        break;
                    case 2:
                        vect = pop_list[i].get().individuals[pos_animal].get_breeding_value(); // vector of bv
                        if ( val == -1.0f ) // calculate bv mean
                        {
                            t_var1 = std::accumulate(vect.begin(), vect.end(), 0.0f);
                            t_var2 = t_var1 / (float)vect.size();
                        }
                        else // use values of specific bv
                        {
                            if ( val >= (float)vect.size() )
                                throw std::string("The expeted last input parameter's value for the category " + category + " should be in the range [0, num_of_traits-1].");
                            t_var2 = vect[ (size_t)val ];
                        }
                        std::get<0>(cand) = t_var2;
                        std::get<1>(cand) = pos_in_activelist;
                        empty_candidate = false;
                        break;
                    case 3:
                        t_var1 = (float)pop_list[i].get().individuals[pos_animal].get_age();
                        if ( val == -1.0f || _mode == 1 ) // select all available if val is default or the mode set to random
                        {
                            std::get<0>(cand) = t_var1;
                            std::get<1>(cand) = pos_in_activelist;
                            empty_candidate = false;
                        }
                        else
                        {
                            //std::cout<<"here 1, _mode "<<_mode<<" t_var1 = "<<t_var1<<" val = "<<val<<'\n';
                            if ( _mode == 2 && t_var1 >= val ) // select above the threshold value defined by val
                            {
                                //std::cout<<"here 2"<<'\n';
                                std::get<0>(cand) = t_var1;
                                std::get<1>(cand) = pos_in_activelist;
                                empty_candidate = false;
                            }
                            if ( _mode == 3 && t_var1 <= val ) // select bellow the threshold value defined by val
                            {
                                std::get<0>(cand) = t_var1;
                                std::get<1>(cand) = pos_in_activelist;
                                empty_candidate = false;
                            }
                        }
                        break;
                    case 4:
                        t_var1 = (float)pop_list[i].get().individuals[pos_animal].get_sex();
                        if ( t_var1 == val )
                        {
                            std::get<0>(cand) = t_var1;
                            std::get<1>(cand) = pos_in_activelist;
                            empty_candidate = false;
                            //std::cout<<"id: "<<pop_list[i].get().individuals[pos_animal].get_id()<<" sex: "<<pop_list[i].get().individuals[pos_animal].get_sex()<<" pos: "<<pos_in_activelist<<'\n';
                        }
                        break;
                    case 5:
                        t_var1 = (float)pop_list[i].get().individuals[pos_animal].get_alive();
                        if ( t_var1 == val )
                        {
                            std::get<0>(cand) = t_var1;
                            std::get<1>(cand) = pos_in_activelist;
                            empty_candidate = false;
                        }
                        break;
                    case 6:
                        t_var1 = (float)pop_list[i].get().individuals[pos_animal].get_active();
                        if ( t_var1 == val )
                        {
                            std::get<0>(cand) = t_var1;
                            std::get<1>(cand) = pos_in_activelist;
                            empty_candidate = false;
                        }
                        break;
                    case 7:
                        t_var1 = (float)pop_list[i].get().individuals[pos_animal].get_isgenotyped();
                        if ( t_var1 == val )
                        {
                            std::get<0>(cand) = t_var1;
                            std::get<1>(cand) = pos_in_activelist;
                            empty_candidate = false;
                        }
                        break;
                    case 8:
                        std::get<0>(cand) = (float)pop_list[i].get().individuals[pos_animal].get_id();
                        std::get<1>(cand) = pos_in_activelist;
                        empty_candidate = false;
                        break;
                    default:
                        throw std::string("Used incorect selection categoty. Allowed categories: phen, bv, age, sex, alive, active, genot.");
                        break;
                    }

                    if ( !empty_candidate )
                        candidates.push_back(cand);
                }

                if ( candidates.empty() )
                    continue;
                
                auto rng = std::default_random_engine {};

                switch (_mode)
                {
                case 1:
                    std::shuffle(std::begin(candidates), std::end(candidates), rng);
                    break;
                case 2:
                    std::sort(candidates.begin(), candidates.end(), sortdesc); // use descending order: best comes first
                    break;
                case 3:
                    std::sort(candidates.begin(), candidates.end()); // use the (default) ascending order: worst comes first
                    break;
                default:
                    throw std::string("Used incorect selection mode. Allowed modes: rand, best, worst.");
                    break;
                }
                
                int num_selected = std::min(num, (int)candidates.size());

                for (size_t l = 0; l < (size_t)num_selected; l++)
                    selected.add( pop_list[i], std::get<1>(candidates[l]) );
                
                for (size_t l = (size_t)num_selected; l < candidates.size(); l++)
                    notselected.add( pop_list[i], std::get<1>(candidates[l]) );
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::select_basic( Group &, Group &, const std::string, const std::string, int, float)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::select_basic( Group &, Group &, const std::string, const std::string, num, float)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::select_basic( Group &, Group &, const std::string, const std::string, int, float)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Group::select( const std::string category, const std::string mode, int num, float val )
    {
        try
        {
            Group selected_grp;
            Group not_selected_grp;
            select_basic( selected_grp, not_selected_grp, category, mode, num, val );
            clear();
            add( selected_grp );
            selected_grp.clear();
            not_selected_grp.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::select( const std::string, const std::string, int, float)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::select( const std::string, const std::string, int, float)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::select( const std::string, const std::string, int, float)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Group::select_into_group( Group &grp, const std::string category, const std::string mode, int num, float val )
    {
        try
        {
            Group selected_grp;
            Group not_selected_grp;
            select_basic( selected_grp, not_selected_grp, category, mode, num, val );
            clear();
            grp.add( selected_grp );
            add( not_selected_grp );
            selected_grp.clear();
            not_selected_grp.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::select_into_group( Group &, const std::string, const std::string, int, float)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::select_into_group( Group &, const std::string, const std::string, int, float)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::select_into_group( Group &, const std::string, const std::string, int, float)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Group::aging(int delta_t)
    {
        // Adds delta_t to the age of every individual in the group
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
                for (size_t j = 0; j < individuals_list[i].size(); j++)
                    pop_list[i].get().isgenotyped_at( individuals_list[i][j], true );
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
        // Disables all individuals in the calling group.
        // This method changes the alive status of an individual to FALSE.
        try
        {
            for (size_t i = 0; i < pop_list.size(); i++)
                for (size_t j = 0; j < individuals_list[i].size(); j++)
                    pop_list[i].get().alive_at( individuals_list[i][j], false );
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

#ifdef PYBIND

    void Group::make_observation(Trait &trt, pybind11::array_t<float> py_env)
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

                //---------------------------
                std::vector<float> env;

                pybind11::buffer_info buf = py_env.request();

                if (buf.ndim != 1)
                    throw std::runtime_error("Number of dimensions must be one");

                float *ptr = static_cast<float *>(buf.ptr);

                for (pybind11::ssize_t i = 0; i < buf.shape[0]; i++)
                    env.push_back(ptr[i]);
                //---------------------------

                if (env.size() != n_traits)
                    throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

                if (trt.base_genome_structure != in_pop.get_genome_structure())
                    throw std::string("The genome structure of the base population is not the same as in the population being observed!");

                trt.realloc_traits(n_individuals, n_traits);

                trt.calculate_trait(in_pop, individuals_list, env, n_traits);

                evolm::matrix<float> t(n_individuals, n_traits);

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
                    std::vector<float> obs;
                    for (size_t j = 0; j < n_traits; j++)
                    {
                        obs.push_back( t(i,j) );
                    }
                    in_pop.phenotype_at_cpp(individuals_list[i],obs);
                }

                trt.ta.clear();
                trt.te.clear();
                t.clear();
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, pybind11::array_t<float>)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, pybind11::array_t<float>)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, pybind11::array_t<float>)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::make_observation(Trait &trt, pybind11::array_t<float> py_env, const std::string &out_trvalues)
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

                //---------------------------
                std::vector<float> env;

                pybind11::buffer_info buf = py_env.request();

                if (buf.ndim != 1)
                    throw std::runtime_error("Number of dimensions must be one");

                float *ptr = static_cast<float *>(buf.ptr);

                for (pybind11::ssize_t i = 0; i < buf.shape[0]; i++)
                    env.push_back(ptr[i]);
                //---------------------------

                if (env.size() != n_traits)
                    throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

                if (trt.base_genome_structure != in_pop.get_genome_structure())
                    throw std::string("The genome structure of the base population is not the same as in the population being observed!");

                trt.realloc_traits(n_individuals, n_traits);

                trt.calculate_trait(in_pop, individuals_list, env, n_traits);

                evolm::matrix<float> t(n_individuals, n_traits);

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
                    std::vector<float> obs;
                    for (size_t j = 0; j < n_traits; j++)
                    {
                        obs.push_back( t(i,j) );
                    }
                    in_pop.phenotype_at_cpp(individuals_list[i],obs);
                }

                t.printf(out_trvalues, true); // in append mode

                trt.ta.clear();
                trt.te.clear();
                t.clear();
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, pybind11::array_t<float>, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, pybind11::array_t<float>, const std::string &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, pybind11::array_t<float>, const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::make_observation(Trait &trt, pybind11::array_t<float> py_env, const std::string &out_trvalues, const std::string &out_genotypes)
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

                //---------------------------
                std::vector<float> env;

                pybind11::buffer_info buf = py_env.request();

                if (buf.ndim != 1)
                    throw std::runtime_error("Number of dimensions must be one");

                float *ptr = static_cast<float *>(buf.ptr);

                for (pybind11::ssize_t i = 0; i < buf.shape[0]; i++)
                    env.push_back(ptr[i]);
                //---------------------------

                if (env.size() != n_traits)
                    throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

                if (trt.base_genome_structure != in_pop.get_genome_structure())
                    throw std::string("The genome structure of the base population is not the same as in the population being observed!");

                trt.realloc_traits(n_individuals, n_traits);

                trt.calculate_trait(in_pop, individuals_list, env, n_traits);

                evolm::matrix<float> t(n_individuals, n_traits);

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
                    std::vector<float> obs;
                    for (size_t j = 0; j < n_traits; j++)
                    {
                        obs.push_back( t(i,j) );
                    }
                    in_pop.phenotype_at_cpp(individuals_list[i],obs);
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
            std::cerr << "Exception in Group::make_observation(Trait &, pybind11::array_t<float>, const std::string &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, pybind11::array_t<float>, const std::string &, const std::string &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, pybind11::array_t<float>, const std::string &, const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::make_observation(Trait &trt, pybind11::array_t<float> py_env, pybind11::array_t<float> py_trvalues)
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

                //---------------------------
                std::vector<float> env;

                pybind11::buffer_info buf1 = py_env.request();

                if (buf1.ndim != 1)
                    throw std::runtime_error("Number of dimensions must be one");

                float *ptr1 = static_cast<float *>(buf1.ptr);

                for (pybind11::ssize_t i = 0; i < buf1.shape[0]; i++)
                    env.push_back(ptr1[i]);
                //---------------------------

                if (env.size() != n_traits)
                    throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

                if (trt.base_genome_structure != in_pop.get_genome_structure())
                    throw std::string("The genome structure of the base population is not the same as in the population being observed!");

                trt.realloc_traits(n_individuals, n_traits);

                trt.calculate_trait(in_pop, individuals_list, env, n_traits);

                evolm::matrix<float> t(n_individuals, n_traits);

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
                    std::vector<float> obs;
                    for (size_t j = 0; j < n_traits; j++)
                    {
                        obs.push_back( t(i,j) );
                    }
                    in_pop.phenotype_at_cpp(individuals_list[i],obs);
                }

                //------------------------------
                std::vector<std::vector<float>> out_trvalues;

                t.to_vector2d(out_trvalues);

                size_t N = out_trvalues.size();
                size_t M = out_trvalues[0].size();

                // reshape array to match input shape
                py_trvalues.resize( {N,M}, false );

                //  allocate the buffer
                //pybind11::array_t<float> py_trvalues2 = pybind11::array_t<float>( N * M );
                //int dim1 = py_trvalues.shape()[0], dim2 = py_trvalues.shape()[1];

                pybind11::buffer_info buf2 = py_trvalues.request();

                if (buf2.ndim != 2)
                    throw std::runtime_error("Number of dimensions must be two");

                float *ptr2 = (float *) buf2.ptr;

                for (size_t i = 0; i < N; i++) {
                    for (size_t j = 0; j < M; j++) {
                        ptr2[ i * M + j ] = out_trvalues[i][j];
                    }
                }
                // --------------------------------------------------

                trt.ta.clear();
                trt.te.clear();
                t.clear();
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, pybind11::array_t<float>, pybind11::array_t<float>)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, pybind11::array_t<float>, pybind11::array_t<float>)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, pybind11::array_t<float>, pybind11::array_t<float>)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::make_observation(Trait &trt, pybind11::array_t<float> py_env, pybind11::array_t<float> py_trvalues, pybind11::array_t<int> py_genotypes)
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

                //---------------------------
                std::vector<float> env;

                pybind11::buffer_info buf1 = py_env.request();

                if (buf1.ndim != 1)
                    throw std::runtime_error("Number of dimensions must be one");

                float *ptr1 = static_cast<float *>(buf1.ptr);

                for (pybind11::ssize_t i = 0; i < buf1.shape[0]; i++)
                    env.push_back(ptr1[i]);
                //---------------------------

                if (env.size() != n_traits)
                    throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

                if (trt.base_genome_structure != in_pop.get_genome_structure())
                    throw std::string("The genome structure of the base population is not the same as in the population being observed!");

                trt.realloc_traits(n_individuals, n_traits);

                trt.calculate_trait(in_pop, individuals_list, env, n_traits);

                evolm::matrix<float> t(n_individuals, n_traits);

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
                    std::vector<float> obs;
                    for (size_t j = 0; j < n_traits; j++)
                    {
                        obs.push_back( t(i,j) );
                    }
                    in_pop.phenotype_at_cpp(individuals_list[i],obs);
                }

                //------------------------------
                std::vector<std::vector<float>> out_trvalues;

                t.to_vector2d(out_trvalues);

                size_t N = out_trvalues.size();
                size_t M = out_trvalues[0].size();

                py_trvalues.resize( {N,M}, false );

                //  allocate the buffer
                pybind11::buffer_info buf2 = py_trvalues.request();

                if (buf2.ndim != 2)
                    throw std::runtime_error("Number of dimensions must be two");

                float *ptr2 = (float *) buf2.ptr;

                for (size_t i = 0; i < N; i++) {
                    for (size_t j = 0; j < M; j++) {
                        ptr2[ i * M + j ] = out_trvalues[i][j];
                    }
                }
                //------------------------------

                trt.ta.clear();
                trt.te.clear();
                t.clear();

                //--------------------------------------
                std::vector<std::vector<short>> out_genotypes;

                in_pop.get_all_genotypes(out_genotypes);

                N = M = 0;

                N = out_genotypes.size();
                M = out_genotypes[0].size();

                py_genotypes.resize( {N,M}, false );

                //std::cout<<"indiv, snps: "<<N<<" "<<M<<"\n";

                //  allocate the buffer
                pybind11::buffer_info buf3 = py_genotypes.request();

                if (buf3.ndim != 2)
                    throw std::runtime_error("Number of dimensions must be two");

                int *ptr3 = (int *) buf3.ptr;

                for (size_t i = 0; i < N; i++) {
                    for (size_t j = 0; j < M; j++) {
                        ptr3[ i * M + j ] = (int)out_genotypes[i][j];
                        //std::cout<<"ptr3: "<<ptr3[ i * M + j ]<<"\n";
                    }
                }
                //--------------------------------------
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<int>)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<int>)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<int>)" << '\n';
            throw;
        }
    }

#else

    //===============================================================================================================

    void Group::make_observation(Trait &trt, std::vector<float> &env)
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

                evolm::matrix<float> t(n_individuals, n_traits);

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
                    std::vector<float> obs;
                    for (size_t j = 0; j < n_traits; j++)
                    {
                        obs.push_back( t(i,j) );
                    }
                    in_pop.phenotype_at_cpp(individuals_list[i],obs);
                }

                trt.ta.clear();
                trt.te.clear();
                t.clear();
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<float> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<float> &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<float> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::make_observation(Trait &trt, std::vector<float> &env, const std::string &out_trvalues)
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

                evolm::matrix<float> t(n_individuals, n_traits);

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
                    std::vector<float> obs;
                    for (size_t j = 0; j < n_traits; j++)
                    {
                        obs.push_back( t(i,j) );
                    }
                    in_pop.phenotype_at_cpp(individuals_list[i],obs);
                }

                t.printf(out_trvalues, true); // in append mode

                trt.ta.clear();
                trt.te.clear();
                t.clear();
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<float> &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<float> &, const std::string &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<float> &, const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::make_observation(Trait &trt, std::vector<float> &env, const std::string &out_trvalues, const std::string &out_genotypes)
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

                evolm::matrix<float> t(n_individuals, n_traits);

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
                    std::vector<float> obs;
                    for (size_t j = 0; j < n_traits; j++)
                    {
                        obs.push_back( t(i,j) );
                    }
                    in_pop.phenotype_at_cpp(individuals_list[i],obs);
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
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<float> &, const std::string &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<float> &, const std::string &, const std::string &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<float> &, const std::string &, const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::make_observation(Trait &trt, std::vector<float> &env, std::vector<std::vector<float>> &out_trvalues)
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

                evolm::matrix<float> t(n_individuals, n_traits);

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
                    std::vector<float> obs;
                    for (size_t j = 0; j < n_traits; j++)
                    {
                        obs.push_back( t(i,j) );
                    }
                    in_pop.phenotype_at_cpp(individuals_list[i],obs);
                }

                t.to_vector2d(out_trvalues);

                trt.ta.clear();
                trt.te.clear();
                t.clear();
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<float> &, std::vector<std::vector<float>> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<float> &, std::vector<std::vector<float>> &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<float> &, std::vector<std::vector<float>> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Group::make_observation(Trait &trt, std::vector<float> &env, std::vector<std::vector<float>> &out_trvalues, std::vector<std::vector<short>> &out_genotypes)
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

                evolm::matrix<float> t(n_individuals, n_traits);

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
                    std::vector<float> obs;
                    for (size_t j = 0; j < n_traits; j++)
                    {
                        obs.push_back( t(i,j) );
                    }
                    in_pop.phenotype_at_cpp(individuals_list[i],obs);
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
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<float> &, std::vector<std::vector<float>> &, std::vector<std::vector<short>> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<float> &, std::vector<std::vector<float>> &, std::vector<std::vector<short>> &)" << '\n';
            std::cerr <<"Reason: "<< e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Group::make_observation(Trait &, std::vector<float> &, std::vector<std::vector<float>> &, std::vector<std::vector<short>> &)" << '\n';
            throw;
        }
    }

#endif
    //===============================================================================================================


} // end of namespace evogen