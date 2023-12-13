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
                        //std::cout<<"copied id: "<<pop_list[i].get().individuals[ pop_list[i].get().active_individuals[ individuals_list[i][j] - 0 ] ].get_id()<<"\n";
                        //std::cout<<"removed id_at: "<<pop_list[i].get().id_at( individuals_list[i][j] - j )<<"\n";
                        
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


} // end of namespace evogen