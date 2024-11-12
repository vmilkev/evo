import numpy as np
from random import randint
import math
import evogen

def select_parents(select_male, from_pop, to_group, n_select, do_random):
    
    sex = 0 # female is the default sex
    if select_male: # if true, select males
        sex = 1
    
    males_list = []
    for i in range( from_pop.size() ): # calculate the total number of males in population
         if from_pop.sex_at(i) == sex:
              males_list.append(i)

    all_males = len(males_list)

    if all_males == 0:
        return

    #if select_male: # if true, select males
    #    print("available males: ", all_males)
    #else:
    #    print("available females: ", all_males)

    if all_males <= n_select: # if there are not enough males, add all of them to the group
         for i in range(all_males):
            #print("not enough, selecting id: ", from_pop.id_at(males_list[i]), " sex: ", from_pop.sex_at(males_list[i]), " no. ", i)
            to_group.add( from_pop, males_list[i] )
         return
    
    if do_random: # perform random selection from the list of males
        up_bound = 100*n_select/all_males
        for i in range(all_males):
            r_num = randint(1, 100)
            if r_num < up_bound:
                #print("random, selecting id: ", from_pop.id_at(males_list[i]), " sex: ", from_pop.sex_at(males_list[i]))
                to_group.add( from_pop, males_list[i] )
    else:
        selected = 0
        for i in range(all_males):
            if selected < n_select:
                #print("non-random, selecting id: ", from_pop.id_at(males_list[i]), " sex: ", from_pop.sex_at(males_list[i]))
                to_group.add( from_pop, males_list[i] )
                selected = selected + 1

def simulate_ancestral_pop(starting_from, proportion_for_reprod, pop_limit, repr_rate, n_gen):
    
    n_offsprings = 2 * repr_rate
    #print("expected offsprings per mmating: ", n_offsprings)

    pop = evogen.Population() #simulated
    pop.set_population(starting_from, "tests/data/discrete_gen/struct_haplotypes_pop1.dat", 1.0, 2); #Needs description in the doc for the haplotype structure
    
    par_gr = evogen.Group() # group for parents
    off_gr = evogen.Group() # group for offspring
    pop_gr = evogen.Group() # group for individuals to be removed at each iteration

    pop.get_ld("first_ld_file");

    for i in range( n_gen ):
        
        print(" ")
        print("==> NEW iteration: ", i)
        print(" ")

        pop_gr.add(pop)

        all_animals = pop.size()

        proportion2 = pop_limit / ( repr_rate * all_animals )

        current_proportion = min(proportion_for_reprod, proportion2)

        #print("current proportion for reproduction: ", current_proportion)

        if proportion_for_reprod <= 1.0: # if the number in the range [0,1] is provided
            to_be_selected = math.floor( all_animals * current_proportion / 2.0 )
        else: # otherwise, proportion_for_reprod should be in the range [1, inf)
            to_be_selected = proportion_for_reprod
        
        print("to_be_selected (of each sex): ", to_be_selected)
        
        select_parents(True, pop, par_gr, to_be_selected, True); # selecting males
        n_males = par_gr.size_at(0)

        select_parents(False, pop, par_gr, to_be_selected, True); # selecting females
        n_females = par_gr.size_at(0) - n_males
        
        print("Pop size: ", all_animals, " selected males: ", n_males, " selected females: ", n_females)
                
        #print("Mating ...")
        par_gr.mate(True, n_offsprings, 1.0) # sexual mating
        
        print("Pop size after mating = ", pop.size())

        #print("Regrouping ...")
        par_gr.regroup_newborn(off_gr)

        print("number of offspring: ", off_gr.size_at(0))

        pop_gr.remove()

        print("Pop size = ", pop.size(), "; current capacity = ", pop.capacity() )
        
        #print("Clearing ...")
        pop_gr.clear()
        par_gr.clear()
        off_gr.clear()

        if i % 10 == 0:
            print("Reshaping ...")
            pop.reshape()
            print("Reshaped pop size = ", pop.size(), "; current capacity = ", pop.capacity() )
    
    pop.get_ld("test_results/last_ld_file");

def main():

    init_num_animals = 18
    generations = 500
    reproduction_rate = 1.5
    proportion_in_reproduction = 1.0
    population_limit = 30

    simulate_ancestral_pop( init_num_animals, proportion_in_reproduction, population_limit, reproduction_rate, generations )

if __name__ == '__main__':
    main()