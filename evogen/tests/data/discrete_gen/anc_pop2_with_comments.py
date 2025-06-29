import numpy as np
import random
import math
import sys
#sys.path.append('/usr/home/qgg/vimi/evo/bin') # cluster
sys.path.append('../../../release') # Mac
import evogen

def simulate_ancestral_pop(starting_from, num_for_reprod, pop_limit, n_gen):
    
    n_offsprings = 4

    pop_1 = evogen.Population() #simulated
    pop_1.set_population(starting_from, "struct_haplotypes_pop1.dat", 1.0, 2); #Needs description in the doc for the haplotype structure
    #pop_1.set_population(starting_from, "genetic_map_chr_1_3.txt", 1.0, 2);

    pop_2 = evogen.Population() #simulated
    pop_2.set_population(starting_from, "struct_haplotypes_pop1.dat", 1.0, 2);
    #pop_2.set_population(starting_from, "genetic_map_chr_1_3.txt", 1.0, 2);

    pop_3 = evogen.Population()
    #pop_3.set_population(10, "struct_haplotypes_pop1.dat", 1.0, 2);
    
    # Groups
    # off_gr = evogen.Group() # group for offspring
    # sel_gr0 = evogen.Group() # group to test select
    # sel_gr1 = evogen.Group() # group to test select
    # sel_gr2 = evogen.Group() # group to test select

    #pop_1.get_ld("first_ld_file");

    for gen in range( n_gen ):

        off_gr = evogen.Group() # group for offspring
        sel_gr0 = evogen.Group() # group to test select
        sel_gr1 = evogen.Group() # group to test select
        sel_gr2 = evogen.Group() # group to test select

        print(" ")
        print("==> Generation: ", gen+1)
        print(" ")

        print("Pop 1:", pop_1.get_popid(), "size = ", pop_1.size())
        for i in range( pop_1.size() ):
            pos = pop_1.get_valid_pos(i)
            print("id:", pop_1.id_at( pos ), "sire:",pop_1.sire_at(pos), "dame:",pop_1.dame_at(pos), "sex:", pop_1.sex_at(pos), "pos:", pos )
        
        print(" ")
        print("Pop 2:", pop_2.get_popid(), "size = ", pop_2.size())
        for i in range( pop_2.size() ):
            pos = pop_2.get_valid_pos(i)
            print("id:", pop_2.id_at(pos), "sire:",pop_2.sire_at(pos), "dame:",pop_2.dame_at(pos), "sex:", pop_2.sex_at(pos), "pos:", pos )
        
        if (gen+1) <= 1:
            print("=> gen: ", gen+1)
            sel_gr1.add(pop_1)
            sel_gr2.add(pop_2)

            sel_gr1.select("sex", "rand", num_for_reprod, 1)
            sel_gr2.select("sex", "rand", num_for_reprod, 0)

            print("individuals (males) in selection group 1 after selection:", sel_gr1.size_at(0))
            print("individuals (females) in selection group 2 after selection:", sel_gr2.size_at(0))
            print(" ")

            sel_gr0.add(sel_gr1) # nothing happens to the group sel_gr1 at this stage
            sel_gr0.add(sel_gr2) # nothing happens to the group sel_gr2 at this stage
            
            print("individuals in selection group 0 after adding group 1 & 2:", sel_gr0.size_at(0)+sel_gr0.size_at(1))

            print("Mating between pop 1 and pop 2 ...")

            sel_gr0.mate(True, n_offsprings, 1.0, 0.2e-4, 2) # sexual mating
            
            print("individuals in selection group 0 after mating:", sel_gr0.size_at(0)+sel_gr0.size_at(1))

            sel_gr0.regroup_newborn(off_gr)

            print("individuals in selection group 0 after regrouping:", sel_gr0.size_at(0)+sel_gr0.size_at(1))

            sel_gr1.clear()
            sel_gr2.clear()
            sel_gr0.clear()

            to_select = math.floor( (off_gr.size_at(0)+off_gr.size_at(1)) /2)

            off_gr.select_into_group(sel_gr1,"id", "rand", to_select, 1)

            # sel_gr1.move(pop_1)
            # off_gr.move(pop_2)
            sel_gr1.move(pop_3)
            off_gr.move(pop_3)

            print(" ")
            print("Pop 3:", pop_3.get_popid(), "size = ", pop_3.size())
            for i in range( pop_3.size() ):
                pos = pop_3.get_valid_pos(i)
                print("id:", pop_3.id_at(pos), "sire:",pop_3.sire_at(pos), "dame:",pop_3.dame_at(pos), "sex:", pop_3.sex_at(pos), "pos:", pos )
            print(" ")
        else:
            print("=> gen: ", gen+1)
            print("Mating in pop 3 ...")
            sel_gr1.add(pop_3)
            sel_gr1.mate(True, n_offsprings, 1.0, 0.2e-4, 2)
        
            print(" ")
            print("Pop 3:", pop_3.get_popid(), "size = ", pop_3.size())
            for i in range( pop_3.size() ):
                pos = pop_3.get_valid_pos(i)
                print("id:", pop_3.id_at(pos), "sire:",pop_3.sire_at(pos), "dame:",pop_3.dame_at(pos), "sex:", pop_3.sex_at(pos), "pos:", pos )
            print(" ")

        if i % 10 == 0:
            print("Reshaping ...")
            pop_1.reshape()
            pop_2.reshape()
            print("Reshaped pop_1 size = ", pop_1.size(), "; current capacity = ", pop_1.capacity() )
            print("Reshaped pop_2 size = ", pop_2.size(), "; current capacity = ", pop_2.capacity() )

            #print("Calculating LD ...")
            #pop_1.get_ld("ld_at_"+str(i), False, 0, 10);
        
        print("... end of loop")

    #pop_1.get_ld("ld_at_"+str(i), False, 0, 10);
    pop_1.get_ancestry("pop1_ancestry.txt");
    pop_2.get_ancestry("pop2_ancestry.txt");
    pop_3.get_ancestry("pop3_ancestry.txt");
    
    pop_1.get_haplotypes("pop1_haplotypes.txt");
    pop_2.get_haplotypes("pop2_haplotypes.txt");
    pop_3.get_haplotypes("pop3_haplotypes.txt");
    print("End of simulation")


def main():

    init_num_animals = 5
    generations = 2
    num_in_reproduction = 1
    population_limit = 10

    simulate_ancestral_pop( init_num_animals, num_in_reproduction, population_limit, generations )

if __name__ == '__main__':
    main()