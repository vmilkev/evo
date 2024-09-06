import sys

sys.path.append('/Users/au383883/Documents/MY/codebase/evo/evogen/release') # cluster

import evogen
import numpy as np

def simple_test_all():
    print("Create LD through simulation of ancestral population:")
    print("100 individuals mating at random for 10 generations")
    
    pop = evogen.Population() #simulated
    pop.set_population(100, "struct_haplotypes_pop1.dat", 0.5, 2); #Needs description in the doc for the haplotype structure
    print("Pop size:", pop.size())

    #Nb of generations
    cycles = 10;
    for i in range( cycles):
        
        print("generation: ", i)
        
        G0=evogen.Group() # this is local group var
        G0.add(pop)
        
        print("Nb of parents:", G0.size_at(0))
                        
        G0.mate(True, 2, 1.0)
        
        G1=evogen.Group()
        G0.regroup_newborn(G1)

        print("Off size:", G1.size_at(0))
        
        print("Pop size = ", pop.size(), "; current capacity = ", pop.capacity() )

        G0.kill() # this just changes the alive status of animals (will not be able too mate, for example), but they remain part of the population (semen, genotypes, etc. could be used)
        
        G0.remove() # this method trully disables animals and from the simulation scenario point of view they are not exist (even their genotypes) 
                    # however, on 'backstages' all related data remains, which needs to be cleared by calling reshape()

        print("Pop size = ", pop.size(), "; capacity = ", pop.capacity() )

        #G1.add(pop) this is not necessary, since the new-born animals in the G1 are already belong to the pop
        
        pop.reshape() # this should be called after calling remove() and move(); check the pop.capacity() before and after the reshape() call

        print("Pop size = ", pop.size(), "; after reshaping capacity = ", pop.capacity() )
    
    print("Done!")
    #

if __name__ == '__main__':      
#def main():
	simple_test_all()
