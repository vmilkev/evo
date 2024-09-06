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
        
        G0=evogen.Group()
        G0.add(pop)
        print("Nb of parents:", G0.size_at(0))
    
        G0.mate(0, 2, 1.0)
        
        G1=evogen.Group()
        G0.regroup_newborn(G1)

        print("Off size:", G1.size_at(0))
    
        #G0.remove()
        G0.kill()
        G1.add(pop)
        pop.reshape()
        print("Pop size:", pop.size())
    
    print("Done!")

    #

if __name__ == '__main__':      
#def main():
	simple_test_all()

