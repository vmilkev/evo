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
    
    G0=evogen.Group() # global var
    G1=evogen.Group() # global var
    
    #Nb of generations
    cycles = 10;
    
    for i in range( cycles):
        
        print("generation: ", i)
        
        G0.add(pop)
        
        print("Nb of parents:", G0.size_at(0))
    
        G0.mate(True, 2, 1.0)
        
        G0.regroup_newborn(G1)

        print("Off size:", G1.size_at(0))
    
        G0.remove() # needs reshaping after this call! Check the capacity() method.

        # if does not call this, the data will be accumulated
        G0.clear() # because G0 is global
        G1.clear() # because G1 is global
        
        print("Pop size = ", pop.size(), "; capacity = ", pop.capacity() )
    
    print("Done!")
    #

if __name__ == '__main__':      
#def main():
	simple_test_all()
