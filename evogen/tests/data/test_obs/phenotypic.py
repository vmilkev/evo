import sys

sys.path.append('/Users/au383883/Documents/MY/codebase/evo/evogen/release') # cluster

import numpy as np

#void doesnt return any values
import evogen

def simple_test_all():
    print("Testing Population:")

    pop = evogen.Population() #DNK you can get access to single individuals    
    pop.set_population("DNK_hap", "temp_map", False); # false: only haplotypes

    # aging increase the age of the population
    # reshape delete unessary animals
    
    print("Testing Trait:")

    tr_mean = [0] # (2) trait means
    qtl_prop = [ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1] # (3) proportion of snps selected as qtls

    var_g = [1.0] # (6) genomic variances
    var_e = [2.0] # (7) residual variances
    cor_g = [ [0.01] ]
    cor_e = [ [0.01] ]
    
    env = [0.0] # (8) enviironment

    # (i) For uniform distribution, model 1
    k_range_U = [ -1.0, 0.2] # range of k parameter
    # (ii) For normal distribution, model 2
    k_range_N = [ 0.0, 0.2 ] # mean & std
    # (iii) For gamma distribution, mdodel 3
    k_range_G = [ 0.05 ] # expected value of k in loci; is '1/b' param. in the distribution

    which_model = 1

    print("Setting trait:")

    T = evogen.Trait(pop, tr_mean, qtl_prop, cor_g, var_g, cor_e,var_e, env, which_model, k_range_U)

    print("Make observations on this single trait:")
    
    T.get_observations(pop, env, "T_trait_pop.dat")


#    print("Creating groups:")

#    G = evogen.Group() # move animals from one group to another
#    G2 = evogen.Group()
    
#    print("G size = ", G.size())

#    print("Adding a to G:")
#    G.add(a)

#    print("G size = ", G.size(), "; size at = ", G.size_at(0))
    
#    print("Make observations on G:")
    
#    G.make_observation(T,env)
    
#    print("a size = ", a.size())
#    print("a capacity = ", a.capacity())

#    print("Select in a ids < 500000000 to G2:")
#    print("group G2 size = ", G2.size_at(0), " ... before selection.")

#    count = 0;
#    for i in range( a.size() ):
#        if a.id_at(i) < 500000000: # specific individuals, get id depending on the parameter
#            count = count + 1
#            print("selecting ids: ", a.id_at(i))
#            G2.add(a, i)
#
#    print("Number of selected ids = ", count)
#    print("group G2 size = ", G2.size_at(0))
#    
#    print("Testing G2.remove()")
#    print("group G size: ", G.size_at(0) )
#    
#    for i in range( a.size() ):
#        print("all pop ids: ", a.id_at(i), "; active: ", a.alive_at(i) )
#
#    print("Create empty population pop2:")
#    
#    pop2 = evogen.Population()
#
#    print("pop2 size = ", pop2.size(), "; capacity = ", pop2.capacity() )
#    
#    print("Move the content of G2 into pop2, moved ids should dissapiar from pop and G2 should be empty:")
#    
#    G2.move(pop2);
#    
#    print("pop2 size = ", pop2.size(), "; capacity = ", pop2.capacity() )
#    print("pop size = ", a.size(), "; capacity = ", a.capacity() )
#    print("group G2 size = ", G2.size(), "; size_at = ", G2.size_at(0) )
#
#    print("The conteent of pop:")        
#    for i in range( a.size() ):
#        print("all pop ids: ", a.id_at(i), "; active: ", a.alive_at(i) )
#
#    print("The conteent of pop2:")        
#    for i in range( pop2.size() ):
#        print("all pop2 ids: ", pop2.id_at(i), "; active: ", pop2.alive_at(i) )
#
#    print("Reshaping the pop:")
#    a.reshape()
#
#    print("The conteent of pop after reshaping:")
#
#    print("pop size = ", a.size(), "; capacity = ", a.capacity() )
#    print("pop2 size = ", pop2.size(), "; capacity = ", pop2.capacity() )
#    
#    print("After reshaping. Calculating trait in pop2:")
#    T.get_observations(pop2, env, "T_trait_pop2.dat");
#    
#    print("After reshaping. Calculating trait in pop:")
#    T.get_observations(a, env, "T_trait_pop_reduced.dat")
#    T.get_observations(a, env)
#
#    G.clear()
#    G2.clear()
#
#    a.aging(3)
#    pop2.aging(5)
#
#    G.add(a)
#    G.add(pop2)
#
#    G.mate()
#
#    print("Relocating new-borns:")
#
#    G_newborn = evogen.Group()
#    G3 = evogen.Group()
#    pop_newborn = evogen.Population()
#    pop3 = evogen.Population()
#    
#    G.regroup_newborn(G_newborn)
#
#    print("new-borns group, size = ", G_newborn.size() )
#    for i in range( G_newborn.size() ):
#        print( "new-borns group, size at = ", i, ", ", G_newborn.size_at(i) )
#
#    G_newborn.aging(4)
#    G_newborn.genotype()
#    G_newborn.kill()
#
#    G_newborn.move(pop_newborn)
#
#    G3.add(a)
#    G3.add(pop2)
#    G3.add(pop_newborn)
#    G3.move(pop3)
#
#    T.clear()
#    T2.clear()

if __name__ == '__main__':      
    #def main():
	simple_test_all()

