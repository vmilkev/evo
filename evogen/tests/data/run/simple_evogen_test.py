import sys
#sys.path.append('/usr/home/qgg/vimi/evo/bin') # cluster
sys.path.append('../../../release') # Mac

import evogen
import numpy as np

def simple_test_all():
    print("Testing Population:")

    a = evogen.Population() #simulated
    b = evogen.Population() #simulated
    
    a.set_population(5, "../struct_haplotypes_pop1.dat", 0.4, 4);
    b.set_population(5, "../struct_haplotypes_pop1.dat", 0.6, 6);

    c = evogen.Population() # using data from files
    d = evogen.Population() # using data from files
    
    c.set_population("../haplotypes_pop1.dat", "../struct_haplotypes_pop1.dat", True); # true: haplotypes with pedigree and sex
    d.set_population("../haplotypes_pop2.dat", "../struct_haplotypes_pop2.dat", False); # false: only haplotypes

    print("Testing Trait:")

    pop = evogen.Population()

    pop.set_population(20, "../struct_haplotypes_pop3.dat", 0.7, 4) # (1) population

    tr_mean = [ 40.0, 5.0, 0.5 ] # (2) trait means
    qtl_prop = [ 0.65, 0.65, 0.65, 0.65 ] # (3) proportion of snps selected as qtls
    qtl_prop = [ 1.0, 1.0, 1.0, 1.0 ]
    cor_g = [ [1.0, 0.5, 0.7],
                [0.5, 1.0, 0.2],
                [0.7, 0.2, 1.0] ] # (4) genomic correlations
    cor_e = [ [1.0, 0.3, 0.5],
                [0.3, 1.0, 0.4],
                [0.5, 0.4, 1.0] ] # (5) rsidual correlations

    var_g = [ 100.0, 10.0, 0.1 ] # (6) genomic variances
    var_e = [ 200.0, 20.0, 0.3 ] # (7) residual variances

    env = [ 0.0, 0.0, 0.0 ] # (8) enviironment

    # (i) For uniform distribution, model 1
    k_range_U = [ -1.0, 1.0 ] # range of k parameter
    # (ii) For normal distribution, model 2
    k_range_N = [ 0.0, 0.5 ] # mean & std
    # (iii) For gamma distribution, mdodel 3
    k_range_G = [ 0.05 ] # expected value of k in loci; is '1/b' param. in the distribution

    which_model = 1

    print("Setting trait:")

    T = evogen.Trait(pop, tr_mean, qtl_prop, cor_g, var_g, cor_e, var_e, env, which_model, k_range_U)
    
    T2 = evogen.Trait()
    T2.set_trait(pop, tr_mean, qtl_prop, cor_g, var_g, cor_e, var_e, env, which_model, k_range_U)

    print("Make observations on pop:")
    
    T.get_observations(pop, env, "T_trait_pop.dat")
    T2.get_observations(pop, env, "T2_trait.dat", "T2_genotypes.dat")

    print("Creating groups:")

    G = evogen.Group()
    G2 = evogen.Group()
    
    print("G size = ", G.size())

    print("Adding pop to G:")
    G.add(pop)

    print("G size = ", G.size(), "; size at = ", G.size_at(0))
    
    print("Make observations on G:")
    
    G.make_observation(T,env)

    #obs = np.zeros( shape=(1,1), dtype=np.float32 ) # working
    obs = np.zeros( shape=(1), dtype=np.float32 ) # working
    G.make_observation(T,env, obs)

    gen = np.zeros( shape=(1), dtype=np.int32 ) # note! it is important to pass a correct typpe!
    G.make_observation(T,env, obs, gen)

    print("collected observations:")
    print(obs)
    print(obs.shape)
    print(gen.shape)
    print(gen)

    print("pop size = ", pop.size())
    print("pop capacity = ", pop.capacity())

    print("Select in pop ids < 500000000 to G2:")
    print("group G2 size = ", G2.size_at(0), " ... before selection.")

    count = 0;
    for i in range( pop.size() ):
        if pop.id_at(i) < 500000000:
            count = count + 1
            print("selecting ids: ", pop.id_at(i))
            G2.add(pop, i)

    print("Number of selected ids = ", count)
    print("group G2 size = ", G2.size_at(0))
    
    print("Testing G2.remove()")
    print("group G size: ", G.size_at(0) )
    
    for i in range( pop.size() ):
        print("all pop ids: ", pop.id_at(i), "; active: ", pop.alive_at(i) )

    print("Create empty population pop2:")
    
    pop2 = evogen.Population()

    print("pop2 size = ", pop2.size(), "; capacity = ", pop2.capacity() )
    
    print("Move the content of G2 into pop2, moved ids should dissapiar from pop and G2 should be empty:")
    
    G2.move(pop2);
    
    print("pop2 size = ", pop2.size(), "; capacity = ", pop2.capacity() )
    print("pop size = ", pop.size(), "; capacity = ", pop.capacity() )
    print("group G2 size = ", G2.size(), "; size_at = ", G2.size_at(0) )

    print("The conteent of pop:")        
    for i in range( pop.size() ):
        print("all pop ids: ", pop.id_at(i), "; active: ", pop.alive_at(i) )

    print("The conteent of pop2:")        
    for i in range( pop2.size() ):
        print("all pop2 ids: ", pop2.id_at(i), "; active: ", pop2.alive_at(i) )

    print("Reshaping the pop:")
    pop.reshape()

    print("The conteent of pop after reshaping:")

    print("pop size = ", pop.size(), "; capacity = ", pop.capacity() )
    print("pop2 size = ", pop2.size(), "; capacity = ", pop2.capacity() )
    
    print("After reshaping. Calculating trait in pop2:")
    T.get_observations(pop2, env, "T_trait_pop2.dat");
    
    print("After reshaping. Calculating trait in pop:")
    T.get_observations(pop, env, "T_trait_pop_reduced.dat")
    T.get_observations(pop, env)

    G.clear()
    G2.clear()

    pop.aging(3)
    pop2.aging(5)

    G.add(pop)
    G.add(pop2)

    G.mate()

    print("Relocating new-borns:")

    G_newborn = evogen.Group()
    G3 = evogen.Group()
    pop_newborn = evogen.Population()
    pop3 = evogen.Population()
    
    G.regroup_newborn(G_newborn)

    print("new-borns group, size = ", G_newborn.size() )
    for i in range( G_newborn.size() ):
        print( "new-borns group, size at = ", i, ", ", G_newborn.size_at(i) )

    G_newborn.aging(4)
    G_newborn.genotype()
    G_newborn.kill()

    G_newborn.move(pop_newborn)

    G3.add(pop)
    G3.add(pop2)
    G3.add(pop_newborn)
    G3.move(pop3)

    T.clear()
    T2.clear()

if __name__ == '__main__':      
#def main():
	simple_test_all()

