import numpy as np
import random
import math
import sys
#sys.path.append('/usr/home/qgg/vimi/evo/bin') # cluster
sys.path.append('../../../release') # Mac
import evogen as eg

# define some constants:
male_sex = 1
female_sex = 0

# ----------------------------------------------------------------------------
def make_admixtured_population_v1( afrpop, easpop, eurpop, mixpop, n_offspr ):
    
    # generating admixtured individuals for the mix_pop: AFR + EAS + EUR -> MIX
    # selecting males and females into different groups
    # in order to allow a well mixtured population where we do not want
    # males and females from the same population contribute offsprings
    # to the mix_pop, only females mate with males from distinct population;
    # therefore, we design the following mating plan:
    # females from AFR populaation mate with EUR + EAS males;
    # females from EAS population mate with AFR + EUR males;
    # females from EUR population mate with AFR + EAS males;

    num_selected = math.ceil( afrpop.size()/(3*n_offspr) ) # note, here we need integer value (not a floating point)

    print(" ")
    print("number of females required to be selected from each population is", num_selected)

    # instantiating required groups:

    gr_afr_males = eg.Group() # group for selected males from afr_pop
    gr_afr_females = eg.Group() # group for selected males from afr_pop
    gr_eas_males = eg.Group() # group for selected males from eas_pop
    gr_eas_females = eg.Group() # group for selected males from eas_pop
    gr_eur_males = eg.Group() # group for selected males from eur_pop
    gr_eur_females = eg.Group() # group for selected males from eur_pop

    gr_afr_females.add(afrpop) # now the group holds the entire population
    gr_afr_females.select("sex", "rand", num_selected, female_sex) # gr_afr_females now left only with required (selected) females 
    gr_afr_males.add(afrpop) # now the group holds the entire population
    gr_afr_males.select("sex", "rand", num_selected, male_sex) # gr_afr_males now left only with required (selected) males 
    
    gr_eas_females.add(easpop)
    gr_eas_females.select("sex", "rand", num_selected, female_sex)
    gr_eas_males.add(easpop)
    gr_eas_males.select("sex", "rand", num_selected, male_sex)
    
    gr_eur_females.add(eurpop)
    gr_eur_females.select("sex", "rand", num_selected, female_sex)
    gr_eur_males.add(eurpop)
    gr_eur_males.select("sex", "rand", num_selected, male_sex)

    print("selected females in afr_pop:", gr_afr_females.size_at(0))
    print("selected males in afr_pop:", gr_afr_males.size_at(0))
    print("selected females in eas_pop:", gr_eas_females.size_at(0))
    print("selected males in eas_pop:", gr_eas_males.size_at(0))
    print("selected females in eur_pop:", gr_eur_females.size_at(0))
    print("selected males in eur_pop:", gr_eur_males.size_at(0))

    gr_afr_males2 = eg.Group() # will hold hlf of selected males from AFR population
    gr_afr_males.select_into_group( gr_afr_males2, "id", "rand", math.floor(gr_afr_males.size_at(0)/2), 1 ) # relocate 1/2 of males in gr_afr_males group into gr_afr_males2 group
    print("males from afr_pop in group 1:", gr_afr_males.size_at(0))
    print("males from afr_pop in group 2:", gr_afr_males2.size_at(0))

    gr_eas_males2 = eg.Group()
    gr_eas_males.select_into_group( gr_eas_males2, "id", "rand", math.floor(gr_eas_males.size_at(0)/2), 1 ) # relocate 1/2 of males in gr_eas_males group into gr_eas_males2 group
    print("males from afr_pop in group 1:", gr_eas_males.size_at(0))
    print("males from afr_pop in group 2:", gr_eas_males2.size_at(0))
    
    gr_eur_males2 = eg.Group()
    gr_eur_males.select_into_group( gr_eur_males2, "id", "rand", math.floor(gr_eur_males.size_at(0)/2), 1 ) # relocate 1/2 of males in gr_eur_males group into gr_eur_males2 group
    print("males from afr_pop in group 1:", gr_eur_males.size_at(0))
    print("males from afr_pop in group 2:", gr_eur_males2.size_at(0))

    # 1-st mating group: AFR fenales + AES males (sub-gr 1) + EUR males (sub-gr 1)
    gr_afr_females.add(gr_eas_males)
    gr_afr_females.add(gr_eur_males)
    
    # checking the actual number of individuals in the gr_afr_females group (all of them will mate randomly within the group):
    n_in_group = 0
    for ipop in range( gr_afr_females.size() ): # loop over registered populations in the group 
        n_in_group += gr_afr_females.size_at(ipop) # access number of individuals from specific population ipop in the group
    print("the actual number of individuals in the gr_afr_females group:", n_in_group)
    
    # 2-nd mating group:  EAS fenales + AFR males (sub-gr 1) + EUR males (sub-gr 2)
    gr_eas_females.add(gr_afr_males)
    gr_eas_females.add(gr_eur_males2)
    
    n_in_group = 0
    for ipop in range( gr_eas_females.size() ):
        n_in_group += gr_eas_females.size_at(ipop)
    print("the actual number of individuals in the gr_eas_females group:", n_in_group)

    # 3-rd mating group:  EUR fenales + AFR males (sub-gr 2) + EAS males (sub-gr 2)
    gr_eur_females.add(gr_afr_males2)
    gr_eur_females.add(gr_eas_males2)
    
    n_in_group = 0
    for ipop in range( gr_eur_females.size() ):
        n_in_group += gr_eur_females.size_at(ipop)
    print("the actual number of individuals in the gr_eur_females group:", n_in_group)
    print(" ")

    # mating:
    gr_afr_females.mate(True, n_offspr, 1.0, 0.2e-4, 2) # sexual mating
    gr_eas_females.mate(True, n_offspr, 1.0, 0.2e-4, 2) # sexual mating
    gr_eur_females.mate(True, n_offspr, 1.0, 0.2e-4, 2) # sexual mating

    off_gr = eg.Group() # creating the group that will hold all offsprings

    gr_afr_females.regroup_newborn(off_gr)
    gr_eas_females.regroup_newborn(off_gr)
    gr_eur_females.regroup_newborn(off_gr)
    
    off_gr.move(mixpop) # finally "physically" moving all offsprings into mix_pop

# ----------------------------------------------------------------------------
def make_admixtured_population_v2( afrpop, easpop, eurpop, mixpop, n_offspr ):
    
    # generating admixtured individuals for the mix_pop: AFR + EAS + EUR -> MIX
    # this function implement a bit simplified schema of v1 (males are not divided into two sub-groups)
    
    num_selected = math.ceil( afrpop.size()/(3*n_offspr) ) # note, here we need integer value (not a floating point)

    print(" ")
    print("number of females required to be selected from each population is", num_selected)

    gr_afr_males = eg.Group() # group for selected males from afr_pop
    gr_afr_females = eg.Group() # group for selected males from afr_pop
    gr_eas_males = eg.Group() # group for selected males from eas_pop
    gr_eas_females = eg.Group() # group for selected males from eas_pop
    gr_eur_males = eg.Group() # group for selected males from eur_pop
    gr_eur_females = eg.Group() # group for selected males from eur_pop

    gr_afr_females.add(afrpop) # now the group holds the entire population
    gr_afr_females.select("sex", "rand", num_selected, female_sex) # gr_afr_females now left only with required (selected) females 
    gr_afr_males.add(afrpop) # now the group holds the entire population
    gr_afr_males.select("sex", "rand", num_selected, male_sex) # gr_afr_males now left only with required (selected) males 
    
    gr_eas_females.add(easpop)
    gr_eas_females.select("sex", "rand", num_selected, female_sex)
    gr_eas_males.add(easpop)
    gr_eas_males.select("sex", "rand", num_selected, male_sex)
    
    gr_eur_females.add(eurpop)
    gr_eur_females.select("sex", "rand", num_selected, female_sex)
    gr_eur_males.add(eurpop)
    gr_eur_males.select("sex", "rand", num_selected, male_sex)

    print("selected females in afr_pop:", gr_afr_females.size_at(0))
    print("selected males in afr_pop:", gr_afr_males.size_at(0))
    print("selected females in eas_pop:", gr_eas_females.size_at(0))
    print("selected males in eas_pop:", gr_eas_males.size_at(0))
    print("selected females in eur_pop:", gr_eur_females.size_at(0))
    print("selected males in eur_pop:", gr_eur_males.size_at(0))

    # 1-st mating group: AFR fenales + AES males (sub-gr 1) + EUR males (sub-gr 1)
    gr_afr_females.add(gr_eas_males)
    gr_afr_females.add(gr_eur_males)
    
    # checking the actual number of individuals in the gr_afr_females group (all of them will mate randomly within the group):
    n_in_group = 0
    for ipop in range( gr_afr_females.size() ): # loop over registered populations in the group 
        n_in_group += gr_afr_females.size_at(ipop) # access number of individuals from specific population ipop in the group
    print("the actual number of individuals in the gr_afr_females group:", n_in_group)
    
    # 2-nd mating group:  EAS fenales + AFR males (sub-gr 1) + EUR males (sub-gr 2)
    gr_eas_females.add(gr_afr_males)
    gr_eas_females.add(gr_eur_males)
    
    n_in_group = 0
    for ipop in range( gr_eas_females.size() ):
        n_in_group += gr_eas_females.size_at(ipop)
    print("the actual number of individuals in the gr_eas_females group:", n_in_group)

    # 3-rd mating group:  EUR fenales + AFR males (sub-gr 2) + EAS males (sub-gr 2)
    gr_eur_females.add(gr_afr_males)
    gr_eur_females.add(gr_eas_males)
    
    n_in_group = 0
    for ipop in range( gr_eur_females.size() ):
        n_in_group += gr_eur_females.size_at(ipop)
    print("the actual number of individuals in the gr_eur_females group:", n_in_group)
    print(" ")

    gr_afr_females.mate(True, n_offspr, 1.0, 0.2e-4, 2) # sexual mating
    gr_eas_females.mate(True, n_offspr, 1.0, 0.2e-4, 2) # sexual mating
    gr_eur_females.mate(True, n_offspr, 1.0, 0.2e-4, 2) # sexual mating

    off_gr = eg.Group() # creating the group that will hold all offsprings

    gr_afr_females.regroup_newborn(off_gr)
    gr_eas_females.regroup_newborn(off_gr)
    gr_eur_females.regroup_newborn(off_gr)
    
    off_gr.move(mixpop) # finally "physically" moving all offsprings into mix_pop

# ----------------------------------------------------------------------------
def make_admixtured_population_v3( afrpop, easpop, eurpop, mixpop, n_offspr ):
    
    # generating admixtured individuals for the mix_pop: AFR + EAS + EUR -> MIX
    # this function implement very simplified schema compare to v1 and v2
    
    num_selected = math.ceil( afrpop.size()/(3*n_offspr) ) # note, here we need integer value (not a floating point)

    print(" ")
    print("number of females required to be selected from each population is", num_selected)

    gr_afr_males = eg.Group() # group for selected males from afr_pop
    gr_afr_females = eg.Group() # group for selected males from afr_pop
    gr_eas_males = eg.Group() # group for selected males from eas_pop
    gr_eas_females = eg.Group() # group for selected males from eas_pop
    gr_eur_males = eg.Group() # group for selected males from eur_pop
    gr_eur_females = eg.Group() # group for selected males from eur_pop

    gr_afr_females.add(afrpop) # now the group holds the entire population
    gr_afr_females.select("sex", "rand", num_selected, female_sex) # gr_afr_females now left only with required (selected) females 
    gr_afr_males.add(afrpop) # now the group holds the entire population
    gr_afr_males.select("sex", "rand", num_selected, male_sex) # gr_afr_males now left only with required (selected) males 
    
    gr_eas_females.add(easpop)
    gr_eas_females.select("sex", "rand", num_selected, female_sex)
    gr_eas_males.add(easpop)
    gr_eas_males.select("sex", "rand", num_selected, male_sex)
    
    gr_eur_females.add(eurpop)
    gr_eur_females.select("sex", "rand", num_selected, female_sex)
    gr_eur_males.add(eurpop)
    gr_eur_males.select("sex", "rand", num_selected, male_sex)

    print("selected females in afr_pop:", gr_afr_females.size_at(0))
    print("selected males in afr_pop:", gr_afr_males.size_at(0))
    print("selected females in eas_pop:", gr_eas_females.size_at(0))
    print("selected males in eas_pop:", gr_eas_males.size_at(0))
    print("selected females in eur_pop:", gr_eur_females.size_at(0))
    print("selected males in eur_pop:", gr_eur_males.size_at(0))

    mating_gr = eg.Group()

    mating_gr.add(gr_afr_males)
    mating_gr.add(gr_afr_females)
    mating_gr.add(gr_eas_males)
    mating_gr.add(gr_eas_females)
    mating_gr.add(gr_eur_males)
    mating_gr.add(gr_eur_females)

    n_in_group = 0
    for ipop in range( mating_gr.size() ): # loop over registered populations in the group 
        n_in_group += mating_gr.size_at(ipop) # access number of individuals from specific population ipop in the group
    print("the actual number of individuals in the mating_gr group:", n_in_group)
    print(" ")
    
    mating_gr.mate(True, n_offspr, 1.0, 0.2e-4, 2) # sexual mating

    off_gr = eg.Group() # creating the group that will hold all offsprings

    mating_gr.regroup_newborn(off_gr)
    
    off_gr.move(mixpop) # finally "physically" moving all offsprings into mix_pop

# -------------------------------------------------------------------------
def simulate_admixture_pop(num_at_start, limit, n_gen, n_offsprings, create_mix_mode):
    
    # -------------------------
    # this is generation 0
    # -------------------------

    # create empty populations:
    
    afr_pop = eg.Population()
    eas_pop = eg.Population()
    eur_pop = eg.Population()
    mix_pop = eg.Population()

    # generate individuals for each population using haplotypes amd genetic map;
    # without pedigree (parents are unknown); sexes and ids are simulated:
    
    afr_pop.set_population(num_at_start, "struct_haplotypes_pop1.dat", 0.50, 2);
    eas_pop.set_population(num_at_start, "struct_haplotypes_pop1.dat", 0.50, 2);
    eur_pop.set_population(num_at_start, "struct_haplotypes_pop1.dat", 0.50, 2);

    # printing some basic populations' info:
    
    print("population afr_pop is assigned id:", afr_pop.get_popid(), ", and has the initial size of", afr_pop.size(), "individuals.")
    print("population eas_pop is assigned id:", eas_pop.get_popid(), ", and has the initial size of", eas_pop.size(), "individuals.")
    print("population eur_pop is assigned id:", eur_pop.get_popid(), ", and has the initial size of", eur_pop.size(), "individuals.")
    print("population mix_pop is assigned id:", mix_pop.get_popid(), ", and has the initial size of", mix_pop.size(), "individuals.")
    
    # making admixtured population
    
    if create_mix_mode == 1:
        make_admixtured_population_v1( afr_pop, eas_pop, eur_pop, mix_pop, n_offsprings )
    elif create_mix_mode == 2:
        make_admixtured_population_v2( afr_pop, eas_pop, eur_pop, mix_pop, n_offsprings )
    elif create_mix_mode == 3:
        make_admixtured_population_v3( afr_pop, eas_pop, eur_pop, mix_pop, n_offsprings )
    else:
        print("wrong create_mix_mode parameter!")
        return

    print("population afr_pop with id:", afr_pop.get_popid(), ", is of size", afr_pop.size(), "individuals.")
    print("population eas_pop with id:", eas_pop.get_popid(), ", is of size", eas_pop.size(), "individuals.")
    print("population eur_pop with id:", eur_pop.get_popid(), ", is of size", eur_pop.size(), "individuals.")
    print("population mix_pop with id:", mix_pop.get_popid(), ", is of size", mix_pop.size(), "individuals.")

    # print("Calculating and writing to files LD ...")
    # afr_pop.get_ld("mix_ld_gen0", False, 0, 1);

    time = 0;
    time_step = 1
    age_of_death = 50.0
    age_of_young = 1.0

    # -------------------------
    # this are generations 1+
    # -------------------------

    for gen in range( 1, n_gen+1 ):

        print(" ")
        print("==> Generation: ", gen)
        print(" ")

        time += time_step # counting "time", if need for some reason

        # making all individuals in the populations getting older:
        afr_pop.aging( time_step )
        eas_pop.aging( time_step )
        eur_pop.aging( time_step )
        mix_pop.aging( time_step )

        # prepare the working groups (one for each population will be enough)
        gr_afr = eg.Group()
        gr_eas = eg.Group()
        gr_eur = eg.Group()
        gr_mix = eg.Group()
        off_gr = eg.Group()

        # ----- Do the job for AFR population --------------------

        # remove old individuals from the populations

        gr_afr.add(afr_pop) # add entire population to the group
        gr_afr.select("age","best",afr_pop.size(),age_of_death) # select all individuals with age above age_of_death
        if gr_afr.size_at(0) > 0:
            print("the number of individuals above the reproduction age in AFR population:", gr_afr.size_at(0))
            #gr_afr.kill() # kill the selected; actually they are not alive but still can be used for reproduction
            gr_afr.remove() # or remove the selected entirely; they will not be accessible at all
        gr_afr.clear() # make the group ready for re-use within the same scope

        gr_afr.add(afr_pop) # add entire population to the group
        gr_afr.select("age","best",afr_pop.size(),age_of_young) # select all individuals with age above age_ofyoung for mating
        if gr_afr.size_at(0) > 0:
            print("the number of individuals at the reproduction age in AFR population:", gr_afr.size_at(0))
            gr_afr.mate(True, n_offsprings, 1.0, 0.2e-4, 2)
        
        # checking population (if needed):
        # for ind in range(afr_pop.size()):
        #     pos = afr_pop.get_valid_pos(ind)
        #     print("id:", afr_pop.id_at(pos), "sex:", afr_pop.sex_at(pos), "age:", afr_pop.age_at(pos))

        # get the number of offsprings
        gr_afr.regroup_newborn(off_gr)        
        off_num = off_gr.size_at(0)
        off_gr.clear() # because it is going to be reused below
        
        print("the total number of individuals in the  AFR population:", afr_pop.size(), "among them the number of newborn:", off_num)

        # control the population limit
        n_over_limit = afr_pop.size() - limit
        if n_over_limit > 0:
            off_gr.add(afr_pop)
            off_gr.select("id", "rand", n_over_limit, -1)
            print("to be removed from AFR population due to control of population limit:", off_gr.size_at(0))
            off_gr.remove()
            off_gr.clear()
        print(" ")

        # ----- Do the job for EAS population --------------------

        gr_eas.add(eas_pop)
        gr_eas.select("age","best",eas_pop.size(),age_of_death)
        if gr_eas.size_at(0) > 0:
            print("the number of individuals above the reproduction age in EAS population:", gr_eas.size_at(0))
            #gr_eas.kill()
            gr_eas.remove()
        gr_eas.clear()

        gr_eas.add(eas_pop)
        gr_eas.select("age","best",eas_pop.size(),age_of_young)
        if gr_eas.size_at(0) > 0:
            print("the number of individuals at the reproduction age in EAS population:", gr_eas.size_at(0))
            gr_eas.mate(True, n_offsprings, 1.0, 0.2e-4, 2)
        
        # for ind in range(eas_pop.size()):
        #     pos = eas_pop.get_valid_pos(ind)
        #     print("id:", eas_pop.id_at(pos), "sex:", eas_pop.sex_at(pos), "age:", eas_pop.age_at(pos))
        
        gr_eas.regroup_newborn(off_gr)
        off_num = off_gr.size_at(0)
        off_gr.clear()

        print("the total number of individuals in the  EAS population:", eas_pop.size(), "among them the number of newborn:", off_num)

        # control the population limit
        n_over_limit = eas_pop.size() - limit
        if n_over_limit > 0:
            off_gr.add(eas_pop)
            off_gr.select("id", "rand", n_over_limit, -1)
            print("to be removed from EAS population due to control of population limit:", off_gr.size_at(0))
            off_gr.remove()
            off_gr.clear()
        print(" ")

        # ----- Do the job for EUR population --------------------

        gr_eur.add(eur_pop)
        gr_eur.select("age","best",eur_pop.size(),age_of_death)
        if gr_eur.size_at(0) > 0:
            print("the number of individuals above the reproduction age in EUR population:", gr_eur.size_at(0))
            #gr_eur.kill()
            gr_eur.remove()
        gr_eur.clear()

        gr_eur.add(eur_pop)
        gr_eur.select("age","best",eur_pop.size(),age_of_young)
        if gr_eur.size_at(0) > 0:
            print("the number of individuals at the reproduction age in EUR population:", gr_eur.size_at(0))
            gr_eur.mate(True, n_offsprings, 1.0, 0.2e-4, 2)
        
        # for ind in range(eur_pop.size()):
        #     pos = eur_pop.get_valid_pos(ind)
        #     print("id:", eur_pop.id_at(pos), "sex:", eur_pop.sex_at(pos), "age:", eur_pop.age_at(pos))
        
        gr_eur.regroup_newborn(off_gr)
        off_num = off_gr.size_at(0)
        off_gr.clear()

        print("the total number of individuals in the  EUR population:", eur_pop.size(), "among them the number of newborn:", off_num)

        # control the population limit
        n_over_limit = eur_pop.size() - limit
        if n_over_limit > 0:
            off_gr.add(eur_pop)
            off_gr.select("id", "rand", n_over_limit, -1)
            print("to be removed from EUR population due to control of population limit:", off_gr.size_at(0))
            off_gr.remove()
            off_gr.clear()
        print(" ")

        # ----- Do the job for MIX population --------------------

        gr_mix.add(mix_pop)
        gr_mix.select("age","best",mix_pop.size(),age_of_death)
        if gr_mix.size_at(0) > 0:
            print("the number of individuals above the reproduction age in MIX population:", gr_mix.size_at(0))
            #gr_mix.kill()
            gr_mix.remove()
        gr_mix.clear()

        gr_mix.add(mix_pop)
        gr_mix.select("age","best",mix_pop.size(),age_of_young)
        if gr_mix.size_at(0) > 0:
            print("the number of individuals at the reproduction age in MIX population:", gr_mix.size_at(0))
            gr_mix.mate(True, n_offsprings, 1.0, 0.2e-4, 2)
        
        # for ind in range(mix_pop.size()):
        #     pos = mix_pop.get_valid_pos(ind)
        #     print("id:", mix_pop.id_at(pos), "sex:", mix_pop.sex_at(pos), "age:", mix_pop.age_at(pos))
        
        gr_mix.regroup_newborn(off_gr)
        off_num = off_gr.size_at(0)
        off_gr.clear()

        print("the total number of individuals in the  MIX population:", mix_pop.size(), "among them the number of newborn:", off_num)

        # control the population limit
        n_over_limit = mix_pop.size() - limit
        if n_over_limit > 0:
            off_gr.add(mix_pop)
            off_gr.select("id", "rand", n_over_limit, -1)
            print("to be removed from MIX population due to control of population limit:", off_gr.size_at(0))
            off_gr.remove()
            off_gr.clear()
        print(" ")

        # setting up birthdays for newborn:
        afr_pop.set_birthday( time )
        eas_pop.set_birthday( time )
        eur_pop.set_birthday( time )
        mix_pop.set_birthday( time )

        if gen % 100 == 0: # do this every 15-th iteration, the reshaping frequency depends ...
            print("Reshaping ...")
            afr_pop.reshape()
            eas_pop.reshape()
            eur_pop.reshape()
            mix_pop.reshape()
            print("Reshaped afr_pop size = ", afr_pop.size(), "; current capacity = ", afr_pop.capacity() )
            print("Reshaped eas_pop size = ", eas_pop.size(), "; current capacity = ", eas_pop.capacity() )
            print("Reshaped eur_pop size = ", eur_pop.size(), "; current capacity = ", eur_pop.capacity() )

        print("... end of generation", gen)

    print(" ")
    print("establish observations as 3 correlated traits ...")

    # establish observations as 3 correlated traits
    tr_mean = [ 40.0, 5.0, 0.5 ] # (1) trait means
    qtl_prop = [ 1.0 ] # (2) proportion of snps selected as qtls
    cor_g = [ [1.0, 0.5, 0.7], [0.5, 1.0, 0.2], [0.7, 0.2, 1.0] ] # (3) genomic correlations
    cor_e = [ [1.0, 0.3, 0.5], [0.3, 1.0, 0.4], [0.5, 0.4, 1.0] ] # (4) rsidual correlations
    var_g = [ 100.0, 10.0, 0.1 ] # (5) genomic variances
    var_e = [ 200.0, 20.0, 0.3 ] # (6) residual variances
    env = [ 0.0, 0.0, 0.0 ] # (7) enviironment
    # (i) For uniform distribution, model 1
    k_range_U = [ -1.0, 1.0 ] # range of k parameter
    # (ii) For normal distribution, model 2
    k_range_N = [ 0.0, 0.5 ] # mean & std
    # (iii) For gamma distribution, mdodel 3
    k_range_G = [ 0.05 ] # expected value of k in loci; is '1/b' param. in the distribution
    mode = 1
    T = eg.Trait(afr_pop,tr_mean, qtl_prop, cor_g, var_g, cor_e, var_e, env, mode, k_range_U)
    #T.get_observations(afr_pop,env)

    print(" ")
    print("getting observations ...")

    tr_gr = eg.Group()
    tr_gr.add(afr_pop)
    env2 = np.array(env) # needs to convert to numpy array ???
    tr_gr.make_observation(T,env2)

    # print("Calculating and writing to files LD ...")
    # afr_pop.get_ld("afr_ld", False, 0, 1);
    # eas_pop.get_ld("eas_ld", False, 0, 2);
    # eur_pop.get_ld("eur_ld", False, 0, 2);
    # mix_pop.get_ld("mix_ld", False, 0, 1);

    # print("Writing ancestry data to files ...")
    # afr_pop.get_ancestry("afr_ancestry.txt");
    # eas_pop.get_ancestry("eas_ancestry.txt");
    # eur_pop.get_ancestry("eur_ancestry.txt");
    # mix_pop.get_ancestry("mix_ancestry.txt");
    
    # print("Writing haplotypes to files ...")
    # afr_pop.get_haplotypes("afr_haplotypes.txt");
    # eas_pop.get_haplotypes("eas_haplotypes.txt");
    # eur_pop.get_haplotypes("eur_haplotypes.txt");
    # mix_pop.get_haplotypes("mix_haplotypes.txt");

    print("Writing genotypes to files ...")
    afr_pop.get_genotypes("afr_genotypes.txt");
    # eas_pop.get_genotypes("eas_genotypes.txt");
    # eur_pop.get_genotypes("eur_genotypes.txt");
    # mix_pop.get_genotypes("mix_genotypes.txt");

    print("Writing pedigree to files ...")
    afr_pop.get_pedigree("afr_pedigree.txt");
    # eas_pop.get_pedigree("eas_pedigree.txt");
    # eur_pop.get_pedigree("eur_pedigree.txt");
    # mix_pop.get_pedigree("mix_pedigree.txt");

    # print("Writing data to files ...")
    afr_pop.get_data("afr_data.txt");
    # eas_pop.get_data("eas_data.txt");
    # eur_pop.get_data("eur_data.txt");
    # mix_pop.get_data("mix_data.txt");

    print(" ")
    print("setting non genotyped and missing ...")

    # need to set (i) some ids (approx. 30 %) as not genotyped
    # and some ids (approx. 10 %) as missing observations
    gr = eg.Group()
    gr.add(afr_pop)
    gr.select("birth", "best", math.floor(afr_pop.size()*0.4), n_gen )
    #gr.select("age", "worst", math.floor(afr_pop.size()*0.4), 0 )
    print("selected by age:", gr.size_at(0))
    gr.set_not_genotyped()
    gr.clear()
    gr.add(afr_pop)
    gr.select("genot", "rand", math.floor(afr_pop.size()), 1 )
    print("all genotyped:", gr.size_at(0), "total in pop:", afr_pop.size())
    gr.clear()

    gr.add(afr_pop)
    gr.select("id", "rand", math.floor(afr_pop.size()*0.1), -1 )
    gr.set_missing_observations()
    gr.clear()

    afr_pop.get_genotypes("afr_genotypes2.txt");
    afr_pop.get_data("afr_data2.txt");

    print("End of simulation")

# -------------------------------------------------------------------------
def main():

    pop_num_at_start = 10
    generations = 15
    num_offspring = 2
    generate_mixtured_pop_mode = 1
    pop_limit = 50000

    simulate_admixture_pop( pop_num_at_start, pop_limit, generations, num_offspring, generate_mixtured_pop_mode )

# -------------------------------------------------------------------------
if __name__ == '__main__':
    main()