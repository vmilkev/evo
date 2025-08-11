import sys
sys.path.append('../../../release')
import evogen

# Constants
mutation_rate =1e-10
littersize = 4

# Instantiate
pop1                = evogen.Population()
pop2                = evogen.Population()
popmix              = evogen.Population()

gr_pop1_males       = evogen.Group() 
gr_pop2_males       = evogen.Group() 
gr_pop1_females     = evogen.Group() 
gr_pop2_females     = evogen.Group() 
gr_mating           = evogen.Group()
gr_offspring        = evogen.Group()

print("passed 1")

# Populate homogenous populations
## Pop1
pop1.set_population (20, 'viktor_map.csv', 0.5, 2)
gr_pop1_males.add(pop1) 
gr_pop1_females.add(pop1) 
gr_pop1_males.select("sex", "rand", 4, 1) 
gr_pop1_females.select("sex", "rand", 4, 0)

print("passed 2")

## Pop2
pop2.set_population (20, 'viktor_map.csv', 0.5, 2)
gr_pop2_males.add(pop2) 
gr_pop2_females.add(pop2) 
gr_pop2_males.select("sex", "rand", 4, 1) 
gr_pop2_females.select("sex", "rand", 4, 0)

print("passed 3")

## Mating group
gr_mating.add(gr_pop1_males)
gr_mating.add(gr_pop1_females)
gr_mating.add(gr_pop2_males)
gr_mating.add(gr_pop2_females)

print("passed 4")

# Mate
gr_mating.mate(True, littersize, 1.0, mutation_rate, 2)

print("passed 5")

gr_mating.regroup_newborn(gr_offspring)
gr_offspring.move(popmix)
gr_offspring.clear()
gr_mating.clear()
    
# Pre-Grow print 
print("--- Before population growth: ---")
print("")
print("Population:")
print("popmix      ", popmix.size())

print("")
print("Groups (n populations):")
print("Males   pop1   ", gr_pop2_males.size())
print("Females pop1   ", gr_pop2_females.size())
print("Males   pop2   ", gr_pop1_males.size())
print("Females pop2   ", gr_pop1_females.size())
print("Mating         ", gr_mating.size())
print("Offspring      ", gr_offspring.size())

print("")
print("Groups (n individuals):")
print("Males   pop1", gr_pop1_males.size_at(0))
print("Females pop1", gr_pop1_females.size_at(0))
print("Males   pop2", gr_pop2_males.size_at(0))
print("Females pop2", gr_pop2_females.size_at(0))
print("")

# Grow
print("Round", "Class\t", "Instance", "Stage\t", "Variable", "Amount", "Note", sep = "\t")
iter = 0
while popmix.size() < 40:
    iter = iter + 1
    ## Males
    gr_mating.add(gr_pop1_males)
    gr_mating.add(gr_pop2_males)
    gr_mating.add(gr_pop1_females)
    gr_mating.add(gr_pop2_females)
    gr_mating.mate(True, littersize, 1.0, mutation_rate, 2)
    gr_mating.regroup_newborn(gr_offspring)

    # Mid print
    print(iter, "group\t", "mating\t", "preselection", "populations", gr_mating.size(), sep = "\t")
    print(iter, "group\t", "mating\t", "preselection", "individuals", gr_mating.size_at(0), "pop1", sep = "\t")
    print(iter, "group\t", "mating\t", "preselection", "individuals", gr_mating.size_at(1), "pop2", sep = "\t")
    print(iter, "group\t", "offspring", "preselection", "populations", gr_offspring.size(), "", sep = "\t")
    print(iter, "group\t", "offspring", "preselection", "individuals", gr_offspring.size_at(0), "pop1", sep = "\t")
    print(iter, "group\t", "offspring", "preselection", "individuals", gr_offspring.size_at(1), "pop2", sep = "\t")
    gr_offspring.select("id", "rand", 1, 0)
    print(iter, "group\t", "offspring", "postselection", "individuals", gr_offspring.size_at(0), "pop1", sep = "\t")
    print(iter, "group\t", "offspring", "postselection", "individuals", gr_offspring.size_at(1), "pop2", sep = "\t")
    print(iter, "population", "popmix\t", "postselection", "individuals", popmix.size(), "", sep = "\t")
    gr_offspring.move(popmix)
    
    # Clear for next round
    gr_offspring.clear()
    gr_mating.clear()
    print("")

