## phenotypic selection
##
#Q: is it possible to set a seed? #report
#Q: How can I print tbv and residual of each animal. I only know how to print observations
#Q:report: sex ratio of offsrping
#Q: How to make a loop for selection procedures
#Q: whether we need to shuffle by ourselves: I double checked - It's ok - the sires are shuffled without actions by users, and dams remain the same order. 
#Q: phenotyping on group

import sys
sys.path.append('../../release')

import numpy as np
import evogen
import os
import random
import pandas as pd

#sys.path.append('~/BOA/Simulation/phenotypic/') # cluster
#sys.path.append('release') # Mac

def simple_test_all():
    print(f"=========================================================================")
    print("1.0: Testing Population:")
    os.listdir()    
    pop = evogen.Population() #Construct an empty Population object. This is a default class constructor with no parameters
                              #15 individuals 
    pop.set_population("DNK_hap", "temp_map", False); # haplotype, structure, pedigree;  false: only haplotypes, true:hap with pedigree and sex    
    print("1.0: Testing Trait:")

    tr_mean = [ 0] # trait means
    qtl_prop = [ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1] # (3) proportion of snps selected as qtls
    var_g = [ 1.0  ] # genomic variances
    var_e = [ 2.0  ] # residual variances
    cor_g = [[0.01]] # doesnt work if it is 0.0
    cor_e = [[0.01]]    
    env = [ 0.0 ] # enviironment

    # (i) For uniform distribution, model 1
    k_range_U = [ -1.0, 0.2] # range of k parameter
    # (ii) For normal distribution, model 2
    k_range_N = [ 0.0, 0.2 ] # mean & std
    # (iii) For gamma distribution, mdodel 3
    k_range_G = [ 0.05 ] # expected value of k in loci; is '1/b' param. in the distribution

    which_model = 1 # uniform distribution

    #Trait values    
    print("3.0: Setting trait:")
    T = evogen.Trait(pop, tr_mean, qtl_prop, cor_g, var_g, cor_e,var_e, env, which_model, k_range_U)
    print(f"=========================================================================")
    print("4.0: Make observations on this single trait:")
    T.get_observations(pop, env, "T_trait_pop.dat") #printed the 15 observations; the observations are written in file "T_trait_pop.dat"

    # Print and write sex, pedigree, and phenotypes
    print("5.0: Writing sex, pedigree, and phenotypes to 'population_data.txt':")
    with open("population_data.txt", "w") as file:
        file.write("Individual\tSex\tSire\tDame\tPhenotype\n")
        for i in range(pop.size()):
            sex = pop.sex_at(i)
            sire = pop.sire_at(i)
            dame = pop.dame_at(i)
            id = pop.id_at(i)
            phenotype=float(pop.phenotype_at(i).item()) 
            #phenotype = pop.phenotype_at(i)#[0]  # Assuming a single trait # start with trying with more than 1 trait without extra brackets
            file.write(f"{id}\t{sex}\t{sire}\t{dame}\t{phenotype}\n")
            
    print(f"=========================================================================")
    #Extract all relevant data into a list
    selected_parents = []
    for i in range(pop.size()):
        id = pop.id_at(i)
        sex = pop.sex_at(i)
        sire = pop.sire_at(i)
        dame = pop.dame_at(i)
        phenotype = float(pop.phenotype_at(i).item())  # Safely convert to scalar
#        phenotype = pop.phenotype_at(i)[0]  # Assuming a single trait
        selected_parents.append([i, id, sex, sire, dame, float(phenotype)]) 
        
    # Convert the data into a pandas DataFrame for easier manipulation
    df = pd.DataFrame(selected_parents, columns=["Index", "ID", "Sex", "Sire", "Dame", "Phenotype"])
    df["Phenotype"] = pd.to_numeric(df["Phenotype"], errors='coerce')
    print(f"Full data:\n{df}")    
    print(f"=========================================================================")    

    # Get top 5 phenotypes for Sex = 0 and Sex = 1
    top_5_sires = df[df["Sex"] == 1].nlargest(5, "Phenotype") #Filter the DataFrame to include only rows where the "Sex" column has a value of 0
    top_5_sires_r = df[df["Sex"] == 1].nlargest(5, "Phenotype").sample(frac=1).reset_index(drop=True)                                                           # Select the top 5 rows (based on the highest values) from the filtered DataFrame, using the "Phenotype" column for ranking.
    top_5_dams = df[df["Sex"] == 0].nlargest(5, "Phenotype")
    top_5_dams_r = df[df["Sex"] == 0].nlargest(5, "Phenotype").sample(frac=1).reset_index(drop=True)
    print(f"=========================================================================")
    # Not needed
    top_5_sires_r_repeated = top_5_sires_r.loc[np.repeat(top_5_sires_r.index, 2)].reset_index(drop=True)
    top_5_dams_r_repeated = top_5_dams_r.loc[np.repeat(top_5_dams_r.index, 2)].reset_index(drop=True)    

    print(f"Top 5 sires:\n{top_5_sires}")
    print(f"=========================================================================")
    print(f"Top 5 dams:\n{top_5_dams}")
    print(f"=========================================================================")
    print(f"Top 5 sires random mating :\n{top_5_sires_r}")
    print(f"=========================================================================")
    print(f"Top 5 dams random mating :\n{top_5_dams_r}")
    print(f"=========================================================================")

    print(f"Top 5 sires random mating :\n{top_5_sires_r_repeated}")
    print(f"=========================================================================")
    print(f"Top 5 dams random mating :\n{top_5_dams_r_repeated}")
    print(f"=========================================================================")
    
    # Combine the indices of the selected individuals
    selected_indices = pd.concat([top_5_sires, top_5_dams])["Index"] # combine the two data into a single data
    print(f"Print selected_indices:\n{selected_indices}")
    print(f"=========================================================================")
    
    # Create G2 by adding selected individuals
    G2 = evogen.Group()  # Question: How to print the data
    
    for idx in selected_indices:
        #Find the row in the DataFrame df where the "Index" column matches the current idx
        #Select the "ID" column from that row.
        #Extract the first value
        print(f"Adding individual {df.loc[df['Index'] == idx, 'ID'].values[0]} to G2")
        print(idx)
        G2.add(pop, idx)
    # Verify the size of G2
    print(f"G2 size: {G2.size_at(0)}")

    pop.aging(2)
    #G2.mate(True, 2, 0.8)
    G2.mate(True, 2, 1)

    G_newborn = evogen.Group()
    G2.regroup_newborn(G_newborn)

    print(f"G_newborn.size: {G_newborn.size_at(0)}")
    G_newborn.aging(1)
    G_newborn.genotype()

    T.get_observations(pop, env) #printed the 15 observations; the observations are written in file "T_trait_pop.dat"
    selected_parents = []
    for i in range(pop.size()):
        id = pop.id_at(i)
        sex = pop.sex_at(i)
        sire = pop.sire_at(i)
        dame = pop.dame_at(i)
        age = pop.age_at(i)
        phenotype = float(pop.phenotype_at(i).item())  # convert to scalar
#        phenotype = pop.phenotype_at(i)[0]  # Assuming a single trait
        selected_parents.append([i, id, sex, sire, dame, age, phenotype])  
        
    # Convert the data into a pandas DataFrame for easier manipulation
    df = pd.DataFrame(selected_parents, columns=["Index", "ID", "Sex", "Sire", "Dam","age" ,"Phenotype"])
    print(df.dtypes)  # Check column data types
    print(df["Phenotype"].head())  # Print first few values
    print(f"After First round of selection:\n{df}")



    print(f"=========================================================================")
    print("Selecting sires and dams for the next generation:")
    print(f"=========================================================================")    

    # Extract newborn data into DataFrame
    # Select top 5 males (sires) and top 5 females (dams) based on highest phenotype
    top_5_sires_r2 = df[(df["Sex"] == 1) & (df["age"] == 1)].nlargest(5, "Phenotype").sample(frac=1).reset_index(drop=True)
    top_5_dams_r2 = df[(df["Sex"] == 0) & (df["age"] == 1)].nlargest(5, "Phenotype").sample(frac=1).reset_index(drop=True)
    
    print(f"Top 5 selected sires for next generation (shuffled):\n{top_5_sires_r2}")
    print(f"Top 5 selected dams for next generation (shuffled):\n{top_5_dams_r2}")

    # Combine the indices of the selected individuals
    selected_indices = pd.concat([top_5_sires_r2, top_5_dams_r2])["Index"] # combine the two data into a single data
    print(f"Print selected_indices:\n{selected_indices}")
    print(f"=========================================================================")
    G3 = evogen.Group()  # Q: How to print the data
    
    for idx in selected_indices:
        print(f"Adding individual {df.loc[df['Index'] == idx, 'ID'].values[0]} to G2")
        print(idx)
        G3.add(pop, idx)
    print(f"G3 size: {G3.size_at(0)}")

    pop.aging(1)
    G3.mate(True, 2, 1)

    G_newborn = evogen.Group()
    G3.regroup_newborn(G_newborn)

    print(f"G_newborn.size: {G_newborn.size_at(0)}")
    G_newborn.aging(1)
    G_newborn.genotype()

    T.get_observations(pop, env) #printed the 15 observations; the observations are written in file "T_trait_pop.dat"
    selected_parents = []
    for i in range(pop.size()):
        id = pop.id_at(i)
        sex = pop.sex_at(i)
        sire = pop.sire_at(i)
        dame = pop.dame_at(i)
        age = pop.age_at(i)
        phenotype = float(pop.phenotype_at(i).item())  
        selected_parents.append([i, id, sex, sire, dame, age, phenotype])          
    df = pd.DataFrame(selected_parents, columns=["Index", "ID", "Sex", "Sire", "Dam","age" ,"Phenotype"])
    print(f"After Second round of selection:\n{df}")
    
if __name__ == '__main__':      
#def main():
	simple_test_all()
    
        
    # aging increase the age of the population
    # reshape delete unessary animals



    
