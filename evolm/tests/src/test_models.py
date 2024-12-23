import evolm
import evoped

def print_basic_info():
    msg = """
    All 'numbered' examples were taken from:
    Linear Models for the Prediction of Animal Breeding Values,
    3rd Edition, by Raphael A. Mrode and Robin Thompson.
    """
    print(msg) 

def example_4_1_p_62(): # section 4.2 Repeatability model
    amat = evoped.Amat()
    model = evolm.lmm()

    amat.make_matrix("tests/data/section_4/pedigree_4_1.dat", True) # making full A(-1)
    amat.save_matrix("iA", "tests/data/section_4/iA_4_1") # saving matrix to .corbin file
    amat.clear()

    model.define( "effect = tests/data/section_4/ref_id.dat; data = tests/data/section_4/data_4_1.dat" );
    model.define( "fat ~ parity + hys + (1|cow2) + (1|cow{ref_id})a" )
    model.define( "var = (a)*A*G + (cow2)*I*Pe + R; A = tests/data/section_4/iA_4_1.corbin; G = [20]; Pe = [12]; R = [28]" )

    model.solve("pcg", 10, 5, "example_4_1.log", "solution_4_1.dat")

def example_4_2_p_67(): # section 4.3 Model with common environmental effects
    amat = evoped.Amat()
    model = evolm.lmm()

    amat.make_matrix("tests/data/section_4/pedigree_4_2.dat", True) # making full A(-1)
    amat.save_matrix("iA", "tests/data/section_4/iA_4_2") # saving matrix to .corbin file
    amat.clear()

    model.define( "effect = tests/data/section_4/ref_id2.dat; data = tests/data/section_4/data_4_2.dat" );
    model.define( "ww ~ sex + (1|piglet{ref_id2}) + (1|dam)" )
    model.define( "var = (piglet)*A*G + (dam)*I*Pe + R; A = tests/data/section_4/iA_4_2.corbin; G = [20]; Pe = [15]; R = [65]]" )

    model.solve("pcg", 10, 5, "example_4_2.log", "solution_4_2.dat")

def example_5_1_p_72(): # section 5.2 Equal design matrices and no missing records
    amat = evoped.Amat()
    model = evolm.lmm()

    amat.make_matrix("tests/data/section_5/pedigree_5_1.dat", True) # making full A(-1)
    amat.save_matrix("iA", "tests/data/section_5/iA_5_1") # saving matrix to .corbin file
    amat.clear()

    model.define( "effect = tests/data/section_5/ref_id.dat; data = tests/data/section_5/data_5_1.dat" );
    model.define( "wwg, pwg ~ sex + (1|calves{ref_id})a" )
    model.define( "var = (a)*A*G + R; A = tests/data/section_5/iA_5_1.corbin; G = [20 18, 18 40]; R = [40 11, 11 30]" )

    model.solve("pcg", 10, 5, "example_5_1.log", "solution_5_1.dat")

def example_5_2_p_78(): # section 5.3 Equal design matrices with missing records
    amat = evoped.Amat()
    model = evolm.lmm()

    amat.make_matrix("tests/data/section_5/pedigree_5_2.dat", True) # making full A(-1)
    amat.save_matrix("iA", "tests/data/section_5/iA_5_2") # saving matrix to .corbin file
    amat.clear()

    model.define( "effect = tests/data/section_5/ref_id2.dat; data = tests/data/section_5/data_5_2.dat; obs_missing_value = [-9.0]" );
    model.define( "wwg ~ sex + (1|calves{ref_id})a" )
    model.define( "pwg ~ sex + (1|calves{ref_id})a" )
    model.define( "var = (a)*A*G + R; A = tests/data/section_5/iA_5_2.corbin; G = [20 18, 18 40]; R = [40 11, 11 30]" )

    model.solve("pcg", 10, 5, "example_5_2.log", "solution_5_2.dat")

def example_5_3_p_80(): # section 5.4 Unequal design matrices
    amat = evoped.Amat()
    model = evolm.lmm()

    amat.make_matrix("tests/data/section_5/pedigree_5_3.dat", True) # making full A(-1)
    amat.save_matrix("iA", "tests/data/section_5/iA_5_3") # saving matrix to .corbin file
    amat.clear()

    model.define( "effect = tests/data/section_5/ref_id.dat; data = tests/data/section_5/data_5_3.dat" );
    model.define( "fat1 ~ hys1 + (1|cow{ref_id})" )
    model.define( "fat2 ~ hys2 + (1|cow{ref_id})" )
    model.define( "var = (cow)*A*G + R; A = tests/data/section_5/iA_5_3.corbin; G = [35 28, 28 30]; R = [65 27, 27 70]" )

    model.solve("pcg", 10, 5, "example_5_3.log", "solution_5_3.dat")

def example_5_4_p_85(): # section 5.5 Multivariate models with no environmental covariance
    amat = evoped.Amat()
    model = evolm.lmm()

    amat.make_matrix("tests/data/section_5/pedigree_5_4.dat", True) # making full A(-1)
    amat.save_matrix("iA", "tests/data/section_5/iA_5_4") # saving matrix to .corbin file
    amat.clear()

    model.define( "effect = tests/data/section_5/ref_id3.dat; data = tests/data/section_5/data_5_4.dat; obs_missing_value = [-999]" );
    model.define( "yw ~ hys + (1|calf{ref_id})" )
    model.define( "fy ~ hys + (1|calf{ref_id})" )
    model.define( "var = (calf)*A*G + R; A = tests/data/section_5/iA_5_4.corbin; G = [43 18, 18 30]; R = [77 0, 0 70]" )

    model.solve("pcg", 10, 5, "example_5_4.log", "solution_5_4.dat")

def example_7_1_p_110(): # section 7.2 Animal model for maternal trait
    
    amat = evoped.Amat()
    model = evolm.lmm()

    amat.make_matrix("tests/data/section_7/pedigree.dat", True) # making full A(-1)
    amat.save_matrix("iA", "tests/data/section_7/ped_corr") # saving matrix to .corbin file
    amat.clear()
    
    model.define( "effect = tests/data/section_7/ref_id.dat; data = tests/data/section_7/model_2_data_reduced.dat; obs_missing_value = [-9.0]" );
    model.define( "weight ~ herd + pen + (1|id{ref_id})a + (1|dam2{ref_id})b + (1|dam)" )
    model.define( "var = (a, b)*A1*G1 + (dam)*I*G2 + R; A1 = tests/data/section_7/ped_corr.corbin; G1 = [150 -40, -40 90]; G2 = [40]; R = [350]" )
        
    model.solve("pcg", 10, 5, "example_7_1.log", "solution_7_1.dat")

def example_8_1_p_125(): # section 8.2 Animal model with social interaction effects
    
    amat = evoped.Amat()
    model = evolm.lmm()

    amat.make_matrix("tests/data/section_8/pedigree_8_1.dat", True) # making full A(-1)
    amat.save_matrix("iA", "tests/data/section_8/iA_8_1") # saving matrix to .corbin file
    amat.clear()
    
    model.define( "effect = tests/data/section_8/ref_id.dat; data = tests/data/section_8/data_8_1.dat" );
    model.define( "S = tests/data/section_8/S.dat" )
    model.define( "gr ~ sex + (1|animal{ref_id}) + (1|pen) + (1|dam) + (1|S)" )
    model.define( "var = (animal, S)*A1*G2 + (dam)*I*G3 + (pen)*I*G4 + R; A1 = tests/data/section_8/iA_8_1.corbin; G4 = [12.12]; G3 = [12.5]; G2 = [25.7 2.25, 2.25 3.6]; R = [48.48]" )
        
    model.solve("pcg", 10, 5, "example_8_1.log", "solution_8_1.dat")

def example_11_1_p_180(): # section 11.4 Fixed effect model for SNP effects
    amat = evoped.Amat()
    model = evolm.lmm()

def free_example():
    model = evolm.lmm()
        
    model.define("data = tests/data/model_2/model_2_data_reduced.dat; obs_missing_value = [-9.0]")
    model.define("G5 = [1 0, 0 1]; G4 = [5 2 2 2, 2 5 2 2, 2 2 5 2, 2 2 2 5]")
    model.define("G6 = [3 2 1, 2 3 1, 1 1 3]")
    model.define("G1 = [150 -40, -40 90]; G2 = [40]; R = [350]")

    model.define("weight ~ 1 + herd*pen + (f|id)c + (pen|dam{id})d")
    model.define("     f ~ 1 + (1|weight:herd)a + (weight|dam2{id})b")
    model.define("    f2 ~ 1 + (1|weight:herd)a + (f|dam2{id})e + (weight|id)k + (f|id)l")
    
    model.define("var = (c)*A1*G2 + (a)*I*G1 + (b,e)*A1*G1 + (k,l)*A1*G1 + G6; A1 = tests/data/model_2/A1.corbin")
    
    model.solve("pcg", 10, 5, "free_example.log", "solution_free.dat")

def is_solution_correct( sol_fname, corr_sol_fname, txt_msg):
    import csv
    curr_sol = []
    true_sol = []
    with open(sol_fname, "r") as f: # reading obtained solution to curr_sol
        reader = csv.reader(f, delimiter=',')
        next(f)
        for row in reader:
            curr_sol.append(row[4])
        f.close()
    with open(corr_sol_fname, "r") as f: # reading true solution to true_sol
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            true_sol.append(row[0])
        f.close()
    result = True
    if len(curr_sol) != len(true_sol): # checking the sizes of solutions lists
        print("len(curr_sol) != len(true_sol)")
        result = False
    for i in range(0,len(curr_sol)): # element-wise checking of solutions lists
        if curr_sol[i] != true_sol[i]:
            print(curr_sol[i], " != ", true_sol[i])
            result = False
    if result:
        print(txt_msg, "PASSED")
    else:
        print(txt_msg, "FAILED")

def main():
    print_basic_info()

    """
    """
    example_4_1_p_62()
    is_solution_correct("solution_4_1.dat", "tests/data/section_4/correct_solution_4_1.dat", "model 4_1")

    example_4_2_p_67()
    is_solution_correct("solution_4_2.dat", "tests/data/section_4/correct_solution_4_2.dat", "model 4_2")

    example_5_1_p_72()
    is_solution_correct("solution_5_1.dat", "tests/data/section_5/correct_solution_5_1.dat", "model 5_1")
    
    example_5_2_p_78()
    is_solution_correct("solution_5_2.dat", "tests/data/section_5/correct_solution_5_2.dat", "model 5_2")

    example_5_3_p_80()
    is_solution_correct("solution_5_3.dat", "tests/data/section_5/correct_solution_5_3.dat", "model 5_3")

    example_5_4_p_85()
    is_solution_correct("solution_5_4.dat", "tests/data/section_5/correct_solution_5_4.dat", "model 5_4")

    example_7_1_p_110()
    is_solution_correct("solution_7_1.dat", "tests/data/section_7/correct_solution_7_1.dat", "model 7_1")
    
    example_8_1_p_125()
    is_solution_correct("solution_8_1.dat", "tests/data/section_8/correct_solution_8_1.dat", "model 8_1")

    example_11_1_p_180()
    #is_solution_correct("solution_11_1.dat", "tests/data/section_11/correct_solution_11_1.dat", "model 11_1")

    #free_example()

if __name__ == '__main__':
    main()