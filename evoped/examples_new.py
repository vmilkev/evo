import evoped
import time

# some data:
snp_file = "../evolm/tests/data/large_data/SNPS/YY_plink";
#snp_file = "YY_plink";
ped_file = "../evolm/tests/data/large_data/DMU/data/dmu_pedigree_yy_20240125.txt";
#ped_file = "dmu_pedigree_yy_20240125.txt";
gid_file = "../evolm/tests/data/large_data/DMU/YY/yy.grm.id";

# ----------------------------------
def example_Amat():

    amat = evoped.Amat()

    #help(amat)
    
    start  = time.time()

    # making specific matrices:

    print("making A22(-1) and reduced A(-1) ...")
    amat.make_matrix_forgenotyped(ped_file, gid_file, True) # making A22(-1) and reduced A(-1)
    amat.save_matrix("iA22", "iA22") # writing A22(-1) to .corbin file
    amat.save_matrix("irA", "irA") # writing reduced A(-1) to .corbin file
    # because we going to make reduced A(-1) one more time, it is better to call this
    amat.clear() # though, this is not absolutelly necessary here

    print("making full A(-1) ...")
    amat.make_matrix(ped_file, True) # making full A(-1)
    amat.save_matrix("iA", "iA") # writing full A(-1) to .corbin file

    print("making (the direct way) reduced A(-1) ...")
    amat.make_matrix(ped_file, gid_file, True) # making reduced A(-1)
    amat.save_matrix("irA", "irA2") # writing reduced A(-1) to .corbin file

    amat.clear()

    end = time.time()

    print("Time of making specific matrices using Amat class: ", end - start, " sec.")

    # making all matrices required for single-step:

    start  = time.time()

    print("making full A(-1), reduced A(-1), iA22 and A22 ...")
    amat.make_all(ped_file, gid_file)

    print("saving all matrices to .corbin files ...")
    amat.save_matrix("iA", "iA_v2")
    amat.save_matrix("irA", "irA_v2")
    amat.save_matrix("iA22", "iA22_v2")
    amat.save_matrix("A22", "A22_v2")

    end = time.time()

    print("Time of making all matrices using Amat class: ", end - start, " sec.")

# ----------------------------------
def example_Gmat():

    gmat = evoped.Gmat()
    amat = evoped.Amat()

    #help(gmat)

    start  = time.time()

    gmat.make_matrix("tests/data/allele.dat")
    gmat.save_ids("g_ids")

    amat.make_all("tests/data/ped_bkg.dat", "g_ids")
    amat.save_matrix("A22", "A22")

    gmat.scale_matrix("A22", 0.25)
    gmat.invert_matrix(True)
    gmat.save_matrix("G")

    end = time.time()
    print("Time of making Gmatrix using Gmat class: ", end - start, " sec.")

# ----------------------------------
def example_Hmat():

    hmat = evoped.Hmat()
    
    #help(hmat)

    g_matr_file = "tests/data/sstep_050/data/gmat_050.gmatrix"
    ped_file2 = "tests/data/sstep_050/data/id4trace.PED"

    start  = time.time()

    hmat.make_matrix(g_matr_file, ped_file2, "h_matrix")
    
    end = time.time()
    print("Time of making Hmatrix using Hmat class: ", end - start, " sec.")

# ----------------------------------    
def main():

    #example_Amat()
    print("completed: example_Amat")

    #example_Gmat()
    print("completed: example_Gmat")

    example_Hmat()
    print("completed: example_Hmat")

# ----------------------------------	
if __name__ == '__main__':
    main()