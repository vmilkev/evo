import evolm
import evoped
# --------------------------------------------------------------------
def print_info():
    msg = """
    1) Testing big SNP-SSTEP model on simulated data.
    2) Testing big GBLUP model on the same data
    """
    print(msg) 
# --------------------------------------------------------------------
def main():
    print_info()
    """
    """
    test_impute_genotypes()

    # test_snp_sstep()
    # get_snp_bv()

    # test_snpblup()
    # get_snp_bv_small()

    # test_gblup()

    # test_sstep()
# --------------------------------------------------------------------
def test_impute_genotypes():
    gmat = evoped.Gmatd()
    amat = evoped.Amatd()
    amat.make_matrix("tests/data/big_snp_sstep/data/afr_pedigree.txt", True) # making full A(-1)
    amat.save_matrix("iA", "tests/data/big_snp_sstep/iA.corbin") # saving matrix to .corbin file
    amat.save_ids("iA", "tests/data/big_snp_sstep/ref_ids.txt") # saving traced ids to a file
    amat.clear()
    # input
    snp_file = "tests/data/big_snp_sstep/data/afr_genotypes2.txt"
    id_file = "tests/data/big_snp_sstep/id_geno.txt"
    iA_file = "tests/data/big_snp_sstep/iA.corbin"
    # output
    impt_snp_file = "tests/data/big_snp_sstep/impt_snp"
    freq_file = "tests/data/big_snp_sstep/freqs_big.txt";

    gmat.impute_genotypes(snp_file, id_file, iA_file, impt_snp_file, freq_file)
    gmat.clear()
# --------------------------------------------------------------------
def test_snp_sstep():
    model = evolm.lmm()
    model.define("data = tests/data/big_snp_sstep/data/afr_data2.txt; obs_missing_value = [-9999]")
    model.define("M = tests/data/big_snp_sstep/impt_snp.dmbin")
    model.define("iA = tests/data/big_snp_sstep/iA.corbin")
    model.define("trait_1 ~ 1 + (1|id&M)a + (1|z1)")
    model.define("var = (a)*I*Gg + (z1)*iA*Ga + R; Gg = [0.045]; Ga = [10.0]; R = [200.0]")
    model.solve("pcg", 10, 10, "tests/data/big_snp_sstep/snp_sstep.log", "tests/data/big_snp_sstep/solution_snp_sstep.dat")
# --------------------------------------------------------------------
def get_snp_bv():
    model = evolm.lmm()
    impt_snp_file = "tests/data/big_snp_sstep/impt_snp.dmbin"
    snp_sol_file = "tests/data/big_snp_sstep/solution_snp_sstep.dat"
    snp_bv_big = "tests/data/big_snp_sstep/bv_snpsstep.dat"
    model.snp_to_bv(impt_snp_file, snp_sol_file, snp_bv_big, 2, 5001)
# --------------------------------------------------------------------
def test_sstep():
    hmat = evoped.Hmat()
    model = evolm.lmm()

    snp_file = "tests/data/big_snp_sstep/data/afr_genotypes2.txt"
    ped_file = "tests/data/big_snp_sstep/data/afr_pedigree.txt"
    iH_file = "tests/data/big_snp_sstep/iH.corbin"
    hmat.make_matrix(snp_file, ped_file, iH_file) # making H
    hmat.clear()

    model.define("data = tests/data/big_snp_sstep/data/afr_data2.txt; obs_missing_value = [-9999]")
    model.define("H = tests/data/big_snp_sstep/iH.corbin")
    model.define("refeff = tests/data/big_snp_sstep/ref_ids.txt")
    model.define("trait_1 ~ 1 + (1|id{ref_ids})")
    model.define("var = (id)*H*G2 + R; G2 = [100]; R = [200.0]")
    model.solve("pcg", 10, 5, "tests/data/big_snp_sstep/sstep.log", "tests/data/big_snp_sstep/solution_sstep.dat")
# --------------------------------------------------------------------
def test_gblup():
    gmat = evoped.Gmat()
    model = evolm.lmm()
    snp_file = "tests/data/big_snp_sstep/data/afr_genotypes.txt"
    #id_file = "tests/data/big_snp_sstep/id_geno.dat"
    iG_file = "tests/data/big_snp_sstep/iG"

    gmat.make_matrix(snp_file) # making G
    gmat.scale_diag(0.01)
    gmat.invert_matrix()
    gmat.save_matrix(iG_file) # saving matrix to G_file ??? => .dmbin file
    gmat.clear()

    model.define("data = tests/data/big_snp_sstep/data/afr_data2.txt; obs_missing_value = [-9999]")
    model.define("refeff = tests/data/big_snp_sstep/ref_ids.txt")
    model.define("gcorr = tests/data/big_snp_sstep/iG.dmbin")
    model.define("trait_1 ~ 1 + (1|id)")
    model.define("var = (id)*gcorr*G2 + R; G2 = [100.0]; R = [200.0]")
    model.solve("pcg", 10, 5, "tests/data/big_snp_sstep/gblup.log", "tests/data/big_snp_sstep/solution_gblup.dat")

def test_snpblup():
    gmat = evoped.Gmat()
    model = evolm.lmm()

    snp_file = "tests/data/big_snp_sstep/data/afr_genotypes.txt"
    #id_file = "tests/data/big_snp_sstep/id_geno.dat"
    Z_file = "tests/data/big_snp_sstep/Z"
    freq_file = "tests/data/big_snp_sstep/freqs.txt";

    #gmat.scale_genotypes(snp_file, id_file, Z_file, freq_file) # making and saving Z matrix to .dmbin file
    gmat.scale_genotypes(snp_file, Z_file, freq_file)
    gmat.clear()

    model.define("snp = tests/data/big_snp_sstep/Z.dmbin")
    model.define("data = tests/data/big_snp_sstep/data/afr_data2.txt; obs_missing_value = [-9999]")
    model.define("refeff = tests/data/big_snp_sstep/ref_ids.dat")
    model.define("trait_1 ~ 1 + (1|id&snp)")
    model.define("var = (snp)*I*G + R; G = [0.045]; R = [200.0]") # 3.5241/4297 => sigma_g
    model.solve("pcg", 10, 5,  "tests/data/big_snp_sstep/snpblup.log", "tests/data/big_snp_sstep/solution_snpblup.dat")

def get_snp_bv_small():
    model = evolm.lmm()
    impt_snp_file = "tests/data/big_snp_sstep/Z.dmbin"
    snp_sol_file = "tests/data/big_snp_sstep/solution_snpblup.dat"
    snp_bv = "tests/data/big_snp_sstep/bv_snpblup.dat"
    model.snp_to_bv(impt_snp_file, snp_sol_file, snp_bv, 2, 5001)

# --------------------------------------------------------------------
if __name__ == '__main__':
    main()