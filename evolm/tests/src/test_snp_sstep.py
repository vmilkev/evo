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
    #test_impute_genotypes()
    test_snp_sstep()
    #test_gblup()
    #test_sstep()
# --------------------------------------------------------------------
def test_impute_genotypes():
    gmat = evoped.Gmat()
    amat = evoped.Amat()
    amat.make_matrix("tests/data/big_snp_sstep/pedigree.dat", True) # making full A(-1)
    amat.save_matrix("iA", "tests/data/big_snp_sstep/iA.corbin") # saving matrix to .corbin file
    amat.clear()

    snp_file = "tests/data/big_snp_sstep/genotype_withoutID.dat"
    id_file = "tests/data/big_snp_sstep/id_geno.dat"
    iA_file = "tests/data/big_snp_sstep/iA.corbin"
    impt_snp_file = "tests/data/big_snp_sstep/impt_snp"

    gmat.impute_genotypes(snp_file, id_file, iA_file, impt_snp_file)
    gmat.clear()
# --------------------------------------------------------------------
def test_snp_sstep():
    model = evolm.lmm()
    model.define("data = tests/data/big_snp_sstep/phenotype.dat; obs_missing_value = [-9999]")
    model.define("M = tests/data/big_snp_sstep/impt_snp.dmbin")
    model.define("iA = tests/data/big_snp_sstep/iA.corbin")
    model.define("obs ~ 1 + (1|id&M)a + (1|z1)")
    model.define("var = (a)*I*Ga + (z1)*iA*Gg + R; Ga = [2.0]; Gg = [1.52]; R = [24.5]")
    model.solve("pcg", 10, 5, "example_big_snp_sstep.log", "solution_big_snp_sstep.dat")
# --------------------------------------------------------------------
def test_sstep():

    hmat = evoped.Hmat()
    model = evolm.lmm()

    snp_file = "tests/data/big_snp_sstep/genotype_withoutID.dat"
    id_file = "tests/data/big_snp_sstep/id_geno.dat"
    ped_file = "tests/data/big_snp_sstep/pedigree.dat"
    iH_file = "tests/data/big_snp_sstep/iH.corbin"

    hmat.make_matrix(snp_file, id_file, ped_file, iH_file) # making H
    hmat.clear()

    model.define("data = tests/data/big_snp_sstep/phenotype.dat; obs_missing_value = [-9999]")
    model.define("H = tests/data/big_snp_sstep/iH.corbin")
    model.define("obs ~ 1 + (1|id)")
    model.define("var = (id)*H*G2 + R; G2 = [3.5241]; R = [24.5]")
    model.solve("pcg", 10, 5, "big_sstep.log", "solution_big_sstep.dat")
# --------------------------------------------------------------------
def test_gblup():
    gmat = evoped.Gmat()
    model = evolm.lmm()

    snp_file = "tests/data/big_snp_sstep/genotype_withoutID.dat"
    id_file = "tests/data/big_snp_sstep/id_geno.dat"
    iG_file = "tests/data/big_snp_sstep/iG"

    gmat.make_matrix(snp_file, id_file) # making G
    gmat.scale_diag(0.01)
    gmat.invert_matrix()
    gmat.save_matrix(iG_file) # saving matrix to G_file ??? => .dmbin file
    gmat.clear()

    model.define("data = tests/data/big_snp_sstep/phenotype.dat; obs_missing_value = [-9999]")
    model.define("refeff = tests/data/big_snp_sstep/ref_id.dat")
    model.define("gcorr = tests/data/big_snp_sstep/iG.dmbin")
    model.define("obs ~ 1 + (1|id{ref_id})")
    model.define("var = (id)*gcorr*G2 + R; G2 = [3.5241]; R = [24.5]")
    model.solve("pcg", 10, 5, "big_gblup.log", "solution_big_gblup.dat")

# --------------------------------------------------------------------
if __name__ == '__main__':
    main()