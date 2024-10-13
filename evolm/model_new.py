import evolm

def example_7_1_pp_110():
    model = evolm.lmm()
    
    model.define("weight ~ herd + pen + (1|id) + (1|dam2{id})b + (1|dam); data = tests/data/model_2/model_2_data.dat; obs_missing_value = [-9.0]")
    model.define("var = (id, b)*A1*G1 + (dam)*I*G2 + R; A1 = tests/data/model_2/A1.corbin; G1 = [150 -40, -40 90]; G2 = [40]; R = [350]")
        
    model.solve("pcg", 10, 5, "example_7_1.log", "solution_7_1.dat")

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

def main():
    example_7_1_pp_110()
    print("completed: example_7_1_pp_110")

    free_example()
    print("completed: free_example")
	
if __name__ == '__main__':
    main()