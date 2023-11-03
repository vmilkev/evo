import sys

sys.path.append('release')

import evolm

def test_model_6():

        # one trait setup
        # EXAMPLE: 11.2. p.183.
        # y = Xb + Za + e

        # corr: [a], identity matrix

        #----------- DATA ----------------------------
        x = [
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1 ]

        z = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              -1, 1, -1, 1, -1, 1, -1, 1, -1, 1,
              1, -1, -1, 1, 1, 1, 1, -1, -1, 1,
              -1, -1, 1, 1, -1, 1, -1, -1, 1, 1,
              1, 1, 1, -1, -1, 1, 1, 1, 1, -1,
              -1, 1, -1, -1, 1, 1, -1, 1, -1, -1,
              1, -1, -1, -1, -1, 1, 1, -1, -1, -1,
              -1, -1, 1, -1, 1, 1, -1, -1, 1, -1,
              1, 1, 1, 1, 1, -1, -1, -1, -1, -1,
              -1, 1, -1, 1, -1, -1, 1, -1, 1, -1,
              1, -1, -1, 1, 1, -1, -1, 1, 1, -1,
              -1, -1, 1, 1, -1, -1, 1, 1, -1, -1,
              1, 1, 1, -1, -1, -1, -1, -1, -1, 1,
              -1, 1, -1, -1, 1, -1, 1, -1, 1, 1,
              1, -1, -1, -1, -1, -1, -1, 1, 1, 1,
              -1, -1, 1, -1, 1, -1, 1, 1, -1, 1
                  ]
        
        corZ = [
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ]

        y = [18.86, 21.43, 20.00, 66.03, 21.43, 21.20, 10.26, 27.26, 29.60, 22.63, 21.40, 14.66, 15.43, 11.90, 7.03, 17.50]

        iR = [0.0041]

        iG1 = [0.057]

        #---------------------------------------------

        model = evolm.Model()
        solver = evolm.Pcg()

        #----- define the model -------

        model.append_residual(iR, 1)

        model.append_observation(y, 16) #obs := 0

        model.append_effect(x, 16, 1)  #eff := 0
        model.append_effect(z, 16, 10) #eff := 1

        corr_eff = [1]

        eff_trate = [0, 1]
        obs_trate = 0

        model.append_corrstruct(iG1, 1, corZ, 10, corr_eff)

        model.append_traitstruct(obs_trate, eff_trate)

        solver.append_model(model)

        #model.print()

        #----------------------------------

        solver.solve()

        sol = solver.get_solution()
        
        solver.get_solution( "py_out_solution_model_6.dat" )

        model.clear()
	
if __name__ == '__main__':	
#def main():
	test_model_6()

