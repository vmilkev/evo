#include "model.hpp"
#include "solver_pcg.hpp"
#include "iointerface.hpp"
#include "parser.hpp"
#include "effects.hpp"

/*#include <stdio.h>
#include <time.h>
#include <string>
#include <functional>
#include <map>*/

int main(void)
{
    try
    {
        // -------------------------
        std::cout << "Parser:"
                  << "\n";

        char expstr[80];

        evolm::Parser<double> p;

        std::cout << "Parseer of the float type. To end the program print the dot."
                  << "\n";

        for (;;)
        {
            std::cout << "Ennter the expression: ";
            std::cin.getline(expstr, 79);
            if (*expstr == '.')
                break;
            std::cout << "The answer: " << p.eval_exp(expstr) << "\n\n";
        }

        std::cout << "Selected file reading:"
                  << "\n";

        /*
            Example:
            var_f1 var_i1 var_f2 var_cat var_str
            12.2   20     51.1   1       aple
            15.5   30     10     2       plum
            21.0   45     562    3       aple
            30.5   50     452    3       plum
            40     61     231    4       tomato
            51.3   71     125    2       tomato
            60.6   80     121    1       plum
            70.001 91     121    1       aple
            82.012 10     110.0  4       tomato
       */
        evolm::IOInterface in;
        in.set_fname("tests/data/diverse_data.dat");

        std::string var_name("var_i1");

        evolm::Effects res;

        in.fgetvar(var_name, res);

        res.print(var_name);

        // -------------------------
        exit(0);

        // Testing model

        evolm::Pcg solver;
        evolm::Model model;

        std::vector<float> iR{0.041}; //, 0.041, 0.041, 0.041};

        std::vector<float> iG1{10.1}; //, 0.1, 0.1, 0.1};

        // size_t dim = 1;

        std::cout << "In main."
                  << "\n";

        model.append_residual(iR, 1);

        model.append_observation("tests/data/model_4/obs_1.dat"); // obs := 0

        size_t type1 = 0;
        float type2 = 0.0f;

        model.append_effect("tests/data/model_4/obs_489_snp_1000.txt", type1); // eff := 0
        model.append_effect("tests/data/model_4/fixed_1.dat", type2);          // eff := 1

        std::vector<int> corr_eff{0};

        std::vector<int> eff_trate{1, 0};
        int obs_trate = 0;

        std::string identity("I");

        model.append_corrstruct(iG1, 1, identity, 1000, corr_eff);

        model.append_traitstruct(obs_trate, eff_trate);
        // model.append_traitstruct(obs_trate, eff_trate);

        solver.append_model(model);

        solver.solve(2);

        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }
}