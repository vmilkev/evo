#include "model.hpp"
#include "sparse_pcg.hpp"
#include "iointerface.hpp"
#include "parser.hpp"
#include "modelparser.hpp"
#include "effects.hpp"

int main(void)
{
    try
    {
        /*
        std::cout << "Parser:"
                  << "\n";

        //evolm::Modelparser m;

        //m.eval_expr( "y1,y2 , another_var, ,= x1+ x2 - 1 +(1+f|x3*x4) - x4 +(-1|x4) + x1+x2:x3 + x4*x5*x6 - x4:x5:x6" );

        // --------------------------
        std::cout << "Selected file reading:"
                  << "\n";

            // Example:
            // id   obs_f var_i1 var_f2 var_cat var_str var_str2
            // 1002 12.2  20     51.1   1       aple    true
            // 1003 15.5  30     10     2       plum    true
            // 1004 21.0  45     562    3       aple    false
            // 1052 30.5  50     452    3       plum    true
            // 1062 40    61     231    4       tomato  false
            // 1072 51.3  71     125    2       tomato  true
            // 1082 -9.0  80     121    1       plum    true
            // 1092 70.01 91     121    1       aple    false
            // 1102 82.12 10     110.0  4       tomato  false

        evolm::IOInterface in;
        evolm::Modelparser m;

        in.set_fname("tests/data/diverse_data.dat");
        in.set_missing(-9.0);

        m.eval_expr("obs_f = -1 + var_cat + var_i1 + var_str*var_str2 + id");

        std::vector<std::string> obs;
        std::vector<std::string> fixed;
        std::vector<std::string> random;

        m.get_modelvars(obs,fixed,random);

         evolm::Effects obs_res;

        in.fgetvar(obs[0], obs_res);
        obs_res.print( obs[0] );

        // here we still need to parse "." (each col elementwise multiplication)
        // and concatenate effect matrices (+);
        // and vector of intercept (1).
        for (size_t i = 0; i < fixed.size(); i++)
        {
             evolm::Effects fix_res;
            in.fgetvar(fixed[i], fix_res);
            fix_res.print( fixed[i] );
        }

        // -------------------------
*/
/*
        // Testing model
        std::vector<float> R{40, 11,
                             11, 30}; // full matrix

        std::vector<float> y1{4.5, 2.9, 3.9, 3.5, 5.0, 4.0};

        std::vector<float> y2_miss{-999.0, 5.0, 6.8, 6.0, 7.5, -999.0};
        std::vector<float> y2{10.0, 5.0, 6.8, 6.0, 7.5, 10.0};

        std::vector<int> x1{1, 0,
                            0, 1,
                            0, 1,
                            1, 0,
                            1, 0,
                            0, 1};

        std::vector<int> x2_miss{
            0, 0,
            0, 1,
            0, 1,
            1, 0,
            1, 0,
            0, 0};

        std::vector<int> z1{0, 0, 0, 1, 0, 0, 0, 0, 0,
                            0, 0, 0, 0, 1, 0, 0, 0, 0,
                            0, 0, 0, 0, 0, 1, 0, 0, 0,
                            0, 0, 0, 0, 0, 0, 1, 0, 0,
                            0, 0, 0, 0, 0, 0, 0, 1, 0,
                            0, 0, 0, 0, 0, 0, 0, 0, 1};

        std::vector<int> z2_miss{0, 0, 0, 0, 0, 0, 0, 0, 0,
                                 0, 0, 0, 0, 1, 0, 0, 0, 0,
                                 0, 0, 0, 0, 0, 1, 0, 0, 0,
                                 0, 0, 0, 0, 0, 0, 1, 0, 0,
                                 0, 0, 0, 0, 0, 0, 0, 1, 0,
                                 0, 0, 0, 0, 0, 0, 0, 0, 0};

        std::vector<float> iA{1.8333, 0.5, 0, -0.66667, 0, -1, 0, 0, 0,
                              0.5, 2, 0.5, 0, -1, -1, 0, 0, 0,
                              0, 0.5, 2, 0, -1, 0.5, 0, -1, 0,
                              -0.66667, 0, 0, 1.8333, 0.5, 0, -1, 0, 0,
                              0, -1, -1, 0.5, 2.5, 0, -1, 0, 0,
                              -1, -1, 0.5, 0, 0, 2.5, 0, -1, 0,
                              0, 0, 0, -1, -1, 0, 2.3333, 0, -0.66667,
                              0, 0, -1, 0, 0, -1, 0, 2, 0,
                              0, 0, 0, 0, 0, 0, -0.66667, 0, 1.3333};

        std::vector<float> G{20, 18,
                             18, 40};

        evolm::model_sparse model;
        evolm::sparse_pcg solver;

        // define the model
        model.append_residual(R, 2);

        model.append_observation(y1, 6);      // obs := 0
        model.append_observation(y2_miss, 6); // obs := 1
std::cout<<"=> appending x1 ..."<<"\n";
        model.append_effect(x1, 6, 2); // eff := 0
std::cout<<"=> appending z1 ..."<<"\n";
        model.append_effect(z1, 6, 9); // eff := 1
std::cout<<"=> appending x2_miss ..."<<"\n";
        model.append_effect(x2_miss, 6, 2); // eff := 2
std::cout<<"=> appending z2_miss ..."<<"\n";
        model.append_effect(z2_miss, 6, 9); // eff := 3

        std::vector<int> eff_trate_1{0, 1};
        int obs_trate_1 = 0;

        std::vector<int> eff_trate_2{2, 3};
        int obs_trate_2 = 1;

        std::vector<int> corr_eff{1, 3};

        model.append_corrstruct(G, 2, iA, 9, corr_eff);
        model.append_traitstruct(obs_trate_1, eff_trate_1);
        model.append_traitstruct(obs_trate_2, eff_trate_2);

        solver.append_model(model);

        solver.solve();

        std::vector<float> sol = solver.get_solution();

        solver.get_solution("sparse_solution_model_mv_2.dat");

        solver.remove_model();

        model.clear();
*/
            // Prepare effects
            size_t n_snps = 141123;
            //---------------------------------
            /*std::vector<std::vector<float>> in;
            std::vector<std::vector<float>> in2;

            evolm::IOInterface datstream;
            datstream.set_fname("tests/data/model_4/obs_489_snp_141123.txt");
            datstream.fgetdata(in);
std::cout<<"        ==> completed reading SNPs..."<<"\n";
            datstream.set_fname("tests/data/model_4/fixed_1.dat");
            datstream.fgetdata(in2);
std::cout<<"        ==> completed reading fixed effects..."<<"\n";

            datstream.clear();

            std::vector<float> eff_snp;
            for (size_t i = 0; i < in.size(); i++)
            {
                for (size_t j = 0; j < in[0].size(); j++)
                    eff_snp.push_back(in[i][j]);
            }
            std::vector<float> eff_fixed;
            for (size_t i = 0; i < in2.size(); i++)
            {
                for (size_t j = 0; j < in2[0].size(); j++)
                    eff_fixed.push_back(in2[i][j]);
            }

            evolm::compact_storage<float> snp(eff_snp, in.size(), in[0].size());
            snp.fwrite("tests/data/model_4/snp_489_141123.bin");
            evolm::compact_storage<float> fixed(eff_fixed, in2.size(), in2[0].size());
            fixed.fwrite("tests/data/model_4/fixed_489.bin");*/
            
            //---------------------------------

            evolm::sparse_pcg solver;
            evolm::model_sparse model;

            model.set_sparsity_threshold(0.0);

std::cout<<"        ==> Appending model 4 on sparse solver..."<<"\n";
auto start = std::chrono::high_resolution_clock::now();

            std::vector<float> iR{0.00001};

            std::vector<float> iG1{0.00003};

            model.append_residual(iR, 1);

            model.append_observation("tests/data/model_4/obs_1.dat"); // obs := 0
            //bool t_var = true;
            //while (t_var)
            //    std::cout<<"waiting ..."<<"\r";

            //model.append_effect(eff_snp, in.size(), in[0].size()); // eff := 0
            //eff_snp.clear(); eff_snp.shrink_to_fit();
            model.append_effect("tests/data/model_4/snp_489_141123.bin"); // eff := 0

            //model.append_effect(eff_fixed, in2.size(), in2[0].size());    // eff := 1
            //eff_fixed.clear(); eff_fixed.shrink_to_fit();
            model.append_effect("tests/data/model_4/fixed_489.bin");    // eff := 1

            std::vector<int> corr_eff{0};

            std::vector<int> eff_trate{1, 0};
            //std::vector<int> eff_trate{0, 1};
            int obs_trate = 0;

            std::string identity("I");

            model.append_corrstruct(iG1, 1, identity, n_snps, corr_eff);

            model.append_traitstruct(obs_trate, eff_trate);
            
            solver.append_model(model);

            //clear initial data structure
            iR.clear();
            iG1.clear();
            /*for (auto &v: in)
                v.clear();
            for (auto &v: in2)
                v.clear();
            in.clear();
            in.shrink_to_fit();
            in2.clear();
            in2.shrink_to_fit();*/
            //----------------------------

            //bool t_var = true;
            //while (t_var)
            //    std::cout<<"waiting ..."<<"\r";

            solver.set_memory_limit(10);
            solver.set_cpu_limit(10);

std::cout<<"            ==> solving ..."<<"\n";
            
            solver.solve();

auto stop = std::chrono::high_resolution_clock::now();
auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
std::cout <<"model 4 on sparse solver (milliseconds): "<< duration.count() << std::endl;

            solver.get_solution("sparse_solution_model_4.dat");

            model.clear();


        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }
}