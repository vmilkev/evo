#include "model_sparse.hpp"
#include "iointerface.hpp"
#include "model_parser.hpp"

void get_effect_from_data(evolm::IOInterface &in_data, std::vector<std::string> &in_var_vect, evolm::effects_storage &out_s, std::string &effect_name);
bool consist_snp_substring(std::string &var_name);

int main(void)
{
    try
    {
        std::cout << "Parser:"<< "\n\n";
        
        evolm::model_sparse model; // model constructor, all results should go here
        model.set_missing(-9.0);

        evolm::model_parser m_parser; // parser constructor

        //-------------------------------------------------
        // Example of data file:
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
        //--------------------------------------------------

        std::string model_expression1("obs_f, var_f2 ~ x1 + (1+var_f2|id)a +(1|var_cat)random2 +(1|x1)random_3; data = tests/data/diverse_data.dat; x1 = tests/data/data_matr.dat");
        std::string model_expression("weight ~ 1 + herd + pen + (1|id) + (1|id)b + (1|dam); data = tests/data/model_2/model_2_data.dat; obs_missing_value = [-9.0]");
        //std::string model_expression("obs_f,var_f2 , , ,~ 1 + (1+var_f2|id) + x1+ x2 +(1+f|x3*x4) - x4 +(1 + x1:x2|x4) + x1+x2:x3 + (1|var_cat) + x4*x5*x6 - x4:x5:x6");

        std::string extra_vars_expression(
                                "x2 = tests/data/x3_9_obs.dat;"
                                "obs_missing_value = [-9.0];"
                                "Z1_snp = tests/data/snp_with_id_9_20.dat;"
                                "Z2_snp = tests/data/snp_no_id_9_20.dat");

        std::string extra_vars_expression2("Z3_snp = tests/data/large_data/SNPS/YY_plink_");

        std::string extra_vars_expression3("G = [1 2 3, 4 5.0 6, 7 8 9, 10 11 12.0 ]");

        std::string variance_expression("var = (id, b)*A1*G1 + (dam)*I*G2 + R; A1 = tests/data/model_2/A1.corbin; G1 = [150 -40, -40 90]; G2 = [40]; R = [350]");

        m_parser.process_expression(extra_vars_expression);
        m_parser.process_expression(model_expression);
        m_parser.process_expression(model_expression);
        m_parser.process_expression(extra_vars_expression3);
        m_parser.process_expression(variance_expression);

        //m_parser.print();

        //m_parser.report();

        std::cout << "Completed model parsing ..."<< "\n";

        //------------------------------------------------

/*
            // Prepare effects
            size_t n_snps = 141123;
            //---------------------------------
            // std::vector<std::vector<float>> in;
            // std::vector<std::vector<float>> in2;

            // evolm::IOInterface datstream;
            // datstream.set_fname("tests/data/model_4/obs_489_snp_141123.txt");
            // datstream.fgetdata(in);

            // datstream.set_fname("tests/data/model_4/fixed_1.dat");
            // datstream.fgetdata(in2);

            // datstream.clear();

            // std::vector<float> eff_snp;
            // for (size_t i = 0; i < in.size(); i++)
            // {
            //     for (size_t j = 0; j < in[0].size(); j++)
            //         eff_snp.push_back(in[i][j]);
            // }
            // std::vector<float> eff_fixed;
            // for (size_t i = 0; i < in2.size(); i++)
            // {
            //     for (size_t j = 0; j < in2[0].size(); j++)
            //         eff_fixed.push_back(in2[i][j]);
            // }

            // evolm::compact_storage<float> snp(eff_snp, in.size(), in[0].size());
            // snp.fwrite("tests/data/model_4/snp_489_141123.bin");
            // evolm::compact_storage<float> fixed(eff_fixed, in2.size(), in2[0].size());
            // fixed.fwrite("tests/data/model_4/fixed_489.bin");
            
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
            // for (auto &v: in)
            //     v.clear();
            // for (auto &v: in2)
            //     v.clear();
            // in.clear();
            // in.shrink_to_fit();
            // in2.clear();
            // in2.shrink_to_fit();
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
*/
//std::cout<<"call return"<<"\n";
        return 0;
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
    }
    catch (const std::string &e)
    {
        std::cerr << e << '\n';
    }
    catch (...)
    {
        std::cerr << "Exception -> exit." << '\n';
    }
}