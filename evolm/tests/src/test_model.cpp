#include "catch_amalgamated.hpp"
#include "model_sparse.hpp"
#include "sparse_solver.hpp"
#include "sparse_pcg.hpp"
#include "iointerface.hpp"
#include <string>
#include <float.h>

TEST_CASE("Testing model set-up")
{
    bool is_ok = true;
    CHECK(is_ok == true);
}

TEST_CASE("Testing on model 1")
{
    // ---------------------------
    // model:
    // effect: 1       2
    // y1 = b1*X1 + a1*Z1 + e1;
    // effect: 3       4
    // y2 = b2*X2 + a2*Z2 + e2;
    // ---------------------------

    // DATA

    std::vector<float> iR{40, 11,
                          11, 30}; // full matrix, not inverse!

    std::vector<float> y1{4.5, 2.9, 3.9, 3.5, 5.0};

    std::vector<float> y2{6.8, 5.0, 6.8, 6.0, 7.5};

    std::vector<int> x1{1, 0,
                        0, 1,
                        0, 1,
                        1, 0,
                        1, 0};

    std::vector<int> x2(x1);

    std::vector<int> val_x1{1, 1, 1, 1, 1};
    std::vector<size_t> row_x1{0, 1, 2, 3, 4};
    std::vector<size_t> col_x1{0, 1, 1, 0, 0};

    std::vector<int> val_x2{1, 1, 1, 1, 1};
    std::vector<size_t> row_x2{0, 1, 2, 3, 4};
    std::vector<size_t> col_x2{0, 1, 1, 0, 0};

    std::vector<int> z1{   0, 0, 0, 1, 0, 0, 0, 0,
                           0, 0, 0, 0, 1, 0, 0, 0,
                           0, 0, 0, 0, 0, 1, 0, 0,
                           0, 0, 0, 0, 0, 0, 1, 0,
                           0, 0, 0, 0, 0, 0, 0, 1};

    std::vector<int> z2(z1);

    std::vector<int> val_z1{1, 1, 1, 1, 1};
    std::vector<size_t> row_z1{0, 1, 2, 3, 4};
    std::vector<size_t> col_z1{3, 4, 5, 6, 7};

    std::vector<int> val_z2{1, 1, 1, 1, 1};
    std::vector<size_t> row_z2{0, 1, 2, 3, 4};
    std::vector<size_t> col_z2{3, 4, 5, 6, 7};

    std::vector<float> iA{1.833, 0.5, 0, -0.667, 0, -1, 0, 0,
                          0.5, 2, 0.5, 0, -1, -1, 0, 0,
                          0, 0.5, 2, 0, -1, 0.5, 0, -1,
                          -0.667, 0, 0, 1.833, 0.5, 0, -1, 0,
                          0, -1, -1, 0.5, 2.5, 0, -1, 0,
                          -1, -1, 0.5, 0, 0, 2.5, 0, -1,
                          0, 0, 0, -1, -1, 0, 2, 0,
                          0, 0, -1, 0, 0, -1, 0, 2}; // full matrix

    std::vector<float> iG1{20, 18,
                           18, 40}; // full matrix, not inverse!

    std::vector<std::vector<float>> true_z{
        {1, 0, 0, 1, 1},
        {0, 1, 1, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0},
        {0, 1, 0, 0, 0},
        {0, 0, 1, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 0, 1}};

    std::vector<float> _rhs{
        0.15449,
        0.068767377,
        0,
        0,
        0,
        0.0557924,
        0.0296570,
        0.03911028,
        0.036144578,
        0.062557,
        0.620018,
        0.36811862,
        0,
        0,
        0,
        0.20620945,
        0.1557924,
        0.212326227,
        0.18674698,
        0.22706209};

    std::vector<float> _sol{
        4.360827,
        3.397239,
        0.150977,
        -0.015388,
        -0.078379,
        -0.010185,
        -0.270317,
        0.275839,
        -0.316081,
        0.243784,
        6.799822,
        5.880252,
        0.279713,
        -0.007601,
        -0.170316,
        -0.01257,
        -0.477801,
        0.517297,
        -0.478915,
        0.392017};

    std::vector<float> dval_true{
        0.08341,
        0.0556,
        0.154,
        0.168,
        0.168,
        0.1818,
        0.2378,
        0.2378,
        0.1958,
        0.1958,
        0.1112,
        0.0741,
        0.077,
        0.084,
        0.084,
        0.114,
        0.1421,
        0.1421,
        0.1211,
        0.1211};

    size_t n_all_levels = 4;

    std::vector<size_t> ordered_random_levels{2, 8, 2, 8};

    std::vector<std::vector<size_t>> rcov_offsets{{0, 0}, {0, 2}, {0, 10}, {0, 12}, {2, 0}, {2, 2}, {2, 10}, {2, 12}, {10, 0}, {10, 2}, {10, 10}, {10, 12}, {12, 0}, {12, 2}, {12, 10}, {12, 12}};

    std::vector<std::vector<float>> CoeffMatrix{{0.08341, 0, 0, 0, 0, 0.0278, 0, 0, 0.0278, 0.0278, -0.030583, 0, 0, 0, 0, -0.0101, 0, 0, -0.0101, -0.0101},
                                                {0, 0.0556, 0, 0, 0, 0, 0.0278, 0.0278, 0, 0, 0, -0.0203, 0, 0, 0, 0, -0.0101, -0.0101, 0, 0},
                                                {0, 0, 0.154, 0.042, 0, -0.056, 0, -0.084, 0, 0, 0, 0, -0.06931, -0.0189, 0, 0.02522, 0, 0.03781, 0, 0},
                                                {0, 0, 0.042, 0.168, 0.042, 0, -0.084, -0.084, 0, 0, 0, 0, -0.0189, -0.07563, -0.0189, 0, 0.03781, 0.03781, 0, 0},
                                                {0, 0, 0, 0.042, 0.168, 0, -0.084, 0.042, 0, -0.084, 0, 0, 0, -0.0189, -0.07563, 0, 0.03781, -0.0189, 0, 0.0378},
                                                {0.0278, 0, -0.056, 0, 0, 0.18183, 0.042, 0, -0.084, 0, -0.0101, 0, 0.02522, 0, 0, -0.0795, -0.0189, 0, 0.03781, 0},
                                                {0, 0.0278, 0, -0.084, -0.084, 0.042, 0.237887, 0, -0.084, 0, 0, -0.0101, 0, 0.03781, 0.03781, -0.0189, -0.10473, 0, 0.03781, 0},
                                                {0, 0.0278, -0.084, -0.084, 0.042, 0, 0, 0.2378, 0, -0.084, 0, -0.0101, 0.0378, 0.0378, -0.0189, 0, 0, -0.10473, 0, 0.0378},
                                                {0.0278, 0, 0, 0, 0, -0.084, -0.084, 0, 0.19587, 0, -0.0101, 0, 0, 0, 0, 0.0378, 0.0378, 0, -0.0858, 0},
                                                {0.0278, 0, 0, 0, -0.08403, 0, 0, -0.084, 0, 0.19587, -0.0101, 0, 0, 0, 0.0378, 0, 0, 0.0378, 0, -0.0858},
                                                {-0.030583, 0, 0, 0, 0, -0.0101, 0, 0, -0.0101, -0.0101, 0.1112, 0, 0, 0, 0, 0.03707, 0, 0, 0.03707, 0.03707},
                                                {0, -0.020389, 0, 0, 0, 0, -0.0101, -0.0101, 0, 0, 0, 0.07414, 0, 0, 0, 0, 0.037, 0.037, 0, 0},
                                                {0, 0, -0.0693, -0.0189, 0, 0.0252, 0, 0.0378, 0, 0, 0, 0, 0.077, 0.021, 0, -0.028, 0, -0.04201, 0, 0},
                                                {0, 0, -0.0189, -0.07563, -0.0189, 0, 0.0378, 0.0378, 0, 0, 0, 0, 0.021, 0.084, 0.021, 0, -0.04201, -0.042, 0, 0},
                                                {0, 0, 0, -0.0189, -0.07563, 0, 0.0378, -0.0189, 0, 0.0378, 0, 0, 0, 0.021, 0.084, 0, -0.042, 0.021, 0, -0.042},
                                                {-0.0101, 0, 0.0252, 0, 0, -0.0795, -0.0189, 0, 0.0378, 0, 0.037, 0, -0.028, 0, 0, 0.114, 0.021, 0, -0.042, 0},
                                                {0, -0.0101, 0, 0.0378, 0.0378, -0.0189, -0.10473, 0, 0.0378, 0, 0, 0.037, 0, -0.042, -0.042, 0.021, 0.1421, 0, -0.042, 0},
                                                {0, -0.0101, 0.0378, 0.0378, -0.0189, 0, 0, -0.10473, 0, 0.0378, 0, 0.037, -0.042, -0.042, 0.021, 0, 0, 0.1421, 0, -0.042},
                                                {-0.0101, 0, 0, 0, 0, 0.0378, 0.0378, 0, -0.0858, 0, 0.037, 0, 0, 0, 0, -0.042, -0.042, 0, 0.121, 0},
                                                {-0.0101, 0, 0, 0, 0.0378, 0, 0, 0.0378, 0, -0.0858, 0.037, 0, 0, 0, -0.042, 0, 0, -0.042, 0, 0.121}};
    // DATA for testing IO interface:

    std::vector<std::vector<int>> allele{
        {2, 0, 1, 1, 0, 0, 0, 2, 1, 2},
        {5, 0, 0, 0, 0, 2, 0, 2, 1, 0},
        {1, 5, 2, 1, 1, 0, 0, 2, 1, 2},
        {0, 0, 2, 1, 0, 1, 0, 2, 2, 1},
        {0, 1, 1, 2, 0, 0, 0, 2, 1, 2},
        {1, 1, 5, 1, 0, 2, 0, 2, 2, 1}};

    std::vector<std::vector<float>> matr{
        {12.2, 20, 51.1},
        {15.5, 30, 10},
        {21.0, 45, 562},
        {30.5, 50, 452},
        {40, 61, 231},
        {51.3, 71, 125},
        {60.6, 80, 121},
        {70.001, 91, 121},
        {82.012, 10, 110.0}};

    std::vector<std::vector<int>> allelebin{
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0},
        {0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};

    std::vector<std::vector<int>> alleledat{
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0},
        {0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0},
        {0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0},
        {0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0},
        {0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0}};

        evolm::compact_storage<int> z(z1,5,8);
        evolm::compact_storage<int> x(x1,5,2);
        z.fwrite("tests/data/z.bin");
        x.fwrite("tests/data/x.bin");
        
        evolm::compact_storage<float> ainv(iA,8,8);
        ainv.fwrite("tests/data/ainv.bin");

    // ========================================================================================

    SECTION("1. Testing sizes of residual and observation")
    {
        try
        {
            evolm::model_sparse model;

            CHECK(model.size_of("res") == 0);
            CHECK(model.size_of("obs") == 0);

            model.append_residual(iR, 2);

            model.append_observation(y1, 5); // obs := 0

            CHECK(model.size_of("res") == 4);

            CHECK(model.size_of("obs") == 5);

            model.append_residual(iR, 2);

            model.append_observation(y2, 5); // obs := 1

            CHECK(model.size_of("res") == 8);
            CHECK(model.size_of("obs") == 10);

            model.clear_residuals();
            model.clear_observations();

            CHECK(model.size_of("res") == 0);
            CHECK(model.size_of("obs") == 0);
        }
        catch (const std::exception &e)
        {
            std::cerr << "1. Testing sizes of residual and observation:" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "1. Testing sizes of residual and observation:" << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "1. Testing sizes of residual and observation." << '\n';
        }
    }
    // ========================================================================================
    SECTION("2. Testing sizes of effects: from full vectors on memory")
    {
        try
        {
            evolm::model_sparse model;

            model.set_sparsity_threshold(0.4);

            CHECK(model.size_of("eff") == 0);

            model.append_effect(x1, 5, 2); // eff:= 0, size_of == 5 because of sparse

            CHECK(model.size_of("eff") == 5);

            model.append_effect(x2, 5, 2); // eff:= 1, size_of == 5 because of sparse

            CHECK(model.size_of("eff") == 10);

            model.append_effect(z1, 5, 8); // eff:= 2, size_of == 5 because of sparse

            CHECK(model.size_of("eff") == 15);

            model.append_effect(z2, 5, 8); // eff:= 3, size_of == 5 because of sparse

            CHECK(model.size_of("eff") == 20);

            model.clear_effects();

            CHECK(model.size_of("eff") == 0);
        }
        catch (const std::exception &e)
        {
            std::cerr << "2. Testing sizes of effects: from full vectors on memory." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "2. Testing sizes of effects: from full vectors on memory." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "2. Testing sizes of effects: from full vectors on memory." << '\n';
        }
    }
    // ========================================================================================
    SECTION("3. Testing sizes of effects: from sparse vectors on memory")
    {
        try
        {
            evolm::model_sparse model;

            model.set_sparsity_threshold(0.4);

            CHECK(model.size_of("eff") == 0);

            model.append_effect(val_x1, row_x1, col_x1, 5, 2); // eff:= 0, size_of == 5 because of sparse

            CHECK(model.size_of("eff") == 5);

            model.append_effect(val_x2, row_x2, col_x2, 5, 2); // eff:= 1, size_of == 5 because of sparse

            CHECK(model.size_of("eff") == 10);

            model.append_effect(val_z1, row_z1, col_z1, 5, 8); // eff:= 2, size_of == 5 because of sparse

            CHECK(model.size_of("eff") == 15);

            model.append_effect(val_z2, row_z2, col_z2, 5, 8); // eff:= 3, size_of == 5 because of sparse

            CHECK(model.size_of("eff") == 20);

            model.clear_effects();

            CHECK(model.size_of("eff") == 0);
        }
        catch (const std::exception &e)
        {
            std::cerr << "3. Testing sizes of effects: from sparse vectors on memory." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "3. Testing sizes of effects: from sparse vectors on memory." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "3. Testing sizes of effects: from sparse vectors on memory." << '\n';
        }
    }
    // ========================================================================================
    SECTION("4. Testing effects: data from files")
    {
        try
        {
            evolm::model_sparse model;

            model.set_sparsity_threshold(0.4);

            evolm::compact_storage<int> eff_x1(x1, 5, 2);
            evolm::compact_storage<int> eff_x2(x2, 5, 2);
            evolm::compact_storage<int> eff_z1(z1, 5, 8);
            evolm::compact_storage<int> eff_z2(z2, 5, 8);

            eff_x1.set_sparsity_threshold(0.4); eff_x1.optimize();
            eff_x2.set_sparsity_threshold(0.4); eff_x2.optimize();
            eff_z1.set_sparsity_threshold(0.4); eff_z1.optimize();
            eff_z2.set_sparsity_threshold(0.4); eff_z2.optimize();

            eff_x1.fwrite("x1.bin");
            eff_x2.fwrite("x2.bin");
            eff_z1.fwrite("z1.bin");
            eff_z2.fwrite("z2.bin");

            model.append_effect("z1.bin"); // eff_0
            model.append_effect("z2.bin"); // eff_1
            model.append_effect("x1.bin"); // eff_2
            model.append_effect("x2.bin"); // eff_3

            eff_x1.fclear("x1.bin");
            eff_x2.fclear("x2.bin");
            eff_z1.fclear("z1.bin");
            eff_z2.fclear("z2.bin");

            evolm::compact_storage<float> eff_0_st = model.test_effects(0);
            evolm::compact_storage<float> eff_1_st = model.test_effects(1);
            evolm::compact_storage<float> eff_2_st = model.test_effects(2);
            evolm::compact_storage<float> eff_3_st = model.test_effects(3);

            std::vector<size_t> eff_0_key;
            std::vector<float> eff_0_val;
            std::vector<size_t> eff_1_key;
            std::vector<float> eff_1_val;
            std::vector<size_t> eff_2_key;
            std::vector<float> eff_2_val;
            std::vector<size_t> eff_3_key;
            std::vector<float> eff_3_val;

            eff_0_st.to_sparse(eff_0_val, eff_0_key);
            eff_1_st.to_sparse(eff_1_val, eff_1_key);
            eff_2_st.to_sparse(eff_2_val, eff_2_key);
            eff_3_st.to_sparse(eff_3_val, eff_3_key);

            size_t counter = 0;
            for (size_t i = 0; i < z1.size(); i++)
            {
                if (z1[i] != 0)
                {
                    CHECK( static_cast<float>(z1[i]) == eff_0_val[counter] );
                    counter++;
                }
            }
            counter = 0;
            for (size_t i = 0; i < z2.size(); i++)
            {
                if (z2[i] != 0)
                {
                    CHECK( (float)(z2[i]) == eff_1_val[counter] );
                    counter++;
                }
            }
            counter = 0;
            for (size_t i = 0; i < x1.size(); i++)
            {
                if (x1[i] != 0)
                {
                    CHECK( (float)x1[i] == eff_2_val[counter] );
                    counter++;
                }
            }
            counter = 0;
            for (size_t i = 0; i < x2.size(); i++)
            {
                if (x2[i] != 0)
                {
                    CHECK( (float)x2[i] == eff_3_val[counter] );
                    counter++;
                }
            }
            counter = 0;

            model.append_observation("tests/data/y1.dat"); // obs_0
            model.append_observation("tests/data/y2.dat"); // obs_1

            std::vector<float> obs_0 = model.test_observations(0);
            std::vector<float> obs_1 = model.test_observations(1);

            for (size_t i = 0; i < obs_0.size(); i++)
                CHECK(y1[i] == obs_0[i]);

            for (size_t i = 0; i < obs_1.size(); i++)
                CHECK(y2[i] == obs_1[i]);

            CHECK(model.size_of("eff") == 20);
            CHECK(model.size_of("obs") == 10);

            model.append_residual("tests/data/iR.dat");

            std::vector<int> corr_eff{0, 1};

            model.append_corrstruct("tests/data/iG.dat", "tests/data/ainv.bin", corr_eff);

            std::vector<float> res_test = model.test_residual(0);
            std::vector<float> var_test = model.test_variance(0);
            std::vector<float> cor_test = model.test_correlation(0);

            CHECK(res_test.size() == iR.size());
            CHECK(var_test.size() == iG1.size());
            CHECK(cor_test.size() == iA.size());

            for (size_t i = 0; i < res_test.size(); i++)
                CHECK(iR[i] == res_test[i]);

            for (size_t i = 0; i < var_test.size(); i++)
                CHECK(iG1[i] == var_test[i]);

            for (size_t i = 0; i < cor_test.size(); i++)
                CHECK(iA[i] == cor_test[i]);

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "4. Testing effects: data from files." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "4. Testing effects: data from files." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "4. Testing effects: data from files." << '\n';
        }
    }
    // ========================================================================================
    SECTION("5. Testing effects: data from memory")
    {
        try
        {
            evolm::model_sparse model;

            model.set_sparsity_threshold(0.4);

            evolm::compact_storage<int> eff_x1(val_x1, row_x1, col_x1, 5, 2);
            evolm::compact_storage<int> eff_x2(val_x2, row_x2, col_x2, 5, 2);
            evolm::compact_storage<int> eff_z1(val_z1, row_z1, col_z1, 5, 8);
            evolm::compact_storage<int> eff_z2(val_z2, row_z2, col_z2, 5, 8);

            eff_x1.set_sparsity_threshold(0.4); eff_x1.optimize();
            eff_x2.set_sparsity_threshold(0.4); eff_x2.optimize();
            eff_z1.set_sparsity_threshold(0.4); eff_z1.optimize();
            eff_z2.set_sparsity_threshold(0.4); eff_z2.optimize();

            model.append_effect(eff_z1); // eff_0
            model.append_effect(eff_z2); // eff_1
            model.append_effect(eff_x1); // eff_2
            model.append_effect(eff_x2); // eff_3

            evolm::compact_storage<float> eff_0_st = model.test_effects(0);
            evolm::compact_storage<float> eff_1_st = model.test_effects(1);
            evolm::compact_storage<float> eff_2_st = model.test_effects(2);
            evolm::compact_storage<float> eff_3_st = model.test_effects(3);

            std::vector<size_t> eff_0_key;
            std::vector<float> eff_0_val;
            std::vector<size_t> eff_1_key;
            std::vector<float> eff_1_val;
            std::vector<size_t> eff_2_key;
            std::vector<float> eff_2_val;
            std::vector<size_t> eff_3_key;
            std::vector<float> eff_3_val;

            eff_0_st.to_sparse(eff_0_val, eff_0_key);
            eff_1_st.to_sparse(eff_1_val, eff_1_key);
            eff_2_st.to_sparse(eff_2_val, eff_2_key);
            eff_3_st.to_sparse(eff_3_val, eff_3_key);

            size_t counter = 0;

            for (size_t i = 0; i < z1.size(); i++)
            {
                if (z1[i] != 0)
                {
                    CHECK( static_cast<float>(z1[i]) == eff_0_val[counter] );
                    counter++;
                }
            }
            counter = 0;
            for (size_t i = 0; i < z2.size(); i++)
            {
                if (z2[i] != 0)
                {
                    CHECK( (float)(z2[i]) == eff_1_val[counter] );
                    counter++;
                }
            }
            counter = 0;
            for (size_t i = 0; i < x1.size(); i++)
            {
                if (x1[i] != 0)
                {
                    CHECK( (float)x1[i] == eff_2_val[counter] );
                    counter++;
                }
            }
            counter = 0;
            for (size_t i = 0; i < x2.size(); i++)
            {
                if (x2[i] != 0)
                {
                    CHECK( (float)x2[i] == eff_3_val[counter] );
                    counter++;
                }
            }
            counter = 0;

            CHECK(model.size_of("eff") == 20);

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "5. Testing effects: data from memory." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "5. Testing effects: data from memory." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "5. Testing effects: data from memory." << '\n';
        }
    }
    // ========================================================================================
    SECTION("6. Testing sizes of structures")
    {
        try
        {
            evolm::model_sparse model;

            CHECK(model.size_of("var") == 0);
            CHECK(model.size_of("cor") == 0);
            CHECK(model.size_of("cor_eff") == 0);

            std::vector<int> corr_eff{2, 3};

            std::vector<int> eff_trate_1{0, 2};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{1, 3};
            int obs_trate_2 = 1;

            model.append_corrstruct(iG1, 2, iA, 8, corr_eff);
            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            CHECK(model.size_of("var") == 4);
            CHECK(model.size_of("cor") == 64);
            CHECK(model.size_of("cor_eff") == 2);
            CHECK(model.size_of("obs_trt") == 2);
            CHECK(model.size_of("eff_trt") == 4);

            model.clear_corrstruct();
            model.clear_traitstruct();

            CHECK(model.size_of("var") == 0);
            CHECK(model.size_of("cor") == 0);
            CHECK(model.size_of("cor_eff") == 0);
            CHECK(model.size_of("obs_trt") == 0);
            CHECK(model.size_of("eff_trt") == 0);
        }
        catch (const std::exception &e)
        {
            std::cerr << "6. Testing sizes of structures." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "6. Testing sizes of structures." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "6. Testing sizes of structures." << '\n';
        }
    }
    // ========================================================================================
    SECTION("7. Testing shapes")
    {
        try
        {
            evolm::model_sparse model;

            std::vector<std::vector<size_t>> shapes;

            std::vector<int> corr_eff{2, 3};

            std::vector<int> eff_trate_1{0, 2};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{1, 3};
            int obs_trate_2 = 1;

            model.append_residual(iR, 2);
            model.append_observation(y1, 5); // obs := 0

            model.append_residual(iR, 2);
            model.append_observation(y2, 5); // obs := 1

            model.append_effect(x1, 5, 2); // eff:= 0

            model.append_effect(x2, 5, 2); // eff:= 1

            model.append_effect(z1, 5, 8); // eff:= 2

            model.append_effect(z2, 5, 8); // eff:= 3

            model.append_corrstruct(iG1, 2, iA, 8, corr_eff);
            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            shapes = model.shape_of("eff");

            CHECK(shapes[0].at(0) == 5);
            CHECK(shapes[0].at(1) == 2);
            CHECK(shapes[1].at(0) == 5);
            CHECK(shapes[1].at(1) == 2);
            CHECK(shapes[2].at(0) == 5);
            CHECK(shapes[2].at(1) == 8);
            CHECK(shapes[3].at(0) == 5);
            CHECK(shapes[3].at(1) == 8);

            shapes.clear();

            shapes = model.shape_of("obs");
            for (auto &e : shapes)
            {
                CHECK(e[0] == 5);
                CHECK(e[1] == 1);
            }
            shapes.clear();

            shapes = model.shape_of("res");
            for (auto &e : shapes)
            {
                CHECK(e[0] == 2);
                CHECK(e[1] == 2);
            }
            shapes.clear();

            shapes = model.shape_of("var");
            for (auto &e : shapes)
            {
                CHECK(e[0] == 2);
                CHECK(e[1] == 2);
            }
            shapes.clear();

            shapes = model.shape_of("cor");
            for (auto &e : shapes)
            {
                CHECK(e[0] == 8);
                CHECK(e[1] == 8);
            }
            shapes.clear();

            shapes = model.shape_of("cor_eff");
            for (auto &e : shapes)
            {
                CHECK(e[0] == 2);
                CHECK(e[1] == 1);
            }
            shapes.clear();

            shapes = model.shape_of("eff_trt");
            for (auto &e : shapes)
            {
                CHECK(e[0] == 2);
                CHECK(e[1] == 1);
            }
            shapes.clear();

            shapes = model.shape_of("obs_trt");
            for (auto &e : shapes)
            {
                CHECK(e[0] == 1);
                CHECK(e[1] == 1);
            }
            shapes.clear();

            model.clear_residuals();
            model.clear_observations();
            model.clear_effects();
            model.clear_corrstruct();
            model.clear_traitstruct();
        }
        catch (const std::exception &e)
        {
            std::cerr << "7. Testing shapes." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "7. Testing shapes." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "7. Testing shapes." << '\n';
        }
    }
    // ========================================================================================
    SECTION("8. Testing vect_z_uni")
    {
        try
        {
            evolm::sparse_solver solver;
            evolm::model_sparse model;

            std::vector<int> eff_trate_1{0, 1};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{2, 3};
            int obs_trate_2 = 1;

            // define the model
            model.append_residual(iR, 2);

            model.append_observation(y1, 5); // obs := 0
            model.append_observation(y2, 5); // obs := 1

            model.append_effect(x1, 5, 2); // eff := 0
            model.append_effect(z1, 5, 8); // eff := 1
            model.append_effect(x2, 5, 2); // eff := 2
            model.append_effect(z2, 5, 8); // eff := 3

            std::vector<int> corr_eff{1, 3};

            model.append_corrstruct(iG1, 2, iA, 8, corr_eff);
            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            evolm::matrix<float> z;

            solver.test_vect_z_uni(0, z);

            //z.print("testing Z1");

            for (size_t i = 0; i < true_z.size(); i++)
            {
                for (size_t j = 0; j < true_z[i].size(); j++)
                {
                    CHECK(z(i, j) == true_z[i][j]);
                }
            }

            z.clear();

            solver.test_vect_z_uni(1, z);

            //z.print("testing Z2");

            for (size_t i = 0; i < true_z.size(); i++)
            {
                for (size_t j = 0; j < true_z[i].size(); j++)
                {
                    CHECK(z(i, j) == true_z[i][j]);
                }
            }

            z.clear();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "8. Testing vect_z_uni." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "8. Testing vect_z_uni." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "8. Testing vect_z_uni." << '\n';
        }
    }
    // ========================================================================================
    SECTION("9. Testing RHS")
    {
        try
        {
            evolm::sparse_pcg solver;
            evolm::model_sparse model;

            std::vector<int> eff_trate_1{0, 1};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{2, 3};
            int obs_trate_2 = 1;

            // define the model
            model.append_residual(iR, 2);

            model.append_observation(y1, 5); // obs := 0
            model.append_observation(y2, 5); // obs := 1

            model.append_effect(x1, 5, 2); // eff := 0
            model.append_effect(z1, 5, 8); // eff := 1
            model.append_effect(x2, 5, 2); // eff := 2
            model.append_effect(z2, 5, 8); // eff := 3

            std::vector<int> corr_eff{1, 3};

            model.append_corrstruct(iG1, 2, iA, 8, corr_eff);
            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            evolm::matrix<float> model_rhs = solver.test_rhs();

            for (size_t i = 0; i < _rhs.size(); i++)
                CHECK((_rhs[i]) == Catch::Approx(model_rhs[i]).margin(0.0001).epsilon(1e-3));

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "9. Testing RHS." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "9. Testing RHS." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "9. Testing RHS." << '\n';
        }
    }
    // ========================================================================================
    SECTION("10. Testing intermediate data structures")
    {
        try
        {
            evolm::sparse_pcg solver;
            evolm::model_sparse model;

            std::vector<int> eff_trate_1{0, 1};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{2, 3};
            int obs_trate_2 = 1;

            // define the model
            model.append_residual(iR, 2);

            model.append_observation(y1, 5); // obs := 0
            model.append_observation(y2, 5); // obs := 1

            model.append_effect(x1, 5, 2); // eff := 0
            model.append_effect(z1, 5, 8); // eff := 1
            model.append_effect(x2, 5, 2); // eff := 2
            model.append_effect(z2, 5, 8); // eff := 3

            std::vector<int> corr_eff{1, 3};

            model.append_corrstruct(iG1, 2, iA, 8, corr_eff);
            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            size_t n_all_levels_checking = solver.test_num_all_levels();

            CHECK(n_all_levels_checking == n_all_levels);

            std::vector<size_t> ordered_random_levels_checking = solver.test_ordered_levels();

            CHECK(ordered_random_levels_checking.size() == ordered_random_levels.size());

            for (size_t i = 0; i < ordered_random_levels_checking.size(); i++)
                CHECK(ordered_random_levels_checking[i] == ordered_random_levels[i]);

            std::vector<std::vector<size_t>> rcov_offsets_checking = solver.test_cov_offsets();

            CHECK(rcov_offsets_checking.size() == rcov_offsets.size());
            CHECK(rcov_offsets_checking[0].size() == rcov_offsets[0].size());

            for (size_t i = 0; i < rcov_offsets.size(); i++)
            {
                for (size_t j = 0; j < rcov_offsets[0].size(); j++)
                {
                    CHECK(rcov_offsets_checking[i][j] == rcov_offsets[i][j]);
                }
            }

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "10. Testing intermediate data structures." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "10. Testing intermediate data structures." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "10. Testing intermediate data structures." << '\n';
        }
    }
    // ========================================================================================
    SECTION("11. Testing dval")
    {
        try
        {
            evolm::sparse_pcg solver;
            evolm::model_sparse model;

            std::vector<int> eff_trate_1{0, 1};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{2, 3};
            int obs_trate_2 = 1;

            // define the model
            model.append_residual(iR, 2);

            model.append_observation(y1, 5); // obs := 0
            model.append_observation(y2, 5); // obs := 1

            model.append_effect(x1, 5, 2); // eff := 0
            model.append_effect(z1, 5, 8); // eff := 1
            model.append_effect(x2, 5, 2); // eff := 2
            model.append_effect(z2, 5, 8); // eff := 3

            std::vector<int> corr_eff{1, 3};

            model.append_corrstruct(iG1, 2, iA, 8, corr_eff);
            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            std::vector<double> dval = solver.test_dval();

            CHECK(dval.size() == dval_true.size());

            for (size_t i = 0; i < dval_true.size(); i++)
                CHECK(dval[i] == Catch::Approx(dval_true[i]).margin(0.0001).epsilon(1e-4));

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "11. Testing dval." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "11. Testing dval." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "11. Testing dval." << '\n';
        }
    }
    // ========================================================================================
    SECTION("12. Testing CoeffMatrix")
    {
        try
        {
            evolm::sparse_pcg solver;
            evolm::model_sparse model;

            std::vector<int> eff_trate_1{0, 1};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{2, 3};
            int obs_trate_2 = 1;

            // define the model
            model.append_residual(iR, 2);

            model.append_observation(y1, 5); // obs := 0
            model.append_observation(y2, 5); // obs := 1

            model.append_effect(x1, 5, 2); // eff := 0
            model.append_effect(z1, 5, 8); // eff := 1
            model.append_effect(x2, 5, 2); // eff := 2
            model.append_effect(z2, 5, 8); // eff := 3

            std::vector<int> corr_eff{1, 3};

            model.append_corrstruct(iG1, 2, iA, 8, corr_eff);
            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            std::vector<std::vector<float>> A = solver.test_A();

            CHECK(A.size() == CoeffMatrix.size());
            CHECK(A[0].size() == CoeffMatrix[0].size());

            for (size_t i = 0; i < A.size(); i++)
            {
                for (size_t j = 0; j < A[0].size(); j++)
                {
                    CHECK(A[i][j] == Catch::Approx(CoeffMatrix[i][j]).margin(0.001).epsilon(1e-3));
                }
            }

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "12. Testing CoeffMatrix." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "12. Testing CoeffMatrix." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "12. Testing CoeffMatrix." << '\n';
        }
    }
    // ========================================================================================
    SECTION("13. Testing solution")
    {
        try
        {
            evolm::sparse_pcg solver;
            evolm::model_sparse model;

            std::vector<int> eff_trate_1{0, 1};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{2, 3};
            int obs_trate_2 = 1;

            // define the model
            model.append_residual(iR, 2);

            model.append_observation(y1, 5); // obs := 0
            model.append_observation(y2, 5); // obs := 1

            model.append_effect(x1, 5, 2); // eff := 0
            model.append_effect(z1, 5, 8); // eff := 1
            model.append_effect(x2, 5, 2); // eff := 2
            model.append_effect(z2, 5, 8); // eff := 3

            std::vector<int> corr_eff{1, 3};

            model.append_corrstruct(iG1, 2, iA, 8, corr_eff);
            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            solver.solve();

            std::vector<float> sol = solver.get_solution();

            solver.get_solution("sparse_solution_model_1.dat");

            for (size_t i = 0; i < sol.size(); i++)
            {
                CHECK((_sol[i]) == Catch::Approx(sol[i]).margin(0.0001).epsilon(1e-3));
            }

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "13. Testing solution." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "13. Testing solution." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "13. Testing solution." << '\n';
        }
    }
    // ========================================================================================
    SECTION("14. The overloaded methods (IO interface): testing solution")
    {
        try
        {
            std::cout<<"   ===> Testing 14 ..."<<"\n";
            evolm::model_sparse model;
            evolm::sparse_pcg solver;

            model.append_effect("tests/data/z.bin"); // eff_0
            model.append_effect("tests/data/x.bin"); // eff_1
            model.append_effect("tests/data/z.bin"); // eff_2
            model.append_effect("tests/data/x.bin"); // eff_3

            model.append_observation("tests/data/y1.dat"); // obs_0
            model.append_observation("tests/data/y2.dat"); // obs_1

            model.append_residual("tests/data/iR.dat");

            std::vector<int> eff_trate_1{1, 2};
            int obs_trate_1 = 0;

            std::vector<int> eff_trate_2{3, 0};
            int obs_trate_2 = 1;

            std::vector<int> corr_eff{1, 3};
            std::vector<int> corr_eff2{0, 2}; // same matrix but diff order of effects

            model.append_corrstruct("tests/data/iG_diff_order.dat", "tests/data/ainv.bin", corr_eff2);

            model.append_traitstruct(obs_trate_1, eff_trate_1);
            model.append_traitstruct(obs_trate_2, eff_trate_2);

            solver.append_model(model);

            solver.set_memory_limit(0.0000005);
            solver.set_cpu_limit(2);

            size_t n_all_levels_checking = solver.test_num_all_levels();

            CHECK(n_all_levels_checking == n_all_levels);

            std::vector<size_t> ordered_random_levels_checking = solver.test_ordered_levels();

            CHECK(ordered_random_levels_checking.size() == ordered_random_levels.size());

            for (size_t i = 0; i < ordered_random_levels_checking.size(); i++)
                CHECK(ordered_random_levels_checking[i] == ordered_random_levels[i]);

            std::vector<std::vector<size_t>> rcov_offsets_checking = solver.test_cov_offsets();

            CHECK(rcov_offsets_checking.size() == rcov_offsets.size());
            CHECK(rcov_offsets_checking[0].size() == rcov_offsets[0].size());

            for (size_t i = 0; i < rcov_offsets.size(); i++)
            {
                for (size_t j = 0; j < rcov_offsets[0].size(); j++)
                {
                    CHECK(rcov_offsets_checking[i][j] == rcov_offsets[i][j]);
                }
            }

            std::vector<std::vector<float>> A = solver.test_A();

            CHECK(A.size() == CoeffMatrix.size());
            CHECK(A[0].size() == CoeffMatrix[0].size());

            for (size_t i = 0; i < A.size(); i++)
            {
                for (size_t j = 0; j < A[0].size(); j++)
                {
                    CHECK(A[i][j] == Catch::Approx(CoeffMatrix[i][j]).margin(0.001).epsilon(1e-3));
                }
            }

            solver.solve();

            std::vector<float> sol = solver.get_solution();

            solver.get_solution("sparse_solution_model_1_v2.dat");

            for (size_t i = 0; i < sol.size(); i++)
            {
                CHECK((_sol[i]) == Catch::Approx(sol[i]).margin(0.0001).epsilon(1e-3));
            }

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "14. The overloaded methods (IO interface): testing solution." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "14. The overloaded methods (IO interface): testing solution." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "14. The overloaded methods (IO interface): testing solution." << '\n';
        }
    }
     // ========================================================================================
}

TEST_CASE("Testing on model 2")
{
    // one trait setup
    // EXAMPLE: 7.1. p.110.
    // y = Xb + Ua + Wm + Sp + e
    // corr1: [a m]
    // corr2: [p]

    // changes to model interpretation:
    // X moves to random effects but without covariance structure
    // corr1: [a m]
    // corr2: [p]
    // corr0 == 0: [b], because it is a fixed effect

    // ----------- DATA ----------------------------
    std::vector<int> x{
        1, 0, 0, 1, 0,
        1, 0, 0, 0, 1,
        1, 0, 0, 0, 1,
        1, 0, 0, 1, 0,
        0, 1, 0, 1, 0,
        0, 1, 0, 0, 1,
        0, 1, 0, 0, 1,
        0, 0, 1, 0, 1,
        0, 0, 1, 1, 0,
        0, 0, 1, 0, 1};

    std::vector<int> s{
        1, 0, 0, 0,
        1, 0, 0, 0,
        0, 0, 1, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        1, 0, 0, 0,
        0, 0, 0, 1,
        0, 0, 0, 1,
        1, 0, 0, 0,
        0, 0, 1, 0};

    std::vector<int> w{
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0};

    std::vector<int> u{
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

    std::vector<float> iA{
        2, 0.5, 0, 0, -1, 0.5, 0, 0, -1, 0, 0, 0, 0, 0,
        0.5, 3, 1, 0, -1, -1, 0, 0, 0.5, -1, 0, 0, -1, 0,
        0, 1, 3.5, 0, 0.5, -0.5, 0.5, -1, 0, -1, -1, 0, 0, -1,
        0, 0, 0, 1.5, 0, 0.5, -1, 0, 0, 0, 0, 0, 0, 0,
        -1, -1, 0.5, 0, 2.5, 0, 0, -1, 0, 0, 0, 0, 0, 0,
        0.5, -1, -0.5, 0.5, 0, 3.5, -1, 0, -1, 0, 0, 0, 0, -1,
        0, 0, 0.5, -1, 0, -1, 3, 0.5, 0, 0, -1, -1, 0, 0,
        0, 0, -1, 0, -1, 0, 0.5, 2.5, 0, 0, 0, -1, 0, 0,
        -1, 0.5, 0, 0, 0, -1, 0, 0, 2.5, 0, 0, 0, -1, 0,
        0, -1, -1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0,
        0, 0, -1, 0, 0, 0, -1, 0, 0, 0, 2, 0, 0, 0,
        0, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 2, 0, 0,
        0, -1, 0, 0, 0, 0, 0, 0, -1, 0, 0, 0, 2, 0,
        0, 0, -1, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 2};

    std::vector<float> corS{
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1};

    std::vector<float> y{35, 20, 25, 40, 42, 22, 35, 34, 20, 40};

    std::vector<float> iR{350};

    std::vector<float> iG1{1};

    std::vector<float> iG2{150, -40, -40, 90};

    std::vector<float> iG3{40};

    std::vector<float> _sol{
        14.492235,
        17.878691,
        15.926014,
        20.047707,
        13.198704,
        0.56389,
        -1.244389,
        1.164936,
        -0.484436,
        0.629533,
        -0.858629,
        -1.155968,
        1.917404,
        -0.553263,
        -1.055081,
        0.385358,
        0.863317,
        -2.979573,
        1.751301,
        0.261565,
        -1.583161,
        0.735732,
        0.585864,
        -0.507471,
        0.841041,
        1.299317,
        -0.157915,
        0.659541,
        -0.152954,
        0.915958,
        0.442008,
        0.093056,
        0.362213,
        -1.70083,
        0.415397,
        0.824915,
        0.460519};
    // ---------------------------------------------
    evolm::sparse_pcg solver;
    evolm::model_sparse model;

    // ----- define the model -------

    model.append_residual(iR, 1);

    model.append_observation(y, 10); // obs := 0

    model.append_effect(x, 10, 5);  // eff := 0
    model.append_effect(s, 10, 4);  // eff := 1
    model.append_effect(w, 10, 14); // eff := 2
    model.append_effect(u, 10, 14); // eff := 3

    std::vector<int> corr_eff_1{3, 2};
    std::vector<int> corr_eff_2{1};

    std::vector<int> eff_trate{0, 3, 2, 1};
    int obs_trate = 0;

    model.append_corrstruct(iG2, 2, iA, 14, corr_eff_1);

    std::string identity("I");
    //model.append_corrstruct(iG3, 1, corS, 4, corr_eff_2);
    model.append_corrstruct(iG3, 1, identity, 4, corr_eff_2);

    model.append_traitstruct(obs_trate, eff_trate);

    solver.append_model(model);
    // ========================================================================================

    SECTION("15. Testing solution")
    {
        try
        {
            solver.solve();

            std::vector<float> sol = solver.get_solution();

            solver.get_solution("sparse_solution_model_2.dat");

            for (size_t i = 0; i < sol.size(); i++)
                CHECK((_sol[i]) == Catch::Approx(sol[i]).margin(0.03).epsilon(1e-3));

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "15. Testing solution." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "15. Testing solution." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "15. Testing solution." << '\n';
        }
    }
    // ========================================================================================
}

TEST_CASE("Testing on model 3")
{
    // one trait setup
    // EXAMPLE: 11.2. p.183.
    // y = Xb + Za + e

    // corr: [a], identity matrix

    // ----------- DATA ----------------------------
    std::vector<int> x{
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1};

    std::vector<float> z{
        1.357, -0.357, 0.286, 0.286, -0.286, -1.214, -0.143, 0.071, -0.143, 1.214,
        0.357, -0.357, -0.714, -0.714, -0.286, 0.786, -0.143, 0.071, -0.143, -0.786,
        0.357, 0.643, 1.286, 0.286, 0.714, -1.214, -0.143, 0.071, -0.143, 1.214,
        -0.643, -0.357, 1.286, 0.286, -0.286, -0.214, -0.143, 0.071, 0.857, 0.214,
        -0.643, 0.643, 0.286, 1.286, -0.286, -1.214, -0.143, 0.071, -0.143, 1.214,
        0.357, 0.643, -0.714, 0.286, -0.286, 0.786, -0.143, 0.071, 0.857, 0.214,
        -0.643, -0.357, 0.286, 0.286, -0.286, 0.786, -0.143, 0.071, 0.857, -0.786,
        -0.643, 0.643, 0.286, -0.714, -0.286, -0.214, -0.143, 0.071, 0.857, -0.786};

    std::vector<float> corZ{
        1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 1};

    std::vector<float> y{9.0, 13.4, 12.7, 15.4, 5.9, 7.7, 10.2, 4.8};

    std::vector<float> iR{245};

    std::vector<float> iG1{10};

    std::vector<float> _sol{
        9.94415,
        0.0873718,
        -0.312064,
        0.263549,
        -0.0805778,
        0.110685,
        0.139673,
        -2.3004e-07,
        -1.62874e-07,
        -0.0609667,
        -0.0158181};

    // ---------------------------------------------

    evolm::sparse_pcg solver;
    evolm::model_sparse model;

    // ----- define the model -------

    model.append_residual(iR, 1);

    model.append_observation(y, 8); // obs := 0

    model.append_effect(x, 8, 1);  // eff := 0
    model.append_effect(z, 8, 10); // eff := 1

    std::vector<int> corr_eff{1};

    std::vector<int> eff_trate{0, 1};
    int obs_trate = 0;

    model.append_corrstruct(iG1, 1, corZ, 10, corr_eff);

    model.append_traitstruct(obs_trate, eff_trate);

    solver.append_model(model);
    // ========================================================================================
    
    SECTION("16. Testing solution")
    {
        try
        {
            solver.solve();

            std::vector<float> sol = solver.get_solution();

            solver.get_solution("sparse_solution_model_3.dat");

            for (size_t i = 0; i < sol.size(); i++)
                CHECK((_sol[i]) == Catch::Approx(sol[i]).margin(0.003).epsilon(1e-3));

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "16. Testing solution." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "16. Testing solution." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "16. Testing solution." << '\n';
        }
    }
    // ========================================================================================
    SECTION("17. Test solution with identity corr structure")
    {
        try
        {
            model.clear();
            solver.remove_model();

            model.append_residual(iR, 1);

            model.append_observation(y, 8); // obs := 0

            model.append_effect(x, 8, 1);  // eff := 0
            model.append_effect(z, 8, 10); // eff := 1

            std::vector<int> _corr_eff{1};

            std::vector<int> _eff_trate{0, 1};
            int _obs_trate = 0;

            std::string identity("I");

            model.append_corrstruct(iG1, 1, identity, 10, _corr_eff);
            //model.append_corrstruct(iG1, 1, corZ, 10, _corr_eff);

            model.append_traitstruct(_obs_trate, _eff_trate);

            solver.append_model(model);

            solver.solve();

            std::vector<float> sol = solver.get_solution();

            solver.get_solution("sparse_solution_model_3(I).dat");

            for (size_t i = 0; i < sol.size(); i++)
                CHECK((_sol[i]) == Catch::Approx(sol[i]).margin(0.003).epsilon(1e-3));

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "17. Test solution with identity corr structure." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "17. Test solution with identity corr structure." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "17. Test solution with identity corr structure." << '\n';
        }
    }
    // ========================================================================================
}

TEST_CASE("Testing on model 4")
{
    SECTION("18. Testing full SNP blup")
    {
        try
        {
            // Prepare effects
            //---------------------------------
            std::vector<std::vector<int>> in;
            std::vector<std::vector<float>> in2;

            evolm::IOInterface datstream;
            datstream.set_fname("tests/data/model_4/obs_489_snp_100.txt");
            datstream.fgetdata(in);
            datstream.set_fname("tests/data/model_4/fixed_1.dat");
            datstream.fgetdata(in2);

            //std::cout<<"dim of in: "<<in.size()<<" "<<in[0].size()<<"\n";
            //std::cout<<"dim of in2: "<<in2.size()<<" "<<in2[0].size()<<"\n";

            std::vector<int> eff_snp;
            for (size_t i = 0; i < in.size(); i++)
            {
                for (size_t j = 0; j < in[0].size(); j++)
                    eff_snp.push_back(in[i][j]);
            }
            std::vector<int> eff_fixed;
            for (size_t i = 0; i < in2.size(); i++)
            {
                for (size_t j = 0; j < in2[0].size(); j++)
                    eff_fixed.push_back(in2[i][j]);
            }
            //---------------------------------

            evolm::sparse_pcg solver;
            evolm::model_sparse model;

            std::vector<float> iR{0.041};

            std::vector<float> iG1{0.1};

            model.append_residual(iR, 1);

            model.append_observation("tests/data/model_4/obs_1.dat"); // obs := 0

            //model.append_effect("tests/data/model_4/obs_489_snp_100.txt"); // eff := 0
            //model.append_effect("tests/data/model_4/fixed_1.dat");         // eff := 1
            model.append_effect(eff_snp, in.size(), in[0].size()); // eff := 0
            model.append_effect(eff_fixed, in2.size(), in2[0].size());         // eff := 1

            std::vector<int> corr_eff{0};

            std::vector<int> eff_trate{1, 0};
            int obs_trate = 0;

            std::string identity("I");

            model.append_corrstruct(iG1, 1, identity, 100, corr_eff);

            model.append_traitstruct(obs_trate, eff_trate);

            solver.append_model(model);

            solver.set_memory_limit(0.005);
            solver.set_cpu_limit(4);

            solver.solve();

            // std::vector<float> sol = solver.get_solution();

            solver.get_solution("sparse_solution_model_4.dat");

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "18. Testing full SNP blup." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "18. Testing full SNP blup." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "18. Testing full SNP blup." << '\n';
        }
    }
}

TEST_CASE("Testing multivariate model with missing values")
{
    // ---------------------------
    // model 5.2 pp. 78:
    // effect: 1       2
    // y1 = b1*X1 + a1*Z1 + e1;
    // effect: 3       4
    // y2 = b2*X2 + a2*Z2 + e2;
    // ---------------------------

    // DATA

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

    std::vector<std::vector<float>> true_z{
        {1, 0, 0, 1, 1},
        {0, 1, 1, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {1, 0, 0, 0, 0},
        {0, 1, 0, 0, 0},
        {0, 0, 1, 0, 0},
        {0, 0, 0, 1, 0},
        {0, 0, 0, 0, 1}};

    std::vector<float> _rhs{
        0.1543399,
        0.068767377,
        0,
        0,
        0,
        0.0557924,
        0.0296570,
        0.03911028,
        0.036144578,
        0.062557,
        0.620018,
        0.36811862,
        0,
        0,
        0,
        0.20620945,
        0.1557924,
        0.212326227,
        0.18674698,
        0.22706209};

    std::vector<float> _sol{
        4.367,
        3.657,
        0.1298,
        -0.0836,
        -0.098,
        0.007,
        -0.343,
        0.1915,
        -0.308,
        0.2006,
        -0.0184,
        6.834,
        6.007,
        0.266,
        -0.0752,
        -0.194,
        0.01557,
        -0.555,
        0.440,
        -0.483,
        0.349,
        -0.11937};

    // ========================================================================================
    SECTION("19. Testing model using different fixed and raandom matrices (some with a missing records) for all traits (observations)")
    {
        try
        {
            evolm::sparse_pcg solver;
            evolm::model_sparse model;

            // define the model
            model.append_residual(R, 2);

            model.append_observation(y1, 6);      // obs := 0
            model.append_observation(y2_miss, 6); // obs := 1

            model.append_effect(x1, 6, 2);      // eff := 0

            model.append_effect(z1, 6, 9);      // eff := 1

            model.append_effect(x2_miss, 6, 2); // eff := 2
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

            // std::vector<std::vector<float>> A = solver.test_A();
            // std::ofstream out_a("sparse_a.txt");
            // for (size_t i = 0; i < A.size(); i++)
            // {
            //     for (size_t j = 0; j < A[0].size(); j++)
            //         out_a << std::setprecision(8) << A[i][j]<<" ";
            //     out_a <<"\n";
            // }
            // out_a.close();

            solver.solve();

            std::vector<float> sol = solver.get_solution();

            solver.get_solution("sparse_solution_model_mv_2.dat");

            for (size_t i = 0; i < sol.size(); i++)
                CHECK((_sol[i]) == Catch::Approx(sol[i]).margin(0.0001).epsilon(1e-3));

            solver.remove_model();

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "19. Testing model using different fixed and raandom matrices (some with a missing records)." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "19. Testing model using different fixed and raandom matrices (some with a missing records)." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "19. Testing model using different fixed and raandom matrices (some with a missing records)." << '\n';
        }
    }
    // ========================================================================================
}

TEST_CASE("Testing on big model 4")
{
    SECTION("20. Testing big full SNP blup")
    {
        try
        {
            // Prepare effects
            //---------------------------------
            std::vector<std::vector<int>> in;
            std::vector<std::vector<float>> in2;

            evolm::IOInterface datstream;
            datstream.set_fname("tests/data/model_4/obs_1000_snp_1000.txt");
            datstream.fgetdata(in);
            datstream.set_fname("tests/data/model_4/fixed_1000.dat");
            datstream.fgetdata(in2);

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
            //---------------------------------

            evolm::sparse_pcg solver;
            evolm::model_sparse model;

            model.set_sparsity_threshold(0.2);

std::cout<<"        ==> Testing model 4 on sparse solver..."<<"\n";
auto start = std::chrono::high_resolution_clock::now();

            std::vector<float> iR{0.01};

            std::vector<float> iG1{0.03};

            model.append_residual(iR, 1);

            model.append_observation("tests/data/model_4/obs_1000.dat"); // obs := 0
std::cout<<"            ==> appending eff_snp ..."<<"\n";
            model.append_effect(eff_snp, in.size(), in[0].size()); // eff := 0
std::cout<<"            ==> appending eff_fixed ..."<<"\n";
            model.append_effect(eff_fixed, in2.size(), in2[0].size());         // eff := 1

            std::vector<int> corr_eff{0};

            std::vector<int> eff_trate{1, 0};
            //std::vector<int> eff_trate{0, 1};
            int obs_trate = 0;

            std::string identity("I");

            model.append_corrstruct(iG1, 1, identity, 1000, corr_eff);

            model.append_traitstruct(obs_trate, eff_trate);

            solver.append_model(model);

            solver.set_memory_limit(0.05);
            solver.set_cpu_limit(4);

            //solver.set_memory_limit(15);

            // std::vector<std::vector<float>> A = solver.test_A();
            // std::ofstream out_a("sparse_a.txt");
            // for (size_t i = 0; i < A.size(); i++)
            // {
            //     for (size_t j = 0; j < A[0].size(); j++)
            //         out_a << std::setprecision(16) << A[i][j]<<" ";
            //     out_a <<"\n";
            // }
            // out_a.close();

std::cout<<"            ==> solving ..."<<"\n";
            solver.solve();

auto stop = std::chrono::high_resolution_clock::now();
auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
std::cout <<"model 4 on sparse solver (milliseconds): "<< duration.count() << std::endl;

            // std::vector<float> sol = solver.get_solution();

            solver.get_solution("sparse_solution_model_4.dat");

            model.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "20. Testing full SNP blup." << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (const std::string &err)
        {
            std::cerr << "20. Testing full SNP blup." << '\n';
            std::cerr << err << '\n';
        }
        catch (...)
        {
            std::cerr << "20. Testing full SNP blup." << '\n';
        }
    }
}
