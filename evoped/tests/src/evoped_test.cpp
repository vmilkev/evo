#include "catch_amalgamated.hpp"
#include "cs_matrix.hpp"
#include "Amat.hpp"
#include "Gmat.hpp"
#include "Hmat.hpp"

TEST_CASE("Testing Amat class")
{
    // ------------------- Data: -----------------------
    std::vector<double> iAcorr{ // full A(-1)
                                1, 1, 2,
                                2, 1, 1,
                                2, 2, 2,
                                13, 1, -1,
                                13, 2, -1,
                                13, 13, 4.64286,
                                14, 1, -1,
                                14, 2, -1,
                                14, 13, 1.5,
                                14, 14, 4.78205,
                                15, 13, 0.142857,
                                15, 14, -1,
                                15, 15, 4.91462,
                                16, 13, -1.14286,
                                16, 14, 0.615385,
                                16, 15, -0.415584,
                                16, 16, 3.62837,
                                17, 15, -1.45455,
                                17, 16, -1.45455,
                                17, 17, 3.57576,
                                18, 13, -1.14286,
                                18, 15, -1.14286,
                                18, 18, 2.28571,
                                19, 13, -1,
                                19, 14, -1,
                                19, 19, 2,
                                20, 14, -0.666667,
                                20, 20, 1.33333,
                                21, 15, -0.727273,
                                21, 17, 0.666667,
                                21, 21, 2.12121,
                                22, 14, -0.666667,
                                22, 22, 1.33333,
                                23, 14, -1.23077,
                                23, 16, -1.23077,
                                23, 23, 2.46154,
                                24, 13, -1,
                                24, 14, -1,
                                24, 24, 2,
                                25, 15, 0.680851,
                                25, 17, -1.33333,
                                25, 21, -1.33333,
                                25, 25, 3.34752,
                                26, 15, -1.3617,
                                26, 25, -1.3617,
                                26, 26, 2.7234
                                };
    
    std::vector<double> irAcorr(iAcorr); // reduced A(-1)
    
    std::vector<double> Acorr{ // full A
                                1, 1, 1,
                                2, 2, 1,
                                13, 1, 0.5,
                                13, 2, 0.5,
                                13, 13, 1,
                                14, 1, 0.5,
                                14, 2, 0.5,
                                14, 13, 0.5,
                                14, 14, 1,
                                15, 1, 0.5,
                                15, 2, 0.5,
                                15, 13, 0.75,
                                15, 14, 0.75,
                                15, 15, 1.25,
                                16, 1, 0.5,
                                16, 2, 0.5,
                                16, 13, 0.875,
                                16, 14, 0.625,
                                16, 15, 1,
                                16, 16, 1.375,
                                17, 1, 0.5,
                                17, 2, 0.5,
                                17, 13, 0.8125,
                                17, 14, 0.6875,
                                17, 15, 1.125,
                                17, 16, 1.1875,
                                17, 17, 1.5,
                                18, 1, 0.5,
                                18, 2, 0.5,
                                18, 13, 0.875,
                                18, 14, 0.625,
                                18, 15, 1,
                                18, 16, 0.9375,
                                18, 17, 0.96875,
                                18, 18, 1.375,
                                19, 1, 0.5,
                                19, 2, 0.5,
                                19, 13, 0.75,
                                19, 14, 0.75,
                                19, 15, 0.75,
                                19, 16, 0.75,
                                19, 17, 0.75,
                                19, 18, 0.75,
                                19, 19, 1.25,
                                20, 1, 0.25,
                                20, 2, 0.25,
                                20, 13, 0.25,
                                20, 14, 0.5,
                                20, 15, 0.375,
                                20, 16, 0.3125,
                                20, 17, 0.34375,
                                20, 18, 0.3125,
                                20, 19, 0.375,
                                20, 20, 1,
                                21, 1, 0.25,
                                21, 2, 0.25,
                                21, 13, 0.375,
                                21, 14, 0.375,
                                21, 15, 0.625,
                                21, 16, 0.5,
                                21, 17, 0.5625,
                                21, 18, 0.5,
                                21, 19, 0.375,
                                21, 20, 0.1875,
                                21, 21, 1,
                                22, 1, 0.25,
                                22, 2, 0.25,
                                22, 13, 0.25,
                                22, 14, 0.5,
                                22, 15, 0.375,
                                22, 16, 0.3125,
                                22, 17, 0.34375,
                                22, 18, 0.3125,
                                22, 19, 0.375,
                                22, 20, 0.25,
                                22, 21, 0.1875,
                                22, 22, 1,
                                23, 1, 0.5,
                                23, 2, 0.5,
                                23, 13, 0.6875,
                                23, 14, 0.8125,
                                23, 15, 0.875,
                                23, 16, 1,
                                23, 17, 0.9375,
                                23, 18, 0.78125,
                                23, 19, 0.75,
                                23, 20, 0.40625,
                                23, 21, 0.4375,
                                23, 22, 0.40625,
                                23, 23, 1.3125,
                                24, 1, 0.5,
                                24, 2, 0.5,
                                24, 13, 0.75,
                                24, 14, 0.75,
                                24, 15, 0.75,
                                24, 16, 0.75,
                                24, 17, 0.75,
                                24, 18, 0.75,
                                24, 19, 0.75,
                                24, 20, 0.375,
                                24, 21, 0.375,
                                24, 22, 0.375,
                                24, 23, 0.75,
                                24, 24, 1.25,
                                25, 1, 0.375,
                                25, 2, 0.375,
                                25, 13, 0.59375,
                                25, 14, 0.53125,
                                25, 15, 0.875,
                                25, 16, 0.84375,
                                25, 17, 1.03125,
                                25, 18, 0.734375,
                                25, 19, 0.5625,
                                25, 20, 0.265625,
                                25, 21, 0.78125,
                                25, 22, 0.265625,
                                25, 23, 0.6875,
                                25, 24, 0.5625,
                                25, 25, 1.28125,
                                26, 1, 0.4375,
                                26, 2, 0.4375,
                                26, 13, 0.671875,
                                26, 14, 0.640625,
                                26, 15, 1.0625,
                                26, 16, 0.921875,
                                26, 17, 1.07812,
                                26, 18, 0.867188,
                                26, 19, 0.65625,
                                26, 20, 0.320312,
                                26, 21, 0.703125,
                                26, 22, 0.320312,
                                26, 23, 0.78125,
                                26, 24, 0.65625,
                                26, 25, 1.07812,
                                26, 26, 1.4375
    };
    
    std::vector<double> rAcorr(Acorr); // reduced A

    std::vector<double> A22corr{
            1.375, 
            0.75, 1.25, 
            0.3125, 0.375, 1, 
            0.5, 0.375, 0.1875, 1, 
            0.3125, 0.375, 0.25, 0.1875, 1, 
            0.78125, 0.75, 0.40625, 0.4375, 0.40625, 1.3125, 
            0.75, 0.75, 0.375, 0.375, 0.375, 0.75, 1.25, 
            0.734375, 0.5625, 0.265625, 0.78125, 0.265625, 0.6875, 0.5625, 1.28125, 
            0.867188, 0.65625, 0.320312, 0.703125, 0.320312, 0.78125, 0.65625, 1.07812, 1.4375
    };

    std::vector<double> iA22corr{
            1.53926, 
            -0.318426, 1.57368, 
            0.0227579, -0.13411, 1.20773, 
            -0.0721947, -0.00765392, -0.0284804, 1.93081, 
            0.0227579, -0.13411, -0.125599, -0.0284804, 1.20773, 
            -0.252518, -0.323534, -0.173773, 0.0441258, -0.173773, 1.61325, 
            -0.318426, -0.426317, -0.13411, -0.00765392, -0.13411, -0.323534, 1.57368, 
            -0.0881696, -0.0260629, 0.0303647, -1.02853, 0.0303647, -0.210556, -0.0260629, 2.73362, 
            -0.409304, -0.0728126, -0.0468092, -0.13377, -0.0468092, -0.215256, -0.0728126, -1.36924, 2.23926
    };
    // -------------------------------------------------

    SECTION("Making full A")
    {
        try
        {
            evoped::Amat ap;
            evolm::matrix<double> A;
            std::vector<std::int64_t> a_id;

            ap.make_matrix("tests/data/ped_bkg.dat", false);

            ap.get_matrix("A", A, a_id, false);

            std::vector<double> Avect;
            for (size_t i = 0; i < a_id.size(); i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    if ( A(i,j) != 0.0 )
                    {
                        Avect.push_back( a_id[i] );
                        Avect.push_back( a_id[j] );
                        Avect.push_back( A(i,j) );
                    }
                }
            }

            CHECK( Acorr.size() == Avect.size() );

            for (size_t i = 0; i < Avect.size();)
            {
                CHECK( Acorr[i] == Avect[i] );
                CHECK( Acorr[i+1] == Avect[i+1] );
                CHECK( Acorr[i+2] == Catch::Approx(Avect[i+2]) );
                i = i + 3;
            }

        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        catch(const std::string& e)
        {
            std::cerr << e << '\n';
        }        
    }

    SECTION("Making redused A")
    {
        try
        {
            evoped::Amat ap;
            evolm::matrix<double> rA;
            std::vector<std::int64_t> a_id;

            ap.make_matrix("tests/data/ped_bkg.dat", "tests/data/typed2", false);

            ap.get_matrix("rA", rA, a_id, false);

            std::vector<double> rAvect;
            for (size_t i = 0; i < a_id.size(); i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    if ( rA(i,j) != 0.0 )
                    {
                        rAvect.push_back( a_id[i] );
                        rAvect.push_back( a_id[j] );
                        rAvect.push_back( rA(i,j) );
                    }
                }
            }

            CHECK( rAcorr.size() == rAvect.size() );

            for (size_t i = 0; i < rAvect.size();)
            {
                CHECK( rAcorr[i] == rAvect[i] );
                CHECK( rAcorr[i+1] == rAvect[i+1] );
                CHECK( rAcorr[i+2] == Catch::Approx(rAvect[i+2]) );
                i = i + 3;
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        catch(const std::string& e)
        {
            std::cerr << e << '\n';
        }        
    }

    SECTION("Making full A(-1)")
    {
        try
        {
            evoped::Amat ap;
            evolm::matrix<double> iA;
            std::vector<std::int64_t> a_id;

            ap.make_matrix("tests/data/ped_bkg.dat", true);

            ap.get_matrix("iA", iA, a_id, false);

            std::vector<double> iAvect;
            for (size_t i = 0; i < a_id.size(); i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    if ( iA(i,j) != 0.0 )
                    {
                        iAvect.push_back( a_id[i] );
                        iAvect.push_back( a_id[j] );
                        iAvect.push_back( iA(i,j) );
                    }
                }
            }

            CHECK( iAcorr.size() == iAvect.size() );

            for (size_t i = 0; i < iAvect.size();)
            {
                CHECK( iAcorr[i] == iAvect[i] );
                CHECK( iAcorr[i+1] == iAvect[i+1] );
                CHECK( iAcorr[i+2] == Catch::Approx(iAvect[i+2]) );
                i = i + 3;
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        catch(const std::string& e)
        {
            std::cerr << e << '\n';
        }        
    }

    SECTION("Making reduced A(-1)")
    {
        try
        {
            evoped::Amat ap;
            evolm::matrix<double> iA;
            std::vector<std::int64_t> a_id;

            ap.make_matrix("tests/data/ped_bkg.dat", "tests/data/typed2", true);

            ap.get_matrix("irA", iA, a_id, false);

            std::vector<double> iAvect;
            for (size_t i = 0; i < a_id.size(); i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    if ( iA(i,j) != 0.0 )
                    {
                        iAvect.push_back( a_id[i] );
                        iAvect.push_back( a_id[j] );
                        iAvect.push_back( iA(i,j) );
                    }
                }
            }

            CHECK( irAcorr.size() == iAvect.size() );

            for (size_t i = 0; i < iAvect.size();)
            {
                CHECK( irAcorr[i] == iAvect[i] );
                CHECK( irAcorr[i+1] == iAvect[i+1] );
                CHECK( irAcorr[i+2] == Catch::Approx(iAvect[i+2]) );
                i = i + 3;
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        catch(const std::string& e)
        {
            std::cerr << e << '\n';
        }        
    }

    SECTION("Making All")
    {
        try
        {
            evoped::Amat ap;
            evolm::matrix<double> iA;
            std::vector<std::int64_t> a_id;
            evolm::matrix<double> irA;
            std::vector<std::int64_t> ra_id;
            evolm::matrix<double> iA22;
            std::vector<std::int64_t> ai22_id;
            evolm::matrix<double> A22;
            std::vector<std::int64_t> a22_id;

            ap.make_all("tests/data/ped_bkg.dat", "tests/data/typed2");

            ap.get_matrix("iA", iA, a_id, false);
            ap.get_matrix("irA", irA, ra_id, false);
            ap.get_matrix("iA22", iA22, ai22_id, false);
            ap.get_matrix("A22", A22, a22_id, false);

            std::vector<double> iAvect;
            for (size_t i = 0; i < a_id.size(); i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    if ( iA(i,j) != 0.0 )
                    {
                        iAvect.push_back( a_id[i] );
                        iAvect.push_back( a_id[j] );
                        iAvect.push_back( iA(i,j) );
                    }
                }
            }

            CHECK( iAcorr.size() == iAvect.size() );

            for (size_t i = 0; i < iAvect.size();)
            {
                CHECK( iAcorr[i] == iAvect[i] );
                CHECK( iAcorr[i+1] == iAvect[i+1] );
                CHECK( iAcorr[i+2] == Catch::Approx(iAvect[i+2]) );
                i = i + 3;
            }

            std::vector<double> irAvect;
            for (size_t i = 0; i < ra_id.size(); i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    if ( irA(i,j) != 0.0 )
                    {
                        irAvect.push_back( ra_id[i] );
                        irAvect.push_back( ra_id[j] );
                        irAvect.push_back( irA(i,j) );
                    }
                }
            }

            CHECK( irAcorr.size() == irAvect.size() );

            for (size_t i = 0; i < irAvect.size();)
            {
                CHECK( irAcorr[i] == irAvect[i] );
                CHECK( irAcorr[i+1] == irAvect[i+1] );
                CHECK( irAcorr[i+2] == Catch::Approx(irAvect[i+2]) );
                i = i +3;
            }

            std::vector<double> A22vect;
            for (size_t i = 0; i < a22_id.size(); i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    A22vect.push_back( A22(i,j) );
                }
            }

            CHECK( A22corr.size() == A22vect.size() );

            for (size_t i = 0; i < A22vect.size(); i++)
                CHECK( A22corr[i] == Catch::Approx(A22vect[i]) );

            std::vector<double> iA22vect;
            for (size_t i = 0; i < ai22_id.size(); i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    iA22vect.push_back( iA22(i,j) );
                }
            }

            CHECK( iA22corr.size() == iA22vect.size() );

            for (size_t i = 0; i < iA22vect.size(); i++)
                CHECK( iA22corr[i] == Catch::Approx(iA22vect[i]) );

        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        catch(const std::string& e)
        {
            std::cerr << e << '\n';
        }        
    }
}

TEST_CASE( "Testing Gmat class" )
{
    // --------------- Data -----------------------
    std::vector<double> G_true{
        18, 18, 1.04474,
        19, 18, -0.538185,
        19, 19, 1.54696,
        20, 18, 0.634713,
        20, 19, -0.976534,
        20, 20, 1.19089,
        21, 18, -0.327422,
        21, 19, -0.284721,
        21, 20, 0.115789,
        21, 21, 0.638775,
        22, 18, 0.269347,
        22, 19, -0.900899,
        22, 20, 0.549643,
        22, 21, 0.0663641,
        22, 22, 1.24772,
        23, 18, -0.393878,
        23, 19, 0.2452,
        23, 20, -0.454686,
        23, 21, 0.019099,
        23, 22, -0.284043,
        23, 23, 0.662528,
        24, 18, -0.838939,
        24, 19, 0.737872,
        24, 20, -0.726399,
        24, 21, 0.419554,
        24, 22, -0.445152,
        24, 23, 0.414502,
        24, 24, 1.07722,
        25, 18, -0.481374,
        25, 19, 0.0627879,
        25, 20, -0.216471,
        25, 21, 0.165394,
        25, 22, -0.0574784,
        25, 23, 0.136427,
        25, 24, 0.165394,
        25, 25, 0.653164,
        26, 18, 0.557582,
        26, 19, 0.195985,
        26, 20, -0.181954,
        26, 21, -0.814581,
        26, 22, -0.510107,
        26, 23, -0.288437,
        26, 24, -0.741507,
        26, 25, -0.376492,
        26, 26, 2.11649
    };

    double beta_true = 0.537181;
    double alpha_true = 0.604285;
    double a_diag_true = 1.21181;
    double a_ofd_true = 0.528429;
    double g_diag_true = 1.13094;
    double g_ofd_true = -0.141211;

    std::vector<double> G_scaled_true{
        1.21787, 
        0.423887, 1.38896, 
        0.787055, 0.153532, 1.18301, 
        0.4463, 0.432254, 0.546738, 0.960567, 
        0.639854, 0.184004, 0.737157, 0.526826, 1.2059, 
        0.489838, 0.739501, 0.37159, 0.570283, 0.440339, 1.04826, 
        0.302717, 0.937991, 0.254308, 0.715996, 0.367618, 0.80771, 1.19971, 
        0.442869, 0.619135, 0.432407, 0.715161, 0.496463, 0.680053, 0.660474, 1.03668, 
        0.894652, 0.696236, 0.459985, 0.300812, 0.327777, 0.532319, 0.318534, 0.571061, 1.66529
    };

    std::vector<double> iG_true{
        2.42388, 
        -0.314841, 2.2437, 
        -0.909926, 0.12813, 2.03803, 
        -0.389615, 0.831076, -0.727092, 3.7043, 
        -0.463894, 0.30563, -0.575403, 0.0156385, 1.60403, 
        -0.368775, -0.0915967, -0.029816, 0.306612, -0.128099, 2.65289, 
        0.417486, -1.83483, 0.28858, -1.96199, -0.2555, -1.30352, 3.88307, 
        0.487369, -0.493326, 0.0440369, -1.67387, -0.375967, -0.768891, 0.438543, 2.87667, 
        -0.886641, -0.465173, 0.0561584, 0.24176, 0.180579, -0.120528, 0.391362, -0.51595, 1.31718
    };
    
    std::vector<double> iG_sparse_true{
        2.51922, 
        -1.07622, 2.04631, 
        -0.481927, -0.54046, 1.60344, 
        0.580117, -0.167238, -0.28004, 3.70845, 
        -0.386788, 0.265574, 0.151647, -0.666562, 1.08568, 
        0.0648192, -0.549682, -0.202345, -1.42132, 0, 2.49298, 
        -0.373436, 0.12984, -0.11783, -1.01098, 0, 0, 1.7878, 
        -0.0465776, 0.113674, -0.110325, -0.800005, 0, 0, 0, 1.29542, 
        -0.936183, 0.172983, 0.254283, -0.428974, 0, 0, 0, 0, 1.15272
    };

    std::vector<double> iH_sparse_true{
            1, 1, 2.000000,
            2, 1, 1.000000,
            2, 2, 2.000000,
            13, 1, -1.000000,
            13, 2, -1.000000,
            13, 13, 4.642857,
            14, 1, -1.000000,
            14, 2, -1.000000,
            14, 13, 1.500000,
            14, 14, 4.782051,
            15, 13, 0.142857,
            15, 14, -1.000000,
            15, 15, 4.914617,
            16, 13, -1.142857,
            16, 14, 0.615385,
            16, 15, -0.415584,
            16, 16, 3.628372,
            17, 15, -1.454545,
            17, 16, -1.454545,
            17, 17, 3.575758,
            18, 13, -1.142857,
            18, 15, -1.142857,
            18, 18, 3.265676,
            19, 13, -1.000000,
            19, 14, -1.000000,
            19, 18, -0.068362,
            19, 19, 1.511997,
            20, 14, -0.666667,
            20, 18, -1.098977,
            20, 19, 0.399684,
            20, 20, 2.171910,
            21, 15, -0.727273,
            21, 17, 0.666667,
            21, 18, 0.137014,
            21, 19, 0.007654,
            21, 20, -0.521202,
            21, 21, 2.683388,
            22, 14, -0.666667,
            22, 18, -0.504685,
            22, 19, 0.285758,
            22, 20, -0.414861,
            22, 21, -0.173864,
            22, 22, 1.729035,
            23, 14, -1.230769,
            23, 16, -1.230769,
            23, 18, -0.120919,
            23, 19, 0.323534,
            23, 20, 0.303613,
            23, 21, -0.044126,
            23, 22, 0.055943,
            23, 23, 2.636088,
            24, 13, -1.000000,
            24, 14, -1.000000,
            24, 18, 0.271849,
            24, 19, 0.426317,
            24, 20, 0.247785,
            24, 21, 0.007654,
            24, 22, 0.023785,
            24, 23, 0.323534,
            24, 24, 1.721740,
            25, 15, 0.680851,
            25, 17, -1.333333,
            25, 18, 0.668286,
            25, 19, -0.640499,
            25, 20, -0.197603,
            25, 21, -1.726125,
            25, 22,-0.310405,
            25, 23, -0.800426,
            25, 24, -0.773942,
            25, 25, 4.322353,
            26, 15, -1.361702,
            26, 18, -0.526879,
            26, 19, 0.072813,
            26, 20, 0.219792,
            26, 21, 0.133770,
            26, 22, 0.301092,
            26, 23, 0.215256,
            26, 24, 0.072813,
            26, 25, -0.421433,
            26, 26, 1.636861
};

    // ------------------------------------

    SECTION( "Making G" )
    {
        try
        {
            evoped::Gmat gmat;
            evolm::matrix<double> G;
            std::vector<std::int64_t> g_id;

            gmat.make_matrix("tests/data/allele.dat");
            gmat.get_matrix(G,g_id);

            std::vector<double> Gvect;
            for (size_t i = 0; i < g_id.size(); i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    if ( G(i,j) != 0.0 )
                    {
                        Gvect.push_back( g_id[i] );
                        Gvect.push_back( g_id[j] );
                        Gvect.push_back( G(i,j) );
                    }
                }
            }

            CHECK( G_true.size() == Gvect.size() );

            for (size_t i = 0; i < Gvect.size();)
            {
                CHECK( G_true[i] == Gvect[i] );
                CHECK( G_true[i+1] == Gvect[i+1] );
                CHECK( G_true[i+2] == Catch::Approx(Gvect[i+2]) );
                i = i + 3;
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        catch(const std::string& e)
        {
            std::cerr << e << '\n';
        }                
    }

    SECTION( "Reading G" )
    {
        try
        {
            evoped::Gmat gmat;
            evolm::matrix<double> G;
            std::vector<std::int64_t> g_id;

            gmat.read_matrix("tests/data/g_mat");
            gmat.get_matrix(G,g_id);

            std::vector<double> Gvect;
            for (size_t i = 0; i < g_id.size(); i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    if ( G(i,j) != 0.0 )
                    {
                        Gvect.push_back( g_id[i] );
                        Gvect.push_back( g_id[j] );
                        Gvect.push_back( G(i,j) );
                    }
                }
            }

            CHECK( G_true.size() == Gvect.size() );

            for (size_t i = 0; i < Gvect.size();)
            {
                CHECK( G_true[i] == Gvect[i] );
                CHECK( G_true[i+1] == Gvect[i+1] );
                CHECK( G_true[i+2] == Catch::Approx(Gvect[i+2]) );
                i = i + 3;
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        catch(const std::string& e)
        {
            std::cerr << e << '\n';
        }                
    }

    SECTION( "Scalling G by A22 and (non-sparse) inversion" )
    {
        try
        {
            evoped::Gmat gmat;
            evolm::matrix<double> G;
            std::vector<std::int64_t> g_id;

            gmat.read_matrix("tests/data/g_mat");

            evoped::Amat ap;
            evolm::matrix<double> A22;
            std::vector<std::int64_t> a22_id;

            ap.make_all("tests/data/ped_bkg.dat", "tests/data/typed2");
            ap.get_matrix("A22", A22, a22_id, false);

            gmat.scale_matrix(A22, 0.25);
            gmat.get_matrix(G,g_id);

            double a, b, a_diag, a_ofd, g_diag, g_ofd;

            gmat.get_alpha_beta(a,b,a_diag,a_ofd,g_diag,g_ofd);

            CHECK( Catch::Approx(a) == alpha_true );
            CHECK( Catch::Approx(b) == beta_true );
            CHECK( Catch::Approx(a_diag) == a_diag_true );
            CHECK( Catch::Approx(a_ofd) == a_ofd_true );
            CHECK( Catch::Approx(g_diag) == g_diag_true );
            CHECK( Catch::Approx(g_ofd) == g_ofd_true );

            CHECK( G_scaled_true.size() == G.size() );

            for (size_t i = 0; i < G.size(); i++)
                CHECK( G_scaled_true[i] == Catch::Approx(G[i]) );

            gmat.invert_matrix(true);
            gmat.get_matrix(G,g_id);

            CHECK( iG_true.size() == G.size() );

            for (size_t i = 0; i < G.size(); i++)
                CHECK( iG_true[i] == Catch::Approx(G[i]) );

        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        catch(const std::string& e)
        {
            std::cerr << e << '\n';
        }                
    }

    SECTION( "Sparse inverse of (read + scaling) G" )
    {
        try
        {
            evoped::Gmat gmat;
            evolm::matrix<double> G;
            std::vector<std::int64_t> g_id;

            gmat.read_matrix("tests/data/g_mat");

            evoped::Amat ap;
            evolm::matrix<double> A22;
            std::vector<std::int64_t> a22_id;

            ap.make_all("tests/data/ped_bkg.dat", "tests/data/typed2");
            ap.get_matrix("A22", A22, a22_id, false);

            gmat.scale_matrix(A22, 0.25);

            std::vector<std::int64_t> core_id{ 18, 20, 22, 25 };

            gmat.invert_matrix(core_id);
            gmat.get_matrix(G,g_id);

            CHECK( iG_sparse_true.size() == G.size() );

            for (size_t i = 0; i < G.size(); i++)
                CHECK( iG_sparse_true[i] == Catch::Approx(G[i]) );

        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        catch(const std::string& e)
        {
            std::cerr << e << '\n';
        }                
    }

    SECTION( "Sparse inverse of (making + scaling) G" )
    {
        try
        {
            evoped::Gmat gmat;
            evolm::matrix<double> G;
            std::vector<std::int64_t> g_id;

            gmat.make_matrix("tests/data/allele.dat");

            evoped::Amat ap;
            evolm::matrix<double> A22;
            std::vector<std::int64_t> a22_id;

            ap.make_all("tests/data/ped_bkg.dat", "tests/data/typed2");
            ap.get_matrix("A22", A22, a22_id, false);

            gmat.scale_matrix(A22, 0.25);

            std::vector<std::int64_t> core_id{ 18, 20, 22, 25 };

            gmat.invert_matrix(core_id);
            gmat.get_matrix(G,g_id);

            CHECK( iG_sparse_true.size() == G.size() );

            for (size_t i = 0; i < G.size(); i++)
                CHECK( iG_sparse_true[i] == Catch::Approx(G[i]) );
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        catch(const std::string& e)
        {
            std::cerr << e << '\n';
        }                
    }

    SECTION( "Sparse inverse of H" )
    {
        try
        {
            evoped::Amat ap;
            evolm::matrix<double> iA;
            std::vector<std::int64_t> a_id;
            evolm::matrix<double> iA22;
            std::vector<std::int64_t> ai22_id;
            evolm::matrix<double> A22;
            std::vector<std::int64_t> a22_id;

            ap.make_all("tests/data/ped_bkg.dat", "tests/data/typed2");
            
            ap.get_matrix("iA", iA, a_id, false);
            ap.get_matrix("iA22", iA22, ai22_id, false);
            ap.get_matrix("A22", A22, a22_id, false);

            evoped::Gmat gmat;
            evolm::matrix<double> G;
            std::vector<std::int64_t> g_id;

            //gmat.make_matrix("tests/data/allele.dat");
            gmat.read_matrix("tests/data/g_mat");
            gmat.scale_matrix(A22, 0.25);

            std::vector<std::int64_t> core_id{ 18, 20, 22, 25 };

            gmat.invert_matrix(core_id);
            gmat.get_matrix(G,g_id);

            evoped::Hmat h;
            evolm::matrix<double> H;
            std::vector<std::int64_t> h_id;
            h.make_matrix(iA, a_id, iA22, ai22_id, G, g_id);
            
            h.get_matrix(H,h_id,false);

            std::vector<double> Hvect;
            for (size_t i = 0; i < h_id.size(); i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    if ( H(i,j) != 0.0 )
                    {
                        Hvect.push_back( h_id[i] );
                        Hvect.push_back( h_id[j] );
                        Hvect.push_back( H(i,j) );
                    }
                }
            }

            CHECK( iH_sparse_true.size() == Hvect.size() );

            for (size_t i = 0; i < Hvect.size();)
            {
                CHECK(iH_sparse_true[i] == Hvect[i] );
                CHECK(iH_sparse_true[i+1] == Hvect[i+1] );
                CHECK(iH_sparse_true[i+2] == Catch::Approx(Hvect[i+2]) );
                i = i + 3;
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
        catch(const std::string& e)
        {
            std::cerr << e << '\n';
        }                
    }

}