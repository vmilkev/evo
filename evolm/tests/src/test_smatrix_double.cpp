#include "catch_amalgamated.hpp"
#include "sparse_matrix.hpp"
#include <chrono>

TEST_CASE("Sparse matrix, checking class constructors, type = double")
{
    // Data:
    evolm::smatrix<double> B(11, 11); // symmetric
    B(1, 1) = 5.0000;
    B(10, 1) = 0.9649;
    B(2, 2) = 5.0000;
    B(10, 2) = 0.9595;
    B(3, 3) = 5.0000;
    B(10, 3) = 0.1712;
    B(4, 4) = 5.0000;
    B(5, 4) = 0.0971;
    B(9, 4) = 0.9502;
    B(10, 4) = 0.0344;
    B(4, 5) = 0.0971;
    B(5, 5) = 5.0000;
    B(10, 5) = 0.7547;
    B(6, 6) = 5.0000;
    B(10, 6) = 0.2238;
    B(7, 7) = 5.0000;
    B(9, 7) = 0.1493;
    B(10, 7) = 0.2575;
    B(8, 8) = 5.0000;
    B(10, 8) = 0.4733;
    B(4, 9) = 0.9502;
    B(7, 9) = 0.1493;
    B(9, 9) = 5.0000;
    B(10, 9) = 0.5678;
    B(1, 10) = 0.9649;
    B(2, 10) = 0.9595;
    B(3, 10) = 0.1712;
    B(4, 10) = 0.0344;
    B(5, 10) = 0.7547;
    B(6, 10) = 0.2238;
    B(7, 10) = 0.2575;
    B(8, 10) = 0.4733;
    B(9, 10) = 0.5678;
    B(10, 10) = 5.3371;

    evolm::smatrix<double> symB(11); // symmetric
    symB(1, 1) = 5.0000;
    symB(10, 1) = 0.9649;
    symB(2, 2) = 5.0000;
    symB(10, 2) = 0.9595;
    symB(3, 3) = 5.0000;
    symB(10, 3) = 0.1712;
    symB(4, 4) = 5.0000;
    symB(5, 4) = 0.0971;
    symB(9, 4) = 0.9502;
    symB(10, 4) = 0.0344;
    symB(5, 5) = 5.0000;
    symB(10, 5) = 0.7547;
    symB(6, 6) = 5.0000;
    symB(10, 6) = 0.2238;
    symB(7, 7) = 5.0000;
    symB(9, 7) = 0.1493;
    symB(10, 7) = 0.2575;
    symB(8, 8) = 5.0000;
    symB(10, 8) = 0.4733;
    symB(9, 9) = 5.0000;
    symB(10, 9) = 0.5678;
    symB(10, 10) = 5.3371;

    evolm::smatrix<double> C(11, 11); // NON symmetric
    evolm::smatrix<double> D(11, 11); // NON symmetric, is the transpose of C

    C(1, 1) = 0.0782;
    C(6, 1) = 0.7749;
    C(7, 1) = 0.8173;
    C(8, 1) = 0.8687;
    C(9, 1) = 0.0844;
    C(1, 4) = 0.1839;
    C(6, 4) = 0.9448;
    C(7, 4) = 0.4909;
    C(8, 4) = 0.4893;
    C(9, 4) = 0.3377;
    C(1, 9) = 0.8176;
    C(6, 9) = 0.5328;
    C(7, 9) = 0.3507;
    C(8, 9) = 0.9390;
    C(9, 9) = 0.8759;
    C(1, 10) = 0.6225;
    C(2, 10) = 0.5870;
    C(3, 10) = 0.2077;
    C(4, 10) = 0.3012;
    C(5, 10) = 0.4709;
    C(6, 10) = 0.2305;
    C(7, 10) = 0.8443;
    C(8, 10) = 0.1948;
    C(9, 10) = 0.2259;
    C(10, 10) = 0.1707;

    D(1, 1) = 0.0782;
    D(4, 1) = 0.1839;
    D(9, 1) = 0.8176;
    D(10, 1) = 0.6225;
    D(10, 2) = 0.5870;
    D(10, 3) = 0.2077;
    D(10, 4) = 0.3012;
    D(10, 5) = 0.4709;
    D(1, 6) = 0.7749;
    D(4, 6) = 0.9448;
    D(9, 6) = 0.5328;
    D(10, 6) = 0.2305;
    D(1, 7) = 0.8173;
    D(4, 7) = 0.4909;
    D(9, 7) = 0.3507;
    D(10, 7) = 0.8443;
    D(1, 8) = 0.8687;
    D(4, 8) = 0.4893;
    D(9, 8) = 0.9390;
    D(10, 8) = 0.1948;
    D(1, 9) = 0.0844;
    D(4, 9) = 0.3377;
    D(9, 9) = 0.8759;
    D(10, 9) = 0.2259;
    D(10, 10) = 0.1707;

    evolm::smatrix<double> C2(10, 10); // NON symmetric
    evolm::smatrix<double> D2(10, 10); // NON symmetric, is the transpose of C
    C2(0, 0) = 0.0782;
    C2(5, 0) = 0.7749;
    C2(6, 0) = 0.8173;
    C2(7, 0) = 0.8687;
    C2(8, 0) = 0.0844;
    C2(0, 3) = 0.1839;
    C2(5, 3) = 0.9448;
    C2(6, 3) = 0.4909;
    C2(7, 3) = 0.4893;
    C2(8, 3) = 0.3377;
    C2(0, 8) = 0.8176;
    C2(5, 8) = 0.5328;
    C2(6, 8) = 0.3507;
    C2(7, 8) = 0.9390;
    C2(8, 8) = 0.8759;
    C2(0, 9) = 0.6225;
    C2(1, 9) = 0.5870;
    C2(2, 9) = 0.2077;
    C2(3, 9) = 0.3012;
    C2(4, 9) = 0.4709;
    C2(5, 9) = 0.2305;
    C2(6, 9) = 0.8443;
    C2(7, 9) = 0.1948;
    C2(8, 9) = 0.2259;
    C2(9, 9) = 0.1707;

    D2(0, 0) = 0.0782;
    D2(3, 0) = 0.1839;
    D2(8, 0) = 0.8176;
    D2(9, 0) = 0.6225;
    D2(9, 1) = 0.5870;
    D2(9, 2) = 0.2077;
    D2(9, 3) = 0.3012;
    D2(9, 4) = 0.4709;
    D2(0, 5) = 0.7749;
    D2(3, 5) = 0.9448;
    D2(8, 5) = 0.5328;
    D2(9, 5) = 0.2305;
    D2(0, 6) = 0.8173;
    D2(3, 6) = 0.4909;
    D2(8, 6) = 0.3507;
    D2(9, 6) = 0.8443;
    D2(0, 7) = 0.8687;
    D2(0, 8) = 0.0844;
    D2(3, 7) = 0.4893;
    D2(8, 7) = 0.9390;
    D2(9, 7) = 0.1948;
    D2(3, 8) = 0.3377;
    D2(8, 8) = 0.8759;
    D2(9, 8) = 0.2259;
    D2(9, 9) = 0.1707;

    evolm::smatrix<double> D_mul_C(10, 10);
    D_mul_C(0, 0) = 2.0363;
    D_mul_C(3, 0) = 1.6012;
    D_mul_C(8, 0) = 1.6531;
    D_mul_C(9, 0) = 1.1056;
    D_mul_C(0, 3) = 1.6012;
    D_mul_C(3, 3) = 1.5208;
    D_mul_C(8, 3) = 1.5812;
    D_mul_C(9, 3) = 0.9183;
    D_mul_C(0, 8) = 1.6531;
    D_mul_C(3, 8) = 1.5812;
    D_mul_C(8, 8) = 2.7244;
    D_mul_C(9, 8) = 1.3087;
    D_mul_C(0, 9) = 1.1056;
    D_mul_C(3, 9) = 0.9183;
    D_mul_C(8, 9) = 1.3087;
    D_mul_C(9, 9) = 1.9719;

    evolm::smatrix<double> C_mul_D(10, 10); // note to change indexes!
    C_mul_D(0, 0) = 1.0959;
    C_mul_D(1, 0) = 0.3654;
    C_mul_D(2, 0) = 0.1293;
    C_mul_D(3, 0) = 0.1875;
    C_mul_D(4, 0) = 0.2931;
    C_mul_D(5, 0) = 0.8135;
    C_mul_D(6, 0) = 0.9665;
    C_mul_D(7, 0) = 1.0469;
    C_mul_D(8, 0) = 0.9255;
    C_mul_D(9, 0) = 0.1063;
    C_mul_D(0, 1) = 0.3654;
    C_mul_D(1, 1) = 0.3446;
    C_mul_D(2, 1) = 0.1220;
    C_mul_D(3, 1) = 0.1768;
    C_mul_D(4, 1) = 0.2765;
    C_mul_D(5, 1) = 0.1353;
    C_mul_D(6, 1) = 0.4956;
    C_mul_D(7, 1) = 0.1143;
    C_mul_D(8, 1) = 0.1326;
    C_mul_D(9, 1) = 0.1002;
    C_mul_D(0, 2) = 0.1293;
    C_mul_D(1, 2) = 0.1220;
    C_mul_D(2, 2) = 0.0432;
    C_mul_D(3, 2) = 0.0626;
    C_mul_D(4, 2) = 0.0978;
    C_mul_D(5, 2) = 0.0479;
    C_mul_D(6, 2) = 0.1754;
    C_mul_D(7, 2) = 0.0405;
    C_mul_D(8, 2) = 0.0469;
    C_mul_D(9, 2) = 0.0355;
    C_mul_D(0, 3) = 0.1875;
    C_mul_D(1, 3) = 0.1768;
    C_mul_D(2, 3) = 0.0626;
    C_mul_D(3, 3) = 0.0907;
    C_mul_D(4, 3) = 0.1419;
    C_mul_D(5, 3) = 0.0694;
    C_mul_D(6, 3) = 0.2543;
    C_mul_D(7, 3) = 0.0587;
    C_mul_D(8, 3) = 0.0681;
    C_mul_D(9, 3) = 0.0514;
    C_mul_D(0, 4) = 0.2931;
    C_mul_D(1, 4) = 0.2765;
    C_mul_D(2, 4) = 0.0978;
    C_mul_D(3, 4) = 0.1419;
    C_mul_D(4, 4) = 0.2218;
    C_mul_D(5, 4) = 0.1085;
    C_mul_D(6, 4) = 0.3976;
    C_mul_D(7, 4) = 0.0917;
    C_mul_D(8, 4) = 0.1064;
    C_mul_D(9, 4) = 0.0804;
    C_mul_D(0, 5) = 0.8135;
    C_mul_D(1, 5) = 0.1353;
    C_mul_D(2, 5) = 0.0479;
    C_mul_D(3, 5) = 0.0694;
    C_mul_D(4, 5) = 0.1085;
    C_mul_D(5, 5) = 1.8301;
    C_mul_D(6, 5) = 1.4786;
    C_mul_D(7, 5) = 1.6806;
    C_mul_D(8, 5) = 0.9033;
    C_mul_D(9, 5) = 0.0393;
    C_mul_D(0, 6) = 0.9665;
    C_mul_D(1, 6) = 0.4956;
    C_mul_D(2, 6) = 0.1754;
    C_mul_D(3, 6) = 0.2543;
    C_mul_D(4, 6) = 0.3976;
    C_mul_D(5, 6) = 1.4786;
    C_mul_D(6, 6) = 1.7448;
    C_mul_D(7, 6) = 1.4439;
    C_mul_D(8, 6) = 0.7327;
    C_mul_D(9, 6) = 0.1441;
    C_mul_D(0, 7) = 1.0469;
    C_mul_D(1, 7) = 0.1143;
    C_mul_D(2, 7) = 0.0405;
    C_mul_D(3, 7) = 0.0587;
    C_mul_D(4, 7) = 0.0917;
    C_mul_D(5, 7) = 1.6806;
    C_mul_D(6, 7) = 1.4439;
    C_mul_D(7, 7) = 1.9137;
    C_mul_D(8, 7) = 1.1051;
    C_mul_D(9, 7) = 0.0332;
    C_mul_D(0, 8) = 0.9255;
    C_mul_D(1, 8) = 0.1326;
    C_mul_D(2, 8) = 0.0469;
    C_mul_D(3, 8) = 0.0681;
    C_mul_D(4, 8) = 0.1064;
    C_mul_D(5, 8) = 0.9033;
    C_mul_D(6, 8) = 0.7327;
    C_mul_D(7, 8) = 1.1051;
    C_mul_D(8, 8) = 0.9395;
    C_mul_D(9, 8) = 0.0386;
    C_mul_D(0, 9) = 0.1063;
    C_mul_D(1, 9) = 0.1002;
    C_mul_D(2, 9) = 0.0355;
    C_mul_D(3, 9) = 0.0514;
    C_mul_D(4, 9) = 0.0804;
    C_mul_D(5, 9) = 0.0393;
    C_mul_D(6, 9) = 0.1441;
    C_mul_D(7, 9) = 0.0332;
    C_mul_D(8, 9) = 0.0386;
    C_mul_D(9, 9) = 0.0291;

    evolm::smatrix<double> D_min_C(11, 11);

    D_min_C(4, 1) = 0.1839;
    D_min_C(6, 1) = -0.7749;
    D_min_C(7, 1) = -0.8173;
    D_min_C(8, 1) = -0.8687;
    D_min_C(9, 1) = 0.7332;
    D_min_C(10, 1) = 0.6225;
    D_min_C(10, 2) = 0.5870;
    D_min_C(10, 3) = 0.2077;
    D_min_C(1, 4) = -0.1839;
    D_min_C(6, 4) = -0.9448;
    D_min_C(7, 4) = -0.4909;
    D_min_C(8, 4) = -0.4893;
    D_min_C(9, 4) = -0.3377;
    D_min_C(10, 4) = 0.3012;
    D_min_C(10, 5) = 0.4709;
    D_min_C(1, 6) = 0.7749;
    D_min_C(4, 6) = 0.9448;
    D_min_C(9, 6) = 0.5328;
    D_min_C(10, 6) = 0.2305;
    D_min_C(1, 7) = 0.8173;
    D_min_C(4, 7) = 0.4909;
    D_min_C(9, 7) = 0.3507;
    D_min_C(10, 7) = 0.8443;
    D_min_C(1, 8) = 0.8687;
    D_min_C(4, 8) = 0.4893;
    D_min_C(9, 8) = 0.9390;
    D_min_C(10, 8) = 0.1948;
    D_min_C(1, 9) = -0.7332;
    D_min_C(4, 9) = 0.3377;
    D_min_C(6, 9) = -0.5328;
    D_min_C(7, 9) = -0.3507;
    D_min_C(8, 9) = -0.9390;
    D_min_C(10, 9) = 0.2259;
    D_min_C(1, 10) = -0.6225;
    D_min_C(2, 10) = -0.5870;
    D_min_C(3, 10) = -0.2077;
    D_min_C(4, 10) = -0.3012;
    D_min_C(5, 10) = -0.4709;
    D_min_C(6, 10) = -0.2305;
    D_min_C(7, 10) = -0.8443;
    D_min_C(8, 10) = -0.1948;
    D_min_C(9, 10) = -0.2259;

    evolm::smatrix<double> D_plus_C(11, 11);

    D_plus_C(1, 1) = 0.1564;
    D_plus_C(4, 1) = 0.1839;
    D_plus_C(6, 1) = 0.7749;
    D_plus_C(7, 1) = 0.8173;
    D_plus_C(8, 1) = 0.8687;
    D_plus_C(9, 1) = 0.902;
    D_plus_C(10, 1) = 0.6225;
    D_plus_C(10, 2) = 0.5870;
    D_plus_C(10, 3) = 0.2077;
    D_plus_C(1, 4) = 0.1839;
    D_plus_C(6, 4) = 0.9448;
    D_plus_C(7, 4) = 0.4909;
    D_plus_C(8, 4) = 0.4893;
    D_plus_C(9, 4) = 0.3377;
    D_plus_C(10, 4) = 0.3012;
    D_plus_C(10, 5) = 0.4709;
    D_plus_C(1, 6) = 0.7749;
    D_plus_C(4, 6) = 0.9448;
    D_plus_C(9, 6) = 0.5328;
    D_plus_C(10, 6) = 0.2305;
    D_plus_C(1, 7) = 0.8173;
    D_plus_C(4, 7) = 0.4909;
    D_plus_C(9, 7) = 0.3507;
    D_plus_C(10, 7) = 0.8443;
    D_plus_C(1, 8) = 0.8687;
    D_plus_C(4, 8) = 0.4893;
    D_plus_C(9, 8) = 0.9390;
    D_plus_C(10, 8) = 0.1948;
    D_plus_C(1, 9) = 0.902;
    D_plus_C(4, 9) = 0.3377;
    D_plus_C(6, 9) = 0.5328;
    D_plus_C(7, 9) = 0.3507;
    D_plus_C(8, 9) = 0.9390;
    D_plus_C(9, 9) = 1.7518;
    D_plus_C(10, 9) = 0.2259;
    D_plus_C(1, 10) = 0.6225;
    D_plus_C(2, 10) = 0.5870;
    D_plus_C(3, 10) = 0.2077;
    D_plus_C(4, 10) = 0.3012;
    D_plus_C(5, 10) = 0.4709;
    D_plus_C(6, 10) = 0.2305;
    D_plus_C(7, 10) = 0.8443;
    D_plus_C(8, 10) = 0.1948;
    D_plus_C(9, 10) = 0.2259;
    D_plus_C(10, 10) = 0.3414;

    SECTION("Default constructor")
    {
        evolm::smatrix<double> M;

        CHECK(M.empty() == true);
        CHECK(M.ncols() == 0);
        CHECK(M.nrows() == 0);
    }

    SECTION("Symmetric & rectangular matrix onstructor")
    {
        evolm::smatrix<double> M(5);

        CHECK(M.empty() == true);
        CHECK(M.ncols() == 5);
        CHECK(M.nrows() == 5);
        CHECK(M.max_key() == (5 * 5 + 5) / 2);

        M.resize(5, 5);

        CHECK(M.empty() == true);
        CHECK(M.ncols() == 5);
        CHECK(M.nrows() == 5);
        CHECK(M.max_key() == 5 * 5);
    }

    SECTION("nonzero() method")
    {
        size_t n = 100;
        evolm::smatrix<double> M(n);

        auto start = std::chrono::high_resolution_clock::now();

        for (size_t i = 0; i < M.nrows(); i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                CHECK(M.nonzero(i, j) == false);
            }
        }

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        // std::cout <<"symmetrical (milliseconds): "<< duration.count() <<" for n elements: "<< (n*n+n)/2 << std::endl;
    }

    SECTION("&operator=(const smatrix &rhs)")
    {
        evolm::smatrix<double> M(5, 4);

        M(0, 0) = 1.0;
        M(2, 0) = 2.0;
        M(2, 2) = 3.0;
        M(4, 1) = 4;
        M(4, 2) = 5.0;

        evolm::smatrix<double> N = M;

        evolm::smatrix<double> L;

        CHECK(N.size() == 5);
        CHECK(M.size() == 5);

        CHECK(N(0, 0) == 1.0);
        CHECK(N(2, 0) == 2.0);
        CHECK(N(2, 2) == 3.0);
        CHECK(N(4, 1) == 4.0);
        CHECK(N(4, 2) == 5.0);

        evolm::smatrix<double> K;

        CHECK(K.size() == 0);
        K = M;

        CHECK(K.size() == 5);
        CHECK(M.size() == 5);

        CHECK(K(0, 0) == 1.0);
        CHECK(K(2, 0) == 2.0);
        CHECK(K(2, 2) == 3.0);
        CHECK(K(4, 1) == 4.0);
        CHECK(K(4, 2) == 5.0);
    }

    SECTION("smatrix(const smatrix &obj)")
    {
        evolm::smatrix<double> M(5, 4);

        M(0, 0) = 1.0;
        M(2, 0) = 2.0;
        M(2, 2) = 3.0;
        M(4, 1) = 4;
        M(4, 2) = 5.0;

        evolm::smatrix<double> N(M);

        CHECK(M.size() == 5);

        CHECK(N(0, 0) == 1.0);
        CHECK(N(2, 0) == 2.0);
        CHECK(N(2, 2) == 3.0);
        CHECK(N(4, 1) == 4.0);
        CHECK(N(4, 2) == 5.0);

        CHECK(N.size() == 5);
        N.clear();
        CHECK(N.size() == 0);
    }

    SECTION("fwrite() / fread()")
    {
        evolm::smatrix<double> M(5, 4);

        M(0, 0) = 1.0;
        M(2, 0) = 2.0;
        M(2, 2) = 3.0;
        M(4, 1) = 4;
        M(4, 2) = 5.0;

        CHECK(M.size() == 5);
        CHECK(M.empty() == false);

        size_t count = 0;
        for (size_t i = 0; i < M.nrows(); i++)
        {
            for (size_t j = 0; j < M.ncols(); j++)
            {
                if (M.nonzero(i, j))
                {
                    CHECK(M(i, j) == ++count);
                }
            }
        }

        CHECK(M.size() == 5);

        CHECK(M.size() == 5);

        M.fwrite();

        CHECK(M.size() == 0);
        CHECK(M.empty() == true);

        M.fread();

        CHECK(M.size() == 5);

        count = 0;
        for (size_t i = 0; i < M.nrows(); i++)
        {
            for (size_t j = 0; j < M.ncols(); j++)
            {
                if (M.nonzero(i, j))
                {
                    CHECK(M(i, j) == ++count);
                }
            }
        }

        M.fclear();
    }

    SECTION("transpose()")
    {
        evolm::smatrix<double> M(5, 4);
        evolm::smatrix<double> tM(4, 5);

        M(0, 0) = 1.0;
        M(2, 0) = 2.0;
        M(2, 2) = 3.0;
        M(4, 1) = 4;
        M(4, 2) = 5.0;

        tM(0, 0) = 1.0;
        tM(0, 2) = 2.0;
        tM(2, 2) = 3.0;
        tM(1, 4) = 4;
        tM(2, 4) = 5.0;

        M.transpose();

        CHECK(M.size() == 5);
        CHECK(M.ncols() == 5);
        CHECK(M.nrows() == 4);

        CHECK(M(0, 0) == 1.0);
        CHECK(M(0, 2) == 2.0);
        CHECK(M(2, 2) == 3.0);
        CHECK(M(1, 4) == 4.0);
        CHECK(M(2, 4) == 5.0);

        M.resize();
        CHECK(M.size() == 0);

        evolm::smatrix<double> tB(B);
    }

    SECTION("transpose() symmetric in full format")
    {
        evolm::smatrix<double> M(B);

        M.transpose();

        CHECK(M.size() == B.size());
        CHECK(M.ncols() == 11);
        CHECK(M.nrows() == 11);

        for (size_t i = 0; i < B.nrows(); i++)
            for (size_t j = 0; j < B.ncols(); j++)
                CHECK(M(i, j) == B(i, j));

        M.resize();
        CHECK(M.size() == 0);
    }

    SECTION("transpose() non-symmetric")
    {
        evolm::smatrix<double> M(C);

        M.transpose();

        CHECK(M.size() == D.size());
        CHECK(M.ncols() == 11);
        CHECK(M.nrows() == 11);

        for (size_t i = 0; i < D.nrows(); i++)
            for (size_t j = 0; j < D.ncols(); j++)
                CHECK(M(i, j) == D(i, j));

        M.resize();
        CHECK(M.size() == 0);
    }

    SECTION("transpose(), big")
    {
        size_t n = 5000;
        evolm::smatrix<double> M(n, n);
        evolm::smatrix<double> tM(n, n);

        for (size_t i = 0; i < n * n - 1; i++)
        {
            if (i % 100 == 0)
                M[i] = i + 1;
        }

        auto start = std::chrono::high_resolution_clock::now();

        M.transpose();

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout << "Transpose takes (milliseconds): " << duration.count() << " for row, col: " << n << " " << n << std::endl;
    }

    SECTION("rectosym() and symtorec()")
    {
        size_t n = 6;
        size_t max_sym = (n * n + n) / 2;
        size_t max_rec = n * n;
        evolm::smatrix<double> M(n);
        evolm::smatrix<double> N(n, n);

        M(1, 1) = 2;
        M(3, 1) = 7;
        M(4, 0) = 10;
        M(5, 2) = 17;
        M(5, 5) = 20;

        CHECK(M.size() == 5);
        CHECK(M.max_key() == max_sym);

        N = M;

        CHECK(N.size() == 5);
        CHECK(N.max_key() == max_sym);

        N.symtorec();

        CHECK(N.size() == 8);
        CHECK(N.max_key() == max_rec);

        evolm::smatrix<double> K(N); // K = (N); // Check copy constructor here !!! K(N) does not work!

        CHECK(K.size() == 8);
        CHECK(K.max_key() == max_rec);

        K.rectosym();

        CHECK(K.size() == 5);
        CHECK(K.max_key() == max_sym);

        for (size_t i = 0; i < K.nrows(); i++)
        {
            for (size_t j = 0; j <= i; j++)
                if (K.nonzero(i, j))
                    CHECK(K(i, j) == M(i, j));
        }
    }

    SECTION("rectosym()")
    {
        evolm::smatrix<double> M(B); // rectangular symmetric

        M.rectosym();

        CHECK(M.size() == symB.size());
        CHECK(M.ncols() == 11);
        CHECK(M.nrows() == 11);

        for (size_t i = 0; i < symB.nrows(); i++)
            for (size_t j = 0; j <= i; j++)
                CHECK(M(i, j) == symB(i, j));

        M.resize();
        CHECK(M.size() == 0);
    }

    SECTION("symtorec()")
    {
        evolm::smatrix<double> M(symB); // symmetric

        M.symtorec();

        CHECK(M.size() == B.size());
        CHECK(M.ncols() == 11);
        CHECK(M.nrows() == 11);

        for (size_t i = 0; i < B.nrows(); i++)
            for (size_t j = 0; j < B.ncols(); j++)
                CHECK(M(i, j) == B(i, j));

        M.resize();
        CHECK(M.size() == 0);
    }

    SECTION("operator*, D2 * C2")
    {
        evolm::smatrix<double> M;

        M = D2 * C2;

        //D2.print("D2");
        //C2.print("C2");
        //M.print("D2*C2");
        //D_mul_C.print("D_mul_C");

        CHECK(M.size() == D_mul_C.size());
        CHECK(M.ncols() == 10);
        CHECK(M.nrows() == 10);

        for (size_t i = 0; i < M.nrows(); i++)
            for (size_t j = 0; j < M.ncols(); j++)
                CHECK(Catch::Approx( M(i, j) ).margin(0.0005).epsilon(1e-12) == D_mul_C(i, j));

        M.resize();
        CHECK(M.size() == 0);
    }

    SECTION("operator*, C2 * D2")
    {
        evolm::smatrix<double> M;

        M = C2 * D2;

        //C2.print("C2");
        //D2.print("D2");
        //M.print("C2*D2");
        //C_mul_D.print("C_mul_D");

        CHECK(M.size() == C_mul_D.size());
        CHECK(M.ncols() == 10);
        CHECK(M.nrows() == 10);

        for (size_t i = 0; i < M.nrows(); i++)
            for (size_t j = 0; j < M.ncols(); j++)
                CHECK(Catch::Approx( M(i, j) ).margin(0.0005).epsilon(1e-12) == C_mul_D(i, j));

        M.resize();
        CHECK(M.size() == 0);
    }

    SECTION("operator*")
    {
        evolm::smatrix<double> k(3, 4);
        evolm::smatrix<double> l(4, 3);
        evolm::smatrix<double> e;

        size_t count = 0;
        for (size_t i = 0; i < k.nrows(); i++)
        {
            for (size_t j = 0; j < k.ncols(); j++)
            {
                k(i, j) = ++count;
            }
        }
        count = 0;
        for (size_t i = 0; i < l.nrows(); i++)
        {
            for (size_t j = 0; j < l.ncols(); j++)
            {
                l(i, j) = ++count;
            }
        }
        e = k * l;

        CHECK(e.size() == 9);
        CHECK(e(0, 0) == 70);
        CHECK(e(0, 1) == 80);
        CHECK(e(0, 2) == 90);
        CHECK(e(1, 0) == 158);
        CHECK(e(1, 1) == 184);
        CHECK(e(1, 2) == 210);
        CHECK(e(2, 0) == 246);
        CHECK(e(2, 1) == 288);
        CHECK(e(2, 2) == 330);

        size_t a = 10, b = 13, c = 15;
        evolm::smatrix<double> m(10, 13);
        evolm::smatrix<double> n(13, 15);
        evolm::smatrix<double> res;

        m(2, 0) = 1;
        m(6, 1) = 2;
        m(9, 1) = 8;
        m(9, 3) = 9;
        m(1, 5) = 3;
        m(9, 6) = 10;
        m(7, 7) = 4;
        m(7, 10) = 5;
        m(5, 12) = 6;

        n(7, 0) = 4;
        n(4, 2) = 1;
        n(4, 5) = 2;
        n(4, 6) = 3;
        n(4, 8) = 4;
        n(9, 11) = 5;
        n(9, 13) = 6;
        n(12, 14) = 7;

        res = m * n;

        CHECK(res.size() == 2);
        CHECK(res.max_key() == a * c);

        CHECK(res(7, 0) == 16);
        CHECK(res(5, 14) == 42);

        // check time
        size_t s = 1000;
        evolm::smatrix<double> M(s, s);
        evolm::smatrix<double> result;

        auto start = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < s * s - 1; i++)
        {
            // if ( i%100 == 0 )
            M[i] = i + 1;
        }
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout << "duration (milliseconds): " << duration.count() << " for n elements: " << s * s << std::endl;

        start = std::chrono::high_resolution_clock::now();
        evolm::smatrix<double> N(M);

        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout << "duration (milliseconds): " << duration.count() << " for n elements: " << s * s << std::endl;

        start = std::chrono::high_resolution_clock::now();

        result = M * N;

        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout << "Duration (milliseconds): " << duration.count() << " for n elements: " << s * s << std::endl;
    }

    SECTION("operator+, on small")
    {
        evolm::smatrix<double> M; // result

        M = D + C;

        CHECK(M.size() == D_plus_C.size());
        CHECK(M.ncols() == 11);
        CHECK(M.nrows() == 11);

        for (size_t i = 0; i < M.nrows(); i++)
            for (size_t j = 0; j < M.ncols(); j++)
                CHECK(M(i, j) == D_plus_C(i, j));

        M.resize();
        CHECK(M.size() == 0);
    }

    SECTION("operator+")
    {
        size_t n = 1000;
        evolm::smatrix<double> k(3, 4);
        evolm::smatrix<double> l(3, 4);
        evolm::smatrix<double> e;

        size_t count = 0;
        for (size_t i = 0; i < k.nrows(); i++)
        {
            for (size_t j = 0; j < k.ncols(); j++)
            {
                k(i, j) = ++count;
            }
        }
        count = 0;
        for (size_t i = 0; i < l.nrows(); i++)
        {
            for (size_t j = 0; j < l.ncols(); j++)
            {
                l(i, j) = ++count * 2;
            }
        }
        e = k + l;

        CHECK(e.size() == 12);
        CHECK(e(0, 0) == 3);
        CHECK(e(0, 1) == 6);
        CHECK(e(0, 2) == 9);
        CHECK(e(0, 3) == 12);
        CHECK(e(1, 0) == 15);
        CHECK(e(1, 1) == 18);
        CHECK(e(1, 2) == 21);
        CHECK(e(1, 3) == 24);
        CHECK(e(2, 0) == 27);
        CHECK(e(2, 1) == 30);
        CHECK(e(2, 2) == 33);
        CHECK(e(2, 3) == 36);

        e.resize();
        l.resize(3, 4);

        count = 0;
        for (size_t i = 0; i < l.nrows() * l.ncols() - 1; i++)
        {
            if (i % 2 == 0)
                l[i] = ++count;
        }

        e = k + l;

        CHECK(e.size() == 12);
        CHECK(e(0, 0) == 2);
        CHECK(e(0, 1) == 2);
        CHECK(e(0, 2) == 5);
        CHECK(e(0, 3) == 4);
        CHECK(e(1, 0) == 8);
        CHECK(e(1, 1) == 6);
        CHECK(e(1, 2) == 11);
        CHECK(e(1, 3) == 8);
        CHECK(e(2, 0) == 14);
        CHECK(e(2, 1) == 10);
        CHECK(e(2, 2) == 17);
        CHECK(e(2, 3) == 12);
    }

    SECTION("operator-, on small")
    {
        evolm::smatrix<double> M; // result

        M = D - C;

        CHECK(M.size() == D_min_C.size());
        CHECK(M.ncols() == 11);
        CHECK(M.nrows() == 11);

        for (size_t i = 0; i < M.nrows(); i++)
            for (size_t j = 0; j < M.ncols(); j++)
                CHECK(M(i, j) == D_min_C(i, j));

        M.resize();
        CHECK(M.size() == 0);
    }

    SECTION("operator-")
    {
        size_t n = 1000;
        evolm::smatrix<double> k(3, 4);
        evolm::smatrix<double> l(3, 4);
        evolm::smatrix<double> e;

        size_t count = 0;
        for (size_t i = 0; i < k.nrows(); i++)
        {
            for (size_t j = 0; j < k.ncols(); j++)
            {
                k(i, j) = ++count;
            }
        }
        count = 0;
        for (size_t i = 0; i < l.nrows(); i++)
        {
            for (size_t j = 0; j < l.ncols(); j++)
            {
                l(i, j) = ++count * 2;
            }
        }

        e = k - l;

        CHECK(e.size() == 12);
        CHECK(e(0, 0) == -1);
        CHECK(e(0, 1) == -2);
        CHECK(e(0, 2) == -3);
        CHECK(e(0, 3) == -4);
        CHECK(e(1, 0) == -5);
        CHECK(e(1, 1) == -6);
        CHECK(e(1, 2) == -7);
        CHECK(e(1, 3) == -8);
        CHECK(e(2, 0) == -9);
        CHECK(e(2, 1) == -10);
        CHECK(e(2, 2) == -11);
        CHECK(e(2, 3) == -12);

        e.resize();
        l.resize(3, 4);

        count = 0;
        for (size_t i = 0; i < l.nrows() * l.ncols() - 1; i++)
        {
            if (i % 2 == 0)
                l[i] = ++count;
        }

        e = k - l;

        CHECK(e.size() == 11);
        CHECK(e(0, 1) == 2);
        CHECK(e(0, 2) == 1);
        CHECK(e(0, 3) == 4);
        CHECK(e(1, 0) == 2);
        CHECK(e(1, 1) == 6);
        CHECK(e(1, 2) == 3);
        CHECK(e(1, 3) == 8);
        CHECK(e(2, 0) == 4);
        CHECK(e(2, 1) == 10);
        CHECK(e(2, 2) == 5);
        CHECK(e(2, 3) == 12);
    }
}