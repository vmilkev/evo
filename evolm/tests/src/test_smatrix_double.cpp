#include "catch_amalgamated.hpp"
#include "sparse_matrix.hpp"
#include <chrono>

TEST_CASE("Sparse matrix, checking class constructors, type = double")
{
    SECTION("Default constructor")
    {
        evolm::smatrix<double> M;

        CHECK( M.empty() == true );
        CHECK( M.ncols() == 0 );
        CHECK( M.nrows() == 0 );
    }

    SECTION("Symmetric & rectangular matrix onstructor")
    {
        evolm::smatrix<double> M(5);

        CHECK( M.empty() == true );
        CHECK( M.ncols() == 5 );
        CHECK( M.nrows() == 5 );
        CHECK( M.max_key() == (5*5+5)/2 );

        M.resize(5,5);

        CHECK( M.empty() == true );
        CHECK( M.ncols() == 5 );
        CHECK( M.nrows() == 5 );
        CHECK( M.max_key() == 5*5 );
    }

    SECTION("EXIST() method")
    {
        size_t n = 100;
        evolm::smatrix<double> M(n);

        auto start = std::chrono::high_resolution_clock::now();

        for (size_t i = 0; i < M.nrows(); i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                CHECK( M.nonzero(i,j) == false );
            }
        }

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        //std::cout <<"symmetrical (milliseconds): "<< duration.count() <<" for n elements: "<< (n*n+n)/2 << std::endl;
    }

    SECTION("&operator=(const smatrix &rhs)")
    {
        evolm::smatrix<double> M(5,4);

        M(0,0) = 1.0;
        M(2,0) = 2.0;
        M(2,2) = 3.0;
        M(4,1) = 4;
        M(4,2) = 5.0;

        evolm::smatrix<double> N = M;

        evolm::smatrix<double> L;

        CHECK( N.size() == 5 );
        CHECK( M.size() == 5 );

        CHECK( N(0,0) == 1.0 );
        CHECK( N(2,0) == 2.0 );
        CHECK( N(2,2) == 3.0 );
        CHECK( N(4,1) == 4.0 );
        CHECK( N(4,2) == 5.0 );

        evolm::smatrix<double> K;

        CHECK( K.size() == 0 );
        K = M;

        CHECK( K.size() == 5 );
        CHECK( M.size() == 5 );

        CHECK( K(0,0) == 1.0 );
        CHECK( K(2,0) == 2.0 );
        CHECK( K(2,2) == 3.0 );
        CHECK( K(4,1) == 4.0 );
        CHECK( K(4,2) == 5.0 );
    }

    SECTION("smatrix(const smatrix &obj)")
    {
        evolm::smatrix<double> M(5,4);

        M(0,0) = 1.0;
        M(2,0) = 2.0;
        M(2,2) = 3.0;
        M(4,1) = 4;
        M(4,2) = 5.0;

        evolm::smatrix<double> N(M);

        CHECK( M.size() == 5 );

        CHECK( N(0,0) == 1.0 );
        CHECK( N(2,0) == 2.0 );
        CHECK( N(2,2) == 3.0 );
        CHECK( N(4,1) == 4.0 );
        CHECK( N(4,2) == 5.0 );

        CHECK( N.size() == 5 );
        N.clear();
        CHECK( N.size() == 0 );
    }

    SECTION("fwrite() / fread()")
    {
        evolm::smatrix<double> M(5,4);

        M(0,0) = 1.0;
        M(2,0) = 2.0;
        M(2,2) = 3.0;
        M(4,1) = 4;
        M(4,2) = 5.0;

        CHECK( M.size() == 5 );
        CHECK( M.empty() == false );

        size_t count = 0;
        for (size_t i = 0; i < M.nrows(); i++)
        {
            for (size_t j = 0; j < M.ncols(); j++)
            {
                if (M.nonzero(i,j))
                {
                    CHECK( M(i,j) == ++count );
                }
            }
        }
        
        CHECK( M.size() == 5 );
        //M.print("M before writing");
        CHECK( M.size() == 5 );

        M.fwrite();

        CHECK( M.size() == 0 );
        CHECK( M.empty() == true );

        M.fread();

        CHECK( M.size() == 5 );

        count = 0;
        for (size_t i = 0; i < M.nrows(); i++)
        {
            for (size_t j = 0; j < M.ncols(); j++)
            {
                if (M.nonzero(i,j))
                {
                    CHECK( M(i,j) == ++count );
                }
            }
        }

        //M.print("M after reading");

        M.fclear();
    }

    SECTION("transpose()")
    {
        evolm::smatrix<double> M(5,4);
        evolm::smatrix<double> tM(4,5);

        M(0,0) = 1.0;
        M(2,0) = 2.0;
        M(2,2) = 3.0;
        M(4,1) = 4;
        M(4,2) = 5.0;

        tM(0,0) = 1.0;
        tM(0,2) = 2.0;
        tM(2,2) = 3.0;
        tM(1,4) = 4;
        tM(2,4) = 5.0;

        //M.print("M before transpose");
        M.transpose();
        //M.print("M after transpose");

        CHECK( M.size() == 5 );
        CHECK( M.ncols() == 5 );
        CHECK( M.nrows() == 4 );

        CHECK( M(0,0) == 1.0 );
        CHECK( M(0,2) == 2.0 );
        CHECK( M(2,2) == 3.0 );
        CHECK( M(1,4) == 4.0 );
        CHECK( M(2,4) == 5.0 );

        M.resize();
        CHECK( M.size() == 0 );
    }

    SECTION("transpose(), big")
    {
        size_t n = 1000;
        evolm::smatrix<double> M(n,n);
        evolm::smatrix<double> tM(n,n);

        for (size_t i = 0; i < n*n-1; i++)
        {
            if ( i%100 == 0 )
                M[i] = i+1;
        }

        auto start = std::chrono::high_resolution_clock::now();

        M.transpose();

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout <<"transpose takes (milliseconds): "<< duration.count() <<" for row, col: "<< n<<" "<<n << std::endl;
    }

    SECTION("rectosym() and symtorec()")
    {
        size_t n = 6;
        size_t max_sym = (n * n + n) / 2;
        size_t max_rec = n * n;
        evolm::smatrix<double> M(n);
        evolm::smatrix<double> N(n,n);

        M(1,1) = 2;
        M(3,1) = 7;
        M(4,0) = 10;
        M(5,2) = 17;
        M(5,5) = 20;

        //M.print("M");

        CHECK( M.size() == 5 );
        CHECK( M.max_key() == max_sym );

        N = M;

        //N.print("N = M");

        CHECK( N.size() == 5 );
        CHECK( N.max_key() == max_sym );

        N.symtorec();

        //N.print("N => rec");

        CHECK( N.size() == 8 );
        CHECK( N.max_key() == max_rec );

        evolm::smatrix<double> K(N);// K = (N); // Check copy constructor here !!! K(N) does not work!

        //K.print("K = N");

        CHECK( K.size() == 8 );
        CHECK( K.max_key() == max_rec );

        K.rectosym();

        //K.print("K => sym");

        CHECK( K.size() == 5 );
        CHECK( K.max_key() == max_sym );

        for (size_t i = 0; i < K.nrows(); i++)
        {
            for (size_t j = 0; j <= i; j++)
                if (K.nonzero(i,j))
                    CHECK( K(i,j) == M(i,j) );
        }
    }

    SECTION("operator*")
    {
        evolm::smatrix<double> k(3,4);
        evolm::smatrix<double> l(4,3);
        evolm::smatrix<double> e;

        size_t count = 0;
        for (size_t i = 0; i < k.nrows(); i++)
        {
            for (size_t j = 0; j < k.ncols(); j++)
            {
                k(i,j) = ++count;
            }
        }
        count = 0;
        for (size_t i = 0; i < l.nrows(); i++)
        {
            for (size_t j = 0; j < l.ncols(); j++)
            {
                l(i,j) = ++count;
            }
        }
        e = k * l;

        CHECK( e.size() == 9 );
        CHECK( e(0,0) == 70 );
        CHECK( e(0,1) == 80 );
        CHECK( e(0,2) == 90 );
        CHECK( e(1,0) == 158 );
        CHECK( e(1,1) == 184 );
        CHECK( e(1,2) == 210 );
        CHECK( e(2,0) == 246 );
        CHECK( e(2,1) == 288 );
        CHECK( e(2,2) == 330 );

        //k.print("k");
        //l.print("l");
        //e.print("k*l");

        size_t a = 10, b = 13, c = 15;
        evolm::smatrix<double> m(10,13);
        evolm::smatrix<double> n(13,15);
        evolm::smatrix<double> res;

        m(2,0) = 1;
        m(6,1) = 2;
        m(9,1) = 8;
        m(9,3) = 9;
        m(1,5) = 3;
        m(9,6) = 10;
        m(7,7) = 4;
        m(7,10) = 5;
        m(5,12) = 6;

        n(7,0) = 4;
        n(4,2) = 1;
        n(4,5) = 2;
        n(4,6) = 3;
        n(4,8) = 4;
        n(9,11) = 5;
        n(9,13) = 6;
        n(12,14) = 7;

        res = m * n;

        //m.print("m");
        //n.print("n");
        //res.print("m * n");

        CHECK( res.size() == 2 );
        CHECK( res.max_key() == a*c );

        CHECK( res(7,0) == 16 );
        CHECK( res(5,14) == 42 );

        // check time
        size_t s = 1000;
        evolm::smatrix<double> M(s,s);        
        evolm::smatrix<double> result;

std::cout<<"assigning M ..."<<"\n";
auto start = std::chrono::high_resolution_clock::now();
        for (size_t i = 0; i < s*s-1; i++)
        {
            //if ( i%100 == 0 )
                M[i] = i+1;
        }
std::cout<<"Done"<<"\n";
auto stop = std::chrono::high_resolution_clock::now();
auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
std::cout <<"duration (milliseconds): "<< duration.count() <<" for n elements: "<< s*s << std::endl;

std::cout<<"assigning N ..."<<"\n";
start = std::chrono::high_resolution_clock::now();
        evolm::smatrix<double> N(M);
std::cout<<"Done."<<"\n";
stop = std::chrono::high_resolution_clock::now();
duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
std::cout <<"duration (milliseconds): "<< duration.count() <<" for n elements: "<< s*s << std::endl;
std::cout<<"Multiplying ..."<<"\n";
        start = std::chrono::high_resolution_clock::now();

        result = M * N;

std::cout<<"Done."<<"\n";
        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout <<"Duration (milliseconds): "<< duration.count() <<" for n elements: "<< s*s << std::endl;
    }

    SECTION("operator+")
    {
        size_t n = 1000;
        evolm::smatrix<double> k(3,4);
        evolm::smatrix<double> l(3,4);
        evolm::smatrix<double> e;

        size_t count = 0;
        for (size_t i = 0; i < k.nrows(); i++)
        {
            for (size_t j = 0; j < k.ncols(); j++)
            {
                k(i,j) = ++count;
            }
        }
        count = 0;
        for (size_t i = 0; i < l.nrows(); i++)
        {
            for (size_t j = 0; j < l.ncols(); j++)
            {
                l(i,j) = ++count * 2;
            }
        }
        e = k + l;
        //e.print("k+l");

        CHECK( e.size() == 12 );
        CHECK( e(0,0) == 3 );
        CHECK( e(0,1) == 6 );
        CHECK( e(0,2) == 9 );
        CHECK( e(0,3) == 12 );
        CHECK( e(1,0) == 15 );
        CHECK( e(1,1) == 18 );
        CHECK( e(1,2) == 21 );
        CHECK( e(1,3) == 24 );
        CHECK( e(2,0) == 27 );
        CHECK( e(2,1) == 30 );
        CHECK( e(2,2) == 33 );
        CHECK( e(2,3) == 36 );

        e.resize();
        l.resize(3,4);

        count = 0;
        for (size_t i = 0; i < l.nrows()*l.ncols()-1; i++)
        {
            if ( i%2 == 0 )
                l[i] = ++count;
        }
        //k.print("k");
        //l.print("new l");
        
        e = k + l;
        //e.print("k+l");

        CHECK( e.size() == 12 );
        CHECK( e(0,0) == 2 );
        CHECK( e(0,1) == 2 );
        CHECK( e(0,2) == 5 );
        CHECK( e(0,3) == 4 );
        CHECK( e(1,0) == 8 );
        CHECK( e(1,1) == 6 );
        CHECK( e(1,2) == 11 );
        CHECK( e(1,3) == 8 );
        CHECK( e(2,0) == 14 );
        CHECK( e(2,1) == 10 );
        CHECK( e(2,2) == 17 );
        CHECK( e(2,3) == 12 );
    }

    SECTION("operator-")
    {
        size_t n = 1000;
        evolm::smatrix<double> k(3,4);
        evolm::smatrix<double> l(3,4);
        evolm::smatrix<double> e;

        size_t count = 0;
        for (size_t i = 0; i < k.nrows(); i++)
        {
            for (size_t j = 0; j < k.ncols(); j++)
            {
                k(i,j) = ++count;
            }
        }
        count = 0;
        for (size_t i = 0; i < l.nrows(); i++)
        {
            for (size_t j = 0; j < l.ncols(); j++)
            {
                l(i,j) = ++count * 2;
            }
        }
        
        //k.print("k");
        //l.print("l");

        e = k - l;

        //e.print("k-l");

        CHECK( e.size() == 12 );
        CHECK( e(0,0) == -1 );
        CHECK( e(0,1) == -2 );
        CHECK( e(0,2) == -3 );
        CHECK( e(0,3) == -4 );
        CHECK( e(1,0) == -5 );
        CHECK( e(1,1) == -6 );
        CHECK( e(1,2) == -7 );
        CHECK( e(1,3) == -8 );
        CHECK( e(2,0) == -9 );
        CHECK( e(2,1) == -10 );
        CHECK( e(2,2) == -11 );
        CHECK( e(2,3) == -12 );

        e.resize();
        l.resize(3,4);

        count = 0;
        for (size_t i = 0; i < l.nrows()*l.ncols()-1; i++)
        {
            if ( i%2 == 0 )
                l[i] = ++count;
        }

        //k.print("k");
        //l.print("new l");
        
        e = k - l;
        
        //e.print("k-l");

        CHECK( e.size() == 11 );
        CHECK( e(0,1) == 2 );
        CHECK( e(0,2) == 1 );
        CHECK( e(0,3) == 4 );
        CHECK( e(1,0) == 2 );
        CHECK( e(1,1) == 6 );
        CHECK( e(1,2) == 3 );
        CHECK( e(1,3) == 8 );
        CHECK( e(2,0) == 4 );
        CHECK( e(2,1) == 10 );
        CHECK( e(2,2) == 5 );
        CHECK( e(2,3) == 12 );
    }

}