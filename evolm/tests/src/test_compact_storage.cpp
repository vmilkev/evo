#include "catch_amalgamated.hpp"
#include "compact_storage.hpp"
#include <chrono>

TEST_CASE("Compact storage, checking class constructors, type = double")
{
    SECTION("1. Default empty constructor")
    {
        evolm::compact_storage<double> M;

        CHECK(M.empty() == true);
        CHECK(M.ncols() == 0);
        CHECK(M.nrows() == 0);
    }

    SECTION("2. Symmetric empty constructor")
    {
        evolm::compact_storage<double> M(5);

        CHECK(M.empty() == true);
        CHECK(M.ncols() == 5);
        CHECK(M.nrows() == 5);
        CHECK(M.max_key()+1 == (5 * 5 + 5) / 2);
    }

    SECTION("3. Rectangullar empty constructor")
    {
        evolm::compact_storage<double> M(5,5);

        CHECK(M.empty() == true);
        CHECK(M.ncols() == 5);
        CHECK(M.nrows() == 5);
        CHECK(M.max_key()+1 == 5 * 5);
    }

    SECTION("4. Default empty constructor with resizimg")
    {
        evolm::compact_storage<double> M;

        M.resize(5);

        CHECK(M.empty() == true);
        CHECK(M.ncols() == 5);
        CHECK(M.nrows() == 5);
        CHECK(M.max_key()+1 == (5 * 5 + 5) / 2);

        M.resize(5, 5);

        CHECK(M.empty() == true);
        CHECK(M.ncols() == 5);
        CHECK(M.nrows() == 5);
        CHECK(M.max_key()+1 == 5 * 5);
    }

    SECTION("5. compact_storage(std::vector<T> &, std::vector<size_t> &, std::vector<size_t> &, size_t, size_t)")
    {
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;

        size_t dim = 10;
        size_t records = (dim*dim);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                rows.push_back(i);
                cols.push_back(j);
                values.push_back(i);
            }
        }

        evolm::compact_storage<double> M(values, rows, cols, dim, dim);

        CHECK(M.empty() == false);
        CHECK(M.ncols() == dim);
        CHECK(M.nrows() == dim);
        CHECK(M.max_key()+1 == records);

        evolm::matrix<double> dense_M;

        M.to_dense(dense_M);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                CHECK( dense_M(i,j) == i );
            }
        }

        evolm::smatrix<double> sparse_M;

        M.to_sparse(sparse_M);

        CHECK( sparse_M.size() == 90 );

        for (size_t i = 0; i < dim; i++)
        {
            if ( i != 0 )
            {
                for (size_t j = 0; j < dim; j++)
                {
                    if ( sparse_M.nonzero(i,j) )
                        CHECK( sparse_M(i,j) == i );
                }
            }
        }
    }

    SECTION("6. compact_storage(std::vector<T> &, std::vector<size_t> &, std::vector<size_t> &, size_t)")
    {
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;

        size_t dim = 10;
        size_t records = (dim*dim+dim)/2;

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                rows.push_back(i);
                cols.push_back(j);
                values.push_back(i);
            }
        }

        evolm::compact_storage<double> M(values, rows, cols, dim);

        CHECK(M.empty() == false);
        CHECK(M.ncols() == dim);
        CHECK(M.nrows() == dim);
        CHECK(M.max_key()+1 == records);

        evolm::matrix<double> dense_M;

        M.to_dense(dense_M);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                CHECK( dense_M(i,j) == i );
            }
        }

        evolm::smatrix<double> sparse_M;

        M.to_sparse(sparse_M);

        CHECK( sparse_M.size() == records-1 );

        for (size_t i = 0; i < dim; i++)
        {
            if ( i != 0 )
            {
                for (size_t j = 0; j <= i; j++)
                {
                    if ( sparse_M.nonzero(i,j) )
                        CHECK( sparse_M(i,j) == i );
                }
            }
        }
    }

    SECTION("7. compact_storage(std::vector<T> &, size_t, size_t)")
    {
        std::vector<double> values;

        size_t dim = 10;
        size_t records = (dim*dim);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                values.push_back(i);
            }
        }

        evolm::compact_storage<double> M(values, dim, dim);

        CHECK(M.empty() == false);
        CHECK(M.ncols() == dim);
        CHECK(M.nrows() == dim);
        CHECK(M.max_key()+1 == records);

        evolm::matrix<double> dense_M;

        M.to_dense(dense_M);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                CHECK( dense_M(i,j) == i );
            }
        }

        evolm::smatrix<double> sparse_M;

        M.to_sparse(sparse_M);

        CHECK( sparse_M.size() == 90 );

        for (size_t i = 0; i < dim; i++)
        {
            if ( i != 0 )
            {
                for (size_t j = 0; j < dim; j++)
                {
                    if ( sparse_M.nonzero(i,j) )
                        CHECK( sparse_M(i,j) == i );
                }
            }
        }
    }

    SECTION("8. compact_storage(std::vector<T> &, size_t)")
    {
        std::vector<double> values;

        size_t dim = 10;
        size_t records = (dim*dim+dim)/2;

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                values.push_back(i);
            }
        }

        evolm::compact_storage<double> M(values, dim);

        CHECK(M.empty() == false);
        CHECK(M.ncols() == dim);
        CHECK(M.nrows() == dim);
        CHECK(M.max_key()+1 == records);

        evolm::matrix<double> dense_M;

        M.to_dense(dense_M);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                CHECK( dense_M(i,j) == i );
            }
        }

        evolm::smatrix<double> sparse_M;

        M.to_sparse(sparse_M);

        CHECK( sparse_M.size() == records-1 );

        for (size_t i = 0; i < dim; i++)
        {
            if ( i != 0 )
            {
                for (size_t j = 0; j <= i; j++)
                {
                    CHECK( sparse_M(i,j) == i );
                }
            }
        }
    }

    SECTION("9. append(T, size_t, size_t)")
    {
        size_t dim = 10;
        size_t records = (dim*dim+dim)/2;

        evolm::compact_storage<double> M(dim);

        std::vector<double> values;

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                M.append(i,i,j);
            }
        }

        M.optimize();

        CHECK(M.empty() == false);
        CHECK(M.ncols() == dim);
        CHECK(M.nrows() == dim);
        CHECK(M.max_key()+1 == records);

        evolm::matrix<double> dense_M;

        M.to_dense(dense_M);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                CHECK( dense_M(i,j) == i );
            }
        }

        evolm::smatrix<double> sparse_M;

        M.to_sparse(sparse_M);

        CHECK( sparse_M.size() == records-1 );

        for (size_t i = 0; i < dim; i++)
        {
            if ( i != 0 )
            {
                for (size_t j = 0; j <= i; j++)
                {
                    CHECK( sparse_M(i,j) == i );
                }
            }
        }
    }

    SECTION("10. append(std::vector<T> &, std::vector<size_t> &, std::vector<size_t> &)")
    {
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;

        size_t dim = 10;
        size_t records = (dim*dim);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                rows.push_back(i);
                cols.push_back(j);
                values.push_back(i);
            }
        }

        evolm::compact_storage<double> M(dim, dim);

        M.append(values, rows, cols);

        CHECK(M.empty() == false);
        CHECK(M.ncols() == dim);
        CHECK(M.nrows() == dim);
        CHECK(M.max_key()+1 == records);

        evolm::matrix<double> dense_M;

        M.to_dense(dense_M);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                CHECK( dense_M(i,j) == i );
            }
        }

        evolm::smatrix<double> sparse_M;

        M.to_sparse(sparse_M);

        CHECK( sparse_M.size() == 90 );

        for (size_t i = 0; i < dim; i++)
        {
            if ( i != 0 )
            {
                for (size_t j = 0; j < dim; j++)
                {
                    CHECK( sparse_M(i,j) == i );
                }
            }
        }
    }

    SECTION("11. append(std::vector<T> &)")
    {
        std::vector<double> values;

        size_t dim = 10;
        size_t records = (dim*dim+dim)/2;

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                values.push_back(i);
            }
        }

        evolm::compact_storage<double> M(dim);
        
        M.append(values);

        CHECK(M.empty() == false);
        CHECK(M.ncols() == dim);
        CHECK(M.nrows() == dim);
        CHECK(M.max_key()+1 == records);

        evolm::matrix<double> dense_M;

        M.to_dense(dense_M);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                CHECK( dense_M(i,j) == i );
            }
        }

        evolm::smatrix<double> sparse_M;

        M.to_sparse(sparse_M);

        CHECK( sparse_M.size() == records-1 );

        for (size_t i = 0; i < dim; i++)
        {
            if ( i != 0 )
            {
                for (size_t j = 0; j <= i; j++)
                {
                    if ( sparse_M.nonzero(i,j) )
                        CHECK( sparse_M(i,j) == i );
                }
            }
        }
    }

    SECTION("12. &operator=(const smatrix &rhs)")
    {
        evolm::smatrix<double> M1(5, 4);

        M1(0, 0) = 1.0;
        M1(2, 0) = 2.0;
        M1(2, 2) = 3.0;
        M1(4, 1) = 4;
        M1(4, 2) = 5.0;

        evolm::smatrix<double> N1 = M1;

        evolm::smatrix<double> L1;

        CHECK(N1.size() == 5);
        CHECK(M1.size() == 5);

        CHECK(N1(0, 0) == 1.0);
        CHECK(N1(2, 0) == 2.0);
        CHECK(N1(2, 2) == 3.0);
        CHECK(N1(4, 1) == 4.0);
        CHECK(N1(4, 2) == 5.0);

        evolm::smatrix<double> K1;

        CHECK(K1.size() == 0);
        K1 = M1;

        CHECK(K1.size() == 5);
        CHECK(M1.size() == 5);

        CHECK(K1(0, 0) == 1.0);
        CHECK(K1(2, 0) == 2.0);
        CHECK(K1(2, 2) == 3.0);
        CHECK(K1(4, 1) == 4.0);
        CHECK(K1(4, 2) == 5.0);

        //------------------------------
        evolm::compact_storage<double> M(5, 4);

        M.append(1.0, 0, 0);
        M.append(2.0, 2, 0);
        M.append(3.0, 2, 2);
        M.append(4.0, 4, 1);
        M.append(5.0, 4, 2);

        evolm::compact_storage<double> N = M;

        evolm::compact_storage<double> L;

        CHECK(N.size() == 5);
        CHECK(M.size() == 5);

        N.to_sparse(N1);

        CHECK(N1.size() == 5);

        CHECK(N1(0, 0) == 1.0);
        CHECK(N1(2, 0) == 2.0);
        CHECK(N1(2, 2) == 3.0);
        CHECK(N1(4, 1) == 4.0);
        CHECK(N1(4, 2) == 5.0);

        evolm::compact_storage<double> K;

        CHECK(K.size() == 0);

        K = M;

        CHECK(M.size() == 5);
        CHECK(K.size() == 5);

        K.to_sparse(K1);

        CHECK(K1.size() == 5);

        CHECK(K1(0, 0) == 1.0);
        CHECK(K1(2, 0) == 2.0);
        CHECK(K1(2, 2) == 3.0);
        CHECK(K1(4, 1) == 4.0);
        CHECK(K1(4, 2) == 5.0);
    }

    SECTION("13. compact_storage(const compact_storage &obj)")
    {
        evolm::smatrix<double> M1(5, 4);

        M1(0, 0) = 1.0;
        M1(2, 0) = 2.0;
        M1(2, 2) = 3.0;
        M1(4, 1) = 4;
        M1(4, 2) = 5.0;

        CHECK(M1.size() == 5);

        //------------------------------
        evolm::compact_storage<double> M(5, 4);

        M.append(1.0, 0, 0);
        M.append(2.0, 2, 0);
        M.append(3.0, 2, 2);
        M.append(4.0, 4, 1);
        M.append(5.0, 4, 2);

        evolm::compact_storage<double> N(M);

        evolm::smatrix<double> N1;

        N.to_sparse(N1);

        CHECK(N1(0, 0) == 1.0);
        CHECK(N1(2, 0) == 2.0);
        CHECK(N1(2, 2) == 3.0);
        CHECK(N1(4, 1) == 4.0);
        CHECK(N1(4, 2) == 5.0);

        CHECK(N1.size() == 5);
        N1.clear();
        CHECK(N1.size() == 0);
    }

    SECTION("14. fwrite() / fread()")
    {        
        evolm::compact_storage<double> M(5, 4);

        M.append(1.0, 0, 0);
        M.append(2.0, 2, 0);
        M.append(3.0, 2, 2);
        M.append(4.0, 4, 1);
        M.append(5.0, 4, 2);

        CHECK(M.size() == 5);

        M.fwrite();

        CHECK(M.size() == 0);

        M.fread();

        CHECK(M.size() == 5);

        evolm::compact_storage<double> N(M);

        CHECK(N.size() == 5);

        evolm::smatrix<double> N1;

        N.to_sparse(N1);

        N.clear();
        M.fclear();

        CHECK(N.size() == 0);

        CHECK(N1.size() == 5);

        CHECK(N1(0, 0) == 1.0);
        CHECK(N1(2, 0) == 2.0);
        CHECK(N1(2, 2) == 3.0);
        CHECK(N1(4, 1) == 4.0);
        CHECK(N1(4, 2) == 5.0);
    }

    SECTION("15. Copy by fwrite() / fread()")
    {        
        evolm::compact_storage<double> M(5, 4);

        M.append(1.0, 0, 0);
        M.append(2.0, 2, 0);
        M.append(3.0, 2, 2);
        M.append(4.0, 4, 1);
        M.append(5.0, 4, 2);

        CHECK(M.size() == 5);

        M.fwrite();

        CHECK(M.size() == 0);

        evolm::compact_storage<double> N(M);

        N.fread();

        CHECK(N.size() == 5);

        evolm::smatrix<double> N1;

        N.to_sparse(N1);

        N.clear();
        M.fclear();

        CHECK(N.size() == 0);

        CHECK(N1.size() == 5);

        CHECK(N1(0, 0) == 1.0);
        CHECK(N1(2, 0) == 2.0);
        CHECK(N1(2, 2) == 3.0);
        CHECK(N1(4, 1) == 4.0);
        CHECK(N1(4, 2) == 5.0);
    }

    SECTION("16. Operator =, by fwrite() / fread()")
    {        
        evolm::compact_storage<double> M(5, 4);

        M.append(1.0, 0, 0);
        M.append(2.0, 2, 0);
        M.append(3.0, 2, 2);
        M.append(4.0, 4, 1);
        M.append(5.0, 4, 2);

        CHECK(M.size() == 5);

        M.fwrite();

        CHECK(M.size() == 0);

        evolm::compact_storage<double> N = M;

        N.fread();

        CHECK(N.size() == 5);

        evolm::smatrix<double> N1;

        N.to_sparse(N1);

        N.clear();
        M.fclear();

        CHECK(N.size() == 0);

        CHECK(N1.size() == 5);

        CHECK(N1(0, 0) == 1.0);
        CHECK(N1(2, 0) == 2.0);
        CHECK(N1(2, 2) == 3.0);
        CHECK(N1(4, 1) == 4.0);
        CHECK(N1(4, 2) == 5.0);
    }

    SECTION("17. fwrite(const std::string &) / fread(const std::string &)")
    {
        evolm::compact_storage<double> M(5, 4);

        M.append(1.0, 0, 0);
        M.append(2.0, 2, 0);
        M.append(3.0, 2, 2);
        M.append(4.0, 4, 1);
        M.append(5.0, 4, 2);

        CHECK(M.size() == 5);

        M.fwrite("M_storage.bin");

        M.clear();

        CHECK(M.size() == 0);

        M.fread("M_storage.bin");

        CHECK(M.size() == 5);

        evolm::compact_storage<double> N(M);

        CHECK(N.size() == 5);

        evolm::smatrix<double> N1;

        N.to_sparse(N1);

        N.clear();
        M.fclear("M_storage.bin");

        CHECK(N.size() == 0);

        CHECK(N1.size() == 5);

        CHECK(N1(0, 0) == 1.0);
        CHECK(N1(2, 0) == 2.0);
        CHECK(N1(2, 2) == 3.0);
        CHECK(N1(4, 1) == 4.0);
        CHECK(N1(4, 2) == 5.0);
    }

    SECTION("18. optimize(): dense")
    {
        // optimize dense storage

        size_t dim = 100;

        std::vector<double> values;

        for (size_t i = 0; i < dim; i++)
        {
            if ( i%4 )
                values.push_back(1.0);
            else
                values.push_back(0.0);
        }

        evolm::compact_storage<double> s(values, dim, 1);
        
        //s.optimize();

        CHECK( s.is_sparse() == false );

        evolm::matrix<double> dense_s;

        s.to_dense(dense_s);

        for (size_t i = 0; i < dim; i++)
        {
            CHECK( dense_s(i,0) == values[i] );
        }

        std::vector<double> values2;
        std::vector<size_t> keys2;

        s.to_vector(values2,keys2);

        for (size_t i = 0; i < dim; i++)
        {
            CHECK(values2[i] == values[i]);
        }

        values2.clear();
        keys2.clear();

        // -----------------------------------

        s.set_sparsity_threshold(0.2);

        s.optimize();

        CHECK( s.is_sparse() == true );

        evolm::smatrix<double> sparse_s;

        s.to_sparse(sparse_s);

        for (size_t i = 0; i < dim; i++)
        {
            if ( i%4 )
                CHECK( sparse_s(i,0) == values[i] );
        }

        s.to_vector(values2, keys2);

        CHECK(values2.size() == sparse_s.size());

        size_t index = 0;
        for (size_t i = 0; i < dim; i++)
        {
            if ( i%4 )
            {
                CHECK( sparse_s(i,0) == values2[index] );
                index++;
            }
        }
    }

    SECTION("19. optimize(): sparse")
    {
        // optimize dense storage

        size_t dim = 100;

        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;

        for (size_t i = 0; i < dim; i++)
        {
            if ( i%4 )
            {
                values.push_back(1.0);
                rows.push_back(i);
                cols.push_back(0);
            }
        }

        evolm::compact_storage<double> s(values, rows, cols, dim, 1);
        
        //s.optimize();

        CHECK( s.is_sparse() == false );

        evolm::matrix<double> dense_s;

        s.to_dense(dense_s);

        for (size_t i = 0; i < dim; i++)
        {
            if ( i%4 )
                CHECK( dense_s(i,0) == values[i] );
            else
                CHECK( dense_s(i,0) == 0.0 );
        }

        // -----------------------------------

        s.set_sparsity_threshold(0.2);

        s.optimize();

        CHECK( s.is_sparse() == true );

        evolm::smatrix<double> sparse_s;

        s.to_sparse(sparse_s);

        for (size_t i = 0; i < dim; i++)
        {
            if ( i%4 )
                CHECK( sparse_s(i,0) == values[i] );
        }
    }

    SECTION("20. optimize(): sparse with some zeros")
    {
        // optimize dense storage

        size_t dim = 100;

        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;

        values.push_back(1.0);
        rows.push_back(5);
        cols.push_back(0);

        values.push_back(1.0);
        rows.push_back(15);
        cols.push_back(0);

        values.push_back(1.0);
        rows.push_back(35);
        cols.push_back(0);

        for (size_t i = 0; i < dim; i++)
        {
            if ( i%2 )
            {
                values.push_back(0.0);
                rows.push_back(i);
                cols.push_back(0);
            }
        }

        evolm::compact_storage<double> s(values, rows, cols, dim, 1);

        CHECK( s.size() == 3 );

        CHECK( s.is_sparse() == true );

        evolm::matrix<double> dense_s;

        s.to_dense(dense_s);

        for (size_t i = 0; i < dim; i++)
        {
            if ( i == 5 || i == 15 || i == 35 )
                CHECK( dense_s(i,0) == 1.0 );
            else
                CHECK( dense_s(i,0) == 0.0 );
        }

        // -----------------------------------

        CHECK( s.is_sparse() == true );

        evolm::smatrix<double> sparse_s;

        s.to_sparse(sparse_s);

        CHECK( sparse_s.size() == 3 );

        CHECK( sparse_s(5,0) == 1 );
        CHECK( sparse_s(15,0) == 1 );
        CHECK( sparse_s(35,0) == 1 );

        std::vector<double> values2;
        std::vector<size_t> keys2;

        s.to_vector(values2, keys2);

        CHECK( values2.size() == 3 );

        CHECK( sparse_s(5,0) == values2[0] );
        CHECK( sparse_s(15,0) == values2[1] );
        CHECK( sparse_s(35,0) == values2[2] );
    }

    SECTION("21. to_sparse(smatrix<T> &, size_t *, size_t *): symmetric and dense")
    {
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;
        size_t dim = 10;
        size_t counts = 0;
        size_t counts2 = 0;
        
        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                values.push_back(i);
                rows.push_back(i);
                cols.push_back(j);
            }
        }

        evolm::compact_storage<double> s(values, rows, cols, dim);

        size_t row_range[2];
        size_t col_range[2];

        row_range[0] = 5;
        row_range[1] = 8;
        col_range[0] = 5;
        col_range[1] = 8;

        evolm::smatrix<double> sparse_s;
        evolm::smatrix<double> sparse_s2;

        s.to_sparse(sparse_s, row_range, col_range);

        s.to_sparse(sparse_s2);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                if ( i >= row_range[0] && i <= row_range[1] && j >= col_range[0] && j <= col_range[1] )
                {
                    if ( sparse_s.nonzero(i,j) )
                    {
                        CHECK( sparse_s[counts] == values[ counts2 ] );
                        counts++;
                    }
                }
                counts2++;
            }
        }
    }

    SECTION("22. to_sparse(smatrix<T> &, size_t *, size_t *): rectangular and dense")
    {
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;
        size_t dim = 10;
        size_t counts = 0;
        size_t counts2 = 0;
        
        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                values.push_back(i);
                rows.push_back(i);
                cols.push_back(j);
            }
        }

        evolm::compact_storage<double> s(values, rows, cols, dim, dim);

        size_t row_range[2];
        size_t col_range[2];

        row_range[0] = 5;
        row_range[1] = 8;
        col_range[0] = 5;
        col_range[1] = 8;

        evolm::smatrix<double> sparse_s;
        evolm::smatrix<double> sparse_s2;

        s.to_sparse(sparse_s, row_range, col_range);

        s.to_sparse(sparse_s2);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                if ( i >= row_range[0] && i <= row_range[1] && j >= col_range[0] && j <= col_range[1] )
                {
                    if ( sparse_s.nonzero(i,j) )
                    {
                        CHECK( sparse_s[counts] == values[ counts2 ] );
                        counts++;
                    }
                }
                counts2++;
            }
        }
    }

    SECTION("23. to_sparse(smatrix<T> &, size_t *, size_t *): symmetric and sparse")
    {
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;
        size_t dim = 10;
        size_t counts = 0;
        size_t counts2 = 0;
        
        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                if ( counts % 2 )
                    values.push_back(counts);
                else
                    values.push_back(0);
                rows.push_back(i);
                cols.push_back(j);
                
                counts++;
            }
        }

        counts = 0;

        evolm::compact_storage<double> s(values, rows, cols, dim);

        s.set_sparsity_threshold(0.1);

        CHECK(s.is_sparse() == false);

        s.optimize();

        CHECK(s.is_sparse() == true);

        size_t row_range[2];
        size_t col_range[2];

        row_range[0] = 5;
        row_range[1] = 8;
        col_range[0] = 5;
        col_range[1] = 8;

        evolm::smatrix<double> sparse_s;
        evolm::smatrix<double> sparse_s2;

        s.to_sparse(sparse_s, row_range, col_range);

        s.to_sparse(sparse_s2);


        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                if ( i >= row_range[0] && i <= row_range[1] && j >= col_range[0] && j <= col_range[1] )
                {
                    if ( sparse_s.nonzero(i,j) )
                    {
                        CHECK( sparse_s[counts] == values[ counts2 ] );
                        counts++;
                    }
                }
                counts2++;
            }
        }
    }

    SECTION("24. to_sparse(smatrix<T> &, size_t *, size_t *): rectangular and sparse")
    {
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;
        size_t dim = 10;
        size_t counts = 0;
        size_t counts2 = 0;
        
        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                if ( counts % 2 )
                    values.push_back(counts);
                else
                    values.push_back(0);
                rows.push_back(i);
                cols.push_back(j);
                
                counts++;
            }
        }

        counts = 0;

        evolm::compact_storage<double> s(values, rows, cols, dim, dim);

        s.set_sparsity_threshold(0.1);

        CHECK(s.is_sparse() == false);

        s.optimize();

        CHECK(s.is_sparse() == true);

        size_t row_range[2];
        size_t col_range[2];

        row_range[0] = 5;
        row_range[1] = 8;
        col_range[0] = 5;
        col_range[1] = 8;

        evolm::smatrix<double> sparse_s;
        evolm::smatrix<double> sparse_s2;

        s.to_sparse(sparse_s, row_range, col_range);

        s.to_sparse(sparse_s2);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                if ( i >= row_range[0] && i <= row_range[1] && j >= col_range[0] && j <= col_range[1] )
                {
                    if ( sparse_s.nonzero(i,j) )
                    {
                        CHECK( sparse_s[counts] == values[ counts2 ] );
                        counts++;
                    }
                }
                counts2++;
            }
        }
    }

    SECTION("25. to_dense(matrix<T> &, size_t *, size_t *): symmetric and dense")
    {
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;
        size_t dim = 10;
        size_t counts = 0;
        size_t counts2 = 0;
        
        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                values.push_back(i);
                rows.push_back(i);
                cols.push_back(j);
            }
        }

        evolm::compact_storage<double> s(values, rows, cols, dim);

        size_t row_range[2];
        size_t col_range[2];

        row_range[0] = 5;
        row_range[1] = 8;
        col_range[0] = 5;
        col_range[1] = 8;

        evolm::matrix<double> dense_s;
        evolm::matrix<double> dense_s2;

        s.to_dense(dense_s, row_range, col_range);

        s.to_dense(dense_s2);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                if ( i >= row_range[0] && i <= row_range[1] && j >= col_range[0] && j <= col_range[1] )
                {
                    CHECK( dense_s[counts] == values[ counts2 ] );
                    counts++;
                }
                counts2++;
            }
        }
    }

    SECTION("26. to_dense(matrix<T> &, size_t *, size_t *): rectangular and dense")
    {
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;
        size_t dim = 10;
        size_t counts = 0;
        size_t counts2 = 0;
        
        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                values.push_back(i);
                rows.push_back(i);
                cols.push_back(j);
            }
        }

        evolm::compact_storage<double> s(values, rows, cols, dim, dim);

        size_t row_range[2];
        size_t col_range[2];

        row_range[0] = 5;
        row_range[1] = 8;
        col_range[0] = 5;
        col_range[1] = 8;

        evolm::matrix<double> dense_s;
        evolm::matrix<double> dense_s2;

        s.to_dense(dense_s, row_range, col_range);

        s.to_dense(dense_s2);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                if ( i >= row_range[0] && i <= row_range[1] && j >= col_range[0] && j <= col_range[1] )
                {
                    CHECK( dense_s[counts] == values[ counts2 ] );
                    counts++;
                }
                counts2++;
            }
        }
    }

    SECTION("27. to_dense(matrix<T> &, size_t *, size_t *): symmetric and sparse")
    {
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;
        size_t dim = 10;
        size_t counts = 0;
        size_t counts2 = 0;
        
        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                if ( counts % 2 )
                    values.push_back(i);
                else
                    values.push_back(0);
                rows.push_back(i);
                cols.push_back(j);

                counts++;
            }
        }
        
        counts = 0;

        evolm::compact_storage<double> s(values, rows, cols, dim);

        s.set_sparsity_threshold(0.1);

        CHECK(s.is_sparse() == false);

        s.optimize();

        CHECK(s.is_sparse() == true);

        size_t row_range[2];
        size_t col_range[2];

        row_range[0] = 5;
        row_range[1] = 8;
        col_range[0] = 5;
        col_range[1] = 8;

        evolm::matrix<double> dense_s;
        evolm::matrix<double> dense_s2;

        s.to_dense(dense_s, row_range, col_range);
        s.to_dense(dense_s2);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                if ( i >= row_range[0] && i <= row_range[1] && j >= col_range[0] && j <= col_range[1] )
                {
                    CHECK( dense_s[counts] == values[ counts2 ] );
                    counts++;
                }
                counts2++;
            }
        }
    }

    SECTION("28. to_dense(matrix<T> &, size_t *, size_t *): rectangular and sparse")
    {
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;
        size_t dim = 10;
        size_t counts = 0;
        size_t counts2 = 0;
        
        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                if ( counts % 2 )
                    values.push_back(i);
                else
                    values.push_back(0);
                rows.push_back(i);
                cols.push_back(j);

                counts++;
            }
            
        }
        
        counts = 0;

        evolm::compact_storage<double> s(values, rows, cols, dim, dim);

        s.set_sparsity_threshold(0.1);

        CHECK(s.is_sparse() == false);

        s.optimize();

        CHECK(s.is_sparse() == true);

        size_t row_range[2];
        size_t col_range[2];

        row_range[0] = 5;
        row_range[1] = 8;
        col_range[0] = 5;
        col_range[1] = 8;

        evolm::matrix<double> dense_s;
        evolm::matrix<double> dense_s2;

        s.to_dense(dense_s, row_range, col_range);
        s.to_dense(dense_s2);

        //dense_s.print("dense S");
        //dense_s2.print("dense S2");

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                if ( i >= row_range[0] && i <= row_range[1] && j >= col_range[0] && j <= col_range[1] )
                {
                    CHECK( dense_s[counts] == values[ counts2 ] );
                    counts++;
                }
                counts2++;
            }
        }
    }

    SECTION("29. to_vector(std::vector<T> &, std::vector<size_t> &, size_t *, size_t *): symmetric and dense")
    {
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;
        size_t dim = 10;
        size_t counts = 0;
        size_t counts2 = 0;
        
        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                values.push_back(i);
                rows.push_back(i);
                cols.push_back(j);
            }
        }

        evolm::compact_storage<double> s(values, rows, cols, dim);

        size_t row_range[2];
        size_t col_range[2];

        row_range[0] = 5;
        row_range[1] = 8;
        col_range[0] = 5;
        col_range[1] = 8;

        std::vector<double> vect_s;
        std::vector<size_t> vect_s2;

        s.to_vector(vect_s, vect_s2, row_range, col_range);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                if ( i >= row_range[0] && i <= row_range[1] && j >= col_range[0] && j <= col_range[1] )
                {
                    CHECK( vect_s[counts] == values[ counts2 ] );
                    CHECK( vect_s2[counts] == counts );
                    counts++;
                }
                counts2++;
            }
        }
    }

    SECTION("30. to_vector(std::vector<T> &, std::vector<size_t> &, size_t *, size_t *): rectangular and dense")
    {
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;
        size_t dim = 10;
        size_t counts = 0;
        size_t counts2 = 0;
        
        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                values.push_back(i);
                rows.push_back(i);
                cols.push_back(j);
            }
        }

        evolm::compact_storage<double> s(values, rows, cols, dim, dim);

        size_t row_range[2];
        size_t col_range[2];

        row_range[0] = 5;
        row_range[1] = 8;
        col_range[0] = 5;
        col_range[1] = 8;

        std::vector<double> vect_s;
        std::vector<size_t> vect_s2;

        s.to_vector(vect_s, vect_s2, row_range, col_range);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                if ( i >= row_range[0] && i <= row_range[1] && j >= col_range[0] && j <= col_range[1] )
                {
                    CHECK( vect_s[counts] == values[ counts2 ] );
                    CHECK( vect_s2[counts] == counts );
                    counts++;
                }
                counts2++;
            }
        }
    }

    SECTION("31. to_vector(std::vector<T> &, std::vector<size_t> &, size_t *, size_t *): symmetric and sparse")
    {
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;
        size_t dim = 10;
        size_t counts = 0;
        size_t counts2 = 0;
        size_t counts3 = 0;
        
        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                if ( counts % 2 )
                    values.push_back(counts);
                else
                    values.push_back(0);
                
                rows.push_back(i);
                cols.push_back(j);

                counts++;
            }
        }
        
        counts = 0;

        evolm::compact_storage<double> s(values, rows, cols, dim);

        s.set_sparsity_threshold(0.1);

        CHECK(s.is_sparse() == false);

        s.optimize();

        CHECK(s.is_sparse() == true);

        size_t row_range[2];
        size_t col_range[2];

        row_range[0] = 5;
        row_range[1] = 8;
        col_range[0] = 5;
        col_range[1] = 8;

        std::vector<double> vect_s;
        std::vector<size_t> vect_s2;

        s.to_vector(vect_s, vect_s2, row_range, col_range);

        evolm::smatrix<double> sp_s;
        evolm::smatrix<double> sp_s2;

        s.to_sparse(sp_s,row_range, col_range);
        s.to_sparse(sp_s2);
        
        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j <= i; j++)
            {
                if ( i >= row_range[0] && i <= row_range[1] && j >= col_range[0] && j <= col_range[1] )
                {
                    if (values[ counts2 ] != (double)0)
                    {
                        CHECK( vect_s[counts] == values[ counts2 ] );
                        CHECK( vect_s2[counts] == counts3 );
                        counts++;
                    }
                    counts3++;
                }
                counts2++;
            }
        }
    }

    SECTION("32. to_vector(std::vector<T> &, std::vector<size_t> &, size_t *, size_t *): rectangular and sparse")
    {
        std::vector<double> values;
        std::vector<size_t> rows;
        std::vector<size_t> cols;
        size_t dim = 10;
        size_t counts = 0;
        size_t counts2 = 0;
        size_t counts3 = 0;
        
        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                if ( counts % 2 )
                    values.push_back(i);
                else
                    values.push_back(0);
                rows.push_back(i);
                cols.push_back(j);

                counts++;
            }            
        }
        
        counts = 0;

        evolm::compact_storage<double> s(values, rows, cols, dim, dim);

        s.set_sparsity_threshold(0.1);

        CHECK(s.is_sparse() == false);

        s.optimize();

        CHECK(s.is_sparse() == true);

        size_t row_range[2];
        size_t col_range[2];

        row_range[0] = 5;
        row_range[1] = 8;
        col_range[0] = 5;
        col_range[1] = 8;

        std::vector<double> vect_s;
        std::vector<size_t> vect_s2;

        s.to_vector(vect_s, vect_s2, row_range, col_range);

        for (size_t i = 0; i < dim; i++)
        {
            for (size_t j = 0; j < dim; j++)
            {
                if ( i >= row_range[0] && i <= row_range[1] && j >= col_range[0] && j <= col_range[1] )
                {
                    if (values[ counts2 ] != (double)0)
                    {
                        CHECK( vect_s[counts] == values[ counts2 ] );
                        CHECK( vect_s2[counts] == counts3 );
                        counts++;
                    }
                    counts3++;
                }
                counts2++;
            }
        }
    }
}