#include "sparse_matrix.hpp"

namespace evolm
{
    //===============================================================================================================

    template <typename T>
    smatrix<T>::smatrix(size_t row, size_t col)
    {
        /*
           Constructor for regular rectangular matrix.
           It just points how the dense matrix should looks like,
           though the map container is empty, no memory is preallocated.
        */

        std::random_device rd;
        srand(rd());

        int iNum = rand() % 100000;
        binFilename_keys = "smatrix_key_" + std::to_string(iNum);
        binFilename_vals = "smatrix_val_" + std::to_string(iNum);

        rectangular = true;
        compact = false;
        symetric = false;
        ondisk = false;

        numRow = row;
        numCol = col;
        max_elements = row * col;
        ondisc_elements = 0;

        work_load_perthread = 1;
        zerro_tolerance = (T)1.0e-12;

        debug_file = "SMATRIX.log";
    }

    template smatrix<float>::smatrix(size_t row, size_t col);
    template smatrix<double>::smatrix(size_t row, size_t col);
    template smatrix<int>::smatrix(size_t row, size_t col);

    //===============================================================================================================

    template <typename T>
    smatrix<T>::smatrix(size_t lda)
    {
        /*
           Constructor for symmetrical matrix in a compact form
        */

        std::random_device rd;
        srand(rd());

        int iNum = rand() % 100000;
        binFilename_keys = "smatrix_key_" + std::to_string(iNum);
        binFilename_vals = "smatrix_val_" + std::to_string(iNum);

        rectangular = false;
        symetric = true;
        compact = true;
        ondisk = false;

        work_load_perthread = 1;
        zerro_tolerance = (T)1.0e-12;;

        numRow = numCol = lda;
        max_elements = (lda * lda + lda) / 2;
        ondisc_elements = 0;

        debug_file = "SMATRIX.log";
    }

    template smatrix<float>::smatrix(size_t lda);
    template smatrix<double>::smatrix(size_t lda);
    template smatrix<int>::smatrix(size_t lda);
    
    //===============================================================================================================

    template <typename T>
    void smatrix<T>::resize(size_t row, size_t col)
    {
        try
        {
            //fclear();
            A.clear();

            rectangular = true;
            compact = false;
            symetric = false;
            ondisk = false;

            numRow = row;
            numCol = col;
            max_elements = row * col;
            ondisc_elements = 0;

            // rename file name:
            binFilename_keys.clear();
            binFilename_vals.clear();
            std::random_device rd;
            srand(rd());
            int iNum = rand() % 100000;
            binFilename_keys = "smatrix_keys_" + std::to_string(iNum);
            binFilename_vals = "smatrix_vals_" + std::to_string(iNum);
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::resize(size_t, size_t)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::resize(size_t, size_t)" << '\n';
        }        
    }

    template void smatrix<float>::resize(size_t row, size_t col);
    template void smatrix<double>::resize(size_t row, size_t col);
    template void smatrix<int>::resize(size_t row, size_t col);

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::resize(size_t lda)
    {
        try
        {
            //fclear();
            A.clear();

            rectangular = false;
            symetric = true;
            compact = true;
            ondisk = false;

            numRow = numCol = lda;
            max_elements = (lda * lda + lda) / 2;
            ondisc_elements = 0;

            // rename file name:
            binFilename_keys.clear();
            binFilename_vals.clear();
            std::random_device rd;
            srand(rd());
            int iNum = rand() % 100000;
            binFilename_keys = "smatrix_keys_" + std::to_string(iNum);
            binFilename_vals = "smatrix_vals_" + std::to_string(iNum);
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::resize(size_t)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::resize(size_t)" << '\n';
        }        
    }

    template void smatrix<float>::resize(size_t rlda);
    template void smatrix<double>::resize(size_t lda);
    template void smatrix<int>::resize(size_t lda);

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::resize()
    {
        try
        {
            //fclear();
            A.clear();

            rectangular = false;
            symetric = true;
            compact = true;
            ondisk = false;

            numRow = numCol = 0;
            max_elements = 0;
            ondisc_elements = 0;

            // rename file name:
            binFilename_keys.clear();
            binFilename_vals.clear();
            std::random_device rd;
            srand(rd());
            int iNum = rand() % 100000;

            binFilename_keys = "smatrix_keys_" + std::to_string(iNum);
            binFilename_vals = "smatrix_vals_" + std::to_string(iNum);
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::resize(size_t)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::resize(size_t)" << '\n';
        }        
    }

    template void smatrix<float>::resize();
    template void smatrix<double>::resize();
    template void smatrix<int>::resize();

    //===============================================================================================================

    template <typename T>
    smatrix<T>::smatrix()
    {
        /*
           Default constructor for rectangular matrix;
           no memmory allocation.
        */

        std::random_device rd;
        srand(rd());
        int iNum = rand() % 100000;

        binFilename_keys = "smatrix_key_" + std::to_string(iNum);
        binFilename_vals = "smatrix_val_" + std::to_string(iNum);

        rectangular = true;
        symetric = false;
        compact = false;
        ondisk = false;
        numRow = numCol = 0;
        max_elements = 0;
        ondisc_elements = 0;

        work_load_perthread = 1;
        zerro_tolerance = (T)1.0e-12;;

        debug_file = "SMATRIX.log";
    }

    template smatrix<float>::smatrix();
    template smatrix<double>::smatrix();
    template smatrix<int>::smatrix();

    //===============================================================================================================

    template <typename T>
    smatrix<T>::~smatrix()
    {
        /*
            Class destructor.
        */
            //fclear(); // should this be here ???
            A.clear();
    }

    template smatrix<float>::~smatrix();
    template smatrix<double>::~smatrix();
    template smatrix<int>::~smatrix();

    //===============================================================================================================

    template <typename T>
    smatrix<T>::smatrix(const smatrix<T> &obj)
    {
        try
        {
            /*
                Copy constructor.
            */
            debug_file = "SMATRIX.log";
            compact = obj.compact;
            rectangular = obj.rectangular;
            symetric = obj.symetric;
            numCol = obj.numCol;
            numRow = obj.numRow;
            max_elements = obj.max_elements;
            binFilename_keys = obj.binFilename_keys;
            binFilename_vals = obj.binFilename_vals;
            ondisc_elements = obj.ondisc_elements;
            ondisk = obj.ondisk;
            work_load_perthread = obj.work_load_perthread;
            A = obj.A;
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::smatrix(const smatrix<T> &)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::smatrix(const smatrix<T> &)" << '\n';
        }        
    }

    template smatrix<float>::smatrix(const smatrix<float> &obj);
    template smatrix<double>::smatrix(const smatrix<double> &obj);
    template smatrix<int>::smatrix(const smatrix<int> &obj);

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::set_thread_load( size_t load )
    {
        // seta a thread load in terms of a minimun allowed non-zero elements in a matrix per thread
        try
        {
            work_load_perthread = load;
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in set_thread_load( size_t )" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::set_thread_load( size_t load )" << '\n';
        }
    }

    template void smatrix<float>::set_thread_load( size_t load );
    template void smatrix<double>::set_thread_load( size_t load );
    template void smatrix<int>::set_thread_load( size_t load );

    //===============================================================================================================

    template <typename T>
    size_t smatrix<T>::nrows()
    {
        // Returns the number of rows in a dense matrix

        try
        {
            return numRow;
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::nrows()" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::nrows()" << '\n';
        }

        return 0;
    }

    template size_t smatrix<float>::nrows();
    template size_t smatrix<double>::nrows();
    template size_t smatrix<int>::nrows();

    //===============================================================================================================

    template <typename T>
    size_t smatrix<T>::ncols()
    {
        // Returns the number of rows in a dense matrix
        
        try
        {
            return numCol;
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::ncols()" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::ncols()" << '\n';
        }

        return 0;
    }

    template size_t smatrix<float>::ncols();
    template size_t smatrix<double>::ncols();
    template size_t smatrix<int>::ncols();

    //===============================================================================================================

    template <typename T>
    size_t smatrix<T>::size()
    {
        // Returns the actual number of elements in A
        
        try
        {
            return A.size();
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::size()" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::size()" << '\n';
        }

        return 0;
    }

    template size_t smatrix<float>::size();
    template size_t smatrix<double>::size();
    template size_t smatrix<int>::size();

    //===============================================================================================================

    template <typename T>
    size_t smatrix<T>::max_key()
    {
        // Returns the actual number of elements in A
        
        try
        {
            return max_elements;
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::max_key()" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::max_key()" << '\n';
        }

        return 0;
    }

    template size_t smatrix<float>::max_key();
    template size_t smatrix<double>::max_key();
    template size_t smatrix<int>::max_key();

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::clear()
    {        
        try
        {
            A.clear();
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::clear()" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::clear()" << '\n';
        }
    }

    template void smatrix<float>::clear();
    template void smatrix<double>::clear();
    template void smatrix<int>::clear();

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::fclear()
    {        
        try
        {
            std::ifstream f(binFilename_keys.c_str());
            if (f.good())
                remove(binFilename_keys.c_str());
            
            std::ifstream f2(binFilename_vals.c_str());
            if (f2.good())
                remove(binFilename_vals.c_str());
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::fclear()" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::fclear()" << '\n';
        }
    }

    template void smatrix<float>::fclear();
    template void smatrix<double>::fclear();
    template void smatrix<int>::fclear();

    //===============================================================================================================

    template <typename T>
    bool smatrix<T>::empty()
    {        
        try
        {
            return A.empty();
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::empty()" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::empty()" << '\n';
        }

        return true;
    }

    template bool smatrix<float>::empty();
    template bool smatrix<double>::empty();
    template bool smatrix<int>::empty();

    //===============================================================================================================

    template <typename T>
    T &smatrix<T>::operator()(size_t atRow, size_t atCol)
    {
        /*
            Fast access operator (array-like).

            Return value: the element at the specified position in the container.

            Example:
                    matrix <double> M(2,2);
                    M(1,1) = 0.8;
                    double val = M(1,1);   // val = 0.8.
        */

#ifdef UTEST // use this only when debugging and testing
        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. smatrix<T>::operator()");
#endif
        if (!compact)
            return A[atRow * numCol + atCol];
        else
        {
#ifdef UTEST // use this only when debugging and testing
            if ( atCol > atRow )
                throw std::string("The column value is greater than the row value. This is not allowed for symmetric matrix in compact (L-store) format!");
#endif         
            return A[atRow * (atRow + 1) / 2 + atCol];
        }
    }

    template float &smatrix<float>::operator()(size_t atRow, size_t atCol);
    template double &smatrix<double>::operator()(size_t atRow, size_t atCol);
    template int &smatrix<int>::operator()(size_t atRow, size_t atCol);

    //===============================================================================================================

    template <typename T>
    T &smatrix<T>::operator[](size_t atIndex)
    {
#ifdef UTEST // use this only when debugging and testing
        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. smatrix<T>::operator()");
        
        if ( atIndex > max_elements - 1 )
            throw std::string("The index is greater than the maximum allowed value for this matrix. This is not allowed for symmetric matrix in compact (L-store) format!");
#endif            
            return A[atIndex];
    }

    template float &smatrix<float>::operator[](size_t atIndex);
    template double &smatrix<double>::operator[](size_t atIndex);
    template int &smatrix<int>::operator[](size_t atIndex);

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::thread_loads_unord(smatrix &in, std::vector<size_t> &out)
    {        
        try
        {
            const auto processor_count = std::thread::hardware_concurrency(); //may return 0 when not able to detect

            size_t n_threads = processor_count;
            size_t map_size = in.size();
            size_t work_load = 0;
            work_load = (size_t)( map_size / n_threads ); // expected work load per thread
            size_t max_load = work_load_perthread; // max load (elements in range) per thread (should be probably changed!)

            if (work_load <= max_load) // correct the number of assigned threads (is more the case for a small data or large max_load)
            {
                n_threads = (size_t)n_threads / 2.0;
                work_load = (size_t)( map_size / n_threads );

                if (work_load <= max_load)
                {
                    n_threads = 1;
                    work_load = map_size;
                }
            }

            size_t last_index = 0; // approx. load

            for (size_t i = 0; i < n_threads - 1; i++)
            {
                last_index += work_load;
                out.push_back(last_index);                    
            }
            out.push_back( map_size-1 );
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::thread_loads_unord(smatrix &, std::vector<size_t> &)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::thread_loads_unord(smatrix &, std::vector<size_t> &)" << '\n';
        }
    }

    template void smatrix<float>::thread_loads_unord(smatrix &in, std::vector<size_t> &out);
    template void smatrix<double>::thread_loads_unord(smatrix &in, std::vector<size_t> &out);
    template void smatrix<int>::thread_loads_unord(smatrix &in, std::vector<size_t> &out);

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::thread_loads_ord(ordstorage &in, std::vector<size_t> &out)
    {        
        try
        {
            const auto processor_count = std::thread::hardware_concurrency(); //may return 0 when not able to detect

            size_t n_threads = processor_count;
            size_t map_size = in.size();
            size_t work_load = 0;
            work_load = (size_t)( map_size / n_threads ); // expected work load per thread
            size_t max_load = work_load_perthread; // max load (elements in range) per thread (should be probably changed!)

            if (work_load <= max_load) // correct the number of assigned threads (is more the case for a small data or large max_load)
            {
                n_threads = (size_t)n_threads / 2.0;
                work_load = (size_t)( map_size / n_threads );

                if (work_load <= max_load)
                {
                    n_threads = 1;
                    work_load = map_size;
                }
            }

            size_t last_index = work_load; // approx. load

            typename std::map<size_t,T>::iterator it = in.A.begin(); // because the data is ordered map

            while ( last_index < map_size )
            {
                std::advance(it, last_index); // move iterator to position of expected last element in the row
                size_t i_row = in.row_inrec( it->first );
                size_t next_row = i_row;

                while (next_row == i_row) // find of a very last element of row close to the expected work_load position
                {
                    it++; // do point iterator to the next element
                    if ( last_index++ >= (map_size - 1) )
                        break;
                    i_row = in.row_inrec( it->first ); // determine its row
                }
                last_index--;
                out.push_back(last_index);
                last_index += work_load;
                it = in.A.begin();
            }

            if ( out.empty() || (last_index - work_load) < (map_size -1))
                out.push_back(map_size-1);
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::thread_loads_ord(smatrix &, std::vector<size_t> &)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::thread_loads_ord(smatrix &, std::vector<size_t> &)" << '\n';
        }
    }

    template void smatrix<float>::thread_loads_ord(ordstorage &in, std::vector<size_t> &out);
    template void smatrix<double>::thread_loads_ord(ordstorage &in, std::vector<size_t> &out);
    template void smatrix<int>::thread_loads_ord(ordstorage &in, std::vector<size_t> &out);

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::values_inrow(ordstorage &in, std::vector<size_t> &out_elements)
    {        
        try
        {
            /*
                Returns a number of non-zero values in each row of matrix in.
                This implementation works on ordered map only!
            */
            out_elements.resize(in.numRow, 0);

            typename std::map<size_t,T>::iterator it = in.A.begin();
            
            size_t row = in.row_inrec( it->first ); // very first row
            size_t i_row = 0;

            size_t n_elements = 0;
            for (size_t i = 0; i < in.size(); i++)
            {
                i_row = in.row_inrec( it->first );

                if ( i_row != row )
                {
                    out_elements[ row ] = n_elements; // end of previous row
                    row = i_row;
                    n_elements = 1;
                }
                else
                    n_elements++;
                it++;
            }
            out_elements[ row ] = n_elements; // account for very last element
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::values_inrow(ordstorage &, std::vector<size_t> &)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::values_inrow(ordstorage &, std::vector<size_t> &)" << '\n';
        }
    }

    template void smatrix<float>::values_inrow(ordstorage &in, std::vector<size_t> &out_elements);
    template void smatrix<double>::values_inrow(ordstorage &in, std::vector<size_t> &out_elements);
    template void smatrix<int>::values_inrow(ordstorage &in, std::vector<size_t> &out_elements);

    //===============================================================================================================

    template <typename T>
    smatrix<T> smatrix<T>::operator*(smatrix<T> &rhs)
    {
        /*
            NOTE: this cannot be used for multiply matrix on itself, because rhs is modified inside !!!

            Matrix dot product:
            C = A * B;
            where A & B are rectangular matrices, C is always rectangular.

            Return value: rectangular matrix.
        */

        if (ondisk || rhs.ondisk)
            throw std::string("One of the matrices is empty. Use fread() to relocate data to memory. smatrix<T>::operator*");

        /* Check if matrices can be multiplied. */
        if (numCol != rhs.numRow)
            throw std::string("Matrices are not consistent for multiplication. smatrix<T>::operator*");
        
        if (compact || rhs.compact)
            throw std::string("Cannot provide multiplication of matrices in the compact format! smatrix<T>::operator*");

        smatrix<T> C; // container for the result of multiplication
        
        if ( rhs.size() < size() ) // here the rhs's matrix density is lower than at lhs, hence: res = lhs * transpose(rhs)
        {
            C.resize(numRow, rhs.numCol);

            smatrix<T> _rhs(rhs); // copy rhs to temp container, this'll allow rhs.transpose() only once!
            
            _rhs.transpose(); // transpose because: res = lhs * transpose(rhs)

            ordstorage ord_rhs(_rhs); // copy rhs to ordered map container

            _rhs.resize(); // we do not need it any more
            
            std::vector<size_t> loads;
            thread_loads_ord(ord_rhs,loads); // get the number of rhs's rows per aavailable threads

            std::vector<std::thread> vec_threads;
            std::vector< smatrix<T> > vec_matr;
            
            for(size_t i = 0; i < loads.size(); i++) // for each thread create temp smatrix container holding the results of specific columns
                vec_matr.emplace_back( numRow, ord_rhs.numRow ); // construct smatrix<T> object directly on vec_matr; use rhs.numRow instead of rhs.numCol because rhs is transposed now

            for(size_t i = 0; i < loads.size(); i++) // launch threads
                vec_threads.emplace_back( &smatrix::dot_operation, this, std::ref(*this), std::ref(ord_rhs), std::ref(vec_matr[i]), std::ref(loads), i );

            for (auto &v : vec_threads) // join all threads
                v.join();

            ord_rhs.clear();

            for(size_t i = 0; i < loads.size(); i++) // gathher the results of all threads into one resulting matrix
            {
                C.A.insert(vec_matr[i].A.begin(), vec_matr[i].A.end());
                vec_matr[i].resize();
            }
        }
        else  // the lhs's matrix density is lower than at rhs, hence: transpose(res) = transpose(rhs) * lhs; need to transpose the results at the end
        {
            C.resize(rhs.numCol, numRow);

            smatrix<T> _rhs(rhs); // copy rhs to temp container, this'll allow rhs.transpose() only once!

            _rhs.transpose(); // transpose because we do transpose(res) = transpose(rhs) * lhs

            ordstorage ord_this(*this); // copy lhs to uordered map container

            std::vector<size_t> loads;
            thread_loads_ord(ord_this,loads);

            std::vector<std::thread> vec_threads;
            std::vector< smatrix<T> > vec_matr;

            for(size_t i = 0; i < loads.size(); i++) // for each thread create temp smatrix container holding the results for specific columns
                vec_matr.emplace_back( _rhs.numRow, numRow ); // we use rhs.numRow instead of rhs.numCol because rhs is transposed now

            for(size_t i = 0; i < loads.size(); i++)
                vec_threads.emplace_back( &smatrix::dot_operation, this, std::ref(_rhs), std::ref(ord_this), std::ref(vec_matr[i]), std::ref(loads), i );

            for (auto &v : vec_threads)
                v.join();

            ord_this.clear();
            _rhs.resize();

            for(size_t i = 0; i < loads.size(); i++)
            {
                C.A.insert(vec_matr[i].A.begin(), vec_matr[i].A.end());
                vec_matr[i].resize();
            }

            C.transpose(); // transpose to get correct result because we do transpose(res) = transpose(rhs) * lhs
        }
        
        return C;
    }

    template smatrix<float> smatrix<float>::operator*(smatrix<float> &rhs);
    template smatrix<double> smatrix<double>::operator*(smatrix<double> &rhs);
    template smatrix<int> smatrix<int>::operator*(smatrix<int> &rhs);

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::dot_operation(smatrix &lhs, ordstorage &rhs, smatrix &out, std::vector<size_t> &range_vect, size_t thr_id)
    {
        try
        {
            std::vector<size_t> col_lst;
            std::vector<T> val_lst;
            size_t i_row = 0;
            size_t row = 0;
            size_t col = 0;

            size_t i_first = 0; // very first element in a threads range the iterator it should point to

            if ( thr_id != 0 )
                i_first = range_vect[thr_id - 1] + 1;

            size_t i_last = range_vect[thr_id]; // very last element in a threads range the iterator it should point to

            typename std::map<size_t,T>::iterator it = rhs.A.begin();
            std::advance(it, i_first);

            i_row = rhs.row_inrec( it->first ); // determine from which row of the rhs matrix we start

            size_t n_values = i_last - i_first + 1; // max number of elements in a threads range

            size_t current_element = 0;

            for(size_t i = 0; i < n_values; i++)
            {
                row = rhs.row_inrec(it->first);      // to which row the key belongs 
                col = rhs.col_inrec(it->first, row); // knowing the row find the col
                T value = it->second;

                current_element++;

                if ( current_element == n_values ) // very last element of the matrix
                {
                    if (row != i_row) // indicates we reached last row with a single element in a row; calculate for the previous row data
                    {
                        for (size_t i = 0; i < lhs.nrows(); i++)
                        {
                            T res = (T)0;
                            T zerro_val = (T)0;
                            T lhs_val = (T)0;
                            for (size_t j = 0; j < val_lst.size(); j++)
                            {
                                lhs_val = lhs.get_nonzero( i,col_lst[j] );
                                if ( lhs_val != zerro_val )
                                    res += lhs_val * val_lst[j];
                            }

                            if ( res != zerro_val )
                                out[ out.key_inrec(i,i_row) ] = res;
                        }
                        col_lst.clear(); // clear the previous row data
                        val_lst.clear(); // clear the previous row data
                    }

                    col_lst.push_back(col);
                    val_lst.push_back(value);

                    i_row = row;

                    for (size_t i = 0; i < lhs.nrows(); i++) // calculate a very last column of the result matrix C
                    {
                        T res = (T)0;
                        T zerro_val = (T)0;
                        T lhs_val = (T)0;
                        for (size_t j = 0; j < val_lst.size(); j++)
                        {
                            lhs_val = lhs.get_nonzero( i,col_lst[j] );
                            if ( lhs_val != zerro_val )
                                res += lhs_val * val_lst[j];
                        }

                        if ( res != zerro_val )
                            out[ out.key_inrec(i,i_row) ] = res;
                    }

                    continue;                    
                }

                if ( i_row == row )
                {
                    col_lst.push_back(col);
                    val_lst.push_back(value);
                    it++;
                }
                else
                {
                    if ( val_lst.empty() )
                    {
                        i_row = row;
                        col_lst.push_back(col); // we are in the new row now, so do not miss very first data
                        val_lst.push_back(value); // we are in the new row now, so do not miss very first data
                        it++;
                        continue;
                    }
                    
                    for (size_t i = 0; i < lhs.nrows(); i++)
                    {
                        T res = (T)0;
                        T zerro_val = (T)0;
                        T lhs_val = (T)0;
                        for (size_t j = 0; j < val_lst.size(); j++)
                        {
                            lhs_val = lhs.get_nonzero( i,col_lst[j] );
                            if ( lhs_val != zerro_val )
                                res += lhs_val * val_lst[j];
                        }

                        if ( res != zerro_val )
                            out[ out.key_inrec(i,i_row) ] = res;
                    }

                    i_row = row; // set the new current row
                    col_lst.clear(); // clear the previous row data
                    val_lst.clear(); // clear the previous row data
                    col_lst.push_back(col); // we are in the new row now, so do not miss very first data
                    val_lst.push_back(value); // we are in the new row now, so do not miss very first data
                    it++;
                }
//if (current_element%100 == 0 && thr_id == 0)
//    std::cout<<"\rcompletes, %: "<<current_element*100/n_values<<std::flush;

            }
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::dot_operation_unordered(smatrix &, ordstorage &, smatrix &, std::vector<size_t> &, size_t)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::dot_operation_unordered(smatrix &, ordstorage &, smatrix &, std::vector<size_t> &, size_t)" << '\n';
        }        
    }

    template void smatrix<double>::dot_operation(smatrix &lhs, ordstorage &rhs, smatrix &out, std::vector<size_t> &range_vect, size_t thr_id);
    template void smatrix<float>::dot_operation(smatrix &lhs, ordstorage &rhs, smatrix &out, std::vector<size_t> &range_vect, size_t thr_id);
    template void smatrix<int>::dot_operation(smatrix &lhs, ordstorage &rhs, smatrix &out, std::vector<size_t> &range_vect, size_t thr_id);

    //===============================================================================================================

    template <typename T>
    smatrix<T> smatrix<T>::operator+(smatrix<T> &rhs)
    {
        /*
            Matrix sum:
            C = A + B;
            where A & B are rectangular matrices, C is always rectangular.

            Return value: rectangular matrix.
        */

        if (ondisk || rhs.ondisk)
            throw std::string("One of the matrices is empty. Use fread() to relocate data to memory. smatrix<T>::operator+");

        // Check if matrices can be added.
        if ( (numCol != rhs.numCol) || (numRow != rhs.numRow) )
            throw std::string("Matrices are not consistent for addition. smatrix<T>::operator+");
        
        if ( (compact && !rhs.compact) || (!compact && rhs.compact) )
            throw std::string("Cannot provide addition due to matrices in different format! smatrix<T>::operator+");
        
        smatrix<T> Result(*this);

        //------------------------------------
        T i_val = (T)0;
        T almost_zerro = zerro_tolerance;
        T rslt = (T)0;
        std::vector<size_t> zeros_keys;

        for (auto & e: rhs.A)
        {
            i_val = Result[e.first];
            rslt = i_val + e.second;
            Result[ e.first ] = rslt; // modify (update) the result matrix
            if ( fabs(rslt) <= almost_zerro ) // keep track of the keys with zero values to be removed from matrix later
                zeros_keys.push_back(e.first); // we do not want to store zeros !
        }
        //-------------------------------------
        /*
        std::vector<size_t> loads;

        thread_loads_unord(rhs, loads);

        std::vector<smatrix> vect_matr;
        std::vector<std::thread> vect_thr;

        std::vector<size_t> zeros_keys;

        if (compact)
        {
            for(size_t i = 0; i < loads.size(); i++)
                vect_matr.emplace_back(numRow);
        }
        else
        {
            for(size_t i = 0; i < loads.size(); i++)
                vect_matr.emplace_back(numRow, numCol);
        }

        for(size_t i = 0; i < loads.size(); i++)
            vect_thr.emplace_back( &smatrix::plus_operation, this, std::ref(rhs), std::ref(Result), std::ref(zeros_keys), std::ref(vect_matr[i]), std::ref(loads), i );

        for(size_t i = 0; i < loads.size(); i++)
            vect_thr[i].join();

        for(size_t i = 0; i < loads.size(); i++)
        {
            Result.A.insert( vect_matr[i].A.begin(), vect_matr[i].A.end() );
            vect_matr[i].resize();
        }
        */
        
        // remove from the matrix values equal to zero
        for (size_t i = 0; i < zeros_keys.size(); i++)
            Result.A.erase( zeros_keys[i] );
        
        return Result;
    }

    template smatrix<float> smatrix<float>::operator+(smatrix<float> &rhs);
    template smatrix<double> smatrix<double>::operator+(smatrix<double> &rhs);
    template smatrix<int> smatrix<int>::operator+(smatrix<int> &rhs);

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::plus_operation(smatrix &in, smatrix &res, std::vector<size_t> &zero_keys, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id)
    {
        /*
            Addition operation to be handled by the specific thread.
            in - RHS matrix;
            res - the result of operation final matrix, here is a copy of LHS;
            zero_keys - empty vector where keys corresponding to zero values in the res matrix will be stored;
            out - local (in terms of a specific thread) result matrix with keys not existing in the res matrix, will be merged with res at the end;
            loads_vect - amount of records in the RHS to be processed by the specific thread;
            thr_id - the specific thread ID.
        */
        try
        {
            size_t row = 0;
            size_t col = 0;
            size_t key = 0;
            T i_val = (T)0;
            T zerro_val = (T)0;
            T almost_zerro = zerro_tolerance;
            T rslt = (T)0;

            size_t i_first = 0; // very first element in a threads range the iterator it should point to

            if ( thr_id != 0 )
                i_first = loads_vect[thr_id - 1] + 1;

            size_t i_last = loads_vect[thr_id]; // very last element in a threads range the iterator it should point to

            typename std::unordered_map<size_t, T>::iterator it = in.A.begin();
            std::advance(it, i_first);

            if (compact)
            {
                for(size_t i = 0; i < ( i_last - i_first + 1 ); i++) // work on in.A
                {
                    row = in.row_insym(it->first);      // to which row the key belongs 
                    col = in.col_insym(it->first, row); // knowing the row find the col
                    key = in.key_insym(row,col);        // getting the key

                    i_val = res.get_nonzero(key); // this is lhs's value !

                    if ( i_val != zerro_val ) // if the key in res exists (hence, a non-zero value is)
                    {
                        rslt = i_val + it->second;
                        
                        res[ key ] = rslt; // modify (update) the result matrix

                        if ( fabs(rslt) <= almost_zerro ) // keep track of the keys with zero values to be removed from matrix later
                            zero_keys.push_back(key); // we do not want to store zeros !
                    }
                    else
                        out[ key ] = it->second; // technically new map with keys not exist in res mmatrix, to be merged with  res later

                    it++;
                }
            }
            else
            {
                for(size_t i = 0; i < ( i_last - i_first + 1 ); i++) // work on in.A
                {
                    row = in.row_inrec(it->first);      // to which row the key belongs 
                    col = in.col_inrec(it->first, row); // knowing the row find the col
                    key = in.key_inrec(row,col);        // getting the key

                    i_val = res.get_nonzero(key); // this is lhs's value !

                    if ( i_val != zerro_val ) // if the key in res exists (hence, a non-zero value is)
                    {
                        rslt = i_val + it->second;
                        
                        res[ key ] = rslt; // modify (update) the result matrix

                        if ( fabs(rslt) <= almost_zerro ) // keep track of the keys with zero values to be removed from matrix later
                            zero_keys.push_back(key); // we do not want to store zeros !
                    }
                    else
                        out[ key ] = it->second; // technically new map with keys not exist in res mmatrix, to be merged with  res later

                    it++;
                }
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::plus_operation(smatrix &, smatrix &, std::vector<size_t> &, smatrix &, std::vector<size_t> &, size_t)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::plus_operation(smatrix &, smatrix &, std::vector<size_t> &, smatrix &, std::vector<size_t> &, size_t)" << '\n';
        }        
    }

    template void smatrix<double>::plus_operation(smatrix &in, smatrix &res, std::vector<size_t> &zero_keys, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);
    template void smatrix<float>::plus_operation(smatrix &in, smatrix &res, std::vector<size_t> &zero_keys, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);
    template void smatrix<int>::plus_operation(smatrix &in, smatrix &res, std::vector<size_t> &zero_keys, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);

    //===============================================================================================================

    template <typename T>
    smatrix<T> smatrix<T>::operator-(smatrix<T> &rhs)
    {
        /*
            Matrix substitution:
            C = A - B;
            where A & B are rectangular matrices, C is always rectangular.

            Return value: rectangular matrix.
        */

        if (ondisk || rhs.ondisk)
            throw std::string("One of the matrices is empty. Use fread() to relocate data to memory. smatrix<T>::operator-");

        /* Check if matrices can be added. */
        if ( (numCol != rhs.numCol) || (numRow != rhs.numRow) )
            throw std::string("Matrices are not consistent for substitution. smatrix<T>::operator-");
        
        if ( (compact && !rhs.compact) || (!compact && rhs.compact) )
            throw std::string("Cannot provide substitution due to matrices in different format! smatrix<T>::operator-");

        smatrix<T> Result(*this);
        //------------------------------------
        T i_val = (T)0;
        T almost_zerro = zerro_tolerance;
        T rslt = (T)0;
        std::vector<size_t> zeros_keys;

        for (auto & e: rhs.A)
        {
            i_val = Result[e.first];
            rslt = i_val - e.second;
            Result[ e.first ] = rslt; // modify (update) the result matrix
            if ( fabs(rslt) <= almost_zerro ) // keep track of the keys with zero values to be removed from matrix later
                zeros_keys.push_back(e.first); // we do not want to store zeros !
        }
        //-------------------------------------
        /*
        std::vector<size_t> loads;

        thread_loads_unord(rhs, loads);


        std::vector<smatrix> vect_matr;
        std::vector<std::thread> vect_thr;

        std::vector<size_t> zeros_keys;

        if (compact)
        {
            for(size_t i = 0; i < loads.size(); i++)
                vect_matr.emplace_back(numRow);
        }
        else
        {
            for(size_t i = 0; i < loads.size(); i++)
                vect_matr.emplace_back(numRow, numCol);
        }

        for(size_t i = 0; i < loads.size(); i++)
            vect_thr.emplace_back( &smatrix::minus_operation, this, std::ref(rhs), std::ref(Result), std::ref(zeros_keys), std::ref(vect_matr[i]), std::ref(loads), i );

        for(size_t i = 0; i < loads.size(); i++)
            vect_thr[i].join();


        for(size_t i = 0; i < loads.size(); i++)
        {
            Result.A.insert( vect_matr[i].A.begin(), vect_matr[i].A.end() );
            vect_matr[i].resize();
        }
        */            
        
        // remove from the matrix values equal to zero
        for (size_t i = 0; i < zeros_keys.size(); i++)
            Result.A.erase( zeros_keys[i] );

        return Result;
    }

    template smatrix<float> smatrix<float>::operator-(smatrix<float> &rhs);
    template smatrix<double> smatrix<double>::operator-(smatrix<double> &rhs);
    template smatrix<int> smatrix<int>::operator-(smatrix<int> &rhs);

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::minus_operation(smatrix &in, smatrix &res, std::vector<size_t> &zero_keys, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id)
    {
        /*
            Substraction operation to be handled by the specific thread.
            in - RHS matrix;
            res - the result of operation final matrix, here is a copy of LHS;
            zero_keys - empty vector where keys corresponding to zero values in the res matrix will be stored;
            out - local (in terms of a specific thread) result matrix with keys not existing in the res matrix, will be merged with res at the end;
            loads_vect - amount of records in the RHS to be processed by the specific thread;
            thr_id - the specific thread ID.
        */
        try
        {
            size_t row = 0;
            size_t col = 0;
            size_t key = 0;
            T i_val = (T)0;
            T zerro_val = (T)0;
            T almost_zerro = zerro_tolerance;
            T rslt = (T)0;

            size_t i_first = 0; // very first element in a threads range the iterator it should point to

            if ( thr_id != 0 )
                i_first = loads_vect[thr_id - 1] + 1;

            size_t i_last = loads_vect[thr_id]; // very last element in a threads range the iterator it should point to

            typename std::unordered_map<size_t, T>::iterator it = in.A.begin();
            std::advance(it, i_first);

            if (compact)
            {
                for(size_t i = 0; i < ( i_last - i_first + 1 ); i++) // work on in.A
                {
                    row = in.row_insym(it->first);      // to which row the key belongs 
                    col = in.col_insym(it->first, row); // knowing the row find the col
                    key = in.key_insym(row,col);        // getting the key
                    
                    i_val = res.get_nonzero(key); // this is lhs's value !

                    if ( i_val != zerro_val ) // if the key in res exists (hence, a non-zero value is)
                    {
                        rslt = i_val - it->second;
                        
                        res[ key ] = rslt; // modify (update) the result matrix

                        if ( fabs(rslt) <= almost_zerro ) // keep track of the keys with zero values to be removed from matrix later
                            zero_keys.push_back(key); // we do not want to store zeros !
                    }
                    else
                        out[ key ] = -it->second; // technically new map with keys not exist in res mmatrix, to be merged with  res later
                    
                    it++;
                }
            }
            else
            {
                for(size_t i = 0; i < ( i_last - i_first + 1 ); i++) // work on in.A
                {
                    row = in.row_inrec(it->first);      // to which row the key belongs 
                    col = in.col_inrec(it->first, row); // knowing the row find the col
                    key = in.key_inrec(row,col);        // getting the key
                    
                    i_val = res.get_nonzero(key); // this is lhs's value !

                    if ( i_val != zerro_val ) // if the key in res exists (hence, a non-zero value is)
                    {
                        rslt = i_val - it->second;
                        
                        res[ key ] = rslt; // modify (update) the result matrix

                        if ( fabs(rslt) <= almost_zerro ) // keep track of the keys with zero values to be removed from matrix later
                            zero_keys.push_back(key); // we do not want to store zeros !
                    }
                    else
                        out[ key ] = -it->second; // technically new map with keys not exist in res mmatrix, to be merged with  res later
                    
                    it++;
                }
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::minus_operation(smatrix &, smatrix &, std::vector<size_t> &, smatrix &, std::vector<size_t> &, size_t)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::minus_operation(smatrix &, smatrix &, std::vector<size_t> &, smatrix &, std::vector<size_t> &, size_t)" << '\n';
        }        
    }

    template void smatrix<double>::minus_operation(smatrix &in, smatrix &res, std::vector<size_t> &zero_keys, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);
    template void smatrix<float>::minus_operation(smatrix &in, smatrix &res, std::vector<size_t> &zero_keys, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);
    template void smatrix<int>::minus_operation(smatrix &in, smatrix &res, std::vector<size_t> &zero_keys, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);

    //===============================================================================================================

    template <typename T>
    bool smatrix<T>::nonzero(size_t atRow, size_t atCol)
    {
        try
        {
            /* Checks if the specific matrix element is exists in the container (has the value different from 0).

                NOTE: do not call this method if ...
            */

#ifdef UTEST // use this only when debugging and testing
            if (ondisk)
                throw std::string("Matrix is empty. Use fread() to relocate data to memory. smatrix<T>::operator()");
#endif
            //typename std::map<size_t, T>::iterator it;
            typename std::unordered_map<size_t, T>::iterator it;
            size_t key = 0;

            if (!compact)
                key = atRow * numCol + atCol;
            else
            {
#ifdef UTEST // use this only when debugging and testing
                if ( atCol > atRow )
                    throw std::string("The column value is greater than the row value. This is not allowed for symmetric matrix in compact (L-store) format!");
#endif                
                key = atRow * (atRow + 1) / 2 + atCol;
            }

            it = A.find( key );

            if (it != A.end())
                return true;
            else
                return false;
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::nonzero(size_t, size_t)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(const std::string &e)
        {
            std::cerr << "Exception in smatrix<T>::nonzero(size_t, size_t)" << '\n';
            std::cerr << e << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::nonzero(size_t, size_t)" << '\n';
        }

        return false;      
    }

    template bool smatrix<float>::nonzero(size_t atRow, size_t atCol);
    template bool smatrix<double>::nonzero(size_t atRow, size_t atCol);
    template bool smatrix<int>::nonzero(size_t atRow, size_t atCol);

    //===============================================================================================================

    template <typename T>
    T smatrix<T>::get_nonzero(size_t atRow, size_t atCol)
    {
        try
        {
            /* Checks if the specific matrix element is exists in the container (has the value different from 0).
               If exists, return value, otherwise returns 0.
            */

#ifdef UTEST // use this only when debugging and testing
            if (ondisk)
                throw std::string("Matrix is empty. Use fread() to relocate data to memory. smatrix<T>::operator()");
#endif
            //typename std::map<size_t, T>::iterator it;
            typename std::unordered_map<size_t, T>::iterator it;

            size_t key = 0;

            if (!compact)
                key = atRow * numCol + atCol;
            else
            {
#ifdef UTEST // use this only when debugging and testing
                if ( atCol > atRow )
                    throw std::string("The column value is greater than the row value. This is not allowed for symmetric matrix in compact (L-store) format!");
#endif                
                key = atRow * (atRow + 1) / 2 + atCol;
            }

            it = A.find( key );

            if (it != A.end())
                return it->second;
            else
                return (T)0;
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::get_nonzero(size_t, size_t)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(const std::string &e)
        {
            std::cerr << "Exception in smatrix<T>::get_nonzero(size_t, size_t)" << '\n';
            std::cerr << e << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::get_nonzero(size_t, size_t)" << '\n';
        }

        return (T)0;      
    }

    template float smatrix<float>::get_nonzero(size_t atRow, size_t atCol);
    template double smatrix<double>::get_nonzero(size_t atRow, size_t atCol);
    template int smatrix<int>::get_nonzero(size_t atRow, size_t atCol);

    //===============================================================================================================

    template <typename T>
    bool smatrix<T>::nonzero(size_t key)
    {
        try
        {
            /* Checks if the specific matrix element is exists in the container (has the value different from 0).
            */

#ifdef UTEST // use this only when debugging and testing
            if (ondisk)
                throw std::string("Matrix is empty. Use fread() to relocate data to memory. smatrix<T>::operator()");
#endif
            //typename std::map<size_t, T>::iterator it;
            typename std::unordered_map<size_t, T>::iterator it;

            it = A.find( key );

            if (it != A.end())
                return true;
            else
                return false;
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::nonzero(size_t)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(const std::string &e)
        {
            std::cerr << "Exception in smatrix<T>::nonzero(size_t)" << '\n';
            std::cerr << e << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::nonzero(size_t)" << '\n';
        }

        return false;      
    }

    template bool smatrix<float>::nonzero(size_t atRow);
    template bool smatrix<double>::nonzero(size_t atRow);
    template bool smatrix<int>::nonzero(size_t atRow);

    //===============================================================================================================

    template <typename T>
    T smatrix<T>::get_nonzero(size_t key)
    {
        try
        {
            /* Checks if the specific matrix element is exists in the container (has the value different from 0).
               If exists, return value, otherwise, return 0.
            */

#ifdef UTEST // use this only when debugging and testing
            if (ondisk)
                throw std::string("Matrix is empty. Use fread() to relocate data to memory. smatrix<T>::operator()");
#endif
            //typename std::map<size_t, T>::iterator it;
            typename std::unordered_map<size_t, T>::iterator it;

            it = A.find( key );

            if (it != A.end())
                return it->second;
            else
                return (T)0;
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::get_nonzero(size_t)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(const std::string &e)
        {
            std::cerr << "Exception in smatrix<T>::get_nonzero(size_t)" << '\n';
            std::cerr << e << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::get_nonzero(size_t)" << '\n';
        }

        return (T)0;      
    }

    template float smatrix<float>::get_nonzero(size_t atRow);
    template double smatrix<double>::get_nonzero(size_t atRow);
    template int smatrix<int>::get_nonzero(size_t atRow);

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::fwrite()
    {
        try
        {
            /*
                Moves matrix to a DISK and clears memory allocated for the container A.

                Return value: none.
            */
            fA_keys.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            fA_keys.open(binFilename_keys, fA_keys.binary | fA_keys.trunc | fA_keys.out);

            if (!fA_keys.is_open())
                throw std::string("Error while opening a keys binary file. smatrix<T>::fwrite()");

            fA_vals.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            fA_vals.open(binFilename_vals, fA_vals.binary | fA_vals.trunc | fA_vals.out);

            if (!fA_vals.is_open())
                throw std::string("Error while opening a values binary file. smatrix<T>::fwrite()");

            for (auto const& v : A)
            {
                fA_keys.write(reinterpret_cast<const char *>(&v.first), sizeof(size_t));
                fA_vals.write(reinterpret_cast<const char *>(&v.second), sizeof(T));
            }

            fA_keys.close();
            fA_vals.close();

            ondisk = true;
            ondisc_elements = A.size();

            A.clear();
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::fwrite()" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::fwrite()" << '\n';
        }            
    }

    template void smatrix<float>::fwrite();
    template void smatrix<double>::fwrite();
    template void smatrix<int>::fwrite();

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::fread()
    {
        try
        {
            /*
                Reads matrix from DISK back to the memory.

                Return value: none.
            */

            fA_keys.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            fA_keys.open(binFilename_keys, fA_keys.binary | fA_keys.in);

            if (!fA_keys.is_open())
                throw std::string("Error while opening a keys' binary file. smatrix<T>::fread()");

            fA_vals.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            fA_vals.open(binFilename_vals, fA_vals.binary | fA_vals.in);

            if (!fA_vals.is_open())
                throw std::string("Error while opening a values' binary file. smatrix<T>::fread()");

            for (size_t i = 0; i < ondisc_elements; i++)
            {
                size_t key;
                T val;
                fA_keys.read(reinterpret_cast<char *>(&key), sizeof(size_t));
                fA_vals.read(reinterpret_cast<char *>(&val), sizeof(T));
                A[key] = val;
            }

            fA_keys.close();
            fA_vals.close();

            ondisk = false;
            ondisc_elements = 0;
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::fread()" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::fread()" << '\n';
        }            
    }

    template void smatrix<float>::fread();
    template void smatrix<double>::fread();
    template void smatrix<int>::fread();

    //===============================================================================================================

    template <typename T>
    smatrix<T> &smatrix<T>::operator=(const smatrix<T> &rhs)
    {
        try
        {
            /*
                Overloaded assignment operator.

                Procedure: (1) Copy rhs to tmpObj; (2) Swaps lhs with tmpObj.
            */

            smatrix<T> tmpObj(rhs);

            std::swap(compact, tmpObj.compact);
            std::swap(ondisk, tmpObj.ondisk);
            std::swap(rectangular, tmpObj.rectangular);
            std::swap(symetric, tmpObj.symetric);
            std::swap(numCol, tmpObj.numCol);
            std::swap(numRow, tmpObj.numRow);
            std::swap(max_elements, tmpObj.max_elements);
            std::swap(ondisc_elements, tmpObj.ondisc_elements);
            std::swap(work_load_perthread, tmpObj.work_load_perthread);
            
            if ( ondisk ) // exchange file names only if rhs is not on memory!
            {
                std::swap(binFilename_keys, tmpObj.binFilename_keys);
                std::swap(binFilename_vals, tmpObj.binFilename_vals);
            }
            
            A.swap( tmpObj.A );
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::operator=(const smatrix<T> &)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::operator=(const smatrix<T> &)" << '\n';
        }            
        
        return *this;
    }

    template smatrix<float> &smatrix<float>::operator=(const smatrix<float> &rhs);
    template smatrix<double> &smatrix<double>::operator=(const smatrix<double> &rhs);
    template smatrix<int> &smatrix<int>::operator=(const smatrix<int> &rhs);

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::transpose()
    {
        try
        {
            if (ondisk)
                throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::transpose()");

            if (!compact)
            {
                std::vector<size_t> loads;
                thread_loads_unord(*this, loads);

                std::vector<smatrix> vect_matr;
                std::vector<std::thread> vect_thr;

                for(size_t i = 0; i < loads.size(); i++)
                    vect_matr.emplace_back(numCol, numRow);

                for(size_t i = 0; i < loads.size(); i++)
                    vect_thr.emplace_back( &smatrix::transpose_operation, this, std::ref(*this), std::ref(vect_matr[i]), std::ref(loads), i );

                for(size_t i = 0; i < loads.size(); i++)
                    vect_thr[i].join();

                size_t new_row = numCol;
                size_t new_col = numRow;

                numRow = new_row;
                numCol = new_col;

                A.clear();

                for(size_t i = 0; i < loads.size(); i++)
                {
                    A.insert( vect_matr[i].A.begin(), vect_matr[i].A.end() );
                    vect_matr[i].resize();
                }
            }
            else
            {
                /*We also have to accept symmetric matrices in compact form.
                Due to the reasons described in the upper comment, we return the same matrix
                but in restored (not compact) format.*/

                //symtorec();
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::transpose()" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::transpose()" << '\n';
        }                       
    }

    template void smatrix<float>::transpose();
    template void smatrix<double>::transpose();
    template void smatrix<int>::transpose();

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::transpose_operation(smatrix &in, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id)
    {
        try
        {
            size_t row = 0;
            size_t col = 0;

            size_t i_first = 0; // very first element in a threads range the iterator it should point to

            if ( thr_id != 0 )
                i_first = loads_vect[thr_id - 1] + 1;

            size_t i_last = loads_vect[thr_id]; // very last element in a threads range the iterator it should point to

            typename std::unordered_map<size_t, T>::iterator it = in.A.begin();

            std::advance(it, i_first);

            for(size_t i = 0; i < ( i_last - i_first + 1 ); i++)
            {
                row = in.row_inrec(it->first);      // to which row the key belongs 
                col = in.col_inrec(it->first, row); // knowing the row find the col
                out[ out.key_inrec(col,row) ] = it->second;
                it++;
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::transpose_operation(smatrix &, smatrix &, std::vector<size_t> &, size_t)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::transpose_operation(smatrix &, smatrix &, std::vector<size_t> &, size_t)" << '\n';
        }        
    }

    template void smatrix<double>::transpose_operation(smatrix &in, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);
    template void smatrix<float>::transpose_operation(smatrix &in, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);
    template void smatrix<int>::transpose_operation(smatrix &in, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::rectosym()
    {
        try
        {
            if (ondisk)
                throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::transpose()");

            if (compact)
                return;

            std::vector<size_t> loads;
            thread_loads_unord(*this, loads);

            std::vector<smatrix> vect_matr;
            std::vector<std::thread> vect_thr;

            for(size_t i = 0; i < loads.size(); i++)
                vect_matr.emplace_back(numRow);

            for(size_t i = 0; i < loads.size(); i++)
                vect_thr.emplace_back( &smatrix::rectosym_operation, this, std::ref(*this), std::ref(vect_matr[i]), std::ref(loads), i );

            for(size_t i = 0; i < loads.size(); i++)
                vect_thr[i].join();

            resize(numRow); // resizing the host matrix

            for(size_t i = 0; i < loads.size(); i++)
            {
                A.insert( vect_matr[i].A.begin(), vect_matr[i].A.end() );
                vect_matr[i].resize();
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::rectosym()" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::rectosym()" << '\n';
        }                       
    }

    template void smatrix<float>::rectosym();
    template void smatrix<double>::rectosym();
    template void smatrix<int>::rectosym();

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::rectosym_operation(smatrix &in, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id)
    {
        try
        {
            size_t row = 0;
            size_t col = 0;

            size_t i_first = 0; // very first element in a threads range the iterator it should point to

            if ( thr_id != 0 )
                i_first = loads_vect[thr_id - 1] + 1;

            size_t i_last = loads_vect[thr_id]; // very last element in a threads range the iterator it should point to

            typename std::unordered_map<size_t, T>::iterator it = in.A.begin();

            std::advance(it, i_first);

            for(size_t i = 0; i < ( i_last - i_first + 1 ); i++)
            {
                row = in.row_inrec(it->first);      // to which row the key belongs 
                col = in.col_inrec(it->first, row); // knowing the row find the col
                if (col <= row)
                out[ out.key_insym(row,col) ] = it->second;
                it++;
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::rectosym_operation(smatrix &, smatrix &, std::vector<size_t> &, size_t)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::rectosym_operation(smatrix &, smatrix &, std::vector<size_t> &, size_t)" << '\n';
        }        
    }

    template void smatrix<double>::rectosym_operation(smatrix &in, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);
    template void smatrix<float>::rectosym_operation(smatrix &in, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);
    template void smatrix<int>::rectosym_operation(smatrix &in, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::symtorec()
    {
        try
        {
            if (ondisk)
                throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::transpose()");

            if (!compact)
                return;

            std::vector<size_t> loads;
            thread_loads_unord(*this, loads);

            std::vector<smatrix> vect_matr;
            std::vector<std::thread> vect_thr;

            for(size_t i = 0; i < loads.size(); i++)
                vect_matr.emplace_back(numRow, numRow);

            for(size_t i = 0; i < loads.size(); i++)
                vect_thr.emplace_back( &smatrix::symtorec_operation, this, std::ref(*this), std::ref(vect_matr[i]), std::ref(loads), i );

            for(size_t i = 0; i < loads.size(); i++)
                vect_thr[i].join();

            resize(numRow, numRow); // resizing the host matrix

            for(size_t i = 0; i < loads.size(); i++)
            {
                A.insert( vect_matr[i].A.begin(), vect_matr[i].A.end() );
                vect_matr[i].resize();
            }

        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::symtorec()" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::symtorec()" << '\n';
        }                       
    }

    template void smatrix<float>::symtorec();
    template void smatrix<double>::symtorec();
    template void smatrix<int>::symtorec();

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::symtorec_operation(smatrix &in, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id)
    {
        try
        {
            size_t row = 0;
            size_t col = 0;

            size_t i_first = 0; // very first element in a threads range the iterator it should point to

            if ( thr_id != 0 )
                i_first = loads_vect[thr_id - 1] + 1;

            size_t i_last = loads_vect[thr_id]; // very last element in a threads range the iterator it should point to

            typename std::unordered_map<size_t, T>::iterator it = in.A.begin();

            std::advance(it, i_first);

            for(size_t i = 0; i < ( i_last - i_first + 1 ); i++)
            {
                row = in.row_insym(it->first);      // to which row the key belongs 
                col = in.col_insym(it->first, row); // knowing the row find the col
                out[ out.key_inrec(row,col) ] = it->second;
                out[ out.key_inrec(col,row) ] = it->second;
                it++;
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::symtorec_operation(smatrix &, smatrix &, std::vector<size_t> &, size_t)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::symtorec_operation(smatrix &, smatrix &, std::vector<size_t> &, size_t)" << '\n';
        }        
    }

    template void smatrix<double>::symtorec_operation(smatrix &in, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);
    template void smatrix<float>::symtorec_operation(smatrix &in, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);
    template void smatrix<int>::symtorec_operation(smatrix &in, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);

    //===============================================================================================================

    template <typename T>
    size_t smatrix<T>::find_inrange(std::vector<size_t> &where, size_t what)
    {
        if ( what <= where[0] )
            return 0;
        
        for (size_t i = 1; i < where.size(); i++)
        {
            if ( (what > where[i-1]) && (what <= where[i]) )
                return i;
        }

        return -1;
    }

    template size_t smatrix<float>::find_inrange(std::vector<size_t> &where, size_t what);
    template size_t smatrix<double>::find_inrange(std::vector<size_t> &where, size_t what);
    template size_t smatrix<int>::find_inrange(std::vector<size_t> &where, size_t what);

    //===============================================================================================================

    template <typename T>
    size_t smatrix<T>::find_invect(std::vector<size_t> &where, size_t what)
    {
        std::vector<size_t>::iterator it;
        
        it = std::find(where.begin(), where.end(), what);
        
        if (it != where.end()) 
            return it - where.begin();
        else
            throw std::string("smatrix<T>::find_invect(std::vector<size_t> &, size_t) => The index in the sparse matrix cannot be find in the precalculated vector of keys for the same matrix!");
    }

    template size_t smatrix<float>::find_invect(std::vector<size_t> &where, size_t what);
    template size_t smatrix<double>::find_invect(std::vector<size_t> &where, size_t what);
    template size_t smatrix<int>::find_invect(std::vector<size_t> &where, size_t what);

    //===============================================================================================================

    template <typename T>
    size_t smatrix<T>::key_insym(size_t row, size_t col)
    {
#ifdef UTEST // use this only when debugging and testing
        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::transpose()");

        if (!compact)
            throw std::string("The key_insym(size_t key, size_t row) should be called only on symmetric matrix!");
#endif
        return (size_t)( col + row * (row + 1) / 2.0 );
    }

    template size_t smatrix<float>::key_insym(size_t key, size_t row);
    template size_t smatrix<double>::key_insym(size_t key, size_t row);
    template size_t smatrix<int>::key_insym(size_t key, size_t row);

    //===============================================================================================================

    template <typename T>
    size_t smatrix<T>::key_inrec(size_t row, size_t col)
    {
#ifdef UTEST // use this only when debugging and testing
        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::transpose()");

        if (compact)
            throw std::string("The key_inrec(size_t key, size_t row) should be called only on rectangular matrix!");
#endif
        return (size_t)( col + row * numCol );
    }

    template size_t smatrix<float>::key_inrec(size_t key, size_t row);
    template size_t smatrix<double>::key_inrec(size_t key, size_t row);
    template size_t smatrix<int>::key_inrec(size_t key, size_t row);

    //===============================================================================================================

    template <typename T>
    size_t smatrix<T>::col_insym(size_t key, size_t row)
    {
#ifdef UTEST // use this only when debugging and testing
        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::transpose()");

        if (!compact)
            throw std::string("The col_insym(size_t key, size_t row) should be called only on symmetric matrix!");
#endif
        return (size_t)( key - row * (row + 1) / 2.0 );
    }

    template size_t smatrix<float>::col_insym(size_t key, size_t row);
    template size_t smatrix<double>::col_insym(size_t key, size_t row);
    template size_t smatrix<int>::col_insym(size_t key, size_t row);

    //===============================================================================================================

    template <typename T>
    size_t smatrix<T>::row_insym(size_t key, size_t col)
    {
#ifdef UTEST // use this only when debugging and testing
        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::transpose()");

        if (!compact)
            throw std::string("The row_insym(size_t key, size_t col) should be called only on symmetric matrix!");
#endif
        return round( sqrt( 0.25 + 2.0 * (key - col) ) - 0.5 );
    }

    template size_t smatrix<float>::row_insym(size_t key, size_t col);
    template size_t smatrix<double>::row_insym(size_t key, size_t col);
    template size_t smatrix<int>::row_insym(size_t key, size_t col);

    //===============================================================================================================

    template <typename T>
    size_t smatrix<T>::row_insym(size_t key)
    {
#ifdef UTEST // use this only when debugging and testing
        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::transpose()");

        if (!compact)
            throw std::string("The row_insym(size_t key, size_t col) should be called only on symmetric matrix!");
#endif
        return (size_t)( sqrt( 0.25 + (2 * key + 1.0) ) - 0.5 );
    }

    template size_t smatrix<float>::row_insym(size_t key);
    template size_t smatrix<double>::row_insym(size_t key);
    template size_t smatrix<int>::row_insym(size_t key);

    //===============================================================================================================

    template <typename T>
    size_t smatrix<T>::col_inrec(size_t key, size_t row)
    {
#ifdef UTEST // use this only when debugging and testing
        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::transpose()");

        if (compact)
            throw std::string("The col_inrec(size_t key, size_t row) should be called only on rectangular matrix!");
#endif
        return (size_t) ( key - row * numCol );
    }

    template size_t smatrix<float>::col_inrec(size_t key, size_t row);
    template size_t smatrix<double>::col_inrec(size_t key, size_t row);
    template size_t smatrix<int>::col_inrec(size_t key, size_t row);

    //===============================================================================================================

    template <typename T>
    size_t smatrix<T>::row_inrec(size_t key, size_t col)
    {
#ifdef UTEST // use this only when debugging and testing
        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::transpose()");

        if (compact)
            throw std::string("The row_inrec(size_t key, size_t col) should be called only on rectangular matrix!");
#endif
        return (size_t) (key - col) / numCol;
    }

    template size_t smatrix<float>::row_inrec(size_t key, size_t col);
    template size_t smatrix<double>::row_inrec(size_t key, size_t col);
    template size_t smatrix<int>::row_inrec(size_t key, size_t col);

    //===============================================================================================================

    template <typename T>
    size_t smatrix<T>::row_inrec(size_t key)
    {
#ifdef UTEST // use this only when debugging and testing
        if (ondisk)
            throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::transpose()");

        if (compact)
            throw std::string("The row_inrec(size_t key, size_t col) should be called only on rectangular matrix!");
#endif
        return (size_t)( key / numCol );
    }

    template size_t smatrix<float>::row_inrec(size_t key);
    template size_t smatrix<double>::row_inrec(size_t key);
    template size_t smatrix<int>::row_inrec(size_t key);

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::print(std::string whiichMatrix)
    {
        try
        {
            /*
                Prints part of a matrix into a LOG file.

                Return value: none.
            */

            if (ondisk)
                throw std::string("Matrix is empty. Use fread() to relocate data to memory. smatrix<T>::print(std::string)");

            int integer;
            size_t linteger;
            const std::type_info &ti1 = typeid(integer);
            const std::type_info &ti2 = typeid(linteger);
            //const std::type_info &ti3 = typeid(A[0]);
            const std::type_info &ti3 = typeid(this->A.begin()->second);

            bool isInt = false;

            if (ti3 == ti1 || ti3 == ti2)
                isInt = true;

            FILE *dbgFile;
            dbgFile = fopen(debug_file.c_str(), "a");

            size_t maxRows = 20;

            fprintf(dbgFile, "%s", whiichMatrix.c_str());

            if (rectangular)
            {
                fprintf(dbgFile, "%s%s", ", Rectangular matrix, of type ", typeid(A[0]).name());
                fprintf(dbgFile, "\n\n");
                for (size_t i = 0; i < _min(maxRows, numRow); i++)
                {
                    for (size_t j = 0; j < _min(maxRows, numCol); j++)
                    {
                        if ( nonzero(i,j) )
                        {
                            if (isInt)
                                fprintf(dbgFile, "%12d", (int)A[i * numCol + j]);
                            else
                                fprintf(dbgFile, "%12.5G", (double)A[i * numCol + j]);
                        }
                        else
                        {
                            if (isInt)
                                fprintf(dbgFile, "%12d", 0);
                            else
                                fprintf(dbgFile, "%12.5G", 0.0);
                        }
                    }
                    fprintf(dbgFile, "\n");
                }
            }
            else if (symetric)
            {
                fprintf(dbgFile, "%s%s", ", symetric matrix, of type ", typeid(A[0]).name());
                fprintf(dbgFile, "\n\n");
                for (size_t i = 0; i < _min(maxRows, numRow); i++)
                {
                    for (size_t j = 0; j <= i; j++)
                    {
                        if ( nonzero(i,j) )
                        {
                            if (isInt)
                                fprintf(dbgFile, "%12d", (int)A[i * (i + 1) / 2 + j]);
                            else
                                fprintf(dbgFile, "%12.5G", (double)A[i * (i + 1) / 2 + j]);
                        }
                        else
                        {
                            if (isInt)
                                fprintf(dbgFile, "%12d", 0);
                            else
                                fprintf(dbgFile, "%12.5G", 0.0);
                        }
                    }
                    fprintf(dbgFile, "\n");
                }
            }
            fprintf(dbgFile, "\n");
            fprintf(dbgFile, "\n");

            fclose(dbgFile);
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::print(std::string)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::print(std::string)" << '\n';
        }                 
    }

    template void smatrix<float>::print(std::string whiichMatrix);
    template void smatrix<double>::print(std::string whiichMatrix);
    template void smatrix<int>::print(std::string whiichMatrix);

    //===============================================================================================================

    template <typename T>
    void smatrix<T>::printf(std::string whiichMatrix, bool append)
    {
        try
        {
            /*
                Prints all the matrix into the 'whiichMatrix' file.

                Return value: none.
            */

            if (ondisk)
                throw std::string("Matrix is empty. Use fread() to relocate data to memory. matrix<T>::printf(std::string)");

            int integer;
            size_t linteger;
            const std::type_info &ti1 = typeid(integer);
            const std::type_info &ti2 = typeid(linteger);
            //const std::type_info &ti3 = typeid(A[0]);
            const std::type_info &ti3 = typeid(this->A.begin()->second);
            bool isInt = false;

            if (ti3 == ti1 || ti3 == ti2)
                isInt = true;

            FILE *dbgFile;
            if (append)
                dbgFile = fopen(whiichMatrix.c_str(), "a");
            else
                dbgFile = fopen(whiichMatrix.c_str(), "w");

            if (dbgFile == NULL)
                throw std::string("There is problem with opening file in void matrix<T>::printf(std::string)!");

            if (rectangular)
            {
                for (size_t i = 0; i < numRow; i++)
                {
                    for (size_t j = 0; j < numCol; j++)
                    {
                        if ( nonzero(i,j) )
                        {
                            if (isInt)
                                fprintf(dbgFile, "%12d", (int)A[i * numCol + j]);
                            else
                                fprintf(dbgFile, "%12.5G", (double)A[i * numCol + j]);
                        }
                        else
                        {
                            if (isInt)
                                fprintf(dbgFile, "%12d", 0);
                            else
                                fprintf(dbgFile, "%12.5G", 0.0);
                        }
                    }
                    fprintf(dbgFile, "\n");
                }
            }
            else if (symetric)
            {
                for (size_t i = 0; i < numRow; i++)
                {
                    for (size_t j = 0; j <= i; j++)
                    {
                        if ( nonzero(i,j) )
                        {
                            if (isInt)
                                fprintf(dbgFile, "%12d", (int)A[i * (i + 1) / 2 + j]);
                            else
                                fprintf(dbgFile, "%12.5G", (double)A[i * (i + 1) / 2 + j]);
                        }
                        else
                        {
                            if (isInt)
                                fprintf(dbgFile, "%12d", 0);
                            else
                                fprintf(dbgFile, "%12.5G", 0.0);
                        }
                    }
                    fprintf(dbgFile, "\n");
                }
            }

            fclose(dbgFile);
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in smatrix<T>::printf(std::string, bool)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch(...)
        {
            std::cerr << "Exception in smatrix<T>::printf(std::string, bool)" << '\n';
        }        
    }

    template void smatrix<float>::printf(std::string whiichMatrix, bool append);
    template void smatrix<double>::printf(std::string whiichMatrix, bool append);
    template void smatrix<int>::printf(std::string whiichMatrix, bool append);

    //===============================================================================================================
}