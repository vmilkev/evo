#ifndef compact_storage_hpp__
#define compact_storage_hpp__

/*#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <random>*/

#include "dense_matrix.hpp"
#include "sparse_matrix.hpp"

namespace evolm
{
    template <typename T>
    class compact_storage
    {
    private:

        std::vector<T> vals;
        std::vector<size_t> keys;
        std::vector<std::int64_t> row_first_key; // optional container (on demand) holding indexes of a first kay in each row
        
        size_t nRows;
        size_t nCols;
        bool symmetric;
        bool ondisk;
        std::string binFilename;
        std::fstream fA;

        const size_t data_type = 3003;

        double sparsity_threshold = 0.7;

        size_t col_insym(size_t key, size_t row); /* returns col of symmetric matrix knowing the index (key) and the specific row */
        size_t row_insym(size_t key, size_t col); /* ... the same but for row */
        size_t row_insym(size_t key);
        size_t col_inrec(size_t key, size_t row); /* ... the same as previous but for rectangular matrix */
        size_t row_inrec(size_t key, size_t col);
        size_t row_inrec(size_t key);
        size_t key_insym(size_t row, size_t col); /* position (index, key) of an value in a symmetric matrix based on the specific raw and col */
        size_t key_inrec(size_t row, size_t col); /* position (index, key) of an value in a rectangular matrix based on the specific raw and col */

        void frename();

        void convert_to_sparse();
        void convert_to_dense();
        void convert_to_sparse(smatrix<T> &out);
        void convert_to_dense(matrix<T> &out);

        void sort_vectors();

    public:

        compact_storage(); /* default empty storage */
        compact_storage(size_t nrow, size_t ncol); /* empty rectangular storage */
        compact_storage(size_t nrow); /* empty symmetric (lower-triangular) storage */
        compact_storage(std::vector<T> &value, std::vector<size_t> &row, std::vector<size_t> &col, size_t n_rows, size_t n_cols); /* filled sparse rectangular storage */
        compact_storage(std::vector<T> &value, std::vector<size_t> &row, std::vector<size_t> &col, size_t n_rows); /* filled sparse symmetric storage */
        compact_storage(std::vector<T> &value, size_t n_rows, size_t n_cols); /* filled dense rectangular storage */
        compact_storage(std::vector<T> &value, size_t n_rows); /* filled dense symmetric storage */
        compact_storage(const compact_storage &obj); /* copy constructor */
        ~compact_storage();

        compact_storage &operator=(const compact_storage &rhs);           /* Overloaded assignment '=' operator. */

        bool is_sparse();
        bool empty();
        void clear();
        void remove_data();
        void fclear();
        void fclear(const std::string &fname);
        size_t size();
        size_t ncols();
        size_t nrows();
        size_t max_key();
        double sparsity();
        void set_sparsity_threshold(double threshold);

        void append(T value, size_t irow, size_t icol); /* append single value */
        void append(std::vector<T> &value, std::vector<size_t> &row, std::vector<size_t> &col); /* append (sparse data), even if the storage is not empty */
        void append(std::vector<T> &value); /* append (dense data), even if the storage is not empty */
        void append_with_keys(std::vector<T> &value, std::vector<size_t> &key); /* append values and keys (without rows & cols), even if the storage is not empty */
        
        void fwrite();
        void fread();
        void fwrite(const std::string &fname);
        void fread(const std::string &fname);
        void fwrite(std::fstream &external_fA);
        void fread(std::fstream &external_fA);
        
        void optimize();
        void transpose();

        void resize();
        void resize(size_t row);
        void resize(size_t row, size_t col);

        // OUT interface
        void to_sparse(smatrix<T> &out);
        void to_sparse(std::vector<T> &vals_out, std::vector<size_t> &keys_out);
        void to_dense(matrix<T> &out);
        void to_dense(std::vector<T> &vals_out);
        
        // Methods for getting a sub-set of data ...
        void to_sparse(smatrix<T> &out, size_t *row_range, size_t *col_range);
        void to_dense(matrix<T> &out, size_t *row_range, size_t *col_range);
        void to_dense(T **out, size_t *row_range, size_t *col_range);
        void cast_to_fsparse(smatrix<float> &out, size_t *row_range, size_t *col_range);
        void cast_to_fdense(matrix<float> &out, size_t *row_range, size_t *col_range);
        void cast_to_fdense(float **out, size_t *row_range, size_t *col_range);

        void to_sparse(std::vector<T> &vals_out, std::vector<size_t> &keys_out, size_t *row_range, size_t *col_range);
        void cast_to_fsparse(std::vector<float> &vals_out, std::vector<size_t> &keys_out, size_t *row_range, size_t *col_range);

        void row_dot_float_vect(float **in_vect, size_t irow, size_t icol, float &out_res);
        void row_dot_vect(T **in_vect, size_t irow, size_t icol, T &out_res);
        void add_row_to_vector(std::vector<T> &vect_to_add, size_t vect_first_index, float variance, size_t irow, size_t *icol);

        T value_at(size_t irow, size_t icol);
        void vect_dot_vect(std::vector<double> &in_vect, double &out_res);

        // Type converters
        void to_int(compact_storage<int> &out);
        void to_float(compact_storage<float> &out);
        void to_double(compact_storage<double> &out);

        void make_rows_list();
        void clear_rows_list();
    };

    //===========================================================================================    
    template <typename T>
    compact_storage<T>::compact_storage()
    {
        nRows = 0;
        nCols = 0;
        symmetric = false;
        ondisk = false;

        std::random_device rd;
        srand(rd());
        int iNum = rand() % 100000;

        binFilename = "storage_" + std::to_string(iNum);
    }
    //===========================================================================================    
    template <typename T>
    compact_storage<T>::compact_storage(size_t nrow, size_t ncol)
    {
        nRows = nrow;
        nCols = ncol;
        symmetric = false;
        ondisk = false;

        std::random_device rd;
        srand(rd());
        int iNum = rand() % 100000;

        binFilename = "storage_" + std::to_string(iNum);
    }
    //===========================================================================================    
    template <typename T>
    compact_storage<T>::compact_storage(size_t nrow)
    {
        nRows = nrow;
        nCols = nrow;
        symmetric = true;
        ondisk = false;

        std::random_device rd;
        srand(rd());
        int iNum = rand() % 100000;

        binFilename = "storage_" + std::to_string(iNum);
    }
    //===========================================================================================
    template <typename T>
    compact_storage<T>::compact_storage(std::vector<T> &value, std::vector<size_t> &row, std::vector<size_t> &col, size_t n_rows, size_t n_cols)
    {
        nRows = n_rows;
        nCols = n_cols;
        symmetric = false;
        ondisk = false;

        std::random_device rd;
        srand(rd());
        int iNum = rand() % 100000;

        binFilename = "storage_" + std::to_string(iNum);

        if ( (value.size() != row.size()) || (value.size() != col.size()) )
            throw std::string("compact_storage<T>::compact_storage(std::vector<T> &, std::vector<size_t> &, std::vector<size_t> &, size_t, size_t): Provided arrays are not of the equal sizes!");

        size_t max_row = *max_element(std::begin(row), std::end(row));
        size_t max_col = *max_element(std::begin(col), std::end(col));

        if ( (max_row > nRows-1) || (max_col > nCols-1) )
            throw std::string("compact_storage<T>::compact_storage(std::vector<T> &, std::vector<size_t> &, std::vector<size_t> &, size_t, size_t): Provided row (column) is greater than allowed storage dimension!");

        if ( (n_rows * n_cols) < value.size() )
            throw std::string("compact_storage<T>::compact_storage(std::vector<T> &, std::vector<size_t> &, std::vector<size_t> &, size_t, size_t): Too many elements in the provided array for a givem storage size!");

        vals = value;

        size_t ikey = 0;

        for ( size_t i = 0; i < row.size(); i++ )
        {
            ikey = key_inrec( row[i], col[i] );
            keys.push_back(ikey);
        }

        optimize();
    }
    //===========================================================================================
    template <typename T>
    compact_storage<T>::compact_storage(std::vector<T> &value, std::vector<size_t> &row, std::vector<size_t> &col, size_t n_rows)
    {
        nRows = n_rows;
        nCols = n_rows;
        symmetric = true;
        ondisk = false;

        std::random_device rd;
        srand(rd());
        int iNum = rand() % 100000;

        binFilename = "storage_" + std::to_string(iNum);

        if ( (value.size() != row.size()) || (value.size() != col.size()) )
            throw std::string("compact_storage<T>::compact_storage(std::vector<T> &, std::vector<size_t> &, std::vector<size_t> &, size_t): Provided arrays are not of the equal sizes!");

        size_t max_row = *max_element(std::begin(row), std::end(row));
        size_t max_col = *max_element(std::begin(col), std::end(col));

        if ( (max_row > nRows-1) || (max_col > nCols-1) )
            throw std::string("compact_storage<T>::compact_storage(std::vector<T> &, std::vector<size_t> &, std::vector<size_t> &, size_t): Provided row (column) is greater than allowed storage dimension!");

        if ( (n_rows*n_rows+n_rows)/2 < value.size() )
            throw std::string("compact_storage<T>::compact_storage(std::vector<T> &, std::vector<size_t> &, std::vector<size_t> &, size_t): Too many elements in the provided array for a givem storage size!");

        vals = value;

        size_t ikey = 0;
        size_t irow = 0;
        size_t icol = 0;

        for ( size_t i = 0; i < row.size(); i++ )
        {
            irow = row[i];
            icol = col[i];

            ikey = key_insym(irow,icol);

            if (icol>irow)
                ikey = key_insym(icol,irow);
            
            keys.push_back(ikey);
        }

        optimize();
    }
    //===========================================================================================
    template <typename T>
    compact_storage<T>::compact_storage(std::vector<T> &value, size_t n_rows, size_t n_cols)
    {
        nRows = n_rows;
        nCols = n_cols;
        symmetric = false;
        ondisk = false;

        std::random_device rd;
        srand(rd());
        int iNum = rand() % 100000;

        binFilename = "storage_" + std::to_string(iNum);
        
        if ( (n_rows*n_cols) < value.size() )
            throw std::string("compact_storage<T>::compact_storage(std::vector<T> &, size_t, size_t): Too many elements in the provided array for a givem storage size!");

        if ( value.size() != (n_rows*n_cols) )
            throw std::string("compact_storage<T>::compact_storage(std::vector<T> &, size_t, size_t): The provided data expected to be in a complete dense format, but some records are missing!");

        vals = value;

        optimize();
    }
    //===========================================================================================
    template <typename T>
    compact_storage<T>::compact_storage(std::vector<T> &value, size_t n_rows)
    {
        nRows = n_rows;
        nCols = n_rows;
        symmetric = true;
        ondisk = false;

        std::random_device rd;
        srand(rd());
        int iNum = rand() % 100000;

        binFilename = "storage_" + std::to_string(iNum);
        
        if ( (n_rows*n_rows+n_rows)/2 < value.size() )
            throw std::string("compact_storage<T>::compact_storage(std::vector<T> &, size_t): Too many elements in the provided array for a givem storage size!");

        if ( value.size() != (n_rows*n_rows+n_rows)/2 )
            throw std::string("compact_storage<T>::compact_storage(std::vector<T> &, size_t, size_t): The provided data expected to be in a complete dense format, but some records are missing!");

        vals = value;

        optimize();
    }
    //===========================================================================================    
    template <typename T>
    compact_storage<T>::compact_storage(const compact_storage &obj)
    {
        symmetric = obj.symmetric;
        ondisk = obj.ondisk;
        nRows = obj.nRows;
        nCols = obj.nCols;
        binFilename = obj.binFilename;
        vals = obj.vals;
        keys = obj.keys;
    }
    //===========================================================================================    
    template <typename T>
    compact_storage<T> &compact_storage<T>::operator = (const compact_storage<T> &rhs)
    {
        compact_storage<T> obj(rhs);

        std::swap(symmetric, obj.symmetric);
        std::swap(ondisk, obj.ondisk);
        std::swap(nRows, obj.nRows);
        std::swap(nCols, obj.nCols);

        if ( ondisk ) // exchange file names only if rhs is not on memory!
            std::swap(binFilename, obj.binFilename);

        std::swap(vals, obj.vals);
        std::swap(keys, obj.keys);

        return *this;
    }
    //=========================================================================================== 
    template <typename T>
    compact_storage<T>::~compact_storage()
    {
        //fclear(); // check this !!!
        clear();
    }
    //===========================================================================================
    template <typename T>
    bool compact_storage<T>::is_sparse()
    {
        if ( keys.empty() && !vals.empty() )
            return false;
        return true;
    }
    //===========================================================================================
    template <typename T>
    size_t compact_storage<T>::max_key()
    {
        if ( symmetric )
            return ( (nRows * nRows + nRows) / 2.0 ) - 1.0;
        return nRows * nCols - 1;
    }
    //===========================================================================================
    template <typename T>
    bool compact_storage<T>::empty()
    {
        if ( keys.empty() && vals.empty() )
            return true;
        return false;
    }
    //===========================================================================================
    template <typename T>
    size_t compact_storage<T>::size()
    {
        return vals.size();
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::clear()
    {
        vals.clear();
        vals.shrink_to_fit();
        keys.clear();
        keys.shrink_to_fit();
        nRows = 0;
        nCols = 0;
        symmetric = false;
        ondisk = false;
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::remove_data()
    {
        vals.clear();
        vals.shrink_to_fit();
        keys.clear();
        keys.shrink_to_fit();
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::append(T value, size_t irow, size_t icol)
    {
        if ( (irow > nRows-1) || (icol > nCols-1) )
            throw std::string("compact_storage<T>::append(T, size_t, size_t): Provided row (column) is greater than allowed storage dimension!");
        
        vals.push_back(value);

        size_t ikey = 0;

        if (symmetric)
        {
            ikey = key_insym(irow,icol);
            if (icol>irow)
                ikey = key_insym(icol,irow);
        }
        else
            ikey = key_inrec(irow,icol);
        
        keys.push_back(ikey);
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::append(std::vector<T> &value, std::vector<size_t> &row, std::vector<size_t> &col)
    {
        if ( (value.size() != row.size()) || (value.size() != col.size()) )
            throw std::string("compact_storage<T>::append(std::vector<T> &, std::vector<size_t> &, std::vector<size_t> &): Provided arrays are not of the equal sizes!");

        size_t max_row = *max_element(std::begin(row), std::end(row));
        size_t max_col = *max_element(std::begin(col), std::end(col));

        if ( (max_row > nRows-1) || (max_col > nCols-1) )
            throw std::string("compact_storage<T>::append(std::vector<T> &, std::vector<size_t> &, std::vector<size_t> &): Provided row (column) is greater than allowed storage dimension!");

        if ( (max_key() + 1) < ( value.size() + vals.size() ) )
            throw std::string("compact_storage<T>::append(std::vector<T> &, std::vector<size_t> &, std::vector<size_t> &): Too many elements in the provided array for a givem storage size!");

        if (empty())
            vals = value;
        else
        {
            for ( size_t i = 0; i < value.size(); i++ )
                vals.push_back( value[i] );
        }

        size_t ikey = 0;
        size_t irow = 0;
        size_t icol = 0;

        if (symmetric)
        {
            for ( size_t i = 0; i < row.size(); i++ )
            {
                irow = row[i];
                icol = col[i];

                ikey = key_insym(irow,icol);

                if (icol>irow)
                    ikey = key_insym(icol,irow);
                
                keys.push_back(ikey);
            }            
        }
        else
        {
            for ( size_t i = 0; i < row.size(); i++ )
            {
                ikey = key_inrec( row[i], col[i] );
                keys.push_back(ikey);
            }
        }

        optimize();
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::append(std::vector<T> &value)
    {
        if ( !empty() && is_sparse() )
            throw std::string("compact_storage<T>::append(std::vector<T> &): Trying to append a dense data to a sparse storage!");
        
        if ( (max_key() + 1) < ( value.size() + vals.size() ) )
            throw std::string("compact_storage<T>::append(std::vector<T> &): Too many elements in the provided array for a givem storage size!");
        
        if (empty())
            vals = value;
        else
        {
            for ( size_t i = 0; i < value.size(); i++ )
                vals.push_back( value[i] );
        }

        optimize();
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::append_with_keys(std::vector<T> &value, std::vector<size_t> &key)
    {
        if ( !empty() && is_sparse() )
            throw std::string("compact_storage<T>::append(std::vector<T> &): Trying to append a dense data to a sparse storage!");
        
        if ( (max_key() + 1) < ( value.size() + vals.size() ) )
            throw std::string("compact_storage<T>::append(std::vector<T> &): Too many elements in the provided array for a givem storage size!");
        
        if ( !key.empty() && value.size() != key.size() )
            throw std::string("The sizes of vectors (value and key) are not equal!");
        
        if (empty())
        {
            vals = value;
            keys = key;
        }
        else
        {
            for ( size_t i = 0; i < value.size(); i++ )
            {
                vals.push_back( value[i] );
                keys.push_back( key[i] );
            }
        }

        optimize();
    }
    //===========================================================================================
    template <typename T>
    size_t compact_storage<T>::key_insym(size_t row, size_t col)
    {
        return (size_t)( col + row * (row + 1) * 0.5 );
    }
    //===========================================================================================
    template <typename T>
    size_t compact_storage<T>::key_inrec(size_t row, size_t col)
    {
        return (size_t)( col + row * nCols );
    }
    //===========================================================================================
    template <typename T>
    size_t compact_storage<T>::col_insym(size_t key, size_t row)
    {
        return (size_t)( key - row * (row + 1) / 2.0 );
    }
    //===========================================================================================
    template <typename T>
    size_t compact_storage<T>::row_insym(size_t key, size_t col)
    {
        return round( sqrt( 0.25 + 2.0 * (key - col) ) - 0.5 );
    }
    //===========================================================================================
    template <typename T>
    size_t compact_storage<T>::row_insym(size_t key)
    {
        return (size_t)( sqrt( 0.25 + (2 * key + 1.0) ) - 0.5 );
    }
    //===========================================================================================
    template <typename T>
    size_t compact_storage<T>::col_inrec(size_t key, size_t row)
    {
        return (size_t) ( key - row * nCols );
    }
    //===========================================================================================
    template <typename T>
    size_t compact_storage<T>::row_inrec(size_t key, size_t col)
    {
        return (size_t) (key - col) / nCols;
    }
    //===========================================================================================
    template <typename T>
    size_t compact_storage<T>::row_inrec(size_t key)
    {
        return (size_t)( key / nCols );
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::fwrite()
    {
        fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        fA.open(binFilename, fA.binary | fA.trunc | fA.out);

        if (!fA.is_open())
            throw std::string("compact_storage<T>::fwrite(): Error while opening a binary file.");

        size_t vals_size = vals.size();
        size_t keys_size = keys.size();

        fA.write( reinterpret_cast<const char *>( &vals_size ), sizeof(vals_size) ); // 1. writing size of values
        fA.write( reinterpret_cast<const char *>( vals.data() ), vals_size * sizeof(T) ); // 2. writing all values

        fA.write( reinterpret_cast<const char *>( &keys_size ), sizeof(keys_size) ); // 3. writing size of keys
        fA.write( reinterpret_cast<const char *>( keys.data() ), keys_size * sizeof(keys_size) ); // 4. writing all keys

        fA.close();

        ondisk = true;

        vals.clear();
        keys.clear();
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::fwrite(std::fstream &external_fA)
    {
        size_t vals_size = vals.size();
        size_t keys_size = keys.size();

        external_fA.write( reinterpret_cast<const char *>( &vals_size ), sizeof(vals_size) ); // 1. writing size of values
        external_fA.write( reinterpret_cast<const char *>( vals.data() ), vals_size * sizeof(T) ); // 2. writing all values

        external_fA.write( reinterpret_cast<const char *>( &keys_size ), sizeof(keys_size) ); // 3. writing size of keys
        external_fA.write( reinterpret_cast<const char *>( keys.data() ), keys_size * sizeof(keys_size) ); // 4. writing all keys

        /*ondisk = true;

        vals.clear();
        keys.clear();*/
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::fread()
    {
        if (!ondisk)
            throw std::string("The storage status is not ondisk!");
        
        fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        fA.open(binFilename, fA.binary | fA.in);

        if (!fA.is_open())
            throw std::string("compact_storage<T>::fread(): Error while opening a binary file.");

        size_t vals_size = 0;
        size_t keys_size = 0;

        if ( !empty() )
        {
            vals.clear();
            keys.clear();
        }

        fA.read( reinterpret_cast<char *>( &vals_size ), sizeof(vals_size) ); // 1. reading size of values        
        vals.resize(vals_size);
        fA.read( reinterpret_cast<char *>( vals.data() ), vals_size * sizeof(T) ); // 2. reading all values

        fA.read( reinterpret_cast<char *>( &keys_size ), sizeof(keys_size) ); // 3. reading size of keys
        keys.resize(keys_size);
        fA.read( reinterpret_cast<char *>( keys.data() ), keys_size * sizeof(keys_size) ); // 4. reading all keys

        fA.close();

        ondisk = false;
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::fread(std::fstream &external_fA)
    {
        /*if (!ondisk)
            throw std::string("The storage status is not ondisk!");*/
        
        size_t vals_size = 0;
        size_t keys_size = 0;

        if ( !empty() )
        {
            vals.clear();
            keys.clear();
        }

        external_fA.read( reinterpret_cast<char *>( &vals_size ), sizeof(vals_size) ); // 1. reading size of values        
        vals.resize(vals_size);
        external_fA.read( reinterpret_cast<char *>( vals.data() ), vals_size * sizeof(T) ); // 2. reading all values

        external_fA.read( reinterpret_cast<char *>( &keys_size ), sizeof(keys_size) ); // 3. reading size of keys
        keys.resize(keys_size);
        external_fA.read( reinterpret_cast<char *>( keys.data() ), keys_size * sizeof(keys_size) ); // 4. reading all keys

        //ondisk = false;
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::fwrite(const std::string &fname)
    {
        int integer_var;
        float float_var;
        double double_var;

        const std::type_info &type_1 = typeid(integer_var);
        const std::type_info &type_2 = typeid(float_var);
        const std::type_info &type_3 = typeid(double_var);
        const std::type_info &type_T = typeid(vals[0]);
        
        size_t var_type = 0;

        if (type_T == type_1)
            var_type = 1;
        else if (type_T == type_2)
            var_type = 2;
        else if (type_T == type_3)
            var_type = 3;
        else
            throw std::string("Cannot determine the type of values vector.");

        fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        fA.open(fname, fA.binary | fA.trunc | fA.out);

        if (!fA.is_open())
            throw std::string("compact_storage<T>::fwrite(const std::string &): Error while opening a binary file.");

        size_t B[6];

        B[0] = data_type;
        B[1] = nRows;
        B[2] = nCols;
        B[3] = (size_t)symmetric;
        B[4] = (size_t)sizeof(T);
        B[5] = var_type;

        size_t vals_size = vals.size();
        size_t keys_size = keys.size();

        fA.write( reinterpret_cast<const char *>(&B), 6 * sizeof(size_t) ); // 0. writing a storage data

        fA.write( reinterpret_cast<const char *>( &vals_size ), sizeof(vals_size) ); // 1. writing size of values
        fA.write( reinterpret_cast<const char *>( vals.data() ), vals_size * sizeof(T) ); // 2. writing all values

        fA.write( reinterpret_cast<const char *>( &keys_size ), sizeof(keys_size) ); // 3. writing size of keys
        fA.write( reinterpret_cast<const char *>( keys.data() ), keys_size * sizeof(keys_size) ); // 4. writing all keys

        fA.close();
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::fread(const std::string &fname)
    {
        fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        fA.open(fname, fA.binary | fA.in);

        if (!fA.is_open())
            throw std::string("compact_storage<T>::fread(const std::string &): Error while opening a binary file.");

        size_t B[6];

        fA.read( reinterpret_cast<char *>(&B), 6 * sizeof(size_t) ); // 0. reading a storage data

        if ( data_type != B[0] )
            throw std::string("compact_storage<T>::fread(const std::string &): Trying to read a wrong kind of storage from file into the compact_storage format!");

        nRows = B[1];
        nCols = B[2];
        symmetric = (bool)B[3];
        size_t var_inbytes = B[4];

        if ( !empty() )
        {
            vals.clear();
            keys.clear();
        }

        size_t vals_size;
        size_t keys_size;

        fA.read( reinterpret_cast<char *>( &vals_size ), sizeof(vals_size) ); // 1. reading size of values
        vals.resize(vals_size);
        fA.read( reinterpret_cast<char *>( vals.data() ), vals_size * var_inbytes ); // 2. reading all values

        fA.read( reinterpret_cast<char *>( &keys_size ), sizeof(keys_size) ); // 3. reading size of keys
        keys.resize(keys_size);
        fA.read( reinterpret_cast<char *>( keys.data() ), keys_size * sizeof(keys_size) ); // 4. reading all keys

        fA.close();
    }
    //===========================================================================================    
    template <typename T>
    void compact_storage<T>::optimize()
    {
        if ( is_sparse() )
        {
            std::vector<size_t> zeros_position;
            for (size_t i = 0; i < vals.size(); i++) // Check if there are still some zeros in vals
            {
                if ( vals[i] == (T)0 )
                    zeros_position.push_back(i);
            }

            double local_sparsity = sparsity();

            if ( !zeros_position.empty() )
                local_sparsity = 1.0 - (double)(vals.size()-zeros_position.size()) / ( (double)max_key() + 1.0 );

            if ( local_sparsity > sparsity_threshold ) // If TRUE, keep data sparse
            {
                if ( !zeros_position.empty() ) // if some zerros exist, remove them
                    convert_to_sparse();

                return; // OK, completed
            }

            convert_to_dense(); // otherwise, convert to dense representation
        }
        else
        {
            if ( vals.size() != (max_key()+1) ) // Check if the container is in a correct represesntation
                throw std::string("compact_storage<T>::compact_storage(std::vector<T> &, size_t, size_t): The provided data expected to be in a complete dense format, but some records are missing!");

            if ( sparsity() > sparsity_threshold ) // If TRUE, convert to sparse representation
                convert_to_sparse();
        }
    }
    //===========================================================================================
    template <typename T>
    double compact_storage<T>::sparsity()
    {
        if ( !is_sparse() )
        {
            double num_zeros = 0.0;            
            for (auto const &v: vals)
            {
                if (v == (T)0)
                    num_zeros = num_zeros + 1.0;
            }

            return num_zeros / ( (double)max_key() + 1.0 );
        }

        return 1.0 - (double)vals.size() / ( (double)max_key() + 1.0 );
    }
    //===========================================================================================
    template <typename T>
    size_t compact_storage<T>::ncols()
    {
        return nCols;
    }
    //===========================================================================================
    template <typename T>
    size_t compact_storage<T>::nrows()
    {
        return nRows;
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::resize()
    {
        clear();
        frename();
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::resize(size_t row)
    {
        clear();
        frename();

        nRows = nCols = row;
        symmetric = true;
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::resize(size_t row,size_t col)
    {
        clear();
        frename();

        nRows = row;
        nCols = col;
        symmetric = false;
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::transpose()
    {
        if (symmetric)
        {
            return;
        }

        if (is_sparse())
        {
            std::vector<size_t> tmp_keys(keys);

            size_t i_row = 0;
            size_t i_col = 0;

            for (size_t i = 0; i < vals.size(); i++)
            {
                i_row = row_inrec( tmp_keys[i] );
                i_col = col_inrec( tmp_keys[i], i_row );
                keys[i] = i_col * nRows + i_row;
            }

            tmp_keys.clear();

            size_t t_row = nRows;
            nRows = nCols;
            nCols = t_row;

            sort_vectors();
        }
        else
        {
            std::vector<T> tmp_vals(vals);

            size_t i_row = 0;
            size_t i_col = 0;

            for (size_t i = 0; i < vals.size(); i++)
            {
                i_row = row_inrec(i);
                i_col = col_inrec(i, i_row);
                vals[i_col * nRows + i_row] = tmp_vals[i];
            }

            tmp_vals.clear();

            size_t t_row = nRows;
            nRows = nCols;
            nCols = t_row;
        }
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::fclear()
    {        
        std::ifstream f(binFilename.c_str());
        if (f.good())
            remove(binFilename.c_str());
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::fclear(const std::string &fname)
    {        
        std::ifstream f(fname.c_str());
        if (f.good())
            remove(fname.c_str());
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::frename()
    {        
        // rename file name:
        binFilename.clear();
        std::random_device rd;
        srand(rd());
        int iNum = rand() % 100000;
        binFilename = "smatrix_" + std::to_string(iNum);
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::convert_to_dense(matrix<T> &out)
    {
        if (symmetric)
            out.resize(nRows);
        else
            out.resize(nRows, nCols);
        
        if ( !is_sparse() )
        {
            for (size_t i = 0; i < vals.size(); i++)
                out[i] = vals[i];
        }
        else
        {
            for (size_t i = 0; i < vals.size(); i++)
                out[ keys[i] ] = vals[i];
        }
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::convert_to_dense()
    {
        std::vector<T> second( max_key()+1, (T)0 );
        
        if ( !is_sparse() )
        {
            return;
        }
        else
        {
            for (size_t i = 0; i < vals.size(); i++)
                second[ keys[i] ] = vals[i];
        }

        vals.clear();
        keys.clear();
        vals.shrink_to_fit();
        keys.shrink_to_fit();

        vals = second;

        second.clear();
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::convert_to_sparse(smatrix<T> &out)
    {
        if (symmetric)
            out.resize(nRows);
        else
            out.resize(nRows, nCols);
        
        if ( !is_sparse() ) // dense
        {
            for (size_t i = 0; i < vals.size(); i++)
            {
                if ( vals[i] != (T)0 )
                    out[i] = vals[i];
            }
        }
        else // sparse
        {
            for (size_t i = 0; i < vals.size(); i++)
                out[ keys[i] ] = vals[i];
        }
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::convert_to_sparse()
    {
        std::vector<T> second;
        std::vector<size_t> third;
        
        if ( !is_sparse() )
        {
            for (size_t i = 0; i < vals.size(); i++)
            {
                if ( vals[i] != (T)0 )
                {
                    second.push_back(vals[i]);
                    third.push_back(i);
                }
            }
        }
        else
        {
            for (size_t i = 0; i < vals.size(); i++)
            {
                if ( vals[i] != (T)0 )
                {
                    second.push_back(vals[i]);
                    third.push_back(keys[i]);
                }
            }
        }

        vals.clear();
        keys.clear();
        vals.shrink_to_fit();
        keys.shrink_to_fit();

        vals = second;
        keys = third;

        second.clear();
        third.clear();
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::to_sparse(smatrix<T> &out, size_t *row_range, size_t *col_range)
    {
        // row_range[0] - first row
        // row_range[1] - last row
        // col_range[0] - first col
        // col_range[1] - last col

        size_t first_key = 0;
        size_t last_key = 0;

        size_t i_col = 0;
        size_t i_row = 0;

        size_t n_rows = row_range[1] - row_range[0] + 1;
        size_t n_cols = col_range[1] - col_range[0] + 1;

        if ( symmetric )
            throw std::string("to_dense(smatrix<T> &, size_t *, size_t *): The method isnot allowed on symmetric storage!");

        out.resize(n_rows, n_cols);

        first_key = key_inrec(row_range[0], col_range[0]);
        last_key = key_inrec(row_range[1], col_range[1]);
        
        if ( is_sparse() )
        {
            std::int64_t first_index = 0;
            size_t last_index = keys.size();

            if ( !row_first_key.empty() )
            {
                std::int64_t index = row_first_key[ row_range[0] ];
                size_t next_index = row_range[0];
                while (index == -1 && next_index <= row_range[1])
                {
                    next_index++;
                    index = row_first_key[ next_index ];
                }
                first_index = index; // very first key in the row row_range[0]
            }

            if (first_index == -1) // empty data
                return;

            for (size_t i = (size_t)first_index; i < last_index; i++)
            {
                if (keys[i] > last_key)
                    return;
                
                if ( keys[i] >= first_key ) // assume sorted!
                {
                    i_row = row_inrec(keys[i]);
                    i_col = col_inrec(keys[i], i_row);

                    if ( i_col >= col_range[0] && i_col <= col_range[1] )
                        out( i_row - row_range[0], i_col - col_range[0] ) = vals[i];
                }
            }
        }
        else
        {
            for (size_t i = first_key; i <= last_key; i++) // assume sorted!
            {   
                if ( vals[i] == (T)0 )
                    continue;
                
                i_row = row_inrec(i);
                i_col = col_inrec(i, i_row);

                if ( i_col >= col_range[0] && i_col <= col_range[1] )
                    out( i_row - row_range[0], i_col - col_range[0] ) = vals[i];
            }
        }
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::cast_to_fsparse(smatrix<float> &out, size_t *row_range, size_t *col_range)
    {
        // row_range[0] - first row
        // row_range[1] - last row
        // col_range[0] - first col
        // col_range[1] - last col

        size_t first_key = 0;
        size_t last_key = 0;

        size_t i_col = 0;
        size_t i_row = 0;

        size_t n_rows = row_range[1] - row_range[0] + 1;
        size_t n_cols = col_range[1] - col_range[0] + 1;

        if ( symmetric )
            throw std::string("to_dense(smatrix<float> &, size_t *, size_t *): The method isnot allowed on symmetric storage!");

        out.resize(n_rows, n_cols);

        first_key = key_inrec(row_range[0], col_range[0]);
        last_key = key_inrec(row_range[1], col_range[1]);
        
        if ( is_sparse() )
        {
            std::int64_t first_index = 0;
            size_t last_index = keys.size();

            if ( !row_first_key.empty() )
            {
                std::int64_t index = row_first_key[ row_range[0] ];
                size_t next_index = row_range[0];
                while (index == -1 && next_index <= row_range[1])
                {
                    next_index++;
                    index = row_first_key[ next_index ];
                }
                first_index = index; // very first key in the row row_range[0]
            }

            if (first_index == -1) // empty data
                return;

            for (size_t i = (size_t)first_index; i < last_index; i++)
            {
                if (keys[i] > last_key)
                    return;
                
                if ( keys[i] >= first_key ) // assume sorted!
                {
                    i_row = row_inrec(keys[i]);
                    i_col = col_inrec(keys[i], i_row);

                    if ( i_col >= col_range[0] && i_col <= col_range[1] )
                        out( i_row - row_range[0], i_col - col_range[0] ) = static_cast<float>(vals[i]);
                }
            }
        }
        else
        {
            for (size_t i = first_key; i <= last_key; i++) // assume sorted!
            {   
                if ( vals[i] == (T)0 )
                    continue;
                
                i_row = row_inrec(i);
                i_col = col_inrec(i, i_row);

                if ( i_col >= col_range[0] && i_col <= col_range[1] )
                    out( i_row - row_range[0], i_col - col_range[0] ) = static_cast<float>(vals[i]);
            }
        }
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::to_dense(matrix<T> &out, size_t *row_range, size_t *col_range)
    {
        // row_range[0] - first row
        // row_range[1] - last row
        // col_range[0] - first col
        // col_range[1] - last col

        size_t first_key = 0;
        size_t last_key = 0;

        size_t i_col = 0;
        size_t i_row = 0;

        size_t n_rows = row_range[1] - row_range[0] + 1;
        size_t n_cols = col_range[1] - col_range[0] + 1;

        if ( symmetric )
            throw std::string("to_dense(matrix<T> &, size_t *, size_t *): The method isnot allowed on symmetric storage!");

        out.resize(n_rows, n_cols);

        first_key = key_inrec(row_range[0], col_range[0]);
        last_key = key_inrec(row_range[1], col_range[1]);
        
        if ( is_sparse() )
        {
            std::int64_t first_index = 0;
            size_t last_index = keys.size();

            if ( !row_first_key.empty() )
            {
                std::int64_t index = row_first_key[ row_range[0] ];
                size_t next_index = row_range[0];
                while (index == -1 && next_index <= row_range[1])
                {
                    next_index++;
                    index = row_first_key[ next_index ];
                }
                first_index = index; // very first key in the row row_range[0]
            }

            if (first_index == -1) // empty data
                return;

            for (size_t i = (size_t)first_index; i < last_index; i++)
            {
                if (keys[i] > last_key)
                    return;
                
                if ( keys[i] >= first_key ) // assume sorted!
                {
                    i_row = row_inrec(keys[i]);
                    i_col = col_inrec(keys[i], i_row);

                    if ( i_col >= col_range[0] && i_col <= col_range[1] )
                        out( i_row - row_range[0], i_col - col_range[0] ) = vals[i];
                }
            }
        }
        else
        {
            for (size_t i = first_key; i <= last_key; i++) // assume sorted!
            {   
                if ( vals[i] == (T)0 )
                    continue;
                
                i_row = row_inrec(i);
                i_col = col_inrec(i, i_row);

                if ( i_col >= col_range[0] && i_col <= col_range[1] )
                    out( i_row - row_range[0], i_col - col_range[0] ) = vals[i];
            }
        }
    }
        //===========================================================================================
    template <typename T>
    void compact_storage<T>::cast_to_fdense(matrix<float> &out, size_t *row_range, size_t *col_range)
    {
        // row_range[0] - first row
        // row_range[1] - last row
        // col_range[0] - first col
        // col_range[1] - last col

        size_t first_key = 0;
        size_t last_key = 0;

        size_t i_col = 0;
        size_t i_row = 0;

        size_t n_rows = row_range[1] - row_range[0] + 1;
        size_t n_cols = col_range[1] - col_range[0] + 1;

        if ( symmetric )
            throw std::string("fcast_to_dense(matrix<float> &, size_t *, size_t *): The method isnot allowed on symmetric storage!");

        out.resize(n_rows, n_cols);

        first_key = key_inrec(row_range[0], col_range[0]);
        last_key = key_inrec(row_range[1], col_range[1]);
        
        if ( is_sparse() )
        {
            std::int64_t first_index = 0;
            size_t last_index = keys.size();

            if ( !row_first_key.empty() )
            {
                std::int64_t index = row_first_key[ row_range[0] ];
                size_t next_index = row_range[0];
                while (index == -1 && next_index <= row_range[1])
                {
                    next_index++;
                    index = row_first_key[ next_index ];
                }
                first_index = index; // very first key in the row row_range[0]
            }

            if (first_index == -1) // empty data
                return;

            for (size_t i = (size_t)first_index; i < last_index; i++)
            {
                if (keys[i] > last_key)
                    return;
                
                if ( keys[i] >= first_key ) // assume sorted!
                {
                    i_row = row_inrec(keys[i]);
                    i_col = col_inrec(keys[i], i_row);

                    if ( i_col >= col_range[0] && i_col <= col_range[1] )
                        out( i_row - row_range[0], i_col - col_range[0] ) = static_cast<float>(vals[i]);
                }
            }
        }
        else
        {
            for (size_t i = first_key; i <= last_key; i++) // assume sorted!
            {   
                if ( vals[i] == (T)0 )
                    continue;
                
                i_row = row_inrec(i);
                i_col = col_inrec(i, i_row);

                if ( i_col >= col_range[0] && i_col <= col_range[1] )
                    out( i_row - row_range[0], i_col - col_range[0] ) = static_cast<float>(vals[i]);
            }
        }
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::to_dense(T **out, size_t *row_range, size_t *col_range)
    {
        // row_range[0] - first row
        // row_range[1] - last row
        // col_range[0] - first col
        // col_range[1] - last col

        size_t first_key = 0;
        size_t last_key = 0;

        size_t i_col = 0;
        size_t i_row = 0;

        //size_t n_rows = row_range[1] - row_range[0] + 1;
        //size_t n_cols = col_range[1] - col_range[0] + 1;

        if ( symmetric )
            throw std::string("to_dense(T **, size_t *, size_t *): The method isnot allowed on symmetric storage!");

        first_key = key_inrec(row_range[0], col_range[0]);
        last_key = key_inrec(row_range[1], col_range[1]);
        
        if ( is_sparse() )
        {
            for (size_t i = row_range[0]; i <= row_range[1]; i++)
            {
                for (size_t j = col_range[0]; j <= col_range[1]; j++)
                    out[i - row_range[0]][j - col_range[0]] = (T)0;
            }

            std::int64_t first_index = 0;
            size_t last_index = keys.size();

            if ( !row_first_key.empty() )
            {
                std::int64_t index = row_first_key[ row_range[0] ];
                size_t next_index = row_range[0];
                while (index == -1 && next_index <= row_range[1])
                {
                    next_index++;
                    index = row_first_key[ next_index ];
                }
                first_index = index; // very first key in the row row_range[0]
            }

            if (first_index == -1) // empty data
                return;

            for (size_t i = (size_t)first_index; i < last_index; i++)
            {
                if (keys[i] > last_key)
                    return;
                
                if ( keys[i] >= first_key ) // assume sorted!
                {
                    i_row = row_inrec(keys[i]);
                    i_col = col_inrec(keys[i], i_row);

                    if ( i_col >= col_range[0] && i_col <= col_range[1] )
                        out[ i_row - row_range[0] ][ i_col - col_range[0] ] = vals[i];
                }
            }
        }
        else
        {
            for (size_t i = first_key; i <= last_key; i++) // assume sorted!
            {                       
                i_row = row_inrec(i);
                i_col = col_inrec(i, i_row);

                if ( i_col >= col_range[0] && i_col <= col_range[1] )
                    out[ i_row - row_range[0] ][ i_col - col_range[0] ] = vals[i];
            }
        }
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::cast_to_fdense(float **out, size_t *row_range, size_t *col_range)
    {
        // row_range[0] - first row
        // row_range[1] - last row
        // col_range[0] - first col
        // col_range[1] - last col

        size_t first_key = 0;
        size_t last_key = 0;

        size_t i_col = 0;
        size_t i_row = 0;

        //size_t n_rows = row_range[1] - row_range[0] + 1;
        //size_t n_cols = col_range[1] - col_range[0] + 1;

        if ( symmetric )
            throw std::string("fcast_to_dense(float **, size_t *, size_t *): The method isnot allowed on symmetric storage!");

        first_key = key_inrec(row_range[0], col_range[0]);
        last_key = key_inrec(row_range[1], col_range[1]);
        
        if ( is_sparse() )
        {
            for (size_t i = row_range[0]; i <= row_range[1]; i++)
            {
                for (size_t j = col_range[0]; j <= col_range[1]; j++)
                    out[i - row_range[0]][j - col_range[0]] = (float)0;
            }

            std::int64_t first_index = 0;
            size_t last_index = keys.size();

            if ( !row_first_key.empty() )
            {
                std::int64_t index = row_first_key[ row_range[0] ];
                size_t next_index = row_range[0];
                while (index == -1 && next_index <= row_range[1])
                {
                    next_index++;
                    index = row_first_key[ next_index ];
                }
                first_index = index; // very first key in the row row_range[0]
            }

            if (first_index == -1) // empty data
                return;

            for (size_t i = (size_t)first_index; i < last_index; i++)
            {
                if (keys[i] > last_key)
                    return;
                
                if ( keys[i] >= first_key ) // assume sorted!
                {
                    i_row = row_inrec(keys[i]);
                    i_col = col_inrec(keys[i], i_row);

                    if ( i_col >= col_range[0] && i_col <= col_range[1] )
                        out[ i_row - row_range[0] ][ i_col - col_range[0] ] = static_cast<float>(vals[i]);
                }
            }
        }
        else
        {
            for (size_t i = first_key; i <= last_key; i++) // assume sorted!
            {                       
                i_row = row_inrec(i);
                i_col = col_inrec(i, i_row);

                if ( i_col >= col_range[0] && i_col <= col_range[1] )
                    out[ i_row - row_range[0] ][ i_col - col_range[0] ] = static_cast<float>(vals[i]);
            }
        }
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::to_sparse(std::vector<T> &vals_out, std::vector<size_t> &keys_out, size_t *row_range, size_t *col_range)
    {
        // row_range[0] - first row
        // row_range[1] - last row
        // col_range[0] - first col
        // col_range[1] - last col

        size_t first_key = 0;
        size_t last_key = 0;

        size_t i_col = 0;
        size_t i_row = 0;

        //size_t n_rows = row_range[1] - row_range[0] + 1;
        size_t n_cols = col_range[1] - col_range[0] + 1;

        size_t new_key = 0;

        if ( symmetric )
            throw std::string("to_sparse(std::vector<T> &, std::vector<size_t> &, size_t *, size_t *): The method isnot allowed on symmetric storage!");

        first_key = key_inrec(row_range[0], col_range[0]);
        last_key = key_inrec(row_range[1], col_range[1]);
        
        if ( is_sparse() )
        {
            std::int64_t first_index = 0;
            size_t last_index = keys.size();

            if ( !row_first_key.empty() )
            {
                std::int64_t index = row_first_key[ row_range[0] ];
                size_t next_index = row_range[0];
                while (index == -1 && next_index <= row_range[1])
                {
                    next_index++;
                    index = row_first_key[ next_index ];
                }
                first_index = index; // very first key in the row row_range[0]
            }

            if (first_index == -1) // empty data
                return;

            for (size_t i = (size_t)first_index; i < last_index; i++)
            {
                if (keys[i] > last_key)
                    return;
                
                if ( keys[i] >= first_key ) // assume sorted!
                {
                    i_row = row_inrec(keys[i]);
                    i_col = col_inrec(keys[i], i_row);

                    if ( i_col >= col_range[0] && i_col <= col_range[1] )
                    {
                        new_key = (i_row - row_range[0]) * n_cols + (i_col - col_range[0]);
                        vals_out.push_back(vals[i]);
                        keys_out.push_back( new_key );
                    }
                }
            }
        }
        else
        {
            for (size_t i = first_key; i <= last_key; i++) // assume sorted!
            {
                if ( vals[i] == (T)0 )
                    continue;

                i_row = row_inrec(i);
                i_col = col_inrec(i, i_row);

                if ( i_col >= col_range[0] && i_col <= col_range[1] )
                {
                    new_key = (i_row - row_range[0]) * n_cols + (i_col - col_range[0]);
                    vals_out.push_back(vals[i]);
                    keys_out.push_back( new_key );
                }
            }
        }
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::cast_to_fsparse(std::vector<float> &vals_out, std::vector<size_t> &keys_out, size_t *row_range, size_t *col_range)
    {
        // row_range[0] - first row
        // row_range[1] - last row
        // col_range[0] - first col
        // col_range[1] - last col

        size_t first_key = 0;
        size_t last_key = 0;

        size_t i_col = 0;
        size_t i_row = 0;

        //size_t n_rows = row_range[1] - row_range[0] + 1;
        size_t n_cols = col_range[1] - col_range[0] + 1;

        size_t new_key = 0;

        if ( symmetric )
            throw std::string("cast_to_fsparse(std::vector<T> &, std::vector<size_t> &, size_t *, size_t *): The method isnot allowed on symmetric storage!");

        first_key = key_inrec(row_range[0], col_range[0]);
        last_key = key_inrec(row_range[1], col_range[1]);
        
        if ( is_sparse() )
        {
            std::int64_t first_index = 0;
            size_t last_index = keys.size();

            if ( !row_first_key.empty() )
            {
                std::int64_t index = row_first_key[ row_range[0] ];
                size_t next_index = row_range[0];
                while (index == -1 && next_index <= row_range[1])
                {
                    next_index++;
                    index = row_first_key[ next_index ];
                }
                first_index = index; // very first key in the row row_range[0]
            }

            if (first_index == -1) // empty data
                return;

            for (size_t i = (size_t)first_index; i < last_index; i++)
            {
                if (keys[i] > last_key)
                    return;
                
                if ( keys[i] >= first_key ) // assume sorted!
                {
                    i_row = row_inrec(keys[i]);
                    i_col = col_inrec(keys[i], i_row);

                    if ( i_col >= col_range[0] && i_col <= col_range[1] )
                    {
                        new_key = (i_row - row_range[0]) * n_cols + (i_col - col_range[0]);
                        vals_out.push_back( static_cast<float>(vals[i]) );
                        keys_out.push_back( new_key );
                    }
                }
            }
        }
        else
        {
            for (size_t i = first_key; i <= last_key; i++) // assume sorted!
            {
                if ( vals[i] == (float)0 )
                    continue;

                i_row = row_inrec(i);
                i_col = col_inrec(i, i_row);

                if ( i_col >= col_range[0] && i_col <= col_range[1] )
                {
                    new_key = (i_row - row_range[0]) * n_cols + (i_col - col_range[0]);
                    vals_out.push_back( static_cast<float>(vals[i]) );
                    keys_out.push_back( new_key );
                }
            }
        }
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::to_dense(matrix<T> &out)
    {
        convert_to_dense(out);
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::to_dense(std::vector<T> &vals_out)
    {
        if ( is_sparse() )
            convert_to_dense();

        vals_out = vals;
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::to_sparse(smatrix<T> &out)
    {
        convert_to_sparse(out);
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::to_sparse(std::vector<T> &vals_out, std::vector<size_t> &keys_out)
    {
        convert_to_sparse();

        vals_out = vals;
        keys_out = keys;
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::set_sparsity_threshold(double threshold)
    {
        sparsity_threshold = threshold;
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::sort_vectors()
    {
        std::vector< std::pair <size_t, T> > pair_vect;

        if ( !is_sparse() )
            return;
        
        for (size_t i = 0; i < vals.size(); i++)
        {
            pair_vect.push_back( std::make_pair(keys[i],vals[i]) );
        }

        std::sort(pair_vect.begin(), pair_vect.end());

        vals.clear();
        keys.clear();
        for (size_t i = 0; i < pair_vect.size(); i++)
        {
            keys.push_back(pair_vect[i].first);
            vals.push_back(pair_vect[i].second);    
        }

        pair_vect.clear();
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::to_int(compact_storage<int> &out)
    {
        std::vector<int> t_vect;

        for (size_t i = 0; i < vals.size(); i++)
            t_vect.push_back( static_cast<int>(vals[i]) );

        if (symmetric)
            out.resize(nRows);
        else
            out.resize(nRows,nCols);
        
        out.append_with_keys(t_vect, keys);
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::to_float(compact_storage<float> &out)
    {
        std::vector<float> t_vect;
        
        for (size_t i = 0; i < vals.size(); i++)
            t_vect.push_back( static_cast<float>(vals[i]) );

        if (symmetric)
            out.resize(nRows);
        else
            out.resize(nRows,nCols);
        
        out.append_with_keys(t_vect, keys);
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::to_double(compact_storage<double> &out)
    {
        std::vector<double> t_vect;
        
        for (size_t i = 0; i < vals.size(); i++)
            t_vect.push_back( static_cast<double>(vals[i]) );

        if (symmetric)
            out.resize(nRows);
        else
            out.resize(nRows,nCols);
        
        out.append_with_keys(t_vect, keys);
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::row_dot_float_vect(float **in_vect, size_t irow, size_t icol, float &out_res)
    {
        out_res = 0.0f;

        size_t first_key = key_inrec(irow, 0);
        size_t last_key = key_inrec(irow, icol);
        
        if( is_sparse() )
        {
            std::int64_t first_index = 0;
            size_t last_index = keys.size();

            if ( !row_first_key.empty() )
            {
                std::int64_t index = row_first_key[ irow ];
                if (index == -1) // empty row
                    return;
                first_index = index; // very first key in the row row_range[0]
            }

            for (size_t i = (size_t)first_index; i < last_index; i++)
            {
                if (keys[i] > last_key)
                    break;
                if ( keys[i] >= first_key )
                    out_res = out_res + in_vect[0][ keys[i]-first_key ] * static_cast<float>(vals[i]);
            }
        }
        else
        {
            for (size_t i = first_key; i <= last_key; i++)
            {
                if ( vals[i] == 0 )
                    continue;
                out_res = out_res + in_vect[0][ i-first_key ] * static_cast<float>(vals[i]);
            }
        }
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::row_dot_vect(T **in_vect, size_t irow, size_t icol, T &out_res)
    {
        out_res = (T)0;

        size_t first_key = key_inrec(irow, 0);
        size_t last_key = key_inrec(irow, icol);
        
        if( is_sparse() )
        {
            std::int64_t first_index = 0;
            size_t last_index = keys.size();

            if ( !row_first_key.empty() )
            {
                std::int64_t index = row_first_key[ irow ];
                if (index == -1) // empty row
                    return;
                first_index = index; // very first key in the row row_range[0]
            }

            for (size_t i = (size_t)first_index; i < last_index; i++)
            {
                if (keys[i] > last_key)
                    break;
                if ( keys[i] >= first_key )
                    out_res = out_res + in_vect[0][ keys[i]-first_key ] * vals[i];
            }
        }
        else
        {
            for (size_t i = first_key; i <= last_key; i++)
            {
                if ( vals[i] == 0 )
                    continue;
                out_res = out_res + in_vect[0][ i-first_key ] * vals[i];
            }
        }
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::vect_dot_vect(std::vector<double> &in_vect, double &out_res)
    {
        out_res = 0.0;

        if( is_sparse() )
        {
            for (size_t i = 0; i < keys.size(); i++)
                out_res = out_res + in_vect[ keys[i] ] * static_cast<double>(vals[i]);
        }
        else
        {
            for (size_t i = 0; i < vals.size(); i++)
                out_res = out_res + in_vect[i] * static_cast<double>(vals[i]);
        }
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::add_row_to_vector(std::vector<T> &vect_to_add, size_t vect_first_index, float variance, size_t irow, size_t *icol)
    {
        size_t first_key = key_inrec(irow, icol[0]);
        size_t last_key = key_inrec(irow, icol[1]);

        if( is_sparse() )
        {
            std::int64_t first_index = 0;
            size_t last_index = keys.size();

            if ( !row_first_key.empty() )
            {
                std::int64_t index = row_first_key[ irow ];
                if (index == -1) // empty row
                    return;
                first_index = index; // very first key in the row row_range[0]
            }

            for (size_t i = (size_t)first_index; i < last_index; i++)
            {
                if (keys[i] > last_key)
                    break;
                if ( keys[i] >= first_key )
                    vect_to_add[ vect_first_index+(keys[i]-first_key) ] = vect_to_add[ vect_first_index+(keys[i]-first_key) ] + vals[i] * variance;
            }
        }
        else
        {
            for (size_t i = first_key; i <= last_key; i++)
            {
                if ( vals[i] == 0 )
                    continue;
                vect_to_add[ vect_first_index+(i-first_key) ] = vect_to_add[ vect_first_index+(i-first_key) ] + vals[i] * variance;
            }
        }
    }
    //===========================================================================================
    template <typename T>
    T compact_storage<T>::value_at(size_t irow, size_t icol)
    {
        size_t first_key = key_inrec(irow, icol);
        T out_value = (T)0;

        if( is_sparse() )
        {
            std::int64_t first_index = 0;
            size_t last_index = keys.size();

            if ( !row_first_key.empty() )
            {
                std::int64_t index = row_first_key[ irow ];
                if (index == -1) // empty row
                    return (T)0;
                first_index = index; // very first key in the row row_range[0]
            }

            for (size_t i = (size_t)first_index; i < last_index; i++)
                if (keys[i] == first_key)
                    out_value = vals[i];
        }
        else
            out_value = vals[first_key];
        
        return out_value;
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::make_rows_list()
    {
        if( !row_first_key.empty() )
            row_first_key.clear();
        
        if ( !is_sparse() )
            return;
        
        row_first_key.resize(nRows, -1);

        size_t i_key = 0;
        size_t next_row = 0;

        size_t i_row = row_inrec( keys[i_key] ); // get the row number associated with a very first key
        row_first_key[i_row] = i_key; // key with index 0 associated with row number irow

        while ( i_key < keys.size() )
        {
            next_row = i_row;

            while ( next_row == i_row )
            {
                i_key++;
                i_row = row_inrec( keys[i_key] );
            }

            if( i_key < keys.size() )
                row_first_key[i_row] = i_key;
        }
    }
    //===========================================================================================
    template <typename T>
    void compact_storage<T>::clear_rows_list()
    {
        row_first_key.clear();
        row_first_key.shrink_to_fit();
    }
    //===========================================================================================
    
} // end of the namespace evolm

#endif // compact_storage_hpp__