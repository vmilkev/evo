#ifndef sparse_matrix_hpp__
#define sparse_matrix_hpp__

#include <iostream>
#include <fstream>
#include <string>
#include <random>
#include <thread>
#include <algorithm>
#include <map>
#include <unordered_map>

#ifndef _min
#define _min(x, y) (((x) < (y)) ? (x) : (y))
#endif

#ifndef _max
#define _max(x, y) (((x) > (y)) ? (x) : (y))
#endif

namespace evolm
{
    template <typename T>
    class smatrix
    {
    private:

        struct ustorage // !!! DELETE THIS
        {
            std::unordered_map <size_t, T> A;

            ustorage(const smatrix &obj)
            {
                compact = obj.compact;
                numCol = obj.numCol;
                numRow = obj.numRow;
                A.insert(obj.A.begin(), obj.A.end()); //A = obj.A;
            };

            ~ustorage()
            {
                A.clear();
            };

            T get_nonzero(size_t atRow, size_t atCol)
            {
                typename std::unordered_map<size_t, T>::iterator it;
                size_t key = 0;
                if (!compact)
                    key = atRow * numCol + atCol;
                else
                    key = atRow * (atRow + 1) / 2 + atCol;

                it = A.find( key );

                if (it != A.end())
                    return it->second;
                else
                    return (T)0;
            };

            size_t nrows()
            {
                return numRow;
            };

            size_t ncols()
            {
                return numCol;
            };

            void clear()
            {
                A.clear();
            };

            size_t numRow;
            size_t numCol;
            bool compact;
        };

        struct ordstorage
        {
            std::map <size_t, T> A;

            ordstorage(const smatrix &obj)
            {
                compact = obj.compact;
                numCol = obj.numCol;
                numRow = obj.numRow;
                A.insert(obj.A.begin(), obj.A.end()); //A = obj.A;
            };

            ~ordstorage()
            {
                A.clear();
            };

            T get_nonzero(size_t atRow, size_t atCol)
            {
                typename std::map<size_t, T>::iterator it;
                size_t key = 0;
                if (!compact)
                    key = atRow * numCol + atCol;
                else
                    key = atRow * (atRow + 1) / 2 + atCol;

                it = A.find( key );

                if (it != A.end())
                    return it->second;
                else
                    return (T)0;
            };

            size_t col_inrec(size_t key, size_t row)
            {
#ifdef UTEST // use this only when debugging and testing
                if (compact)
                    throw std::string("The row_inrec(size_t key, size_t col) should be called only on rectangular matrix!");
#endif
                return (size_t) ( key - row * numCol );
            };

            size_t row_inrec(size_t key, size_t col)
            {
#ifdef UTEST // use this only when debugging and testing
                if (compact)
                    throw std::string("The row_inrec(size_t key, size_t col) should be called only on rectangular matrix!");
#endif
                return (size_t) (key - col) / numCol;
            };

            size_t row_inrec(size_t key)
            {
#ifdef UTEST // use this only when debugging and testing
                if (compact)
                    throw std::string("The row_inrec(size_t key, size_t col) should be called only on rectangular matrix!");
#endif
                return (size_t)( key / numCol );
            };

            size_t size()
            {
                return A.size();
            };

            size_t nrows()
            {
                return numRow;
            };

            size_t ncols()
            {
                return numCol;
            };

            void clear()
            {
                A.clear();
            };

            size_t numRow;
            size_t numCol;
            bool compact;
        };

        //std::map <size_t, T> A; /* The main matrix container. */
        std::unordered_map <size_t, T> A; /* The main matrix container. */

        size_t numRow;          /* Number of rows in dense matrix. */
        size_t numCol;          /* Number of columns in dense matrix. */
        size_t max_elements;    /* Number of max possible elements in dense matrix. */
        size_t ondisc_elements; /* Actual number of elements in A writen into disk. */
        bool rectangular;       /* TRUE if the matrix is rectangular. */
        bool symetric;          /* TRUE if the matrix is symmetric. */
        bool compact;           /* TRUE if for the symmetric matrix only a lower triangular part is stored. */
        bool ondisk;
        size_t work_load_perthread;
        T zerro_tolerance;
        std::string debug_file;
        std::string binFilename_keys; /* Name of binary file to store A on disck. */
        std::string binFilename_vals; /* Name of binary file to store A on disck. */
        std::fstream fA_keys;
        std::fstream fA_vals;

        size_t col_insym(size_t key, size_t row); /* returns col of symmetric matrix knowing the index (key) and the specific row */
        size_t row_insym(size_t key, size_t col); /* ... the same but for row */
        size_t row_insym(size_t key);
        size_t col_inrec(size_t key, size_t row); /* ... the same as previous but for rectangular matrix */
        size_t row_inrec(size_t key, size_t col);
        size_t row_inrec(size_t key);
        size_t key_insym(size_t row, size_t col); /* position (index, key) of an value in a symmetric matrix based on the specific raw and col */
        size_t key_inrec(size_t row, size_t col); /* position (index, key) of an value in a rectangular matrix based on the specific raw and col */
        
        size_t find_invect(std::vector<size_t> &where, size_t what); /* finds the specific value in a vector */
        size_t find_inrange(std::vector<size_t> &where, size_t what); /* find the position in vector of the values less then what */
        
        // for multithreading:
        void dot_operation(smatrix &lhs, ordstorage &rhs, smatrix &out, std::vector<size_t> &range_vect, size_t thr_id);
        void plus_operation(smatrix &in, smatrix &res, std::vector<size_t> &zero_keys, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);
        void minus_operation(smatrix &in, smatrix &res, std::vector<size_t> &zero_keys, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);
        void transpose_operation(smatrix &in, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);
        void rectosym_operation(smatrix &in, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);
        void symtorec_operation(smatrix &in, smatrix &out, std::vector<size_t> &loads_vect, size_t thr_id);
        void thread_loads_unord(smatrix &in, std::vector<size_t> &out);
        void thread_loads_ord(ordstorage &in, std::vector<size_t> &out);
        void values_inrow(ordstorage &in, std::vector<size_t> &out_elements);
        
    public:
        smatrix(size_t row, size_t col); /* Constructor for rectangular matrix. */
        smatrix(size_t lda);             /* Constructor for square symmetric half-store (compact) matrix. */
        smatrix();                       /* Default constructor: rectangular matrix, no memmory allocated. */
        smatrix(const smatrix &obj);      /* Copy constructor. */
        ~smatrix();                      /* Destructor. */

        T &operator()(size_t atRow, size_t atCol);      /* Access of any element of a matrix, using raw and col. */
        T &operator[](size_t atIndex);      /* Access of any element of a matrix, using index. */
        smatrix &operator=(const smatrix &rhs);           /* Overloaded assignment '=' operator. */
        smatrix operator*(smatrix &rhs);  /* Overloaded '*' operator to multiply two matrix objects. */
        smatrix operator+(smatrix &rhs);
        smatrix operator-(smatrix &rhs);
        //smatrix operator*(const T &val);

        void transpose();                             /* Transpose matrix. */
        void rectosym();                              /* transform matrix from full storage to compact (lower triangular part), assuming the matrix is symmetric */
        void symtorec();                              /* transform matrix from compact storage (lower triangular part, assuming the matrix is symmetric) to full (rectangular) storage */
        void fwrite();                                /* Move matrix to the disk and clear memory. */
        void fread();                                 /* Restore matrix from the disk into the memory. */
        void fclear();
        void clear();
        bool empty();
        void resize(size_t row, size_t col);
        void resize(size_t lda);
        void resize(); // just empty container, rows = cols = 0
        bool nonzero(size_t atRow, size_t atCol); /* Checks if the specific matrix element is exists in the container (has the value different from 0).*/
        T get_nonzero(size_t atRow, size_t atCol); /* Checks if the specific matrix element is exists in the container (has the value different from 0).*/
        bool nonzero(size_t key);
        T get_nonzero(size_t key);
        size_t nrows(); /* Returns the number of rows in the dense matrix. */
        size_t ncols(); /* Returns the number of cols in the dense matrix. */
        size_t size(); /* Returns the actual size of A. */
        size_t max_key(); /* Returns the maximum index value for the very last element in the dense matrix. */
        void print(std::string whiichMatrix); /* Prints part of a matrix into a LOG file. */
        void printf(std::string whiichMatrix,
                    bool append);                     /* Prints all of a matrix into a specified file file name. */
        void set_thread_load( size_t load );

    };
        
} // end of namespace evolm

#endif // sparse_matrix_hpp__