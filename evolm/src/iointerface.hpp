#ifndef IOInterface_hpp__
#define IOInterface_hpp__

#include <fstream>
#include <vector>
#include <math.h>
#include <iostream>
#include <map>
#include <cstring>
#include <regex>

#include "effects.hpp"

namespace evolm
{
    class IOInterface
    {
    public:
        void set_fname(std::string file);                                 // Setting IO file name
        void fgetdata(size_t samples,
                      size_t variants,
                      std::vector<std::vector<int>> &out);                // Reads a binary SNP data stored in the PLINK bed format
        void fgetdata(bool includes_id,
                      std::vector<std::vector<int>> &out);                // Reads ASCII SNP data with extra info (the very first two columns -> outputs in the additional data structure)
        void fgetdata(std::vector<std::vector<int>> &out);                // Reads ASCII SNP data
        template <typename T>
        void fgetdata(std::vector<std::vector<T>> &out);                  // Reads general ASCII formated data
        void fgetvar(std::string &var_name,
                     Effects &out_var);                                   // Extract a data for a specific variable accessed by the name 'var_name'

    private:
        template <typename T>
        int find_value(std::vector<T> &where,
                       T what);                                           // fins the position of the string 'what' in the vector 'where'
        void str_parse(std::string &snpStr,
                       std::vector<int> &markers);                        // Parsing a SNP string
        int get_datatype(std::string &str_token);                         // accepts data as a string token and return its type
        int define_vartype(std::vector<int> &types_vect);                 // accepts types vector and defines a future variable data type
        void str_tofloat(std::vector<std::string> &data_str,
                         std::vector<float> &fvalues);                    // converts vector of strings to a vector of floating point numbers
        void str_toint(std::vector<std::string> &data_str,
                       std::vector<int> &fvalues);                        // converts vector of strings to a vector of integer numbers
        void catvect_toeffmatrix(std::vector<int> &cvalues,
                                 std::vector<std::vector<int>> &ematrix); // get integers vector and returns an effect matrix
        template <typename T>
        size_t dim_ofeffmatrix(std::vector<T> &ivalues,
                               std::vector<int> &cat_values);             // get integers vector (categorical data) and returns a number of columns of a corresponding effect matrix, and vector of unique values

        std::string io_file;                                              // IO file name
        std::map<size_t, size_t> snp_id;                                  // KEY: the consecutive index; the VALUE: observation ID
    };

} // end of namespace evolm

#endif // IOInterface_hpp__
