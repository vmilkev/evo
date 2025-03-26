#ifndef IOInterface_hpp__
#define IOInterface_hpp__

#include <fstream>
#include <vector>
#include <math.h>
#include <iostream>
#include <map>
#include <cstring>
#include <regex>
#include <sstream>

#include "effects_storage.hpp"
#include <plinkio/plinkio.h>

namespace evolm
{
    class IOInterface
    {
    public:
    
        IOInterface();
        void set_fname(std::string file);                                 // Setting IO file name
        void fgetdata(size_t samples,
                      size_t variants,
                      std::vector<std::vector<int>> &out);                // Reads a binary SNP data stored in the PLINK bed format
        void fgetdata(bool includes_id,
                      std::vector<std::vector<int>> &out);                // Reads ASCII SNP data with extra info (the very first two columns -> outputs in the additional data structure)
        void fgetdata(std::vector<std::vector<int>> &out);                // Reads ASCII SNP data
        
        template <typename T>
        void fgetdata(std::vector<std::vector<T>> &out);                  // Reads general ASCII formated data
        
        void fgetvar(const std::string &var_name,                         // Extract data for a specific variable accessed by the name 'var_name'
                     float miss_constant,
                     std::vector<bool> &missing_vect,
                     effects_storage &out_var,
                     std::vector<std::string> &unique_levels,
                     std::string &messege);         // Reference variable name (observations var name) which used to track missing records
        void fgetvar(const std::string &var_name,
                     std::vector<int> &ref_vals_int,
                     std::vector<std::string> &ref_vals_str,
                     float miss_constant,
                     std::vector<bool> &missing_vect,
                     effects_storage &out_var,
                     std::vector<std::string> &unique_levels,
                     std::string &messege);
        void fget_var_levels(const std::string &var_name,
                             float miss_constant,
                             std::vector<int> &out_int,
                             std::vector<std::string> &out_str,
                             std::string &messege);
        bool is_var_in_header(const std::string &var_name);
        
        void clear();

        void scale_genotypes(std::vector<float> &values, size_t &nrows, size_t &ncols, std::vector<std::string> &snp_names);

        // temporaly in public!
        bool is_plink_file(const std::string &fname);
        
        template <typename T>
        int find_value(std::vector<T> &where,
                       T what);                                           // fins the position of the string 'what' in the vector 'where'
        
        template <typename T>
        void fread_matrix(const std::string &fname, std::vector<T> &vals, std::vector<size_t> &keys, std::vector<std::int64_t> &ids);
        template <typename T>
        void fread_matrix(const std::string &fname, std::vector<T> &vals, std::vector<size_t> &keys, std::vector<std::string> &ids);
        template <typename T>
        void fread_matrix(const std::string &fname, evolm::matrix<T> &matr, std::vector<std::int64_t> &ids);
        template <typename T>
        void fread_matrix(const std::string &fname, evolm::matrix<T> &matr, std::vector<std::string> &ids);
        void fread_matrix_info(const std::string &fname, size_t *info);


    private:
        float missing_constant = 0.0f;
        int missing_int = 999999999;
        //template <typename T>
        //int find_value(std::vector<T> &where,
        //               T what);                                           // fins the position of the string 'what' in the vector 'where'
        void str_parse(std::string &snpStr,
                       std::vector<int> &markers);                        // Parsing a SNP string
        int get_datatype(std::string &str_token, const std::string &var_name);                         // accepts data as a string token and return its type
        int define_vartype(std::vector<int> &types_vect, const std::string &var_name, std::string &messege);                 // accepts types vector and defines a future variable data type
        
        void str_to_float(std::vector<std::string> &data_str,
                         std::vector<float> &fvalues);                    // converts vector of strings to a vector of floating point numbers
        
        void str_to_int(std::vector<std::string> &data_str,
                        std::vector<bool> &missing_pos,
                        std::vector<int> &fvalues);                        // converts vector of strings to a vector of integer numbers
        void cat_to_effect(std::vector<int> &cvalues,
                           std::vector<bool> &missing_vect,
                           compact_storage<int> &ematrix); // get integers vector and returns an effect matrix
        //template <typename T>
        size_t int_to_cat(std::vector<int> &ivalues,
                          std::vector<int> &cat_values,
                          std::vector<std::string> &str_levels);             // get integers vector (categorical data) and returns a number of columns of a corresponding effect matrix, and vector of unique values
        size_t int_to_cat(std::vector<int> &ivalues,
                          std::vector<int> &reference_ivalues,
                          std::vector<int> &cat_values,
                          std::vector<std::string> &str_levels);             // get integers vector (categorical data) and returns a number of columns of a corresponding effect matrix, and vector of unique values
        size_t int_to_cat(std::vector<std::string> &ivalues,
                          std::vector<int> &cat_values,
                          std::vector<std::string> &str_levels);             // get integers vector (categorical data) and returns a number of columns of a corresponding effect matrix, and vector of unique values
        size_t int_to_cat(std::vector<std::string> &ivalues,
                          std::vector<std::string> &reference_ivalues,
                          std::vector<int> &cat_values,
                          std::vector<std::string> &str_levels);             // get integers vector (categorical data) and returns a number of columns of a corresponding effect matrix, and vector of unique values

        template <typename T>
        void get_unique_levels(std::vector<T> &ivalues, std::vector<T> &unique_ivalues);

        std::string io_file;                                              // IO file name
        std::map<size_t, size_t> snp_id;                                  // KEY: the consecutive index; the VALUE: observation ID

        //bool is_plink_file(const std::string &fname);
        void get_m_matrix_plink(evolm::matrix<int> &M, std::vector<std::string> &snp_names);
        void get_m_matrix_ascii(int data_format, evolm::matrix<int> &out, std::vector<std::string> &snp_names);
        void make_zmatrix( evolm::matrix<int> &M, evolm::matrix<float> &Z );

        int detect_data_format_in_snp_txt_file(const std::string &fname);

        bool is_number(const std::string& s);
    };

} // end of namespace evolm

#endif // IOInterface_hpp__
