#ifndef lmm_hpp__
#define lmm_hpp__

#include "model_parser.hpp"
#include "model_sparse.hpp"
#include "sparse_pcg.hpp"

namespace evolm
{
    class lmm
    {
    private:
        model_parser parser; // handler and storage for model's definition and data
        std::vector<std::vector<std::string>> sol_table;
        std::vector<std::vector<std::string>> sol_structure;
        void set_model( model_sparse &model, const std::string &log_file ); // uses parser data to build model
        void read_model_from_file( const std::string &fname, std::vector<std::string> &out_expr );
        void append_solution_table(std::string &var_name, std::vector<std::string> &levels, std::string &trait);
        void append_solution_structure(std::string &var_name, std::vector<std::string> &levels, std::string &trait, size_t record);
        void process_solution(std::vector<double> &sol_vect, const std::string &sol_file);
        std::string to_string_with_precision(const double value, const int n);
        void message_to_log( const std::string &log_file, std::string msg );
    
    public:
        lmm();
        ~lmm();
        void define( const std::string &expression );
        void define_infile( const std::string &fname );
        void solve( const std::string &use_method, int available_memory, int available_cpu, const std::string &log_file, const std::string &sol_file );
        void solve( const std::string &use_method, int available_memory, int available_cpu, const std::string &log_file, const std::string &sol_file, double err_tol );
        void solve( const std::string &use_method, int available_memory, int available_cpu, const std::string &log_file, const std::string &sol_file, double err_tol, size_t max_iter );
        void snp_to_bv( const std::string &fname_zmatr, const std::string &fname_solsnp, const std::string &fname_outres );
        void clear();
    };
} // end of namespace evolm

#endif // lmm_hpp__