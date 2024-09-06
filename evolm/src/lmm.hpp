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
        void set_model( model_sparse &model ); // uses parser data to build model
        void read_model_from_file( const std::string &fname, std::vector<std::string> &out_expr );
        void clear();
    public:
        lmm();
        ~lmm();
        void define( const std::string &expression );
        void define_infile( const std::string &fname );
        void solve( const std::string &use_method, const std::string &sol_file );
    };    
} // end of namespace evolm

#endif // lmm_hpp__