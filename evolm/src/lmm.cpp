#include "lmm.hpp"

namespace evolm
{
    //===============================================================================================================
    lmm::lmm()
    {
    }
    //===============================================================================================================
    lmm::~lmm()
    {
        clear();
    }  
    //===============================================================================================================
    void lmm::define(const std::string &expression)
    {
        try
        {
            parser.process_expression(expression);
        }
        catch(const std::exception& e)
        {
            std::cerr << "define(const std::string &): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "define(const std::string &): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "define(const std::string &): " << "Unknown exception." << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    void lmm::define_infile(const std::string &fname)
    {
        try
        {
            std::vector<std::string> expressions; // read lines from the fname file and store here
            
            read_model_from_file( fname, expressions ); // getting vector of expressions from a file

            for (size_t i = 0; i < expressions.size(); i++)
                parser.process_expression(expressions[i]);
        }
        catch(const std::exception& e)
        {
            std::cerr << "define_infile(const std::string &): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "define_infile(const std::string &): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "define_infile(const std::string &): " << "Unknown exception." << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    void lmm::solve(const std::string &use_method, const std::string &sol_file)
    {
        try
        {
            sparse_pcg solver;
            model_sparse model;

            set_model(model);
            
            solver.append_model(model);
            solver.solve();
            solver.get_solution(sol_file);

            // need to clean solver and model here !!!
        }
        catch(const std::exception& e)
        {
            std::cerr << "solve(const std::string &, const std::string &): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "solve(const std::string &, const std::string &): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "solve(const std::string &, const std::string &): " << "Unknown exception." << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    void lmm::set_model(model_sparse &model)
    {
        try
        {
            // implementation !
        }
        catch(const std::exception& e)
        {
            std::cerr << "set_model(model_sparse &): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "set_model(model_sparse &): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "set_model(model_sparse &): " << "Unknown exception." << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    void lmm::read_model_from_file(const std::string &fname, std::vector<std::string> &out_expr)
    {
        try
        {
            std::ifstream file(fname);

            if ( !file.is_open() )
                throw std::string( "Failed to open file: " + fname );

            std::string line;
            while ( getline(file, line) )
                out_expr.push_back( line );

            file.close();
        }
        catch(const std::exception& e)
        {
            std::cerr << "read_model_from_file(const std::string &, std::vector<std::string> &): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "read_model_from_file(const std::string &, std::vector<std::string> &): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "read_model_from_file(const std::string &, std::vector<std::string> &): " << "Unknown exception." << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    void lmm::clear()
    {
        try
        {
            parser.clear();

            // ... and other data structures ?
        }
        catch(const std::exception& e)
        {
            std::cerr << "clear(): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "clear(): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "clear(): " << "Unknown exception." << '\n';
            throw;
        }        
    }
    //===============================================================================================================
} // end of namespace evolm