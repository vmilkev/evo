#include "lmm.hpp"

namespace evolm
{
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
        // defining model
        try
        {
            // (1) ------- Residual ------------
            compact_storage<float> R_as_storage;
            std::vector<float> R;
            size_t nrows = 0;

            parser.extra_effects[ parser.corr_matr_index_in_extra_vars[parser.corr_matr_index_in_extra_vars.size()-1][0] ].get(R_as_storage);
            R_as_storage.to_dense(R);
            nrows = R_as_storage.nrows();

            model.append_residual(R, nrows);

            R.clear(); R.shrink_to_fit();
            R_as_storage.clear();

            // (2) --------- Observations ----------
            for (auto  const &m: parser.model_definition)
            {
                compact_storage<float> y_as_storage;
                std::vector<float> y;
                size_t nrows = 0;

                size_t obs = m.first;
                int pos_in_extra = parser.obs_in_extra_storage[obs];

                if ( pos_in_extra == -1 )
                    parser.observations[obs].get(y_as_storage);                    
                else
                    parser.extra_effects[ (size_t)pos_in_extra ].get(y_as_storage);
                
                y_as_storage.to_dense(y);
                nrows = y_as_storage.nrows();

                model.append_observation(y, nrows);

                y.clear(); y.shrink_to_fit();
                y_as_storage.clear();
            }

            // (3) ----------- Effects ----------------
            for (auto  const &m: parser.model_definition)
            {
                compact_storage<float> eff_as_storage;

                std::vector<size_t> effects_list = m.second;

                for (auto const &v: effects_list)
                {
                    int pos_in_extra = parser.random_and_fixed_in_extra_storage[v];

                    if ( pos_in_extra == -1 )
                        parser.random_and_fixed_effects[v].get(eff_as_storage);                    
                    else
                        parser.extra_effects[ (size_t)pos_in_extra ].get(eff_as_storage);
                    
                    model.append_effect(eff_as_storage);

                    eff_as_storage.clear();
                }
            }

            // (4) ------------ Trait model --------------
            for (auto  const &m: parser.model_definition)
            {
                int obs = m.first;
                std::vector<int> effects(m.second.begin(), m.second.end());
                model.append_traitstruct(obs, effects);
            }

            // (5) ----------- Correlations --------------
            size_t count_correlations = 0; // which correlation group
            for (auto  const &m: parser.correlated_effects)
            {
                std::vector<int> effects(m.begin(), m.end()); // list of correlated effects in submission order
                compact_storage<float> correlation; // correlation matrix
                std::string identity; // identity correlation matrix

                if (parser.corr_matr_index_in_extra_vars[count_correlations][0] >= 0) // corr matrix
                    parser.extra_effects[ parser.corr_matr_index_in_extra_vars[count_correlations][0] ].get(correlation);
                else
                    identity = "I";

                compact_storage<float> var_as_storage;
                std::vector<float> var;
                size_t nrows = 0;

                parser.extra_effects[ parser.corr_matr_index_in_extra_vars[count_correlations][1] ].get(var_as_storage); // variance matrix
                var_as_storage.to_dense(var);
                nrows = var_as_storage.nrows();

                if ( identity.empty() )
                    model.append_corrstruct(var, nrows, correlation, effects);
                else
                    model.append_corrstruct(var, nrows, identity, effects);

                effects.clear(); effects.shrink_to_fit();
                correlation.clear();
                var.clear(); var.shrink_to_fit();
                var_as_storage.clear();

                count_correlations++;
            }

            // (4) Missing values
            model.set_missing( parser.obs_missing_value );
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