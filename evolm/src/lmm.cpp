#include "lmm.hpp"
#include <sstream>

namespace evolm
{
    //===============================================================================================================
    lmm::lmm()
    {
        std::vector<std::string> sol_header(5);
        sol_header[0] = "Trait";
        sol_header[1] = "Group";
        sol_header[2] = "Name";
        sol_header[3] = "Levels (name x group)";
        sol_header[4] = "Estimate";
        sol_table.push_back(sol_header);
        
        std::vector<std::string> sol_header2(8);
        sol_header2[0] = "Trait";
        sol_header2[1] = "Group";
        sol_header2[2] = "Name";
        sol_header2[3] = "First Level";
        sol_header2[4] = "Last Level";
        sol_header2[5] = "Num Levels";
        sol_header2[6] = "Starts at Row";
        sol_header2[7] = "Ends at Row";
        sol_structure.push_back(sol_header2);
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
            std::cerr << "lmm::define(const std::string &): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "lmm::define(const std::string &): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "lmm::define(const std::string &): " << "Unknown exception." << '\n';
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
            std::cerr << "lmm::define_infile(const std::string &): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "lmm::define_infile(const std::string &): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "lmm::define_infile(const std::string &): " << "Unknown exception." << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    void lmm::solve(const std::string &use_method, int available_memory, int available_cpu, const std::string &log_file, const std::string &sol_file)
    {
        try
        {
            parser.report(log_file);

            if ( available_memory <= (int)0 )
                throw std::string("The provided memory limit (available memory) should be greater than zerro!");

            if ( available_cpu <= (int)0 )
                throw std::string("The provided cpu limit (available cpu) should be greater than zerro!");

            model_sparse model;

            message_to_log( log_file, "Setting the parsed model ..." );

            set_model(model, log_file);

            message_to_log( log_file, "Passing the model to the solver for processing ..." );

            sparse_pcg solver(log_file);
            solver.append_model(model);
            solver.set_memory_limit( double(available_memory) );
            solver.set_cpu_limit( available_cpu );
            solver.solve();

            std::vector<double> sol_vect;
            
            message_to_log( log_file, "Collecting the model solution from the solver ..." );

            solver.get_solution(sol_vect);

            message_to_log( log_file, "Processing the solution ..." );

            process_solution(sol_vect, sol_file);

            message_to_log( log_file, " " );
            message_to_log( log_file, "Completed successfully." );
        }
        catch(const std::exception& e)
        {
            std::cerr << "lmm::solve(const std::string &, int, int, const std::string &, const std::string &): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "lmm::solve(const std::string &, int, int, const std::string &, const std::string &): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "lmm::solve(const std::string &, int, int, const std::string &, const std::string &): " << "Unknown exception." << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    void lmm::solve(const std::string &use_method, int available_memory, int available_cpu, const std::string &log_file, const std::string &sol_file, double err_tol)
    {
        try
        {
            parser.report(log_file);

            if ( available_memory <= (int)0 )
                throw std::string("The provided memory limit (available memory) should be greater than zerro!");

            if ( available_cpu <= (int)0 )
                throw std::string("The provided cpu limit (available cpu) should be greater than zerro!");

            model_sparse model;

            message_to_log( log_file, "Setting the parsed model ..." );

            set_model(model, log_file);

            message_to_log( log_file, "Passing the model to the solver for processing it ..." );

            sparse_pcg solver(log_file);
            solver.append_model(model);
            solver.set_memory_limit( double(available_memory) );
            solver.set_cpu_limit( available_cpu );
            solver.set_tolerance(err_tol); // 1e-6 is the default value in the solver
            solver.solve();

            std::vector<double> sol_vect;

            message_to_log( log_file, "Collecting the model solution from the solver ..." );
            
            solver.get_solution(sol_vect);

            message_to_log( log_file, "Processing the solution ..." );

            process_solution(sol_vect, sol_file);

            message_to_log( log_file, " " );
            message_to_log( log_file, "Completed successfully." );
        }
        catch(const std::exception& e)
        {
            std::cerr << "lmm::solve(const std::string &, int, int, const std::string &, const std::string &, double): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "lmm::solve(const std::string &, int, int, const std::string &, const std::string &, double): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "lmm::solve(const std::string &, int, int, const std::string &, const std::string &, double): " << "Unknown exception." << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    void lmm::solve(const std::string &use_method, int available_memory, int available_cpu, const std::string &log_file, const std::string &sol_file, double err_tol, size_t max_iter)
    {
        try
        {
            parser.report(log_file);

            if ( available_memory <= (int)0 )
                throw std::string("The provided memory limit (available memory) should be greater than zerro!");

            if ( available_cpu <= (int)0 )
                throw std::string("The provided cpu limit (available cpu) should be greater than zerro!");

            model_sparse model;

            message_to_log( log_file, "Setting the parsed model ..." );

            set_model(model, log_file);

            message_to_log( log_file, "Passing the model to the solver for processing it ..." );

            sparse_pcg solver(log_file);
            solver.append_model(model);
            solver.set_memory_limit( double(available_memory) );
            solver.set_cpu_limit( available_cpu );
            solver.set_tolerance(err_tol); // 1e-6 is the default value in the solver
            solver.set_maxiter(max_iter); // default value depends on the size of RHS
            solver.solve();

            std::vector<double> sol_vect;

            message_to_log( log_file, "Collecting the model solution from the solver ..." );
            
            solver.get_solution(sol_vect);

            message_to_log( log_file, "Processing the solution ..." );

            process_solution(sol_vect, sol_file);

            message_to_log( log_file, " " );
            message_to_log( log_file, "Completed successfully." );
        }
        catch(const std::exception& e)
        {
            std::cerr << "lmm::solve(const std::string &, int, int, const std::string &, const std::string &, double, size_t): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "lmm::solve(const std::string &, int, int, const std::string &, const std::string &, double, size_t): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "lmm::solve(const std::string &, int, int, const std::string &, const std::string &, double, size_t): " << "Unknown exception." << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    void lmm::set_model(model_sparse &model, const std::string &log_file)
    {
        // defining model
        try
        {
            message_to_log( log_file, "  ==> Appending Residual ..." );
            
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
            
            message_to_log( log_file, "  ==> Appending Observations ..." );

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

            message_to_log( log_file, "  ==> Appending Effects ..." );

            // (3) ----------- Effects ----------------
            size_t num_records_in_solution = 1;
            for (auto  const &m: parser.model_definition)
            {
                std::vector<size_t> effects_list = m.second;
                size_t obs = m.first;

                for (auto const &v: effects_list)
                {
                    compact_storage<float> eff_as_storage;

                    int pos_in_extra = parser.random_and_fixed_in_extra_storage[v];

                    if ( pos_in_extra == -1 )
                        parser.random_and_fixed_effects[v].get(eff_as_storage);
                    else
                        parser.extra_effects[ (size_t)pos_in_extra ].get(eff_as_storage);

                    std::string trait_name = parser.observations_names[obs];
                    append_solution_table(parser.random_and_fixed_effects_names[v], parser.random_and_fixed_effects_levels[v], trait_name);
                    append_solution_structure(parser.random_and_fixed_effects_names[v], parser.random_and_fixed_effects_levels[v], trait_name, num_records_in_solution);
                    num_records_in_solution = num_records_in_solution + parser.random_and_fixed_effects_levels[v].size();
                    model.append_effect(eff_as_storage);
                }
            }

            message_to_log( log_file, "  ==> Appending Trait structure ..." );

            // (4) ------------ Trait model --------------
            for (auto  const &m: parser.model_definition)
            {
                int obs = m.first;
                std::vector<int> effects(m.second.begin(), m.second.end());
                model.append_traitstruct(obs, effects);
            }

            message_to_log( log_file, "  ==> Appending Correlations ..." );

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
                {
                    model.append_corrstruct(var, nrows, correlation, effects);
                }
                else
                {
                    compact_storage<float> eff_as_storage;
                    int pos_in_extra = parser.random_and_fixed_in_extra_storage[ effects[0] ];
                    if ( pos_in_extra == -1 )
                        parser.random_and_fixed_effects[ effects[0] ].get(eff_as_storage);                    
                    else
                        parser.extra_effects[ (size_t)pos_in_extra ].get(eff_as_storage);                    
                    size_t n_cols = eff_as_storage.ncols();
                    eff_as_storage.clear();

                    model.append_corrstruct(var, nrows, identity, n_cols, effects);
                }

                effects.clear(); effects.shrink_to_fit();
                correlation.clear();
                var.clear(); var.shrink_to_fit();
                var_as_storage.clear();

                count_correlations++;
            }

            message_to_log( log_file, "  ==> Appending Missing values ..." );

            // (4) Missing values
            model.set_missing( parser.obs_missing_value );
        }
        catch(const std::exception& e)
        {
            std::cerr << "lmm::set_model(model_sparse &): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "lmm::set_model(model_sparse &): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "lmm::set_model(model_sparse &): " << "Unknown exception." << '\n';
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
            std::cerr << "lmm::read_model_from_file(const std::string &, std::vector<std::string> &): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "lmm::read_model_from_file(const std::string &, std::vector<std::string> &): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "lmm::read_model_from_file(const std::string &, std::vector<std::string> &): " << "Unknown exception." << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    void lmm::clear()
    {
        try
        {
            parser.clear();
            sol_table.clear();
            sol_table.shrink_to_fit();
        }
        catch(const std::exception& e)
        {
            std::cerr << "lmm::clear(): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "lmm::clear(): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "lmm::clear(): " << "Unknown exception." << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    std::string lmm::to_string_with_precision(const double value, const int n = 6)
    {
        try
        {
            std::ostringstream out;
            out.precision(n);
            out << std::fixed << value;

            return std::move(out).str();
        }
        catch(const std::exception& e)
        {
            std::cerr << "lmm::to_string_with_precision(const double, const int n = 6): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "lmm::to_string_with_precision(const double, const int n = 6): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "lmm::to_string_with_precision(const double, const int n = 6): " << "Unknown exception." << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    void lmm::process_solution(std::vector<double> &sol_vect, const std::string &sol_file)
    {
        try
        {
            int use_precision = 6;

            if ( sol_vect.size() != (sol_table.size() - 1) )
                throw std::string("The size of solution vector does not correspond to the expected size: sol_vect.size() != (sol_table.size() - 1)!");
            
            std::ofstream solution(sol_file);

            if (solution.is_open())
            {
                // do not use header in this file
                //solution << sol_table[0][0] << "," << sol_table[0][1] << "," << sol_table[0][2] << "," << sol_table[0][3] << "," << sol_table[0][4] << "\n";

                for (size_t i = 1; i < sol_table.size(); i++)
                    solution << sol_table[i][0] << "," << sol_table[i][1] << "," << sol_table[i][2] << "," << sol_table[i][3] << "," << to_string_with_precision( sol_vect[i-1], use_precision ) << "\n";

                solution.close();
            }
            else
                throw std::string("Unable to open solution file " + sol_file + " for output!");

            std::string structure_file = "sol_struct_"+sol_file;
            solution.open(structure_file);

            if (solution.is_open())
            {
                solution << "The names of columns in the " << sol_file << " file" << "\n";
                solution << "\n";
                solution << "Column 1: Trait" << "\n";
                solution << "Column 2: Group" << "\n";
                solution << "Column 3: Name" << "\n";
                solution << "Column 4: Levels (name x group)" << "\n";
                solution << "Column 5: Estimate" << "\n";
                solution << "\n";
                
                solution << "Structure of estimates in the " << sol_file << " file:" << "\n";
                solution << "(Access the specific estimates using 'Starts at Row' and 'Ends at Row')"<< "\n";
                solution << "\n";
                int dist = 23;
                int dist2 = 14;
                for (size_t i = 0; i < sol_structure.size(); i++)
                    solution <<std::setw(dist2)<<std::left << sol_structure[i][0] << " "<<std::setw(dist)<<std::left << sol_structure[i][1] << " " <<std::setw(dist)<<std::left << sol_structure[i][2] << " " <<std::setw(dist)<<std::left << sol_structure[i][3] << " " <<std::setw(dist)<<std::left << sol_structure[i][4] << " " <<std::setw(dist2)<<std::left << sol_structure[i][5] << " " <<std::setw(dist2)<<std::left << sol_structure[i][6] << " " <<std::setw(dist2)<<std::left << sol_structure[i][7] << "\n";

                solution.close();
            }
            else
                throw std::string("Unable to open solution structure file " + structure_file + " for output!");    
        }
        catch(const std::exception& e)
        {
            std::cerr << "lmm::process_solution(std::vector<double> &, const std::string &): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "lmm::process_solution(std::vector<double> &, const std::string &): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "lmm::process_solution(std::vector<double> &, const std::string &): " << "Unknown exception." << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    void lmm::append_solution_table(std::string &var_name, std::vector<std::string> &levels, std::string &trait)
    {
        try
        {
            std::string group;
            std::string name;
            std::string bar("|");
            std::string lhs;
            std::string rhs;
            size_t pos = 0;
            pos = var_name.find(bar);
            if (pos != std::string::npos)
            {
                lhs = var_name.substr(0, pos);
                rhs = var_name.substr(pos+1, std::string::npos);

                if (lhs == "1")
                    name = "Intercept";
                else
                    name = lhs;
                
                group = rhs;
            }
            else
            {
                group = "none";
                name = var_name;
            }
            for (size_t i = 0; i < levels.size(); i++)
            {
                std::vector<std::string> table(5);
                table[0] = trait;
                table[1] = group;
                table[2] = name;
                table[3] = levels[i];
                sol_table.push_back(table);
            }
        }
        catch(const std::exception& e)
        {
            std::cerr << "lmm::append_solution_table(std::string &, std::vector<std::string> &, std::string &): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "lmm::append_solution_table(std::string &, std::vector<std::string> &, std::string &): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "lmm::append_solution_table(std::string &, std::vector<std::string> &, std::string &): " << "Unknown exception." << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    void lmm::append_solution_structure(std::string &var_name, std::vector<std::string> &levels, std::string &trait, size_t record)
    {
        try
        {
            std::string group;
            std::string name;
            std::string bar("|");
            std::string lhs;
            std::string rhs;
            size_t pos = 0;
            pos = var_name.find(bar);
            if (pos != std::string::npos)
            {
                lhs = var_name.substr(0, pos);
                rhs = var_name.substr(pos+1, std::string::npos);

                if (lhs == "1")
                    name = "Intercept";
                else
                    name = lhs;
                
                group = rhs;
            }
            else
            {
                group = "none";
                name = var_name;
            }

            std::vector<std::string> table(8);
            table[0] = trait;
            table[1] = group;
            table[2] = name;
            table[3] = levels[0]; // the first level
            table[4] = levels[levels.size()-1]; // the last level
            table[5] = std::to_string(levels.size()); // total num of level
            table[6] = std::to_string(record); // very first row in Estimates file
            table[7] = std::to_string(record + levels.size() - 1); // very last row in Estimates file
            
            sol_structure.push_back(table);
        }
        catch(const std::exception& e)
        {
            std::cerr << "lmm::append_solution_structure(std::string &, std::vector<std::string> &, std::string &, size_t): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "lmm::append_solution_structure(std::string &, std::vector<std::string> &, std::string &, size_t): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "lmm::append_solution_structure(std::string &, std::vector<std::string> &, std::string &, size_t): " << "Unknown exception." << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    void lmm::message_to_log( const std::string &log_file, std::string msg )
    {
        try
        {
            std::ofstream summary(log_file, std::ofstream::out | std::ofstream::app);

            if ( !summary.is_open() )
                throw "Unable to open log file for output!";

            summary<<" "<<msg<<'\n';
            summary.close();
        }
        catch(const std::exception& e)
        {
            std::cerr << "lmm::message_to_log(const std::string &, std::string): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "lmm::message_to_log(const std::string &, std::string): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "lmm::message_to_log(const std::string &, std::string): " << "Unknown exception." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void lmm::snp_to_bv( const std::string &fname_zmatr, const std::string &fname_solsnp, const std::string &fname_outres )
    {
        try
        {
            parser.binmatr_by_txtvect(fname_zmatr, fname_solsnp, fname_outres);
        }
        catch(const std::exception& e)
        {
            std::cerr << "lmm::snp_to_bv(const std::string &, const std::string): " << e.what() << '\n';
            throw e;
        }
        catch(const std::string & e)
        {
            std::cerr << "lmm::snp_to_bv(const std::string &, const std::string): " << e <<'\n';
            throw e;
        }
        catch(...)
        {
            std::cerr << "lmm::snp_to_bv(const std::string &, const std::string): " << "Unknown exception." << '\n';
            throw;
        }
    }
    //===============================================================================================================
} // end of namespace evolm