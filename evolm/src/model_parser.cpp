#include "model_parser.hpp"

namespace evolm
{
    // -------------------------------------------------------------
    model_parser::model_parser()
    {
    }
    // -------------------------------------------------------------
    model_parser::~model_parser()
    {
        clear();
    }
    // -------------------------------------------------------------
    void model_parser::clear()
    {
        special_corr_vars.clear();
        special_corr_vars.shrink_to_fit();

        extra_effects.clear();
        extra_effects.shrink_to_fit();

        extra_effects_names.clear();
        extra_effects_names.shrink_to_fit();

        extra_effects_levels.clear();
        extra_effects_levels.shrink_to_fit();

        observations.clear();
        observations.shrink_to_fit();

        observations_names.clear();
        observations_names.shrink_to_fit();

        obs_in_extra_storage.clear();
        obs_in_extra_storage.shrink_to_fit();

        random_and_fixed_effects.clear();
        random_and_fixed_effects.shrink_to_fit();

        random_and_fixed_effects_names.clear();
        random_and_fixed_effects_names.shrink_to_fit();

        random_and_fixed_in_extra_storage.clear();
        random_and_fixed_in_extra_storage.shrink_to_fit();

        random_and_fixed_effects_levels.clear();
        random_and_fixed_effects_levels.shrink_to_fit();

        corr_vars_index_in_effects.clear();
        corr_vars_index_in_effects.shrink_to_fit();

        corr_matr_index_in_extra_vars.clear();
        corr_matr_index_in_extra_vars.shrink_to_fit();

        correlated_effects.clear();
        correlated_effects.shrink_to_fit();

        identifier_to_eff_name.clear();

        borroow_levels_from.clear();

        model_definition.clear();

        next_subm_eff = 0;
        next_subm_obs = 0;
    }
    // -------------------------------------------------------------
    void model_parser::process_expression(const std::string &expr)
    {
        try
        {
            // detect multiple expressions from the string separated by ";"
            std::vector<std::string> all_expressions;
            std::vector<size_t> type_1_expressions; // obs = model_expression
            std::vector<size_t> type_2_expressions; // var_name = file_value
            std::vector<size_t> type_3_expressions; // var_name = matrix_value
            std::vector<size_t> type_4_expressions; // var = variance_expression

            split_str2(";", expr, all_expressions); // keep spaces for a while, needed to parse matrix_value expression

            for (size_t i = 0; i < all_expressions.size(); i++) // separate by type in order to process types 2 & 3 before the type 1
            {    
                int the_type = detect_expression_type(all_expressions[i]);
                
                switch (the_type)
                {
                case 1:
                    type_1_expressions.push_back(i);
                    remove_space( all_expressions[i] );
                    break;
                case 2:
                    type_2_expressions.push_back(i);
                    remove_space( all_expressions[i] );
                    break;
                case 3:
                    type_3_expressions.push_back(i);
                    break;
                case 4:
                    type_4_expressions.push_back(i);
                    remove_space( all_expressions[i] );
                    break;
                default:
                    break;
                }
            }

            for (size_t i = 0; i < type_3_expressions.size(); i++)
                eval_expr(all_expressions[ type_3_expressions[i] ]);

            for (size_t i = 0; i < type_2_expressions.size(); i++)
                eval_expr(all_expressions[ type_2_expressions[i] ]);

            for (size_t i = 0; i < type_1_expressions.size(); i++)
            {
                eval_expr(all_expressions[ type_1_expressions[i] ]);
                get_list_of_subm_eff(); // get the list of submitted effects
            }

            for (size_t i1 = 0; i1 < type_4_expressions.size(); i1++)
            {
                eval_expr(all_expressions[ type_4_expressions[i1] ]);
                def_corr_struct(); // define the correlation structure
            }
            
            all_expressions.clear(); all_expressions.shrink_to_fit();
            type_1_expressions.clear(); type_1_expressions.shrink_to_fit();
            type_2_expressions.clear(); type_2_expressions.shrink_to_fit();
            type_3_expressions.clear(); type_3_expressions.shrink_to_fit();
            type_4_expressions.clear(); type_4_expressions.shrink_to_fit();
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::process_expression(std::string &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::process_expression(std::string &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::process_expression(std::string &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::get_list_of_subm_eff()
    {
        try
        {
            std::vector<size_t> sub_eff;

            for (size_t i = next_subm_eff; i < random_and_fixed_effects.size(); i++)
                sub_eff.push_back(i);
            
            next_subm_eff = random_and_fixed_effects.size();
            
            // loop over submitted observations to define the model:
            for (size_t i = next_subm_obs; i < observations.size(); i++)
                model_definition[i] = sub_eff;
            
            next_subm_obs = observations.size();
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::get_list_of_subm_eff()"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::get_list_of_subm_eff()"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::get_list_of_subm_eff()"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::def_corr_struct()
    {
        try
        {
            for (size_t i = 0; i < corr_vars_index_in_effects.size(); i++) // loop over all provided model correlations
            {
                std::vector<size_t> corr_effects_subm_order; // consecutive order across all models for the selected correlation                                                                
                
                for (size_t j = 0; j < corr_vars_index_in_effects[i].size(); j++) // for this specific correlation loop over correleted effects and find if they are in the model
                {
                    size_t order = 0; // consecutive order of submitted effects
                    std::string eff_name = random_and_fixed_effects_names[ corr_vars_index_in_effects[i][j] ]; // index of correlated effect we are looking for in the model definition

                    for (auto const &m: model_definition) // loop over map's values of model_definition (the submited effects indeces)
                    {
                        std::vector<size_t> effects_positions = m.second; // submited effects indeces

                        for (size_t l = 0; l < effects_positions.size(); l++)
                        {
                            if ( random_and_fixed_effects_names[ effects_positions[l] ] == eff_name )
                                corr_effects_subm_order.push_back(order);
                            order++;
                        }
                    }
                }
                std::sort( corr_effects_subm_order.begin(), corr_effects_subm_order.end() );
                correlated_effects.push_back(corr_effects_subm_order);
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::def_corr_struct()"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::def_corr_struct()"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::def_corr_struct()"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::eval_expr(const std::string &expression)
    {
        try
        {
            // for type 2: name-value
            std::vector<std::string> name_value_pair; // container for type 2 expression variables
            
            // these four for type 1: model
            std::vector<std::string> obs_vars; // observations variables
            std::vector<std::vector<std::string>> fixedvars; // processed fixed variables
            std::vector<std::vector<std::vector<std::string>>> lhs_randomvars; // processed random variables
            std::vector<std::vector<std::vector<std::string>>> rhs_randomvars; // processed random variables
            
            // vector of identifiers (optional) for random variables
            std::vector<std::string> random_vars_identifiers;

            int expression_type = detect_expression_type(expression);

            switch (expression_type)
            {
            case 1: // process model equation, for example: obs ~ 1 + x1 + (1|z)
                eval_model_expr(expression, obs_vars, lhs_randomvars, rhs_randomvars, fixedvars, random_vars_identifiers); // parsing the model expression to specific vaariables
                process_observation_vars(obs_vars); // append observations container, and var names to observations_names
                obs_vars.clear(); obs_vars.shrink_to_fit();

                get_missing_obs(); // update missing data structure

                process_fixed_vars(fixedvars); // append random_and_fixed_effects container, and var names to random_and_fixed_effects_names
                fixedvars.clear(); fixedvars.shrink_to_fit();
                process_random_vars(lhs_randomvars, rhs_randomvars, random_vars_identifiers); // append random_and_fixed_effects container, and var names to random_and_fixed_effects_names
                lhs_randomvars.clear(); lhs_randomvars.shrink_to_fit();
                rhs_randomvars.clear(); rhs_randomvars.shrink_to_fit();
                random_vars_identifiers.clear(); random_vars_identifiers.shrink_to_fit();
                
                break;
            case 2: // process name-value expression, for example: var1 = var1_file.dat
                separate_name_file_values(expression, name_value_pair);
                process_name_value_pair(name_value_pair); // append processed data matrices to extra_effects container, corresponding var names to extra_effects_names                
                name_value_pair.clear();
                name_value_pair.shrink_to_fit();
                break;
            case 3: // process numerical name-value expression, for example: var1 = [2 0.5, 0.1 3]
                separate_name_file_values(expression, name_value_pair);
                process_name_matrix_pair(name_value_pair);
                name_value_pair.clear();
                name_value_pair.shrink_to_fit();
                break;
            case 4:
                separate_name_file_values(expression, name_value_pair);
                process_model_variance_expression(name_value_pair[1]);
                break;
            default:
                break;
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::eval_expr(std::string &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::eval_expr(std::string &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::eval_expr(std::string &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::process_name_value_pair(std::vector<std::string> &name_value_pair)
    {
        try
        {
            // handle the keywords: data, effect, obs_missing_value, and suffix .corbin

            if (name_value_pair[0] == "data")
            {
                data_file_name = name_value_pair[1];
                return;
            }

            if (name_value_pair[0] == "effect")
            {
                effect_file_name = name_value_pair[1];
                return;
            }

            if (name_value_pair[0] == "obs_missing_value")
            {
                obs_missing_value = std::stof(name_value_pair[1]);
                return;
            }

            if ( is_str_in_expr(name_value_pair[1], ".corbin") ) // check if there is suffix .corbin in name_value_pair[1]
            {
                std::vector<std::string> cor_var;
                cor_var.push_back(name_value_pair[0]); // var name
                cor_var.push_back(name_value_pair[1]); // .corbin file name
                special_corr_vars.push_back(cor_var);
                return;
            }

            // handle variables

            evolm::IOInterface in;
            in.set_fname(name_value_pair[1]);

            if ( in.is_plink_file(name_value_pair[1]) || consist_snp_substring(name_value_pair[0]) ) // Check if the data is SNP variants
            {
                std::vector<float> values;
                size_t nrows = 0;
                size_t ncols = 0;
                std::vector<std::string> snp_names;

                in.scale_genotypes(values, nrows, ncols, snp_names);

                evolm::compact_storage<float> cs(nrows, ncols);

                cs.append(values);
                
                values.clear();
                values.shrink_to_fit();

                evolm::effects_storage eff;

                eff.set(cs);

                cs.clear();

                extra_effects.push_back(eff);
                extra_effects_names.push_back(name_value_pair[0]);
                extra_effects_levels.push_back(snp_names);

                eff.clear();
            }
            else // handle as ASCII text files; general data file, considered as a float-type matrix
            {
                std::vector<std::vector<float>> var_data;
                std::vector<std::string> levels;

                in.fgetdata(var_data);

                evolm::compact_storage<float> cs(var_data.size(), var_data[0].size());

                for (size_t i = 0; i < var_data.size(); i++)
                    for (size_t j = 0; j < var_data[0].size(); j++)
                        cs.append(var_data[i][j], i, j);                
                
                for (size_t i = 0; i < var_data[0].size(); i++)
                    levels.push_back(std::to_string(i+1));

                cs.optimize();

                var_data.clear();
                var_data.shrink_to_fit();

                evolm::effects_storage eff;
                eff.set(cs);

                cs.clear();

                extra_effects.push_back(eff);
                extra_effects_names.push_back(name_value_pair[0]);
                extra_effects_levels.push_back(levels);

                eff.clear();
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::process_name_value_pair(std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::process_name_value_pair(std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::process_name_value_pair(std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::process_name_matrix_pair(std::vector<std::string> &name_value_pair)
    {
        try
        {
            std::vector<std::string> matrix_str; // stores 'matr' out from '[matr]' 
            std::vector<std::string> vectors_in_matr_str; // stores vectors out from 'matr': vect1, vect2, ...
            
            extract_expr_btw_brackets(name_value_pair[1], "[", "]", matrix_str); // extracts 'matr' from '[matr]'

            if ( matrix_str.size() != 1 )
                throw std::string("Incorrect numerical expresion: " + matrix_str[0]);

            split_str2(",", matrix_str[0], vectors_in_matr_str); // extracts vectors out from 'matr': vect1, vect2, ...

            size_t num_rows = vectors_in_matr_str.size();
            size_t num_cols = 0;

            std::vector<std::vector<float>> var_data;

            for (size_t i = 0; i < num_rows; i++)
            {
                std::vector<std::string> vector_str;
                split_str2(" ", vectors_in_matr_str[i], vector_str);
                std::vector<float> vector_num;
                for (auto const &v: vector_str)
                    vector_num.push_back( std::stof(v) );
                var_data.push_back(vector_num);
            }

            num_cols = var_data[0].size();

            for (size_t i = 1; i < num_rows; i++)
                if ( num_cols != var_data[i].size() )
                    throw std::string("Incorrect numerical expresion " + matrix_str[0] + ": inconsistent number of columns in the matrix expression!");            

            matrix_str.clear();
            matrix_str.shrink_to_fit();
            vectors_in_matr_str.clear();
            vectors_in_matr_str.shrink_to_fit();

            remove_space( name_value_pair[0] ); // because, unlike the other types, the spaces were not removed at the beggining

            if ( name_value_pair[0] == "obs_missing_value" )
            {
                if ( var_data.size() != 1 )
                    throw std::string("The missing value constant should be a floating-point scalar but not a vector or matrix!");
                
                if ( var_data[0].size() != 1 )
                    throw std::string("The missing value constant should be a floating-point scalar but not a vector or matrix!");
                
                obs_missing_value = var_data[0][0];

                var_data.clear();
                var_data.shrink_to_fit();

                return;
            }

            evolm::compact_storage<float> cs(var_data.size(), var_data[0].size());

            for (size_t i = 0; i < var_data.size(); i++)
                for (size_t j = 0; j < var_data[0].size(); j++)
                    cs.append(var_data[i][j], i, j);                
            
            cs.optimize();

            var_data.clear();
            var_data.shrink_to_fit();

            evolm::effects_storage eff;
            eff.set(cs);

            cs.clear();

            extra_effects.push_back(eff);

            //remove_space( name_value_pair[0] ); // because, unlike the other types, the spaces were not removed at the beggining

            extra_effects_names.push_back(name_value_pair[0]);

            eff.clear();
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::process_name_matrix_pair(std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::process_name_matrix_pair(std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::process_name_matrix_pair(std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    bool model_parser::consist_snp_substring(std::string &var_name)
    {
        std::string v1("snp");
        std::string v2("SNP");
        std::string v3("Snp");
        
        size_t consists_v1 = var_name.find(v1);
        size_t consists_v2 = var_name.find(v2);
        size_t consists_v3 = var_name.find(v3);

        if ( consists_v1 != std::string::npos )
            return true;
        else if ( consists_v2 != std::string::npos )
            return true;
        else if ( consists_v3 != std::string::npos )
            return true;
        else
            return false;
    }
    // -------------------------------------------------------------
    void model_parser::process_model_variance_expression(const std::string &expr)
    {
        try
        {
            // var = (cor_var_1, cor_var_2) * cor_matrix * var_matr + (cor_var_3) * I * var_matr2 + R

            std::string mdf_expr = expr; // copy to be able to modify
            
            std::vector<std::string> corr_vars; // temporal, to be cleared soon
            std::vector<std::vector<std::string>> correlated_vars; // list of correlated variables/effects: cor_var_1 & cor_var_2; cor_var_3

            extract_expr_btw_brackets(mdf_expr, "(", ")", corr_vars); // extracts: cor_var_1 & cor_var_2; cor_var_3

            for (size_t i = 0; i < corr_vars.size(); i++) // split corr vars for each random effect
            {
                std::vector<std::string> t_vars;
                split_str(",", corr_vars[i], t_vars);
                correlated_vars.push_back(t_vars); // correlated variables
            }

            corr_vars.clear();
            corr_vars.shrink_to_fit();

            std::vector<std::string> corr_groups; // temporal, to be cleared soon
            std::vector<std::vector<std::string>> correlated_groups; // list of pairs: corr matrix variable, var matrix variable. Group1: cor_matrix, var_matr. Group2: I, var_matr2. Group3: R.

            split_str("+", mdf_expr, corr_groups); // split and extract: cor_matrix * var_matr; I * var_matr2; R

            for (size_t i = 0; i < corr_groups.size(); i++) // split pairs: cor_matrix & var_matr; I & var_matr2; R
            {
                std::vector<std::string> t_vars;
                split_str("*", corr_groups[i], t_vars);
                correlated_groups.push_back(t_vars);
            }

            corr_groups.clear();
            corr_groups.shrink_to_fit();
            
            mdf_expr.clear();

            if ( correlated_vars.size() != correlated_groups.size()-1 )
                throw std::string("Incorrect expression for the model variance: correlated_vars.size() != correlated_groups.size()-1");

            find_corr_vars_in_effects(correlated_vars); // for each correlated var find index (position) in effects vector, writing to the class's corr_vars_index_in_effects container

            correlated_vars.clear();
            correlated_vars.shrink_to_fit();

            evolm::IOInterface in;

            // check if the correlated_groups matrices are added
            size_t n_terms = 0;
            for (size_t i = 0; i < correlated_groups.size(); i++) // loop over all groups (number of random effects and residual)
            {
                std::vector<int> indeces;

                for (size_t j = 0; j < correlated_groups[i].size(); j++) // loop over variables (correlations and variances) in a group
                {
                    n_terms++;

                    int index = in.find_value(extra_effects_names, correlated_groups[i][j]); // look at extra_effects_names for already submitted matrix

                    if ( index == -1 && correlated_groups[i][j] != "I" ) // variable is not in extra. vars container and not identity matrix
                    {                        
                        for (size_t l = 0; l < special_corr_vars.size(); l++) // then, need to check in the special_corr_vars for .corbin file
                        {
                            std::vector<std::string> t_str;
                            t_str.push_back(special_corr_vars[l][0]); // this is just to be able to call the find_value(...) method
                            
                            index = in.find_value(t_str, correlated_groups[i][j]);
                            
                            if (index != -1) // found, variable have .corbin file
                            {
                                // here we need to read correlation matrix, process it and add it to the extra vars

                                std::string cor_matr_name = correlated_groups[i][j];
                                std::string cor_matr_file = special_corr_vars[l][1];
                                
                                size_t pos_in_effects = corr_vars_index_in_effects[i][0]; // index for correlated variable in effects list

                                std::string cor_variable = random_and_fixed_effects_names[pos_in_effects]; // var name to look at a data file for num of levels
                                std::vector<std::string> levels = random_and_fixed_effects_levels[pos_in_effects];
                                                                
                                compact_storage<float> corr_storage;

                                load_corbin_file(cor_matr_file, cor_variable, levels, corr_storage);

                                // add corr_storage to extra_vars storage

                                evolm::effects_storage eff;
                                eff.set(corr_storage);

                                corr_storage.clear();

                                extra_effects.push_back(eff);
                                extra_effects_names.push_back(cor_matr_name);

                                eff.clear();

                                index = extra_effects.size() - 1;

                                break;
                            }

                            t_str.clear();
                        }

                        if (index == -1) // Cannot find in special_corr_vars
                            throw std::string("Cannot find the variable " + correlated_groups[i][j] + " variable.");
                    }

                    indeces.push_back(index); // interpretation: positive - index in extra vars; -1 - identity matrix.
                }

                corr_matr_index_in_extra_vars.push_back(indeces);
            }

            if ( n_terms != (2 * correlated_groups.size() - 1) )
                throw std::string("Incorrect number of terms in the model variance expression!");

            correlated_groups.clear();
            correlated_groups.shrink_to_fit();
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::process_model_variance_expression(std::string &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::process_model_variance_expression(std::string &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::process_model_variance_expression(std::string &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::find_corr_vars_in_effects(std::vector<std::vector<std::string>> &corr_vars)
    {        
        try
        {
            evolm::IOInterface in;

            for (size_t i = 0; i < corr_vars.size(); i++) // for each random (correlated) variable find index (position) in effects vector
            {
                std::vector<int> indeces;

                for (size_t j = 0; j < corr_vars[i].size(); j++)
                {
                    std::string t_var_name = corr_vars[i][j];

                    // (1) check the aliases map, return the original name if alias is used
                    std::unordered_map<std::string,std::string>::const_iterator it = identifier_to_eff_name.find(t_var_name);                                        
                    if ( it != identifier_to_eff_name.end() )
                        t_var_name = it->second;
                    
                    // (2) if there is no special identifier provided, check the name as it is
                    int index = in.find_value(random_and_fixed_effects_names, t_var_name);
                    
                    if ( index == -1 ) // if not found
                    {
                        std::string t_var = "1|" + t_var_name;
                        index = in.find_value(random_and_fixed_effects_names, t_var); // check var with '1|' notation
                        if ( index == -1 ) // variable not found, report error
                            throw std::string("Cannot find correlated variable " + t_var_name);
                    }
                    
                    indeces.push_back(index);
                }

                corr_vars_index_in_effects.push_back(indeces); // push indeces to the global scope container
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::find_corr_vars_in_effects(std::vector<std::vector<std::string>> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::find_corr_vars_in_effects(std::vector<std::vector<std::string>> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::find_corr_vars_in_effects(std::vector<std::vector<std::string>> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::load_corbin_file(std::string &cor_matr_file, std::string &cor_variable, std::vector<std::string> &var_levels, compact_storage<float> &corr_storage)
    {
        try
        {
            if( data_file_name.empty() )
                throw std::string("The data file name is empty!");

            IOInterface in;
            
            std::vector<int64_t> int_levels;
            std::vector<std::string> str_levels( var_levels.size() );

            std::string interaction_str(" x ");
            size_t pos = var_levels[0].rfind(interaction_str);

            if ( pos != std::string::npos ) // checking for possible interaction bbetween variables
            {
                size_t start = pos + interaction_str.size();
#pragma omp parallel for
                for (size_t i = 0; i < var_levels.size(); i++) // extract only the levels for grouping variable (very last one): only for .corbin files!
                    str_levels[i] = var_levels[i].substr(start, std::string::npos);
            }
            else
            {
#pragma omp parallel for
                for (size_t i = 0; i < var_levels.size(); i++)
                    str_levels[i] = var_levels[i];
            }

            size_t cor_matr_info[3]; // info storage for a data in .cprbin file
            in.fread_matrix_info( cor_matr_file, cor_matr_info ); // get info

            std::vector<float> f_vals;
            std::vector<double> d_vals;
            std::vector<std::int64_t> i_ids;
            std::vector<std::string> s_ids;
            std::vector<size_t> keys;

            std::vector<size_t> pos_map; // permutation vector
            std::vector<std::int64_t> pos_map_red; // permutation vector, in the case of corr matrix needs to be reduced
            
            if( cor_matr_info[2] == 5 ) // expected type int64 of ids
            {       
                int_levels.resize(str_levels.size());

#pragma omp parallel for
                for (size_t i = 0; i < str_levels.size(); i++) // conversion from std::string to int64
                    int_levels[i] = std::stol( str_levels[i] );

                if ( cor_matr_info[1] == 2 ) // get cor matrix of type float and int64 ids
                {
                    in.fread_matrix(cor_matr_file, f_vals, keys, i_ids);
                }
                else if ( cor_matr_info[1] == 3 ) // get cor matrix of type double and int64 ids
                {
                    in.fread_matrix(cor_matr_file, d_vals, keys, i_ids);
                    
                    f_vals.resize(d_vals.size());
                    std::copy( d_vals.begin(), d_vals.end(), f_vals.begin() ); // need convertion to float
                    d_vals.clear();
                    d_vals.shrink_to_fit();
                }
                else
                    throw std::string("Undefined type of correlation matrix. Expected float or double.");
                
                if ( i_ids.size() < int_levels.size() )
                    throw std::string("The number of rows in the correlation matrix " + cor_variable + " stored in the file " + cor_matr_file + " is lowe than the number of levels in the effect!");

                if ( i_ids.size() > int_levels.size() ) // require permutation and reduction
                {
                    size_t num_of_found = 0;
                    for (size_t i = 0; i < i_ids.size(); i++) // make permutation vector
                    {
                        int found = in.find_value(int_levels, i_ids[i]);
                        pos_map_red.push_back(found);
                        if ( found != -1 )
                            num_of_found++;                                            
                    }
                    if ( num_of_found != int_levels.size() )
                        throw std::string("There are IDs in the correlation matrix which are not in Effect matrix!");
                }
                else // require only permutation
                {
                    for (size_t i = 0; i < i_ids.size(); i++) // make permutation vector
                    {
                        int found = in.find_value(int_levels, i_ids[i]);
                        if ( found != -1 )
                            pos_map.push_back(found);
                        else
                            throw std::string("The following ID " + std::to_string(i_ids[i]) + " from correlation matrix is not in Effect matrix!");
                    }
                }

                corr_storage.resize( i_ids.size() ); // expected symmeric matrix in .corbin

                i_ids.clear();
                i_ids.shrink_to_fit();
            }

            if ( cor_matr_info[2] == 4 ) // expected type std::string of ids
            {                
                if ( cor_matr_info[1] == 2 )
                {
                    in.fread_matrix(cor_matr_file, f_vals, keys, s_ids);
                }
                else if ( cor_matr_info[1] == 3 )
                {
                    in.fread_matrix(cor_matr_file, d_vals, keys, s_ids);

                    f_vals.resize(d_vals.size());
                    std::copy( d_vals.begin(), d_vals.end(), f_vals.begin() ); // need convertion to float
                    d_vals.clear();
                    d_vals.shrink_to_fit();
                }
                else
                    throw std::string("Undefined type of correlation matrix. Expected float or double.");

                if ( s_ids.size() < str_levels.size() )
                    throw std::string("The number of rows in the correlation matrix " + cor_variable + " stored in the file " + cor_matr_file + " is lowe than the number of levels in the effect!");

                if ( s_ids.size() > str_levels.size() ) // require permutation and reduction
                {
                    size_t num_of_found = 0;
                    for (size_t i = 0; i < s_ids.size(); i++) // make permutation vector
                    {
                        int found = in.find_value(str_levels, s_ids[i]);
                        pos_map_red.push_back(found);
                        if ( found != -1 )
                            num_of_found++;                                            
                    }
                    if ( num_of_found != str_levels.size() )
                        throw std::string("There are IDs in the correlation matrix which are not in Effect matrix!");
                }
                else // require only permutation
                {
                    for (size_t i = 0; i < s_ids.size(); i++) // make permutation vector
                    {
                        int found = in.find_value(str_levels, s_ids[i]);
                        if ( found != -1 )
                            pos_map.push_back(found);
                        else
                            throw std::string("The following ID " + s_ids[i] + " from correlation matrix is not in Effect matrix!");
                    }
                }

                corr_storage.resize( s_ids.size() );

                s_ids.clear();
                s_ids.shrink_to_fit();
            }

            corr_storage.append_with_keys(f_vals, keys);

            f_vals.clear();
            f_vals.shrink_to_fit();
            keys.clear();
            keys.shrink_to_fit();

            if ( !pos_map.empty() )
                corr_storage.permute(pos_map); // rearrange matrix according to effect ids list
            else
                corr_storage.permute_and_reduce(pos_map_red); // rearrange and reduce matrix according to effect ids list
            
            pos_map.clear();
            pos_map.shrink_to_fit();
            pos_map_red.clear();
            pos_map_red.shrink_to_fit();
            
            corr_storage.sym_to_rec(); // make it rectangular
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::load_corbin_file(std::string &, std::string &, std::vector<std::string> &, compact_storage<float> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::load_corbin_file(std::string &, std::string &, std::vector<std::string> &, compact_storage<float> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::load_corbin_file(std::string &, std::string &, std::vector<std::string> &, compact_storage<float> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::process_observation_vars(std::vector<std::string> &observ_vars)
    {
        try
        {
            evolm::IOInterface in;

            for (size_t i = 0; i < observ_vars.size(); i++)
            {
                evolm::effects_storage s;

                // first, look at the extra_effects storage if the variable observ_vars[i] is already processed
                int var_find_res = in.find_value(extra_effects_names, observ_vars[i]);

                if ( var_find_res != -1 )
                {
                    // variable is in the extra storage
                    obs_in_extra_storage.push_back(var_find_res); // we know now where to find observ_vars[i] variable
                    observations.push_back(s); // push empty effect to keep the vector size consistent
                    observations_names.push_back(observ_vars[i]);
                    continue; // move to the next variable
                }
                
                obs_in_extra_storage.push_back(-1); // this vector should be the same size as observations

                set_fname(in, observ_vars[i]);

                std::vector<std::string> levels; // this is needed here just to pass all required parameters to the function fgetvar !
                in.fgetvar(observ_vars[i], obs_missing_value, missing_obs, s, levels); // convert string obs[i] to vector data observations[i]

                std::vector<size_t> obs_shape = s.shape();

                if ( obs_shape[0] == 0 || obs_shape[1] == 0 )
                    throw std::string( "The shape of '" + observ_vars[i] + "' is 0!" );
                
                if ( obs_shape[1] != 1 )
                    throw std::string( "The number of columns of the variable '" + observ_vars[i] + "' should be 1! Check if the specific data has a floating-point format. Otherwise it cannot be considered as an observation." );

                observations.push_back(s);
                observations_names.push_back(observ_vars[i]);
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::process_observation_vars(std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::process_observation_vars(std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::process_observation_vars(std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::get_missing_obs()
    {
        // missing_obs (class member) vector - holding the method's results of detection missing observations structure;
        // missing records defined as TRUE elements. 
        try
        {
            evolm::effects_storage s;
            evolm::compact_storage<float> c;
            std::vector<float> obs;
            
            int extra_str_index = obs_in_extra_storage.back(); // check where to find very last observation
            
            if (extra_str_index == -1)
                s = observations.back(); // observation is in normal obs storage
            else
                s = extra_effects[extra_str_index]; // observation is in extra eff storage

            s.get(c); // convert effects_storage to compact_storage
            s.clear();
            c.to_dense(obs); // convert compact_storage to std::vector
            c.clear();

            if ( obs.empty() )
                throw std::string("Vector of very last submitted observation is empty.");
            
            std::vector<bool> missing(obs.size(), false);

            for (size_t i = 0; i < obs.size(); i++)
                if ( obs[i] == obs_missing_value )
                    missing[i] = true;

            missing_obs = missing;
        }
        catch(const std::exception& e)
        {
            std::cerr << "Exception in model_parser::get_missing_obs(std::vector<bool> &)"<< "\n";
            std::cerr << e.what() << '\n';
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::get_missing_obs(std::vector<bool> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::process_fixed_vars(std::vector<std::vector<std::string>> &fixed_vars)
    {
        try
        {
            evolm::IOInterface in;

            for (size_t i = 0; i < fixed_vars.size(); i++) // loop over variables in the model (except random) separated by '+' operator
            {
                std::vector<std::string> var_vect = fixed_vars[i]; // get vector of variable(s) needed to be extracted from a data file for i-th effect
                
                evolm::effects_storage s_effect; // container for a parocessed effect
                
                // build full var name
                std::string var_name(var_vect[0]);
                for (size_t j = 1; j < var_vect.size(); j++) // by design, var_vect.size() > 1 assumes '*' operation between variables
                    var_name = var_name + "*" + var_vect[j];

                // look at the extra_effects storage if the variable observ_vars[i] is already processed
                int var_find_res = in.find_value(extra_effects_names, var_name);

                if ( var_find_res != -1 )
                {
                    // variable is in the extra storage
                    random_and_fixed_in_extra_storage.push_back(var_find_res); // we now know where to (not) find observ_vars[i] variable
                    random_and_fixed_effects.push_back(s_effect); // push empty effect to keep the vector size consistent
                    random_and_fixed_effects_names.push_back(var_name);
                    // need to provide some levels description (as strings) for this variable (obtained not from a data file!)
                    std::vector<size_t> shape_ofvar = extra_effects[var_find_res].shape();
                    std::vector<std::string> s_effect_levels(shape_ofvar[1]);
                    for (size_t i = 1; i <= shape_ofvar[1]; i++ )
                        s_effect_levels[i-1] = std::to_string(i);
                    random_and_fixed_effects_levels.push_back(s_effect_levels);
                    continue; // move to the next variable
                }

                random_and_fixed_in_extra_storage.push_back(-1);

                std::string s_effect_name; // parocessed effect name
                std::vector<std::string> s_effect_levels;

                get_effect_from_data(in, var_vect, s_effect, s_effect_name, s_effect_levels); // extract effect from data and do .* if needed

                random_and_fixed_effects.push_back(s_effect); // stack the processed effect
                random_and_fixed_effects_names.push_back(s_effect_name);
                random_and_fixed_effects_levels.push_back(s_effect_levels);

                s_effect.clear();
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::process_fixed_vars(std::vector<std::vector<std::string>> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::process_fixed_vars(std::vector<std::vector<std::string>> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::process_fixed_vars(std::vector<std::vector<std::string>> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::process_random_vars(std::vector<std::vector<std::vector<std::string>>> &lhs_random,
                                           std::vector<std::vector<std::vector<std::string>>> &rhs_random,
                                           std::vector<std::string> &rand_vars_identifiers)
    {
        try
        {
            evolm::IOInterface in;

            if ( lhs_random.size() != rand_vars_identifiers.size() )
                throw std::string("lhs_random.size() != rand_vars_identifiers.size()");

            for (size_t i = 0; i < lhs_random.size(); i++) // loop over the number of random effects in the model
            {
                std::vector<std::string> lhs_levels;
                evolm::effects_storage i_lhs_eff;
                std::string lhs_effect_name;
                bool simple_lhs = false; // indicate that LHS is 1, that means it is not required to call 'BAR' operator

                // check if the random effect provided as a single variable: (1|var)

                if ( lhs_random[i].size() == 1 && rhs_random[i].size() == 1) // if there is no + operator at the LHS and RHS
                {
                    if ( lhs_random[i][0].size() == 1 && rhs_random[i][0].size() == 1) // if there is no .* operator at the LHS and RHS
                    {
                        if ( lhs_random[i][0][0] == "1" ) // finally, if the random effect in the form of (1|var)
                        {
                            simple_lhs = true; // simple case: (1|var)
                            std::string var_name = rhs_random[i][0][0];

                            int var_find_res = in.find_value(extra_effects_names, var_name); // look at the extra_effects storage if the variable observ_vars[i] is already processed

                            if ( var_find_res != -1 ) // variable is in the extra storage
                            {
                                random_and_fixed_in_extra_storage.push_back(var_find_res); // we now know where to (not) find var_name variable
                                lhs_effect_name = lhs_random[i][0][0] + "|" + rhs_random[i][0][0];
                                random_and_fixed_effects.push_back(i_lhs_eff); // stack the processed effect (empty!)
                                random_and_fixed_effects_names.push_back(lhs_effect_name);
                                
                                lhs_levels = extra_effects_levels[var_find_res];
                                random_and_fixed_effects_levels.push_back(lhs_levels);

                                if ( !rand_vars_identifiers[i].empty() )
                                    identifier_to_eff_name[rand_vars_identifiers[i]] = lhs_effect_name;

                                continue; // move to the next random variable
                            }
                        }
                    }
                }

                random_and_fixed_in_extra_storage.push_back(-1);

                if ( data_file_name.empty() )
                    throw std::string("The data file name is empty!");

                // -------------------- LEFT ---------------------
                
                if (!simple_lhs)
                {
                    get_effect_from_data(in, lhs_random[i][0], i_lhs_eff, lhs_effect_name, lhs_levels); // extract very first LHS effect from data and do .* if needed

                    for (size_t j = 1; j < lhs_random[i].size(); j++) // loop over all parsed variables of i-th effect
                    {
                        std::vector<std::string> i_rhs_levels;

                        evolm::effects_storage j_effect; // parocessed effect
                        std::string j_effect_name; // parocessed effect name

                        get_effect_from_data(in, lhs_random[i][j], j_effect, j_effect_name, i_rhs_levels); // extract effect from data and do .* if needed

                        i_lhs_eff.extend_by(j_effect);

                        for ( auto const &v: i_rhs_levels ) // extend_by for vectors of levels
                            lhs_levels.push_back(v);

                        lhs_effect_name = lhs_effect_name + "+" + j_effect_name;

                        j_effect.clear();
                    }
                }
                else
                {
                    lhs_effect_name = "1";
                    lhs_levels.push_back("Intercept");
                }

                // -------------------- RIGHT --------------------

                std::vector<std::string> rhs_levels;

                evolm::effects_storage i_rhs_eff;                
                std::string rhs_effect_name;

                get_effect_from_data(in, rhs_random[i][0], i_rhs_eff, rhs_effect_name, rhs_levels); // extract very first RHS effect from data and do .* if needed

                for (size_t j = 1; j < rhs_random[i].size(); j++) // loop over all parsed variables of i-th effect
                {
                    std::vector<std::string> i_rhs_levels;
                    evolm::effects_storage j_effect; // parocessed effect
                    std::string j_effect_name; // parocessed effect name

                    get_effect_from_data(in, rhs_random[i][j], j_effect, j_effect_name, i_rhs_levels); // extract effect from data and do .* if needed

                    i_rhs_eff.extend_by(j_effect);

                    for ( auto const &v: i_rhs_levels ) // extend_by for vectors of levels
                        rhs_levels.push_back(v);

                    rhs_effect_name = rhs_effect_name + "+" + j_effect_name;

                    j_effect.clear();
                }

                if (!simple_lhs)
                {
                    i_lhs_eff.element_wise_dot(i_rhs_eff); // do final BAR operation
                    element_wise_dot(lhs_levels, rhs_levels);
                }
                else
                {
                    i_lhs_eff = i_rhs_eff;
                    lhs_levels = rhs_levels;
                }

                lhs_effect_name = lhs_effect_name + "|" + rhs_effect_name;
                
                if ( !rand_vars_identifiers[i].empty() ) // we need this because if the same effect is used multiple times they should have a different names
                    lhs_effect_name = lhs_effect_name + "_" + rand_vars_identifiers[i];

                random_and_fixed_effects.push_back(i_lhs_eff); // stack the processed effect
                random_and_fixed_effects_names.push_back(lhs_effect_name);
                random_and_fixed_effects_levels.push_back(lhs_levels);

                if ( !rand_vars_identifiers[i].empty() )
                    identifier_to_eff_name[rand_vars_identifiers[i]] = lhs_effect_name;

                i_lhs_eff.clear();
                i_rhs_eff.clear(); 
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::process_random_vars(std::vector<std::vector<std::vector<std::string>>> &, std::vector<std::vector<std::vector<std::string>>> &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::process_random_vars(std::vector<std::vector<std::vector<std::string>>> &, std::vector<std::vector<std::vector<std::string>>> &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::process_random_vars(std::vector<std::vector<std::vector<std::string>>> &, std::vector<std::vector<std::vector<std::string>>> &, std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::get_effect_from_data(evolm::IOInterface &in_data, std::vector<std::string> &in_var_vect, evolm::effects_storage &out_s, std::string &effect_name, std::vector<std::string> &var_levels)
    {
        try
        {
            std::vector<std::string> lhs_levels;

            // checking if in_var_vect[0] var is in the extra_effects storage
            int var_find_res = in_data.find_value(extra_effects_names, in_var_vect[0]);

            if ( var_find_res != -1 )
            {
                out_s = extra_effects[var_find_res];
                lhs_levels = extra_effects_levels[var_find_res];
            }
            else
            {
                // find in_var_vect[0] in borroow_levels_from,
                // if success - call overloaded fgetvar
                std::unordered_map<std::string, std::string>::const_iterator found = borroow_levels_from.find( in_var_vect[0] );
                if ( found != borroow_levels_from.end() )
                {
                    set_fname(in_data, found->second);

                    std::vector<int> ref_int_var;
                    std::vector<std::string> ref_str_var;

                    in_data.fget_var_levels(found->second, obs_missing_value, ref_int_var, ref_str_var); // get unique ref. (borroow_levels_from) variable levels

                    set_fname(in_data, in_var_vect[0]);

                    in_data.fgetvar(in_var_vect[0], ref_int_var, ref_str_var, obs_missing_value, missing_obs, out_s, lhs_levels); // very first var
                }
                else
                {
                    set_fname(in_data, in_var_vect[0]);

                    in_data.fgetvar(in_var_vect[0], obs_missing_value, missing_obs, out_s, lhs_levels); // very first var
                }
            }

            std::vector<size_t> var_shape = out_s.shape();

            if ( var_shape[0] == 0 || var_shape[1] == 0 )
                throw std::string( "The shape of '" + in_var_vect[0] + "' is 0!" );
            
            std::string var_name = in_var_vect[0];

            for (size_t j = 1; j < in_var_vect.size(); j++) // do .* operation
            {
                std::vector<std::string> rhs_levels;

                evolm::effects_storage s2;

                var_find_res = in_data.find_value(extra_effects_names, in_var_vect[j]);

                if ( var_find_res != -1 )
                {
                    s2 = extra_effects[var_find_res];
                    rhs_levels = extra_effects_levels[var_find_res];
                }
                else
                {
                    // find in_var_vect[0] in borroow_levels_from,
                    // if success - call overloaded fgetvar
                    std::unordered_map<std::string, std::string>::const_iterator found = borroow_levels_from.find( in_var_vect[j] );
                    if ( found != borroow_levels_from.end() )
                    {
                        set_fname(in_data, found->second);

                        std::vector<int> ref_int_var;
                        std::vector<std::string> ref_str_var;
                        in_data.fget_var_levels(found->second, obs_missing_value, ref_int_var, ref_str_var); // get unique ref. (borroow_levels_from) variable levels

                        set_fname(in_data, in_var_vect[j]);

                        in_data.fgetvar(in_var_vect[j], ref_int_var, ref_str_var, obs_missing_value, missing_obs, s2, rhs_levels); // very first var
                    }
                    else
                    {
                        set_fname(in_data, in_var_vect[j]);

                        in_data.fgetvar(in_var_vect[j], obs_missing_value, missing_obs, s2, rhs_levels); // convert var string to data matrix
                    }
                }

                var_shape = s2.shape();

                if ( var_shape[0] == 0 || var_shape[1] == 0 )
                    throw std::string( "The shape of '" + in_var_vect[j] + "' is 0!" );

                out_s.element_wise_dot(s2);

                element_wise_dot(lhs_levels, rhs_levels);

                var_name = var_name + "*" + in_var_vect[j];

                s2.clear();
            }

            effect_name = var_name;
            var_levels = lhs_levels;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::get_effect_from_data(evolm::IOInterface &, std::vector<std::string> &, evolm::effects_storage &, std::string &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::get_effect_from_data(evolm::IOInterface &, std::vector<std::string> &, evolm::effects_storage &, std::string &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::get_effect_from_data(evolm::IOInterface &, std::vector<std::string> &, evolm::effects_storage &, std::string &, std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::eval_model_expr(const std::string &expression,
                                        std::vector<std::string> &observ_vars,
                                        std::vector<std::vector<std::vector<std::string>>> &lhs_randvars,
                                        std::vector<std::vector<std::vector<std::string>>> &rhs_randvars,
                                        std::vector<std::vector<std::string>> &fixed_vars,
                                        std::vector<std::string> &rand_vars_identifiers)
    {
        try
        {
            std::string rhs_exp = expression;

            // process observation variables
            get_obsvars(rhs_exp, observ_vars); // extract and parse LHS to the container obs_vars

            // Checking if the list of obs variables is unique;
            // NOTE: no sorting on obs variables, since it might affect covariance matrix order!!!
            std::vector<std::string> obs_copy(observ_vars);
            make_unique(obs_copy);

            if ( obs_copy.size() != observ_vars.size() )
                throw std::string("The list of provided observation variables is not unique!");

            obs_copy.clear();

            // process random variables
            std::vector<std::string> rand_vars;

            extract_randvars(rhs_exp, rand_vars, rand_vars_identifiers); // extract random effects parts (between brackets) to the container rand_vars

            for (size_t i = 0; i < rand_vars.size(); i++)
            {
                std::vector<std::vector<std::string>> lhs_vars; // processed random variables
                std::vector<std::vector<std::string>> rhs_vars; // processed random variables

                std::string lhs, rhs;
                split_random(rand_vars[i], lhs, rhs);
                parse_expr(lhs, lhs_vars);
                parse_expr(rhs, rhs_vars);


                if (lhs_vars.empty())
                    throw std::string("Incorect syntaxis in the random term: " + rand_vars[i] + " leaing to empty design matrix!");

                if (rhs_vars.empty())
                    throw std::string("Incorect syntaxis in the random term: " + rand_vars[i] + " leaing to empty design matrix!");

                lhs_randvars.push_back(lhs_vars);
                rhs_randvars.push_back(rhs_vars);
            }            
            rand_vars.clear();

            // process fixed variables
            parse_expr(rhs_exp, fixed_vars);
            // -----------------------------
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::eval_model_expr(std::string &, std::vector<std::string> &, std::vector<std::vector<std::vector<std::string>>> &, std::vector<std::vector<std::vector<std::string>>> &, std::vector<std::vector<std::string>> &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::eval_model_expr(std::string &, std::vector<std::string> &, std::vector<std::vector<std::vector<std::string>>> &, std::vector<std::vector<std::vector<std::string>>> &, std::vector<std::vector<std::string>> &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::eval_model_expr(std::string &, std::vector<std::string> &, std::vector<std::vector<std::vector<std::string>>> &, std::vector<std::vector<std::vector<std::string>>> &, std::vector<std::vector<std::string>> &, std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::parse_expr(std::string &in_exp, std::vector<std::vector<std::string>> &out_container)
    {
        try
        {
            std::vector<std::string> minus_fixedvars; // fixed effect variables with minus sign
            std::vector<std::string> plus_fixedvars; // fixed effect variables with plus sign

            get_minusvars(in_exp, minus_fixedvars); // extract fixed effects variables with minus sign into the container minus_fixedvars
            split_colon(minus_fixedvars);
            split_star(minus_fixedvars);
            make_unique(minus_fixedvars);

            get_plusvars(in_exp, plus_fixedvars); // extract fixed effects variables with plus sign into the container minus_fixedvars

            split_colon(plus_fixedvars);

            split_star(plus_fixedvars);

            make_unique(plus_fixedvars);

            exclude_vars(plus_fixedvars, minus_fixedvars);

            check_borrowed_levels(plus_fixedvars);

            split_dot(plus_fixedvars, out_container);

            minus_fixedvars.clear();
            plus_fixedvars.clear();
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::parse_expr(std::string &, std::vector<std::vector<std::string>> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::parse_expr(std::string &, std::vector<std::vector<std::string>> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::parse_expr(std::string &, std::vector<std::vector<std::string>> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::get_obsvars(std::string &expr, std::vector<std::string> &out_vars)
    {
        try
        {
            std::string delimiter("~");

            size_t pos = 0;
            int n_reads = 0;
            std::string lhs_exp;

            while ( (pos = expr.find(delimiter)) != std::string::npos ) // Extract LHS part of the model equation
            {
                lhs_exp = expr.substr(0, pos);
                expr.erase(0, pos + delimiter.length());
                n_reads++;
            }

            if (n_reads > 1)
                throw std::string("There are more then one ~ operators in the model expression!");
            
            if (n_reads < 1)
                throw std::string("There is no required ~ operator in the model expression!");

            split_str(",", lhs_exp, out_vars);

        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::get_obsvars(std::string &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::get_obsvars(std::string &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::get_obsvars(std::string &, std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    int model_parser:: detect_expression_type(const std::string &expr)
    {
        try
        {
            int the_type = 0;

            std::string model_delimiter("~"); // type 1
            std::string name_value_delimiter("="); // type 2
            std::string matrix_delimiter("["); // type 3
            std::string var_identifier("var"); // type 4

            bool type_1 = is_str_in_expr(expr, model_delimiter);
            bool type_2 = is_str_in_expr(expr, name_value_delimiter);
            bool type_3 = is_str_in_expr(expr, matrix_delimiter);

            if ( type_1 && type_2)
                throw std::string("Cannot detect the expression type, two operators (~ and =) in the same expression are not allowed!");
            if ( type_1 && type_3)
                throw std::string("Cannot detect the expression type, two operators (~ and []) in the same expression are not allowed!");
            
            if ( type_1 )
                the_type = 1;

            if ( type_3 )
                the_type = 3;
            else if ( type_2 && !type_3 )
                the_type = 2;

            if (the_type == 2) // check if the string is a variance expression
            {
                size_t find_delim = expr.find(name_value_delimiter);
                std::string lhs = expr.substr(0, find_delim);
                find_delim = 0;
                find_delim = lhs.find(var_identifier);
                if ( find_delim != std::string::npos )
                    the_type = 4;
            }
            
            if (the_type == 0)
                throw std::string("Cannot detect the expression type for: " + expr + " !");

            return the_type;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::detect_expression_type(std::string &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::detect_expression_type(std::string &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::detect_expression_type(std::string &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    bool model_parser::is_str_in_expr(const std::string &expr, const std::string &str)
    {
        try
        {
            std::size_t res = expr.find(str);
            if ( res != std::string::npos )
                return true;
            else
                return false;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::is_str_in_expr(std::string &, const std::string &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::is_str_in_expr(std::string &, const std::string &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::is_str_in_expr(std::string &, const std::string &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::separate_name_file_values(const std::string &expression, std::vector<std::string> &name_value_pair)
    {
        try
        {
            std::string expr = expression;

            std::string delimiter("=");

            size_t pos = 0;
            int n_reads = 0;
            std::string lhs_exp;
            std::string rhs_exp;

            while ( (pos = expr.find(delimiter)) != std::string::npos ) // Extract LHS part of the model equation
            {
                lhs_exp = expr.substr(0, pos);
                rhs_exp = expr.substr(pos + 1);
                expr.erase(0, pos + delimiter.length());
                n_reads++;
            }

            if (n_reads > 1)
                throw std::string("There are more then one = operator in the model expression!");
            
            if (n_reads < 1)
                throw std::string("There is no = operator in the model expression!");

            if ( lhs_exp.empty() )
                throw std::string("The left side of the name-value expression " + expression + " is empty!");

            if ( rhs_exp.empty() )
                throw std::string("The right side of the name-value expression " + expression + " is empty!");

            name_value_pair.push_back(lhs_exp);
            name_value_pair.push_back(rhs_exp);
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::separate_name_file_values(const std::string &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::separate_name_file_values(const std::string &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::separate_name_file_values(const std::string &, std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::extract_randvars(std::string &expr, std::vector<std::string> &out_vars, std::vector<std::string> &out_identifiers)
    {
        // Get variables between brackets ( )
        try
        {
            // extract identifier: '... + (expression)identifier + ...'
            // we need to remove identifier from the expression string at the end

            std::string operator1 = "+";
            std::string operator2 = "-";
            std::string delim_closed = ")";
            size_t search_pos = 0;
            size_t pos1 = 0;
            size_t pos2 = 0;
            size_t pos3 = 0;
            std::string t_str;

            while ( (pos1 = expr.find(delim_closed, search_pos)) != std::string::npos ) // find ')' starting from the search_pos
            {
                pos2 = expr.find(operator1, pos1); // find operator starting from the pos1
                pos3 = expr.find(operator2, pos1); // find operator starting from the pos1

                if ( pos2 < pos3 )
                    t_str = expr.substr( pos1+1, pos2-(pos1+1) ); // '+' is the first operator after the identifier
                if ( pos3 < pos2 )
                    t_str = expr.substr( pos1+1, pos3-(pos1+1) ); // '-' is the first operator after the identifier
                if ( pos2 == pos3 )
                    t_str = expr.substr( pos1+1 ); // this is the end ofexpression

                out_identifiers.push_back(t_str);
                expr.erase( pos1+1, pos2-(pos1+1) ); // erase identifier from the expression string
                search_pos = pos1 + 1;                
            }

            extract_expr_btw_brackets(expr, "(", ")", out_vars);
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::extract_randvars(std::string &, std::vector<std::string> &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::extract_randvars(std::string &, std::vector<std::string> &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::extract_randvars(std::string &, std::vector<std::string> &, std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::extract_expr_btw_brackets(std::string &expr, std::string delim_opened, std::string delim_closed, std::vector<std::string> &out_vars)
    {
        // extract variables between brackets in the string expr
        try
        {
            size_t pos1 = 0;
            size_t pos2 = 0;
            std::string t_str;

            char *cstr_opened = delim_opened.data();
            char *cstr_closed = delim_closed.data();

            while ( (pos1 = expr.find(delim_opened)) != std::string::npos ) // separate multiple cases of (var) of the expr
            {
                pos2 = expr.find(delim_closed);
                t_str = expr.substr(pos1, pos2-pos1+1); // extracts the string: (expression)

                if ( !is_space(t_str) ) // skip if there are white spaces between commas in the in_str
                {
                    t_str.erase(remove(t_str.begin(), t_str.end(), cstr_opened[0]), t_str.end());
                    t_str.erase(remove(t_str.begin(), t_str.end(), cstr_closed[0]), t_str.end());
                    out_vars.push_back(t_str);
                }
                expr.erase(pos1, pos2-pos1+2); // erase extracted var including sign infront of it
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::extract_expr_btw_brackets(std::string &, std::string, std::string, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::extract_expr_btw_brackets(std::string &, std::string, std::string, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::extract_expr_btw_brackets(std::string &, std::string, std::string, std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::split_random(std::string &expr, std::string &out_lhs, std::string &out_rhs)
    {
        try
        {
            std::string delim("|");
            size_t found = expr.find(delim);
            if ( found != std::string::npos )
            {
                out_lhs = expr.substr(0,found);
                out_rhs = expr.substr(found+1);
            }
            else
                throw std::string("Incorect syntaxis in the random term: " + expr);
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::split_random(std::string &, std::string &, std::string &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::split_random(std::string &, std::string &, std::string &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::split_random(std::string &, std::string &, std::string &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::make_unique(std::vector<std::string> &in_vars)
    {
        try
        {
            sort( in_vars.begin(), in_vars.end() );
            in_vars.erase( unique( in_vars.begin(), in_vars.end() ), in_vars.end() );
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::make_unique(std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::make_unique(std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::make_unique(std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::element_wise_dot(std::vector<std::string> &lhs, std::vector<std::string> &rhs)
    {
        try
        {
            if( lhs.empty() || rhs.empty() )
                throw std::string("lhs.empty() or rhs.empty()");

            size_t expected_size = lhs.size() * rhs.size();

            std::vector<std::string> res;

            for ( auto const &v1: rhs )
                for ( auto const &v2: lhs )
                    res.push_back( v2 + " x " + v1 );
            
            if ( res.size() != expected_size )
                throw std::string("res.size() != expected_size");
            
            lhs.clear();

            lhs = res;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::element_wise_dot(std::vector<std::string> &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::element_wise_dot(std::vector<std::string> &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::element_wise_dot(std::vector<std::string> &, std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::exclude_vars(std::vector<std::string> &where_vars, std::vector<std::string> &which_vars)
    {
        try
        {
            for (size_t i = 0; i < which_vars.size(); i++)
                where_vars.erase( std::remove( where_vars.begin(), where_vars.end(), which_vars[i] ), where_vars.end() );
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::exclude_vars(std::vector<std::string> &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::exclude_vars(std::vector<std::string> &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::exclude_vars(std::vector<std::string> &, std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::split_star(std::vector<std::string> &in_vars)
    {
        try
        {
            std::vector<std::string> all_res;

            for (size_t i = 0; i < in_vars.size(); i++)
            {
                std::vector<std::string> res;
                std::string s = in_vars[i];
                split_str("*", s, res);

                if (res.size() > 1)
                {
                    for (size_t j = 0; j < res.size(); j++)
                        all_res.push_back(res[j]);
                    
                    for (size_t j1 = 0; j1 < res.size(); j1++)
                    {
                        for (size_t j2 = 0; j2 <= j1; j2++)
                        {
                            if (j1 != j2)
                                all_res.push_back(res[j2]+"."+res[j1]);
                        }
                    }
                    std::string s2 = res[0]+".";
                    for (size_t j = 1; j < res.size()-1; j++)
                        s2 = s2 + res[j]+".";
                    s2 = s2 + res[res.size()-1];
                    all_res.push_back(s2);
                }
                else
                    all_res.push_back(s);
            }

            in_vars.clear();
            in_vars.shrink_to_fit();

            in_vars = all_res;
            all_res.clear();
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::split_star(std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::split_star(std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::split_star(std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::split_colon(std::vector<std::string> &in_vars)
    {
        try
        {
            std::vector<std::string> all_res;

            for (size_t i = 0; i < in_vars.size(); i++)
            {
                std::vector<std::string> res;
                std::string s = in_vars[i];
                split_str(":", s, res);
                
                if (res.size() > 1)
                {
                    std::string s2 = res[0]+".";;
                    for (size_t j = 1; j < res.size()-1; j++)
                        s2 = s2 + res[j]+".";
                    s2 = s2 + res[res.size()-1];
                    all_res.push_back(s2);
                }
                else
                    all_res.push_back(s);
            }

            in_vars.clear();
            in_vars.shrink_to_fit();

            in_vars = all_res;
            all_res.clear();
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::split_colon(std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::split_colon(std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::split_colon(std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::check_borrowed_levels(std::vector<std::string> &in_vars)
    {
        try
        {
            std::vector<std::string> all_res;

            for (size_t i = 0; i < in_vars.size(); i++)
            {
                std::vector<std::string> res;
                std::string s = in_vars[i];

                extract_expr_btw_brackets(s, "{", "}", res);
                
                if ( !res.empty() )
                {
                    if (res.size() > 1)
                        throw std::string("Multiple expressions within {} brackets are not allowed!");
                    
                    borroow_levels_from[s] = res[0];
                }
                
                all_res.push_back(s);
            }

            in_vars.clear();
            in_vars.shrink_to_fit();

            in_vars = all_res;
            all_res.clear();
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::check_borrowed_levels(std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::check_borrowed_levels(std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::check_borrowed_levels(std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::split_dot(std::vector<std::string> &in_vars, std::vector<std::vector<std::string>> &out_vars)
    {
        try
        {
            for (size_t i = 0; i < in_vars.size(); i++)
            {
                std::vector<std::string> all_res;

                std::vector<std::string> res;
                std::string s = in_vars[i];
                split_str(".", s, res);
                
                if (res.size() > 1)
                {
                    for (size_t j = 0; j < res.size(); j++)
                        all_res.push_back(res[j]);
                }
                else
                    all_res.push_back(s); // there is no "." delimiter found
                
                out_vars.push_back(all_res);
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::split_dot(std::vector<std::string> &, std::vector<std::vector<std::string>> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::split_dot(std::vector<std::string> &, std::vector<std::vector<std::string>> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::split_dot(std::vector<std::string> &, std::vector<std::vector<std::string>> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::get_minusvars(std::string &expr, std::vector<std::string> &minus_vars)
    {
        try
        {
            std::string delim_opened = "-";
            std::string delim_closed = "+-(";
            size_t pos1 = 0;
            size_t pos2 = 0;
            std::string t_str;

            while ( (pos1 = expr.find(delim_opened)) != std::string::npos ) // separate random variables of the expr
            {
                pos2 = expr.find_first_of(delim_closed, pos1+1);

                if ( pos2 > expr.length() )
                    pos2 = expr.length();

                t_str = expr.substr(pos1+1, pos2-pos1-1);

                if ( !is_space(t_str) ) // skip if there are white spaces between commas in the in_str
                {
                    t_str.erase(remove_if(t_str.begin(), t_str.end(), isspace), t_str.end()); // remove white spaces between var and comma in case of they appiar
                    minus_vars.push_back(t_str);
                }

                expr.erase(pos1, pos2-pos1);
            }        
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::get_minusvars(std::string &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::get_minusvars(std::string &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::get_minusvars(std::string &, std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::get_plusvars(std::string &expr, std::vector<std::string> &plus_vars)
    {
        try
        {
            std::string delim_opened = "+";
            std::string delim_closed = "+-(";
            size_t pos1 = 0;
            size_t pos2 = 0;
            std::string t_str;

            while ( (pos1 = expr.find(delim_opened)) != std::string::npos ) // separate random variables of the expr
            {
                pos2 = expr.find_first_of(delim_closed, pos1+1);

                if ( pos2 > expr.length() )
                    pos2 = expr.length();

                t_str = expr.substr(pos1+1, pos2-pos1-1);

                if ( !is_space(t_str) ) // skip if there are white spaces between commas in the in_str
                {
                    t_str.erase(remove_if(t_str.begin(), t_str.end(), isspace), t_str.end()); // remove white spaces between var and comma in case of they appiar
                    plus_vars.push_back(t_str);
                }

                expr.erase(pos1, pos2-pos1);
            }

            // assume this is very last sub-string left in expr!
            if (expr.length() != 0)
            {
                if ( !is_space(expr) ) // skip if there are white spaces between commas in the in_str
                {
                    expr.erase(remove_if(expr.begin(), expr.end(), isspace), expr.end()); // remove white spaces between var and comma in case of they appiar
                    plus_vars.push_back(expr);
                    expr.clear();
                }
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::get_plusvars(std::string &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::get_plusvars(std::string &, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::get_plusvars(std::string &, std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::split_str(std::string delim, std::string in_str, std::vector<std::string> &out_vars)
    {
        try
        {
            size_t pos = 0;
            std::string t_str;

            while ( (pos = in_str.find(delim)) != std::string::npos ) // separate variables of the in_str
            {
                t_str = in_str.substr(0, pos);

                if ( !is_space(t_str) ) // skip if there are white spaces between commas in the in_str
                {
                    t_str.erase(remove_if(t_str.begin(), t_str.end(), isspace), t_str.end()); // remove white spaces between var and comma in case of they appiar
                    out_vars.push_back(t_str);
                }
                in_str.erase(0, pos + delim.length());               
            }

            if (!is_space(in_str)) // very last tokken in the in_str
            {
                in_str.erase(remove_if(in_str.begin(), in_str.end(), isspace), in_str.end());
                out_vars.push_back(in_str);
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::split_str(std::string, std::string, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::split_str(std::string, std::string, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::split_str(std::string, std::string, std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::split_str2(std::string delim, std::string in_str, std::vector<std::string> &out_vars)
    {
        // Keeps white spaces
        try
        {
            size_t pos = 0;
            std::string t_str;

            while ( (pos = in_str.find(delim)) != std::string::npos ) // separate variables of the in_str
            {
                t_str = in_str.substr(0, pos);

                if ( !is_space(t_str) ) // skip if there are white spaces between commas in the in_str
                    out_vars.push_back(t_str);
                
                in_str.erase(0, pos + delim.length());               
            }

            if (!is_space(in_str)) // very last tokken in the in_str
                out_vars.push_back(in_str);
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::split_str2(std::string, std::string, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::split_str2(std::string, std::string, std::vector<std::string> &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::split_str2(std::string, std::string, std::vector<std::string> &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::remove_space(std::string &in_str)
    {
        in_str.erase(remove_if(in_str.begin(), in_str.end(), isspace), in_str.end());
    }
    // -------------------------------------------------------------
    bool model_parser::is_space(std::string s)
    {
        try
        {
            for(size_t i = 0; i < s.length(); i++){
                if(!std::isspace(s[i]))
                    return false;
            }
            return true;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::is_space(std::string)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::is_space(std::string)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::is_space(std::string)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::set_fname(evolm::IOInterface &io_class_instance, const std::string &var)
    {
        try
        {
            // check where is variable, in which file
            bool in_file = false;
            if ( !data_file_name.empty() )
            {
                io_class_instance.set_fname(data_file_name);

                if ( var == "1" )
                    return;
                
                in_file = io_class_instance.is_var_in_header(var);
            }
            
            if ( !in_file ) // if variable is not in main data file or the file is not defined
            {
                if ( !effect_file_name.empty() )
                {
                    io_class_instance.set_fname(effect_file_name);
                    in_file = io_class_instance.is_var_in_header(var);
                }
                else
                    throw std::string("Trying to find variable in supplemental effects file which is not defined. Variable -> " + var);
                
                if ( !in_file )
                    throw std::string("Failed to find variable in provided file(s). Variable -> " + var);
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::set_fname(evolm::IOInterface &, std::string &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::set_fname(evolm::IOInterface &, std::string &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::set_fname(evolm::IOInterface &, std::string &)"<< "\n";
            throw;
        }
    }
    // -------------------------------------------------------------
    void model_parser::print()
    {
        try
        {
            std::vector<int> printed_from_extra_storage;

            // print observations
            for (size_t i = 0; i < observations.size(); i++)
            {
                if ( obs_in_extra_storage[i] == -1 )
                    observations[i].print( observations_names[i] );
                else
                {
                    extra_effects[ obs_in_extra_storage[i] ].print( observations_names[i] );
                    printed_from_extra_storage.push_back( obs_in_extra_storage[i] );
                }
            }

            // print effects
            for (size_t i = 0; i < random_and_fixed_effects.size(); i++)
            {
                if ( random_and_fixed_in_extra_storage[i] == -1 )
                {
                    random_and_fixed_effects[i].print( random_and_fixed_effects_names[i] );
                }
                else
                {
                    extra_effects[ random_and_fixed_in_extra_storage[i] ].print( random_and_fixed_effects_names[i] );
                    printed_from_extra_storage.push_back( random_and_fixed_in_extra_storage[i] );
                }
            }

            // print unprinted extra_effects
            evolm::IOInterface in;
            for (size_t i = 0; i < extra_effects.size(); i++)
            {
                int find_res = in.find_value(printed_from_extra_storage, (int)i);
                if ( find_res == -1 )
                {
                    std::string mesg("extra:" + extra_effects_names[i]);
                    extra_effects[i].print( mesg );
                }
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::print()"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::print()"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::print()"<< "\n";
            throw;
        }
    }

    // -------------------------------------------------------------
    void model_parser::report(const std::string &log_file)
    {
        try
        {
            std::ofstream summary(log_file, std::ofstream::out | std::ofstream::trunc);

            if ( !summary.is_open() )
                throw "Unable to open log file for output!";                

            summary<<'\n';
            summary<<"--------------------------------------"<<'\n';
            summary<<"Model Parser"<<'\n';
            summary<<"--------------------------------------"<<'\n';
            summary<<'\n';
            summary<<"List of variables defined outside the model equation(s):"<<"\n";
            size_t index = 0;
            for (auto const &v: extra_effects_names)
            {
                summary<<"    name & position: "<<v<<" "<<index<<"\n";
                index++;
            }
            summary<<"\n";
            
            summary<<"Model variables: "<<"\n";
            index = 0;
            for (auto const &v: random_and_fixed_effects_names)
            {
                summary<<"    name & position: "<<v<<" "<<index<<"\n";
                index++;
            }
            summary<<"\n";

            summary<<"Observations:"<<"\n";
            for (size_t i = 0; i < observations.size(); i++)
            {
                std::vector<size_t> n_obs;
                if (obs_in_extra_storage[i] == -1) // is not in extra storage
                {
                    n_obs = observations[i].shape();
                    summary<<"    name, position, records: "<<observations_names[i]<<" "<<i<<" "<<n_obs[0]<<"\n";
                }
                else
                {
                    n_obs = observations[obs_in_extra_storage[i]].shape();
                    summary<<"    name, position (in extra), records: "<<observations_names[i]<<" "<<obs_in_extra_storage[i]<<" "<<n_obs[0]<<"\n";
                }
            }
            summary<<"\n";

            summary<<"Effects:"<<"\n";
            for (size_t i = 0; i < random_and_fixed_effects.size(); i++)
                if (random_and_fixed_in_extra_storage[i] == -1) // is not in extra storage
                    summary<<"    name, position, levels: "<<random_and_fixed_effects_names[i]<<" "<<i<<" "<<random_and_fixed_effects_levels[i].size()<<"\n";
                else
                    summary<<"    name, position (in extra), levels: "<<random_and_fixed_effects_names[i]<<" "<<random_and_fixed_in_extra_storage[i]<<" "<<random_and_fixed_effects_levels[i].size()<<"\n";
            summary<<"\n";

            summary<<"Model variance"<<"\n";
            for (size_t i = 0; i < corr_vars_index_in_effects.size(); i++) // loop over all model correlations
            {
                summary<<"correlation "<<i+1<<":"<<'\n';
                for (size_t j = 0; j < corr_vars_index_in_effects[i].size(); j++)
                    summary<<"    variable & position: "<<random_and_fixed_effects_names[ corr_vars_index_in_effects[i][j] ]<<" "<<corr_vars_index_in_effects[i][j]<<"\n";
                summary<<"    matrix & position (in extra): ";
                if (corr_matr_index_in_extra_vars[i][0] >= 0)
                {
                    summary<<extra_effects_names[ corr_matr_index_in_extra_vars[i][0] ]<<" "<<corr_matr_index_in_extra_vars[i][0]<<"\n";
                }
                else
                    summary<<"I "<<corr_matr_index_in_extra_vars[i][0]<<"\n";
                summary<<"    variance & position in extra.: ";
                if (corr_matr_index_in_extra_vars[i][1] >= 0)
                    summary<<extra_effects_names[ corr_matr_index_in_extra_vars[i][1] ]<<" "<<corr_matr_index_in_extra_vars[i][1]<<"\n";
                else
                    summary<<corr_matr_index_in_extra_vars[i][1]<<"\n";
            }
            summary<<"\n";
            summary<<"residual:"<<"\n";
            summary<<"    name & position in extra: "<<extra_effects_names[ corr_matr_index_in_extra_vars[corr_matr_index_in_extra_vars.size()-1][0] ]<<" "<<corr_matr_index_in_extra_vars[corr_matr_index_in_extra_vars.size()-1][0]<<"\n";
            
            summary<<"\n";
            summary<<"Model equation(s), as defined:"<<"\n";
            for (auto const &v: model_definition)
            {
                summary<<"    "<<observations_names[v.first]<<" ~ ";
                std::vector<size_t> eff = v.second;
                summary<<random_and_fixed_effects_names[eff[0]];
                for (size_t i = 1; i < eff.size(); i++)
                    summary<<" + "<<random_and_fixed_effects_names[eff[i]];
                summary<<'\n';
            }

            summary<<"\n";
            summary<<"Model equation(s), as submitted:"<<"\n";
            for (auto  const &m: model_definition)
            {
                size_t obs = m.first;
                std::vector<size_t> model = m.second;
                summary<<"    "<<observations_names[obs]<<" ~ ";
                summary<<random_and_fixed_effects_names[model[0]];
                for (size_t i = 1; i < model.size(); i++)
                    summary<<" + "<<random_and_fixed_effects_names[model[i]];
                summary<<'\n';
            }
            
            summary<<"\n";
            summary<<"Model equation(s), as variables indeces:"<<"\n";
            for (auto  const &m: model_definition)
            {
                size_t obs = m.first;
                std::vector<size_t> model = m.second;
                summary<<"    "<<obs<<" ~ ";
                summary<<model[0];
                for (size_t i = 1; i < model.size(); i++)
                    summary<<" + "<<model[i];
                summary<<'\n';
            }

            summary<<"\n";
            summary<<"Correlations, in variables indeces as submitted:"<<"\n";
            for (auto  const &m: correlated_effects)
            {
                summary<<"    "<<"correlation: (";
                summary<<m[0];
                for (size_t i = 1; i < m.size(); i++)//for (auto const &v: m)
                    summary<<" ,"<<m[i];
                summary<<")"<<'\n';
            }

            summary.close();
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in model_parser::report(const std::string &)"<< "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in model_parser::report(const std::string &)"<< "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in model_parser::report(const std::string &)"<< "\n";
            throw;
        }
    }

    // -------------------------------------------------------------

} // end of namespace evolm
