#ifndef modelparser_hpp__
#define modelparser_hpp__

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "effects_storage.hpp"
#include "iointerface.hpp"

namespace evolm
{
    class model_parser
    {
    private:
        // Data:
        std::string data_file_name;
        float obs_missing_value = -999.0f;
        std::string effect_file_name;
        std::vector<bool> missing_obs;

        std::vector<std::string> warning_messages;

        std::vector<std::vector<std::string>> special_corr_vars; // strings of name = value for corr matrices (.corbin), to be processed separately
        
        std::vector<evolm::effects_storage> extra_effects; // effects/matrices provided by a name-value expression
        std::vector<std::string> extra_effects_names; // var names for the effects/matrices provided by a name-value expression
        std::vector<std::vector<std::string>> extra_effects_levels;
        
        std::vector<evolm::effects_storage> observations; // obs container
        std::vector<std::string> observations_names;
        std::vector<int> obs_in_extra_storage;

        std::vector<evolm::effects_storage> random_and_fixed_effects; // random effects container
        std::vector<std::string> random_and_fixed_effects_names;
        std::vector<int> random_and_fixed_in_extra_storage; // -1 if not in extra storage
        std::vector<std::vector<std::string>> random_and_fixed_effects_levels;

        std::vector<std::vector<int>> corr_vars_index_in_effects; // points to a random effect
        std::vector<std::vector<int>> corr_matr_index_in_extra_vars; // where to find correlation matrix for random effect, -1 means I (identity matrix)

        std::unordered_map<std::string, std::string> identifier_to_eff_name; // aliases for random effects, if used
        std::unordered_map<std::string, std::string> borroow_levels_from; // key <- variable which uses levels from other variable -> value; uses levels but not observations
        std::unordered_map<std::string, std::string> effects_with_logical_and; // holds effect of the form: var1 & var2, where key <- var2 (main external effect), var1 -> value (effect from data relaated to observations)

        size_t next_subm_eff = 0;
        size_t next_subm_obs = 0;
        std::vector<std::vector<size_t>> correlated_effects; // for each correlation defines correlated effects (in consecutive order as submitted to solver)
        std::map<size_t, std::vector<size_t>> model_definition; // observations and explaining factors (effects); key: obs index, value: list of eff indices defining the model        

        // Methods:
        void get_obsvars(std::string &expr, std::vector<std::string> &out_vars);
        
        void extract_randvars(std::string &expr,
                              std::vector<std::string> &out_vars,
                              std::vector<std::string> &out_identifiers); // extracts random effects
        void get_minusvars(std::string &expr, std::vector<std::string> &minus_vars); // extracts variables with negative sign
        void get_plusvars(std::string &expr, std::vector<std::string> &plus_vars); // extracts variables with negative sign
        bool is_space(std::string s);
        void split_str(std::string delim, std::string in_str, std::vector<std::string> &out_vars);
        void split_str2(std::string delim, std::string in_str, std::vector<std::string> &out_vars);
        void remove_space(std::string &in_str);
        void split_colon(std::vector<std::string> &in_vars);
        void split_random(std::string &expr, std::string &out_lhs, std::string &out_rhs);
        void split_star(std::vector<std::string> &in_vars);
        void split_dot(std::vector<std::string> &in_vars, std::vector<std::vector<std::string>> &out_vars);
        void make_unique(std::vector<std::string> &in_vars);
        void exclude_vars(std::vector<std::string> &where_vars, std::vector<std::string> &which_vars);
        void parse_expr(std::string &in_exp, std::vector<std::vector<std::string>> &out_container);

        void separate_name_file_values(const std::string &expression, std::vector<std::string> &name_value_pair);

        int detect_expression_type(const std::string &expr); // type 1 (~): model; type 2 (=): type name-value; type 3 ([]): numerical matrix
        bool is_str_in_expr(const std::string &expr, const std::string &str);

        void process_name_value_pair(std::vector<std::string> &name_value_pair);
        void process_name_matrix_pair(std::vector<std::string> &name_value_pair);

        void process_model_variance_expression(const std::string &expr);
        bool consist_snp_substring(std::string &var_name);

        void process_observation_vars(std::vector<std::string> &observ_vars);
        void process_fixed_vars(std::vector<std::vector<std::string>> &fixed_vars);
        void process_random_vars(std::vector<std::vector<std::vector<std::string>>> &lhs_random,
                                 std::vector<std::vector<std::vector<std::string>>> &rhs_random,
                                 std::vector<std::string> &rand_vars_identifiers);
        void process_external_vars();
        void process_effects_with_logical_and(std::string it_first, std::string it_second);
        
        void get_effect_from_data(evolm::IOInterface &in_data,
                                  std::vector<std::string> &in_var_vect,
                                  evolm::effects_storage &out_s,
                                  std::string &effect_name,
                                  std::vector<std::string> &var_levels);

        void eval_model_expr(const std::string &expression,
                             std::vector<std::string> &observ_vars,
                             std::vector<std::vector<std::vector<std::string>>> &lhs_randvars,
                             std::vector<std::vector<std::vector<std::string>>> &rhs_randvars,
                             std::vector<std::vector<std::string>> &fixed_vars,
                             std::vector<std::string> &rand_vars_identifiers);
        void eval_expr(const std::string &expression);
        void extract_expr_btw_brackets(std::string &expr,
                                       std::string delim_opened,
                                       std::string delim_closed,
                                       std::vector<std::string> &out_vars);
        void find_corr_vars_in_effects(std::vector<std::vector<std::string>> &corr_vars);
        void load_corr_file(std::string &cor_matr_file, std::string &cor_variable, std::vector<std::string> &var_levels, compact_storage<float> &corr_storage);
        void load_eff_file(std::string &eff_matr_file, smatrix<float> &s_matr, matrix<float> &d_matr, std::vector<std::int64_t> &i_ids, std::vector<std::string> &s_ids);

        void get_list_of_subm_eff();
        void def_corr_struct();

        void check_borrowed_levels(std::vector<std::string> &in_vars);
        void check_logical_and(std::vector<std::string> &in_vars);

        void element_wise_dot(std::vector<std::string> &lhs, std::vector<std::string> &rhs);

        void set_fname(evolm::IOInterface &io_class_instance, const std::string &var);
        void get_missing_obs();

    public:
        model_parser();
        ~model_parser();
        void process_expression(const std::string &expr);
        void clear();

        // for debugging
        void print();
        void report(const std::string &log_file);

        friend class lmm;
    };

} // end of namespace evolm

#endif // parser_hpp__
