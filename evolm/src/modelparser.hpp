#ifndef modelparser_hpp__
#define modelparser_hpp__

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace evolm
{
    class Modelparser
    {
    private:
        std::string rhs_exp;
        std::vector<std::string> obs_vars; // observations variables
        std::vector<std::string> rand_vars; // random effect variables
        std::vector<std::string> minus_fixedvars; // fixed effect variables with minus sign
        std::vector<std::string> plus_fixedvars; // fixed effect variables with plus sign

        void get_obsvars();
        void get_randvars(); // extracts random effects
        void get_minusvars(); // extracts variables with negative sign
        void get_plusvars(); // extracts variables with negative sign
        bool is_space(std::string s);
        void split_str(std::string delim, std::string in_str, std::vector<std::string> &out_vars);
        void split_colon(std::vector<std::string> &in_vars);
        void split_star(std::vector<std::string> &in_vars);
        void make_unique(std::vector<std::string> &in_vars);
        void exclude_vars(std::vector<std::string> &where_vars, std::vector<std::string> &which_vars);

    public:
        Modelparser();
        ~Modelparser();

        void eval_expr(const std::string expression);
        void get_modelvars(std::vector<std::string> &obs, std::vector<std::string> &fixed, std::vector<std::string> &random);
    };

} // end of namespace evolm

#endif // parser_hpp__
