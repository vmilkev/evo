#include "modelparser.hpp"

namespace evolm
{
    // -------------------------------------------------------------

    Modelparser::Modelparser()
    {
        try
        {
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Modelparser::Modelparser()"
                      << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Modelparser::Modelparser()"
                      << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Modelparser::Modelparser()"
                      << "\n";
            throw;
        }
    }

    // -------------------------------------------------------------

    Modelparser::~Modelparser()
    {
        try
        {
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Modelparser::~Modelparser()"
                      << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Modelparser::~Modelparser()"
                      << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Modelparser::~Modelparser()"
                      << "\n";
            throw;
        }
    }

    // -------------------------------------------------------------

    void Modelparser::get_modelvars(std::vector<std::string> &obs, std::vector<std::string> &fixed, std::vector<std::string> &random)
    {
        try
        {
            obs = obs_vars;
            fixed = plus_fixedvars;
            random = rand_vars; // this is not parsed completely yet
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Modelparser::get_modelvars(std::vector<std::string> &, std::vector<std::string> &, std::vector<std::string> &)"
                      << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Modelparser::get_modelvars(std::vector<std::string> &, std::vector<std::string> &, std::vector<std::string> &)"
                      << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Modelparser::get_modelvars(std::vector<std::string> &, std::vector<std::string> &, std::vector<std::string> &)"
                      << "\n";
            throw;
        }
    }

    // -------------------------------------------------------------

    void Modelparser::eval_expr(const std::string expression)
    {
        try
        {
            rhs_exp = expression;

std::cout<<"Expression: "<< rhs_exp<<"\n\n";

            get_obsvars(); // extract and parse LHS
            make_unique(obs_vars);
            
            std::cout<<"observations: ";
            for(const auto &v: obs_vars)
                std::cout<<v<<"|"<<v.size()<<"; ";
            std::cout<<"\n\n";

            get_randvars(); // extract random effects parts

            std::cout<<"random effects vars: ";
            for(const auto &v: rand_vars)
                std::cout<<v<<"|"<<v.size()<<"; ";
            std::cout<<"\n\n";

            get_minusvars();
            split_colon(minus_fixedvars);
            split_star(minus_fixedvars);
            make_unique(minus_fixedvars);

            std::cout<<"excluded fixed vars: ";
            for(const auto &v: minus_fixedvars)
                std::cout<<v<<"|"<<v.size()<<"; ";
            std::cout<<"\n\n";

            get_plusvars();
            split_colon(plus_fixedvars);
            split_star(plus_fixedvars);
            make_unique(plus_fixedvars);

            exclude_vars(plus_fixedvars, minus_fixedvars);

            std::cout<<"fixed vars: ";
            for(const auto &v: plus_fixedvars)
                std::cout<<v<<"|"<<v.size()<<"; ";
            std::cout<<"\n\n";

            /*std::cout<<"num of RHS vars: "<<rhs_vars.size()<<"\n";
            std::cout<<"RHS vars & size: ";
            for(const auto &v: rhs_vars)
                std::cout<<v<<"|"<<v.size()<<"; ";
            std::cout<<"\n";
std::cout<<"RHS: "<< rhs_exp<<"\n";*/
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Modelparser::eval_expr(std::string &)"
                      << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Modelparser::eval_expr(std::string &)"
                      << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Modelparser::eval_expr(std::string &)"
                      << "\n";
            throw;
        }
    }

    // -------------------------------------------------------------

    void Modelparser::get_obsvars()
    {
        try
        {
            std::string delimiter = "=";

            size_t pos = 0;
            int n_reads = 0;
            std::string lhs_exp;

            while ( (pos = rhs_exp.find(delimiter)) != std::string::npos ) // Extract LHS part of the model equation
            {
                lhs_exp = rhs_exp.substr(0, pos);
                rhs_exp.erase(0, pos + delimiter.length());
                n_reads++;
            }

            if (n_reads > 1)
                throw std::string("There are more then one = operator in the model expression!");
            
            if (n_reads < 1)
                throw std::string("There is no = operator in the model expression!");

            split_str(",", lhs_exp, obs_vars);

        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Modelparser::get_obsvars()"
                      << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Modelparser::get_obsvars()"
                      << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Modelparser::get_obsvars()"
                      << "\n";
            throw;
        }
    }

    // -------------------------------------------------------------

    void Modelparser::get_randvars()
    {
        try
        {
            std::string delim_opened = "(";
            std::string delim_closed = ")";
            size_t pos1 = 0;
            size_t pos2 = 0;
            std::string t_str;

            while ( (pos1 = rhs_exp.find(delim_opened)) != std::string::npos ) // separate random variables of the rhs_exp
            {
                pos2 = rhs_exp.find(delim_closed);
                t_str = rhs_exp.substr(pos1, pos2-pos1+1);

                if ( !is_space(t_str) ) // skip if there are white spaces between commas in the in_str
                {
                    t_str.erase(remove_if(t_str.begin(), t_str.end(), isspace), t_str.end()); // remove white spaces between var and comma in case of they appiar
                    t_str.erase(remove(t_str.begin(), t_str.end(), '('), t_str.end());
                    t_str.erase(remove(t_str.begin(), t_str.end(), ')'), t_str.end());
                    rand_vars.push_back(t_str);
                }
                rhs_exp.erase(pos1-1, pos2-pos1+2); // erase extracted var including sign infront of it

            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Modelparser::get_randvars()"
                      << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Modelparser::get_randvars()"
                      << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Modelparser::get_randvars()"
                      << "\n";
            throw;
        }
    }

    // -------------------------------------------------------------

    void Modelparser::make_unique(std::vector<std::string> &in_vars)
    {
        try
        {
            sort( in_vars.begin(), in_vars.end() );
            in_vars.erase( unique( in_vars.begin(), in_vars.end() ), in_vars.end() );
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Modelparser::make_unique(std::vector<std::string> &)"
                      << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Modelparser::make_unique(std::vector<std::string> &)"
                      << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Modelparser::make_unique(std::vector<std::string> &)"
                      << "\n";
            throw;
        }
    }

    // -------------------------------------------------------------

    void Modelparser::exclude_vars(std::vector<std::string> &where_vars, std::vector<std::string> &which_vars)
    {
        try
        {
            for (size_t i = 0; i < which_vars.size(); i++)
                where_vars.erase( std::remove( where_vars.begin(), where_vars.end(), which_vars[i] ), where_vars.end() );
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Modelparser::exclude_vars(std::vector<std::string> &, std::vector<std::string> &)"
                      << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Modelparser::exclude_vars(std::vector<std::string> &, std::vector<std::string> &)"
                      << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Modelparser::exclude_vars(std::vector<std::string> &, std::vector<std::string> &)"
                      << "\n";
            throw;
        }
    }

    // -------------------------------------------------------------

    void Modelparser::split_star(std::vector<std::string> &in_vars)
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
            std::cerr << "Exception in Modelparser::split_star(std::vector<std::string> &)"
                      << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Modelparser::split_star(std::vector<std::string> &)"
                      << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Modelparser::split_star(std::vector<std::string> &)"
                      << "\n";
            throw;
        }
    }

    // -------------------------------------------------------------

    void Modelparser::split_colon(std::vector<std::string> &in_vars)
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
            std::cerr << "Exception in Modelparser::split_colon(std::vector<std::string> &)"
                      << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Modelparser::split_colon(std::vector<std::string> &)"
                      << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Modelparser::split_colon(std::vector<std::string> &)"
                      << "\n";
            throw;
        }
    }

    // -------------------------------------------------------------

    void Modelparser::get_minusvars()
    {
        try
        {
            std::string delim_opened = "-";
            std::string delim_closed = "+-(";
            size_t pos1 = 0;
            size_t pos2 = 0;
            std::string t_str;

            while ( (pos1 = rhs_exp.find(delim_opened)) != std::string::npos ) // separate random variables of the rhs_exp
            {
                pos2 = rhs_exp.find_first_of(delim_closed, pos1+1);

                if ( pos2 > rhs_exp.length() )
                    pos2 = rhs_exp.length();

                t_str = rhs_exp.substr(pos1+1, pos2-pos1-1);

                if ( !is_space(t_str) ) // skip if there are white spaces between commas in the in_str
                {
                    t_str.erase(remove_if(t_str.begin(), t_str.end(), isspace), t_str.end()); // remove white spaces between var and comma in case of they appiar
                    minus_fixedvars.push_back(t_str);
                }

                rhs_exp.erase(pos1, pos2-pos1);
            }        
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Modelparser::get_minusvars()"
                      << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Modelparser::get_minusvars()"
                      << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Modelparser::get_minusvars()"
                      << "\n";
            throw;
        }
    }

    // -------------------------------------------------------------

    void Modelparser::get_plusvars()
    {
        try
        {
            std::string delim_opened = "+";
            std::string delim_closed = "+-(";
            size_t pos1 = 0;
            size_t pos2 = 0;
            std::string t_str;

            while ( (pos1 = rhs_exp.find(delim_opened)) != std::string::npos ) // separate random variables of the rhs_exp
            {
                pos2 = rhs_exp.find_first_of(delim_closed, pos1+1);

                if ( pos2 > rhs_exp.length() )
                    pos2 = rhs_exp.length();

                t_str = rhs_exp.substr(pos1+1, pos2-pos1-1);

                if ( !is_space(t_str) ) // skip if there are white spaces between commas in the in_str
                {
                    t_str.erase(remove_if(t_str.begin(), t_str.end(), isspace), t_str.end()); // remove white spaces between var and comma in case of they appiar
                    plus_fixedvars.push_back(t_str);
                }

                rhs_exp.erase(pos1, pos2-pos1);
            }

            // assume this is very last sub-string left in rhs_exp!
            if (rhs_exp.length() != 0)
            {
                if ( !is_space(rhs_exp) ) // skip if there are white spaces between commas in the in_str
                {
                    rhs_exp.erase(remove_if(rhs_exp.begin(), rhs_exp.end(), isspace), rhs_exp.end()); // remove white spaces between var and comma in case of they appiar
                    plus_fixedvars.push_back(rhs_exp);
                    rhs_exp.clear();
                }
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Modelparser::get_plusvars()"
                      << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Modelparser::get_plusvars()"
                      << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Modelparser::get_plusvars()"
                      << "\n";
            throw;
        }
    }

    // -------------------------------------------------------------

    void Modelparser::split_str(std::string delim, std::string in_str, std::vector<std::string> &out_vars)
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
            std::cerr << "Exception in Modelparser::split_str(std::string, std::string, std::vector<std::string> &)"
                      << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Modelparser::split_str(std::string, std::string, std::vector<std::string> &)"
                      << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Modelparser::split_str(std::string, std::string, std::vector<std::string> &)"
                      << "\n";
            throw;
        }
    }

    // -------------------------------------------------------------

    bool Modelparser::is_space(std::string s)
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
            std::cerr << "Exception in Modelparser::is_space(std::string)"
                      << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Modelparser::is_space(std::string)"
                      << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Modelparser::is_space(std::string)"
                      << "\n";
            throw;
        }
    }

    // -------------------------------------------------------------

} // end of namespace evolm
