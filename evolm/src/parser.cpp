#include "parser.hpp"

namespace evolm
{
    // -------------------------------------------------------------

    template <class PType>
    Parser<PType>::Parser()
    {
        try
        {
            int i;
            exp_ptr = NULL;
            for (i = 0; i < NUMVARS; i++)
                vars[i] = (PType)0;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Parser<PType>::sparser()" << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Parser<PType>::sparser()" << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Parser<PType>::sparser()" << "\n";
            throw;
        }
    }

    template class Parser<float>;
    template class Parser<double>;
    template class Parser<int>;

    // -------------------------------------------------------------

    template <class PType>
    PType Parser<PType>::eval_exp(char *exp)
    {
        // Entry point of the Parser
        try
        {
            PType result;

            exp_ptr = exp;

            get_token();
            
            if ( ! *token )
            {
                throw serror(2); // empty expression
                return (PType) 0;
            }

            eval_expp1(result);

            if ( *token )
                throw serror(0); // last token should be the 0 symbol

            return result;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Parser<PType>::eval_exp(char *)" << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Parser<PType>::eval_exp(char *)" << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Parser<PType>::eval_exp(char *)" << "\n";
            throw;
        }
    }

    template float Parser<float>::eval_exp(char *exp);
    template double Parser<double>::eval_exp(char *exp);
    template int Parser<int>::eval_exp(char *exp);

    // -------------------------------------------------------------

    template <class PType>
    void Parser<PType>::eval_expp1(PType &result)
    {
        // Assignment
        try
        {
            int slot;
            char ttok_type;
            char temp_token[80];

            if ( toke_type == VARIABLE )
            {
                // keep the old token
                strcpy(temp_token, token);
                ttok_type = toke_type;

                // get the variable index
                slot = toupper( *token ) - 'A';

                get_token();

                if ( *token != '=' )
                {
                    putback(); // return the current token. Recover old token (assign nothing).
                    strcpy(token, temp_token);
                    toke_type = ttok_type;
                }
                else
                {
                    get_token(); // get the next part of the expression exp
                    eval_expp2(result);
                    vars[ slot ] = result;
                    return;
                }
            }

            eval_expp2(result);
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp1(PType &)" << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp1(PType &)" << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp1(PType &)" << "\n";
            throw;
        }
    }

    template void Parser<float>::eval_expp1(float &result);
    template void Parser<double>::eval_expp1(double &result);
    template void Parser<int>::eval_expp1(int &result);

    // -------------------------------------------------------------

    template <class PType>
    void Parser<PType>::eval_expp2(PType &result)
    {
        // Addition and subtrucction
        try
        {
            char op;
            PType temp;

            eval_expp3(result);

            while ( ( op = *token ) == '+' || op == '-' )
            {
                get_token();
                
                eval_expp3(temp);

                switch (op)
                {
                case '-':
                    result = result - temp;
                    break;
                case '+':
                    result = result + temp;
                    break;
                default:
                    break;
                }
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp2(PType &)" << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp2(PType &)" << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp2(PType &)" << "\n";
            throw;
        }
    }

    template void Parser<float>::eval_expp2(float &result);
    template void Parser<double>::eval_expp2(double &result);
    template void Parser<int>::eval_expp2(int &result);

    // -------------------------------------------------------------

    template <class PType>
    void Parser<PType>::eval_expp3(PType &result)
    {
        // Multiplication and division of two factors
        try
        {
            char op;
            PType temp;

            eval_expp4(result);

            while ( ( op = *token ) == '*' || op == '/' || op == '%' )
            {
                get_token();
                
                eval_expp4(temp);

                switch (op)
                {
                case '*':
                    result = result * temp;
                    break;
                case '/':
                    result = result / temp;
                    break;
                case '%':
                    result = (int) result % (int) temp;
                    break;
                default:
                    break;
                }
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp3(PType &)" << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp3(PType &)" << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp3(PType &)" << "\n";
            throw;
        }
    }

    template void Parser<float>::eval_expp3(float &result);
    template void Parser<double>::eval_expp3(double &result);
    template void Parser<int>::eval_expp3(int &result);

    // -------------------------------------------------------------

    template <class PType>
    void Parser<PType>::eval_expp4(PType &result)
    {
        // Power operation
        try
        {
            int t;
            PType temp;
            PType ex;

            eval_expp5(result);

            if ( *token  == '^' )
            {
                get_token();
                
                eval_expp4(temp);

                ex = result;

                if ( temp == 0.0 )
                {
                    result = (PType) 1;
                    return;
                }

                for ( t = (int)temp - 1; t > 0; --t )
                    result = result * ex;
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp4(PType &)" << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp4(PType &)" << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp4(PType &)" << "\n";
            throw;
        }
    }

    template void Parser<float>::eval_expp4(float &result);
    template void Parser<double>::eval_expp4(double &result);
    template void Parser<int>::eval_expp4(int &result);

    // -------------------------------------------------------------

    template <class PType>
    void Parser<PType>::eval_expp5(PType &result)
    {
        // Unitary operations '+' and '-'
        try
        {
            char op;

            op = 0;

            if ( ( toke_type == DELIMITER ) && ( *token == '+' || *token == '-' ) )
            {
                op = *token;
                get_token();
            }

            eval_expp6(result);

            if ( op  == '-' )
                result = -result;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp5(PType &)" << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp5(PType &)" << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp5(PType &)" << "\n";
            throw;
        }
    }

    template void Parser<float>::eval_expp5(float &result);
    template void Parser<double>::eval_expp5(double &result);
    template void Parser<int>::eval_expp5(int &result);

    // -------------------------------------------------------------

    template <class PType>
    void Parser<PType>::eval_expp6(PType &result)
    {
        // Computes expression consisting brackets
        try
        {
            if ( *token == '(' )
            {
                get_token();
                eval_expp2(result);
                if ( *token != ')' )
                    throw serror(1);
                get_token();
            }
            else
                atom(result);
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp6(PType &)" << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp6(PType &)" << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Parser<PType>::eval_expp6(PType &)" << "\n";
            throw;
        }
    }

    template void Parser<float>::eval_expp6(float &result);
    template void Parser<double>::eval_expp6(double &result);
    template void Parser<int>::eval_expp6(int &result);

    // -------------------------------------------------------------

    template <class PType>
    void Parser<PType>::atom(PType &result)
    {
        // Returns number or variable's value
        try
        {
            switch (toke_type)
            {
            case VARIABLE:
                result = find_var(token);
                get_token();
                return;
            case NUMBER:
                result = (PType) atof(token);
                get_token();
                return;
            default:
                throw serror(0);
            }
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Parser<PType>::atom(PType &)" << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Parser<PType>::atom(PType &)" << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Parser<PType>::atom(PType &)" << "\n";
            throw;
        }
    }

    template void Parser<float>::atom(float &result);
    template void Parser<double>::atom(double &result);
    template void Parser<int>::atom(int &result);

    // -------------------------------------------------------------

    template <class PType>
    void Parser<PType>::putback()
    {
        // Returns token into an input stream
        try
        {
            char *t;

            t = token;

            for (; *t; t++)
                exp_ptr--;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Parser<PType>::putback()" << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Parser<PType>::putback()" << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Parser<PType>::putback()" << "\n";
            throw;
        }
    }

    template void Parser<float>::putback();
    template void Parser<double>::putback();
    template void Parser<int>::putback();

    // -------------------------------------------------------------

    template <class PType>
    std::string Parser<PType>::serror(int error)
    {
        // Error messages
        try
        {
            static std::vector<std::string> msg{
                "Syntax error!",
                "Wrong order of brackets!",
                "The expression is empty"
            };

            return msg[ error ];
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Parser<PType>::serror(int)" << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Parser<PType>::serror(int)" << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Parser<PType>::serror(int)" << "\n";
            throw;
        }
    }

    template std::string Parser<float>::serror(int error);
    template std::string Parser<double>::serror(int error);
    template std::string Parser<int>::serror(int error);

    // -------------------------------------------------------------

    template <class PType>
    void Parser<PType>::get_token()
    {
        // Extract a next token
        try
        {
            char *temp;

            toke_type = 0;
            temp = token;
            *temp = '\0';

            if ( !*exp_ptr )
                return; // The end of expression

            while ( isspace( *exp_ptr ) )
                ++exp_ptr; // Skip delimiter

            if ( strchr( "+-*/%^=()", *exp_ptr ) )
            {
                toke_type = DELIMITER;
                // Move to the next symbol
                *temp++ = *exp_ptr++;
            }
            else if ( isalpha( *exp_ptr ) )
            {
                while ( !isdelim( *exp_ptr ) )
                    *temp++ = *exp_ptr++;
                toke_type = VARIABLE;
            }
            else if ( isdigit( *exp_ptr ) )
            {
                while ( !isdelim( *exp_ptr ) )
                    *temp++ = *exp_ptr++;
                toke_type = NUMBER;
            }

            *temp = '\0';
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Parser<PType>::get_token()" << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Parser<PType>::get_token()" << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Parser<PType>::get_token()" << "\n";
            throw;
        }
    }

    template void Parser<float>::get_token();
    template void Parser<double>::get_token();
    template void Parser<int>::get_token();

    // -------------------------------------------------------------

    template <class PType>
    int Parser<PType>::isdelim(char c)
    {
        // If the parameter 'c' is delimiter, return true
        try
        {
            if ( strchr( " +-*/%^=()", c ) || c == 9 || c == '\r' || c == 0 )
                return 1;
            
            return 0;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Parser<PType>::isdelim(char)" << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Parser<PType>::isdelim(char)" << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Parser<PType>::isdelim(char)" << "\n";
            throw;
        }
    }

    template int Parser<float>::isdelim(char c);
    template int Parser<double>::isdelim(char c);
    template int Parser<int>::isdelim(char c);

    // -------------------------------------------------------------

    template <class PType>
    PType Parser<PType>::find_var(char *s)
    {
        // If the parameter 'c' is delimiter, return true
        try
        {
            if ( !isalpha( *s ) )
            {
                serror(1);
                return (PType) 0;
            }
            
            return vars[ toupper(*token) - 'A' ];
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Parser<PType>::find_var(char *)" << "\n";
            std::cerr << "Reason => " << e << "\n";
            throw e;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Parser<PType>::ifind_var(char *)" << "\n";
            std::cerr << "Reason => " << e.what() << "\n";
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Parser<PType>::find_var(char *)" << "\n";
            throw;
        }
    }

    template float Parser<float>::find_var(char *s);
    template double Parser<double>::find_var(char *s);
    template int Parser<int>::find_var(char *s);

    // -------------------------------------------------------------

} // end of namespace evolm