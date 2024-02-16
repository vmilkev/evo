#ifndef parser_hpp__
#define parser_hpp__

#include <iostream>
#include <cstdlib>
#include <cctype>
#include <cstring>
#include <string>
#include <vector>

namespace evolm
{

    enum types
    {
        DELIMITER = 1,
        VARIABLE,
        NUMBER
    };

    const int NUMVARS = 26; // ???

    template <class PType>
    class Parser
    {
    public:
        Parser();
        PType eval_exp(char *exp); // entry point of the parser

    private:
        char *exp_ptr;       // pointer to expression
        char token[80];      // container for a current token
        char toke_type;      // type of a tokken
        PType vars[NUMVARS]; // container of values of variables

        void eval_expp1(PType &result);
        void eval_expp2(PType &result);
        void eval_expp3(PType &result);
        void eval_expp4(PType &result);
        void eval_expp5(PType &result);
        void eval_expp6(PType &result);
        void atom(PType &result);
        void get_token();
        void putback();
        std::string serror(int error);
        PType find_var(char *s);
        int isdelim(char c);
    };

} // end of namespace evolm

#endif // parser_hpp__
