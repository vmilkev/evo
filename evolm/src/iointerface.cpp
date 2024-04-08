#include "iointerface.hpp"

namespace evolm
{
        //===============================================================================================================

        IOInterface::IOInterface()
        {
                try
                {
                        miss_constant = -999.0;
                }
                catch (const std::exception &e)
                {
                        std::cerr << "Exception in IOInterface::IOInterface()" << '\n';
                        std::cerr << e.what() << '\n';
                        throw e;
                }
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::IOInterface()" << '\n';
                        throw;
                }
        }

        //===============================================================================================================

        void IOInterface::set_missing(float val)
        {
                try
                {
                        miss_constant = val;
                }
                catch (const std::exception &e)
                {
                        std::cerr << "Exception in IOInterface::set_missing(float)" << '\n';
                        std::cerr << e.what() << '\n';
                        throw e;
                }
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::set_missing(float)" << '\n';
                        throw;
                }
        }

        //===============================================================================================================

        void IOInterface::set_fname(std::string file)
        {
                try
                {
                        io_file = file;
                }
                catch (const std::exception &e)
                {
                        std::cerr << "Exception in IOInterface::set_fname(std::string)." << '\n';
                        std::cerr << e.what() << '\n';
                        throw e;
                }
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::set_fname(std::string)." << '\n';
                        throw;
                }
        }

        //===============================================================================================================

        void IOInterface::str_parse(std::string &snpStr, std::vector<int> &markers)
        {
                /*
                        Parsing SNP string.
                */

                try
                {
                        size_t sz = snpStr.length() + 1;
                        char *cstr = new char[snpStr.length() + 1];
                        std::strcpy(cstr, snpStr.c_str());
                        for (size_t i = 0; i < sz; i++)
                        {
                                if (isdigit(cstr[i]))
                                {
                                        markers.push_back((int)cstr[i] - (int)48);
                                }
                        }

                        delete[] cstr;
                }
                catch (const std::exception &e)
                {
                        std::cerr << "Exception in IOInterface::str_parse(std::string &, std::vector<int> &)." << '\n';
                        std::cerr << e.what() << '\n';
                        throw e;
                }
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::str_parse(std::string &, std::vector<int> &)." << '\n';
                        throw;
                }
        }

        //===============================================================================================================

        void IOInterface::fgetdata(size_t samples, size_t variants, std::vector<std::vector<int>> &out)
        {
                /*
                    The actual binary data are the nine blocks of 8 bits (a byte) in the center: the first 3 bytes have a special meaning.
                    The first two are fixed, a 'magic number' that enables PLINK to confirm that a BED file is really a BED file.
                    That is, BED files should always start 01101100 00011011.
                    The third byte indicates whether the BED file is in SNP-major or individual-major mode:
                    a value of 00000001 indicates SNP-major (i.e. list all individuals for first SNP, all individuals for second SNP, etc)
                    whereas a value of 00000000 indicates individual-major
                    (i.e. list all SNPs for the first individual, list all SNPs for the second individual, etc).
                    By default, all BED files are SNP-major mode.

                    NOTE: The code implements SNP-major mode.
                */
                try
                {
                        FILE *in_stream = fopen(io_file.c_str(), "rb");

                        if (in_stream == NULL)
                                throw std::string("Cannot open file. Error in IOInterface::fgetdata(size_t, size_t, std::vector<std::vector<int>> &)");

                        size_t nbytes = std::ceil(double(samples) / 4);

                        unsigned char *buffer = (unsigned char *)malloc(nbytes);

                        unsigned char buf_k; // 8 bits

                        std::vector<int> map(4);
                        map[0] = 2;
                        map[1] = -1 /*NA_INTEGER*/;
                        map[2] = 1;
                        map[3] = 0;

                        out.resize(variants, std::vector<int>(samples, 0));

                        //  00 01 10 11         bit level  corresponds to
                        //  0  1  2  3          xij level  corresponds to
                        //  2  NA  1  0         number of copies of first allele in bim file

                        for (size_t i = 0; i < variants; i++)
                        {
                                long int offset = (i)*nbytes + 3;
                                fseek(in_stream, offset, SEEK_SET);
                                fread(buffer, sizeof(unsigned char), nbytes, in_stream);
                                size_t j = 0;
                                for (size_t k = 0; k < nbytes; k++)
                                {
                                        buf_k = buffer[k];
                                        for (int pos = 0; pos < 4; pos++, j++)
                                        {
                                                if (j < samples)
                                                {
                                                        out[i][j] = map[buf_k & 3];
                                                        buf_k = (unsigned char)(buf_k >> 2);
                                                }
                                        }
                                }
                        }

                        free(buffer);

                        fclose(in_stream);
                }
                catch (const std::exception &e)
                {
                        std::cerr << "Exception in IOInterface::read_bed(size_t, size_t, std::vector<std::vector<int>> &)." << '\n';
                        std::cerr << e.what() << '\n';
                        throw e;
                }
                catch (const std::string &e)
                {
                        std::cerr << "Exception in IOInterface::read_bed(size_t, size_t, std::vector<std::vector<int>> &)." << '\n';
                        std::cerr << "Reason: " << e << '\n';
                        throw e;
                }
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::read_bed(size_t, size_t, std::vector<std::vector<int>> &)." << '\n';
                        throw;
                }
        }

        //===============================================================================================================

        void IOInterface::fgetdata(bool includes_id, std::vector<std::vector<int>> &out)
        {
                /*
                        Reads file format:
                        [observation ID] [devise code] [list of SNPs with " " delimiter]

                        Example:
                        18 1000 2 0 1 1 0 0 0 2 1 2
                        19 1000 5 0 0 0 0 2 0 2 1 0
                        20 1000 1 5 2 1 1 0 0 2 1 2
                        21 1000 0 0 2 1 0 1 0 2 2 1
                */

                try
                {
                        size_t loop_counter = 0;

                        std::ifstream snpF;
                        std::string line;
                        std::vector<std::string> data_list;

                        snpF.open(io_file.c_str(), std::fstream::in);
                        if (!snpF.good())
                                throw 10;

                        while (getline(snpF, line))
                        {

                                std::string delimiter = " ";
                                size_t pos = 0;
                                std::string token;

                                while ((pos = line.find(delimiter)) != std::string::npos)
                                {

                                        if (pos == 0)
                                                token = " ";
                                        else
                                                token = line.substr(0, pos);

                                        line.erase(0, pos + delimiter.length());

                                        if (token.compare(delimiter) == 0)
                                        {
                                                continue;
                                        }

                                        data_list.push_back(token);

                                        if (data_list.size() == 2)
                                                break;
                                }

                                /* get the last element of the string */
                                data_list.push_back(line);

                                std::vector<int> parsedMarkers;
                                std::string strToParse = data_list[2];
                                str_parse(strToParse, parsedMarkers);

                                out.push_back(parsedMarkers);

                                snp_id[loop_counter] = stoi(data_list[0]);

                                data_list.erase(data_list.begin(), data_list.end());

                                loop_counter++;
                        }

                        snpF.close();
                }
                catch (const std::exception &e)
                {
                        std::cerr << "Exception in IOInterface::fgetdata(bool, std::vector<std::vector<int>> &)." << '\n';
                        std::cerr << e.what() << '\n';
                        throw e;
                }
                catch (int err)
                {
                        std::cerr << "Exception in IOInterface::fgetdata(bool, std::vector<std::vector<int>> &)." << '\n';
                        std::cerr << "Error code => " << err << '\n';
                        throw err;
                }
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::fgetdata(bool, std::vector<std::vector<int>> &)." << '\n';
                        throw;
                }
        }

        //===============================================================================================================

        void IOInterface::fgetdata(std::vector<std::vector<int>> &out)
        {
                /*
                        Reads file format:
                        [list of SNPs with " " delimiter]

                        Example:
                        2 0 1 1 0 0 0 2 1 2
                        5 0 0 0 0 2 0 2 1 0
                        1 5 2 1 1 0 0 2 1 2
                        0 0 2 1 0 1 0 2 2 1
                */

                try
                {
                        std::ifstream snpF;
                        std::string line;
                        std::vector<std::string> data_list;

                        snpF.open(io_file.c_str(), std::fstream::in);
                        if (!snpF.good())
                                throw 10;

                        while (getline(snpF, line))
                        {
                                std::vector<int> parsedMarkers;
                                std::string strToParse = line;
                                str_parse(strToParse, parsedMarkers);

                                out.push_back(parsedMarkers);
                        }

                        snpF.close();
                }
                catch (const std::exception &e)
                {
                        std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<int>> &)." << '\n';
                        std::cerr << e.what() << '\n';
                        throw e;
                }
                catch (int err)
                {
                        std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<int>> &)." << '\n';
                        std::cerr << "Error code => " << err << '\n';
                        throw err;
                }
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<int>> &)." << '\n';
                        throw;
                }
        }

        //===============================================================================================================

        void IOInterface::fgetvar( const std::string &var_name, Effects &out_var, const std::string &ref_var )
        {
                /*
                        Extract a data for a specific variable accessed by the name 'var_name';
                        cConvert it to an effects matrix according to a determined type of the data;
                        Variable ref_var (usually observation var) is used as reference to track missing records,
                        the missing values pointed by miss_constant member.

                        File format:

                        [header]
                        [list of data of different types with " " delimiter]

                        Example:
                        var_f1 var_i1 var_f2 var_cat var_str
                        12.2   20     51.1   1       aple
                        15.5   30     10     2       plum
                        21.0   45     562    3       aple
                        30.5   50     452    3       plum
                        40     61     231    4       tomato
                        51.3   71     125    2       tomato
                        60.6   80     121    1       plum
                        70.001 91     121    1       aple
                        82.012 10     110.0  4       tomato
                */

                try
                {
                        std::ifstream snpF;
                        std::string line;
                        std::string delimiter = " ";
                        std::vector<std::string> vars_header;

                        std::vector<std::string> data_str;
                        std::vector<int> var_types;

                        // -------------- Opening file --------------------

                        snpF.open(io_file.c_str(), std::fstream::in);
                        if (!snpF.good())
                                throw std::string("Cannot open file for reading!");

                        // -------------- Read header ---------------------

                        getline(snpF, line);

                        size_t pos = 0;
                        std::string token1;

                        while ((pos = line.find(delimiter)) != std::string::npos)
                        {
                                if (pos == 0)
                                        token1 = " ";
                                else
                                        token1 = line.substr(0, pos);

                                line.erase(0, pos + delimiter.length());

                                if (token1.compare(delimiter) == 0)
                                        continue;

                                vars_header.push_back(token1);
                        }

                        vars_header.push_back(line);

                        // ---------- Detect var_name column number ------

                        int var_col = find_value(vars_header, var_name);

                        if (var_col == -1)
                        {
                                std::string s("The following variable name is not in the data file header: ");
                                s = s + var_name;
                                throw s;
                        }
                        // ----------------- Read data ------------------

                        while (getline(snpF, line))
                        {
                                size_t pos = 0;
                                size_t which_col = 0;
                                std::string token;

                                while ((pos = line.find(delimiter)) != std::string::npos)
                                {
                                        if (pos == 0)
                                                token = " ";
                                        else
                                                token = line.substr(0, pos);

                                        line.erase(0, pos + delimiter.length());

                                        if (token.compare(delimiter) == 0)
                                                continue;

                                        if (which_col == (size_t)var_col) // do something if at right column
                                        {
                                                var_types.push_back(get_datatype(token));
                                                data_str.push_back(token);
                                        }

                                        which_col++;
                                }

                                // for the very last column:
                                if (which_col == (size_t)var_col) // do something if at right column
                                {
                                        var_types.push_back(get_datatype(line));
                                        data_str.push_back(line);
                                }
                        }

                        snpF.close();

                        // ---------------- Detect types ------------------

                        int detected_type = define_vartype(var_types);

                        var_types.clear();
                        var_types.shrink_to_fit();

                        // ----------- Return specific data matrix --------

                        switch (detected_type)
                        {
                                case 2: // for continious variable, returns a vector
                                {
                                        std::vector<float> fvalues;
                                        str_tofloat(data_str, fvalues);

                                        data_str.clear();
                                        data_str.shrink_to_fit();

                                        out_var.set(fvalues);

                                        fvalues.clear();
                                        fvalues.shrink_to_fit();

                                        break;
                                }
                                case 1:
                                case 3:
                                case 4: // for integer- and string-type variables (categorical), returns a matrix
                                {
                                        std::vector<int> c_values;
                                        size_t n_columns = 0;

                                        if (detected_type == 3)
                                        {
                                                std::vector<int> ivalues;     // temporal container for integer-type values converted from string-type data
                                                str_toint(data_str, ivalues); // convert string-type data to integer (which will be processed further as a categorical-type data)

                                                n_columns = dim_ofeffmatrix(ivalues, c_values); // estimate the number of columns in an effect matrix, and convert integers to categorical data

                                                ivalues.clear();
                                                ivalues.shrink_to_fit();
                                        }
                                        else
                                                n_columns = dim_ofeffmatrix(data_str, c_values); // estimate the number of columns in an effect matrix, and convert integers to categorical data

                                        data_str.clear();
                                        data_str.shrink_to_fit();

                                        if (n_columns == 0)
                                                throw std::string("The estimated number of columns of the effect matrix is zero!");

                                        if (c_values.size() == 0)
                                                throw std::string("The estimated number of rows of the effect matrix is zero!");

                                        std::vector<std::vector<int>> ematrix(c_values.size(),
                                                                        std::vector<int>(n_columns));

                                        catvect_toeffmatrix(c_values, ematrix);

                                        out_var.set(ematrix);

                                        c_values.clear();
                                        c_values.shrink_to_fit();

                                        break;
                                }
                                default:
                                        throw std::string("A variable type has not been determined!");
                        }
                }
                catch (const std::exception &e)
                {
                        std::cerr << "Exception in IOInterface::fgetvar(const std::string &, Effects &, const std::string &)" << '\n';
                        std::cerr << "Reason => " << e.what() << '\n';
                        throw e;
                }
                catch (const std::string err)
                {
                        std::cerr << "Exception in IOInterface::fgetvar(const std::string &, Effects &, const std::string &)" << '\n';
                        std::cerr << "Reason => " << err << '\n';
                        throw err;
                }
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::fgetvar(const std::string &, Effects &, const std::string &)" << '\n';
                        throw;
                }
        }

        //===============================================================================================================

        template <typename T>
        size_t IOInterface::dim_ofeffmatrix(std::vector<T> &ivalues, std::vector<int> &cat_values)
        {
                try
                {
                        std::vector<T> in_values(ivalues);

                        sort(in_values.begin(), in_values.end());

                        if (adjacent_find(in_values.begin(), in_values.end()) != in_values.end())             // if not unique
                                in_values.erase(unique(in_values.begin(), in_values.end()), in_values.end()); // make the vector unique

                        // recoding integer values in the 'ivalues' vector to the vcategorical data
                        for (size_t i = 0; i < ivalues.size(); i++)
                        {
                                int i_observation = find_value(in_values, ivalues[i]); // get categorical value (in the range [0,n]) of a specific observation 'ivalues[i]'
                                cat_values.push_back(i_observation);
                        }

                        return in_values.size();
                }
                catch (const std::exception &e)
                {
                        std::cerr << "Exception in IOInterface::dim_ofeffmatrix(std::vector<T> &, std::vector<int> &)" << '\n';
                        std::cerr << e.what() << '\n';
                        throw;
                }
                catch (const std::string &e)
                {
                        std::cerr << "Exception in IOInterface::dim_ofeffmatrix(std::vector<T> &, std::vector<int> &)" << '\n';
                        std::cerr << "Reason: " << e << '\n';
                        throw;
                }
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::dim_ofeffmatrix(std::vector<T> &, std::vector<int> &)" << '\n';
                        throw;
                }
        }

        template size_t IOInterface::dim_ofeffmatrix(std::vector<int> &ivalues, std::vector<int> &cat_values);
        template size_t IOInterface::dim_ofeffmatrix(std::vector<std::string> &ivalues, std::vector<int> &cat_values);

        //===============================================================================================================

        void IOInterface::catvect_toeffmatrix(std::vector<int> &cvalues, std::vector<std::vector<int>> &ematrix)
        {
                try
                {
                        for (size_t i = 0; i < cvalues.size(); i++)
                                ematrix[i][cvalues[i]] = 1;
                }
                catch (const std::exception &e)
                {
                        std::cerr << "Exception in IOInterface::catvect_toeffmatrix(std::vector<int> &, std::vector<std::vector<int>> &)" << '\n';
                        std::cerr << e.what() << '\n';
                        throw;
                }
                catch (const std::string &e)
                {
                        std::cerr << "Exception in IOInterface::catvect_toeffmatrix(std::vector<int> &, std::vector<std::vector<int>> &)" << '\n';
                        std::cerr << "Reason: " << e << '\n';
                        throw;
                }
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::catvect_toeffmatrix(std::vector<int> &, std::vector<std::vector<int>> &)" << '\n';
                        throw;
                }
        }

        //===============================================================================================================

        void IOInterface::str_toint(std::vector<std::string> &data_str, std::vector<int> &fvalues)
        {
                try
                {
                        for (size_t i = 0; i < data_str.size(); i++)
                                fvalues.push_back(std::stoi(data_str[i]));

                        std::vector<int>::iterator result = std::min_element(fvalues.begin(), fvalues.end());
                        int min_value = *result;

                        if (min_value < 0)
                                throw std::string("There is negative value observed for the integer-type (categorical) data in a data file!");
                }
                catch (const std::exception &e)
                {
                        std::cerr << "Exception in IOInterface::str_toint(std::vector<std::string> &, std::vector<int> &)" << '\n';
                        std::cerr << e.what() << '\n';
                        throw;
                }
                catch (const std::string &e)
                {
                        std::cerr << "Exception in IOInterface::str_toint(std::vector<std::string> &, std::vector<int> &)" << '\n';
                        std::cerr << "Reason: " << e << '\n';
                        throw;
                }
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::str_toint(std::vector<std::string> &, std::vector<int> &)" << '\n';
                        throw;
                }
        }

        //===============================================================================================================

        void IOInterface::str_tofloat(std::vector<std::string> &data_str, std::vector<float> &fvalues)
        {
                try
                {
                        for (size_t i = 0; i < data_str.size(); i++)
                                fvalues.push_back(std::stof(data_str[i]));
                }
                catch (const std::exception &e)
                {
                        std::cerr << "Exception in IOInterface::str_tofloat(std::vector<std::string> &, std::vector<float> &)" << '\n';
                        std::cerr << e.what() << '\n';
                        throw;
                }
                catch (const std::string &e)
                {
                        std::cerr << "Exception in IOInterface::str_tofloat(std::vector<std::string> &, std::vector<float> &)" << '\n';
                        std::cerr << "Reason: " << e << '\n';
                        throw;
                }
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::str_tofloat(std::vector<std::string> &, std::vector<float> &)" << '\n';
                        throw;
                }
        }

        //===============================================================================================================

        int IOInterface::define_vartype(std::vector<int> &types_vect)
        {
                try
                {
                        sort(types_vect.begin(), types_vect.end());

                        if (adjacent_find(types_vect.begin(), types_vect.end()) != types_vect.end())              // if not unique
                                types_vect.erase(unique(types_vect.begin(), types_vect.end()), types_vect.end()); // make the vector unique

                        if (types_vect.size() > 2)
                                throw std::string("More than two different data types detected for a variable in a data file!");

                        if (types_vect.size() == 0)
                                throw std::string("Empty data types vector!");

                        // only one clear type is detected
                        if (types_vect.size() == 1)
                                return types_vect[0];

                        // if some of the integers hafe floating point representation, treat data vector as floating point
                        if ((types_vect[0] == 2 && types_vect[1] == 3) || (types_vect[1] == 2 && types_vect[0] == 3))
                                return types_vect[0];

                        // if a string types is mixed with floating point number
                        if ((types_vect[0] == 4) || (types_vect[1] == 4))
                        {
                                if ((types_vect[0] == 2) || (types_vect[1] == 2))
                                        throw std::string("The mixed types (string & floating point number) is detected at the variable's column in the file!");
                        }

                        // if a string types is mixed with floating integer number
                        if ((types_vect[0] == 4) || (types_vect[1] == 4))
                        {
                                if ((types_vect[0] == 3) || (types_vect[1] == 3))
                                        throw std::string("The mixed types (string & integer number) is detected at the variable's column in the file!");
                        }

                        // if a string types is mixed with floating point number
                        if ((types_vect[0] == 1) || (types_vect[1] == 1))
                        {
                                if ((types_vect[0] == 2) || (types_vect[1] == 2))
                                        throw std::string("The mixed types (boolean & floating point number) is detected at the variable's column in the file!");
                        }

                        // if a string types is mixed with floating integer number
                        if ((types_vect[0] == 1) || (types_vect[1] == 1))
                        {
                                if ((types_vect[0] == 3) || (types_vect[1] == 3))
                                        throw std::string("The mixed types (boolean & integer number) is detected at the variable's column in the file!");
                        }

                        // if a string types is mixed with floating integer number
                        if ((types_vect[0] == 1) || (types_vect[1] == 1))
                        {
                                if ((types_vect[0] == 4) || (types_vect[1] == 4))
                                        throw std::string("The mixed types (boolean & string) is detected at the variable's column in the file!");
                        }

                        return 0;
                }
                catch (const std::exception &e)
                {
                        std::cerr << "Exception in IOInterface::define_vartype( std::vector<int> &)" << '\n';
                        std::cerr << e.what() << '\n';
                        throw;
                }
                catch (const std::string &e)
                {
                        std::cerr << "Exception in IOInterface::define_vartype( std::vector<int> &)" << '\n';
                        std::cerr << "Reason: " << e << '\n';
                        throw;
                }
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::define_vartype( std::vector<int> &)" << '\n';
                        throw;
                }
        }

        //===============================================================================================================

        int IOInterface::get_datatype(std::string &str_token)
        {
                try
                {
                        std::regex boolean_expr = std::regex("^false|true$");  // type 1
                        std::regex float_expr = std::regex("^[+-]?([0-9]+([.][0-9]*)|[.][0-9]+)$");//("^\\d+\\.\\d+$");   // type 2
                        std::regex integer_expr = std::regex("^\\d+$");        // type 3
                        std::regex string_expr = std::regex("[a-zA-Z_#0-9]+"); // type 4

                        int datatype = 0;
//std::cout<<"str_token: "<<str_token<<"\n";
                        if (std::regex_match(str_token, boolean_expr))
                                datatype = 1;
                        else if (std::regex_match(str_token, float_expr))
                                datatype = 2;
                        else if (std::regex_match(str_token, integer_expr))
                                datatype = 3;
                        else if (std::regex_match(str_token, string_expr))
                                datatype = 4;
                        else
                                throw std::string("Undetected data type!");

                        return datatype;
                }
                catch (const std::exception &e)
                {
                        std::cerr << "Exception in IOInterface::get_datatype(std::string &)" << '\n';
                        std::cerr << e.what() << '\n';
                        throw;
                }
                catch (const std::string &e)
                {
                        std::cerr << "Exception in IOInterface::get_datatype(std::string &)" << '\n';
                        std::cerr << "Reason: " << e << '\n';
                        throw;
                }
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::get_datatype(std::string &)" << '\n';
                        throw;
                }
        }

        //===============================================================================================================

        template <typename T>
        int IOInterface::find_value(std::vector<T> &where, T what)
        {
                try
                {
                        // std::vector<T>::iterator it;
                        auto it = find(where.begin(), where.end(), what);
                        if (it != where.end())
                                return it - where.begin();
                        else
                                return -1;
                }
                catch (const std::exception &e)
                {
                        std::cerr << "Exception in IOInterface::find_value(std::vector<T> &, T)" << '\n';
                        std::cerr << e.what() << '\n';
                        throw;
                }
                catch (const std::string &e)
                {
                        std::cerr << "Exception in IOInterface::find_value(std::vector<T> &, T)" << '\n';
                        std::cerr << "Reason: " << e << '\n';
                        throw;
                }
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::find_value(std::vector<T> &, T)" << '\n';
                        throw;
                }
        }

        template int IOInterface::find_value(std::vector<std::string> &where, std::string what);
        template int IOInterface::find_value(std::vector<int> &where, int what);
        template int IOInterface::find_value(std::vector<float> &where, float what);
        template int IOInterface::find_value(std::vector<double> &where, double what);

        //===============================================================================================================

        template <typename T>
        void IOInterface::fgetdata(std::vector<std::vector<T>> &out)
        {
                /*
                        Reads file format:
                        [list of numbers with " " delimiter]

                        Example:
                                12.2 20 51.1
                                15.5 30 10
                                21.0 45 562
                                30.5 50 452
                                40 61 231
                */

                try
                {
                        std::ifstream snpF;
                        std::string line;
                        std::vector<T> data_list;

                        snpF.open(io_file.c_str(), std::fstream::in);

                        if (!snpF.good())
                        {
                                /*std::cout << " good()=" << snpF.good();
                                std::cout << " eof()=" << snpF.eof();
                                std::cout << " fail()=" << snpF.fail();
                                std::cout << " bad()=" << snpF.bad();*/
                                throw 10;
                        }

                        while (getline(snpF, line))
                        {
                                std::string delimiter = " ";
                                size_t pos = 0;
                                std::string token;

                                while ((pos = line.find(delimiter)) != std::string::npos)
                                {

                                        if (pos == 0)
                                                token = " ";
                                        else
                                                token = line.substr(0, pos);

                                        line.erase(0, pos + delimiter.length());

                                        if (token.compare(delimiter) == 0)
                                        {
                                                continue;
                                        }

                                        data_list.push_back(std::stof(token));
                                }

                                data_list.push_back(std::stof(line));

                                out.push_back(data_list);

                                data_list.erase(data_list.begin(), data_list.end());
                        }

                        snpF.close();
                }
                catch (const std::exception &e)
                {
                        std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<float>> &)." << '\n';
                        std::cerr << e.what() << '\n';
                        throw e;
                }
                catch (int err)
                {
                        std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<float>> &)." << '\n';
                        std::cerr << "Error code => " << err << '\n';
                        throw err;
                }
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<float>> &)." << '\n';
                        throw;
                }
        }

        template void IOInterface::fgetdata(std::vector<std::vector<double>> &out);
        template void IOInterface::fgetdata(std::vector<std::vector<float>> &out);
        template void IOInterface::fgetdata(std::vector<std::vector<size_t>> &out);
        template void IOInterface::fgetdata(std::vector<std::vector<bool>> &out);

        //===============================================================================================================
}
