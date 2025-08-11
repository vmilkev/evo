#include "iointerface.hpp"

namespace evolm
{
    //===============================================================================================================
    IOInterface::IOInterface()
    {
    }
    //===============================================================================================================
    /*void IOInterface::set_missing(float val)
    {
            miss_constant = val;
    }*/
    //===============================================================================================================
    void IOInterface::clear()
    {
        snp_id.clear();
    }
    //===============================================================================================================
    void IOInterface::set_fname(std::string file)
    {
        io_file = file;
    }
    //===============================================================================================================
    void IOInterface::str_parse(std::string &snpStr, std::vector<int> &markers)
    {
        // Parsing SNP string.
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
                throw std::string("Cannot open file " + io_file + " to read!");

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
        catch (const std::string &err)
        {
            std::cerr << "Exception in IOInterface::fgetdata(bool, std::vector<std::vector<int>> &)." << '\n';
            std::cerr << "Error => " << err << '\n';
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
                throw std::string("Cannot open file " + io_file + " to read!");

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
        catch (const std::string &err)
        {
            std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<int>> &)." << '\n';
            std::cerr << "Error => " << err << '\n';
            throw err;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<int>> &)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void IOInterface::get_m_matrix_ascii(int data_format, evolm::matrix<int> &out, std::vector<std::string> &snp_names)
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
            size_t tokens_to_skip = data_format - 1;
            size_t n_rows = 0;
            size_t n_cols = 0;

            size_t loop_counter = 0;

            std::ifstream snpF;
            std::string line;
            std::vector<std::string> data_list;

            snpF.open(io_file.c_str(), std::fstream::in);
            if (!snpF.good())
                throw std::string("Cannot open file " + io_file + " to read!");

            while ( getline(snpF, line) ) // detect number of rows
                n_rows++;
            
            snpF.clear();
            snpF.seekg (0, std::ios::beg);

            while ( getline(snpF, line) ) // loop once to get the number of variants
            {
                std::string delimiter = " ";
                size_t pos = 0;
                std::string token;

                while ((pos = line.find(delimiter)) != std::string::npos && tokens_to_skip != 0)
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

                    if (data_list.size() == tokens_to_skip)
                        break;
                }

                /* get the last element of the string */
                data_list.push_back(line);

                std::vector<int> parsedMarkers;
                std::string strToParse = data_list[tokens_to_skip];
                str_parse(strToParse, parsedMarkers);

                n_cols = parsedMarkers.size();

                data_list.erase(data_list.begin(), data_list.end());

                loop_counter++;

                if (loop_counter == 1)
                    break;
            }

            snpF.clear();
            snpF.seekg (0, std::ios::beg);

            loop_counter = 0;

            out.resize(n_rows, n_cols);

            size_t i_snp = 0;

            while ( getline(snpF, line) )
            {
                std::string delimiter = " ";
                size_t pos = 0;
                std::string token;

                // get SNPs names
                if ( tokens_to_skip != 0 ) // snp ids should be provided
                {
                    if ( (pos = line.find(delimiter)) != std::string::npos )
                    {
                        std::string id = line.substr(0, pos);
                        snp_names.push_back( id );
                    }
                    else
                        throw std::string("line.find(delimiter)) == std::string::npos, while reading beggining of snp row in a file; where the 'space' delimiter is expected!");
                }
                else // there are no snp ids provided
                    snp_names.push_back( std::to_string(i_snp) );
                
                i_snp++;
                // --------------

                while ((pos = line.find(delimiter)) != std::string::npos && tokens_to_skip != 0) // we found delimiter and need to skip some columns
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

                    if (data_list.size() == tokens_to_skip)
                        break;
                }

                /* get the last element of the string */
                data_list.push_back(line);

                std::vector<int> parsedMarkers;
                std::string strToParse = data_list[tokens_to_skip];
                str_parse(strToParse, parsedMarkers);

                for (size_t j = 0; j < parsedMarkers.size(); j++)
                    out(loop_counter, j) = parsedMarkers[j];

                data_list.erase(data_list.begin(), data_list.end());

                loop_counter++;
            }

            snpF.close();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::fgetdata(bool, std::vector<std::vector<int>> &, std::vector<std::string> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &err)
        {
            std::cerr << "Exception in IOInterface::fgetdata(bool, std::vector<std::vector<int>> &, std::vector<std::string> &)." << '\n';
            std::cerr << "Error => " << err << '\n';
            throw err;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::fgetdata(bool, std::vector<std::vector<int>> &, std::vector<std::string> &)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    int IOInterface::detect_data_format_in_snp_txt_file(const std::string &fname)
    {
        try
        {
            std::ifstream snpF;
            std::string line;

            snpF.open(fname.c_str(), std::fstream::in);
            if (!snpF.good())
                throw std::string("Cannot open file " + fname + " to read!");

            getline(snpF, line);
            std::string strToParse = line;
    
            std::string str_token; // token obtained from the original string        
            std::stringstream ss(strToParse);        
            std::vector<int> str_vect;        
            int counts = 0;
            while ( getline(ss, str_token, ' ') && counts < 3)
            {
                str_vect.push_back( std::stoi(str_token) );
                counts++;
            }

            int format_type = 0;
            int missing_variant = 3;

            if ( ( str_vect[0] >= 0 && str_vect[0] <= missing_variant ) && ( str_vect[1] >= 0 && str_vect[1] <= missing_variant ) && ( str_vect[2] >= 0 && str_vect[2] <= missing_variant ) )
                format_type = 1; // no id and no device code columns
            else if ( ( str_vect[0] < 0 || str_vect[0] > missing_variant ) && ( str_vect[1] >= 0 && str_vect[1] <= missing_variant ) && ( str_vect[2] >= 0 && str_vect[2] <= missing_variant ) )
                format_type = 2; // there is id but no device code columns
            else if ( ( str_vect[0] < 0 || str_vect[0] > missing_variant ) && ( str_vect[1] < 0 || str_vect[1] > missing_variant ) && ( str_vect[2] >= 0 && str_vect[2] <= missing_variant ) )
                format_type = 3; // there is id and device code columns
            else
                throw std::string("Cannot detect the right data format in snp txt file " + fname + " !");
        
            snpF.close();

            return format_type;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::detect_data_format_in_snp_txt_file(const std::string &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &err)
        {
            std::cerr << "Exception in IOInterface::detect_data_format_in_snp_txt_file(const std::string &)." << '\n';
            std::cerr << "Error => " << err << '\n';
            throw err;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::detect_data_format_in_snp_txt_file(const std::string &)." << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void IOInterface::fgetvar(const std::string &var_name, float miss_constant, std::vector<bool> &missing_vect, effects_storage &out_var, std::vector<std::string> &unique_levels, std::string &messege)
    {
        /*
                Extract data for a specific variable accessed by the name 'var_name';
                convert it to effects matrix according to a determined type of the data;
                Variable ref_var (usually observation var) is used as reference to track missing records,
                the missing values pointed by miss_constant member.

                File format:

                [header]
                [list of data of different types with "space" delimiter]

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
            missing_constant = miss_constant;

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
                if (var_name == "1") // in the case of the intercept variable
                {
                    size_t num_rows = 0;
                    while (getline(snpF, line))
                        num_rows++;

                    compact_storage<int> fvalues(num_rows, 1);

                    for (size_t i = 0; i < num_rows; i++)
                        fvalues.append(1, i, 0);

                    fvalues.optimize();

                    out_var.set(fvalues);

                    fvalues.clear();

                    unique_levels.push_back("Intercept");

                    return;
                }

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
                        var_types.push_back(get_datatype(token,var_name));
                        data_str.push_back(token);
                    }

                    which_col++;
                }

                // for the very last column:
                if (which_col == (size_t)var_col) // do something if at right column
                {
                    var_types.push_back(get_datatype(line,var_name));
                    data_str.push_back(line);
                }
            }

            snpF.close();

            // ---------------- Detect types ------------------

            std::vector<bool> where_is_missing(var_types.size(), false);
            for (size_t i = 0; i < var_types.size(); i++)
                if (var_types[i] == 5)
                    where_is_missing[i] = true;
            
            int detected_type = define_vartype(var_types, var_name, messege);

            var_types.clear();
            var_types.shrink_to_fit();

            // ----------- Return specific data matrix --------

            if ( data_str.size() != where_is_missing.size() )
                throw std::string("data_str.size() != where_is_missing.size()");

            switch (detected_type)
            {
            case 2: // for continious variable, returns a vector
            {
                compact_storage<float> fvalues(data_str.size(), 1);

                std::vector<float> vect_fvalues;
                str_to_float(data_str, vect_fvalues);

                data_str.clear();
                data_str.shrink_to_fit();

                fvalues.append(vect_fvalues);

                out_var.set(fvalues);

                fvalues.clear();
                vect_fvalues.clear();

                unique_levels.push_back(var_name);

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
                    std::vector<int> ivalues;      // temporal container for integer-type values converted from string-type data
                    str_to_int(data_str, where_is_missing, ivalues); // convert string-type data to integer (which will be processed further as a categorical-type data)

                    n_columns = int_to_cat(ivalues, c_values, unique_levels); // estimate the number of columns in an effect matrix, and convert integers to categorical data

                    ivalues.clear();
                    ivalues.shrink_to_fit();
                }
                else
                    n_columns = int_to_cat(data_str, c_values, unique_levels); // estimate the number of columns in an effect matrix, and convert integers to categorical data

                data_str.clear();
                data_str.shrink_to_fit();

                if (n_columns == 0)
                    throw std::string("The estimated number of columns of the effect matrix is zero!");

                if (c_values.size() == 0)
                    throw std::string("The estimated number of rows of the effect matrix is zero!");

                compact_storage<int> ematrix(c_values.size(), n_columns);

                if ( !missing_vect.empty() && ( c_values.size() != missing_vect.size() ) )
                    throw std::string("c_values.size() != missing_vect.size() while getting variable -> " + var_name);

                cat_to_effect(c_values, missing_vect, ematrix);

                out_var.set(ematrix);

                c_values.clear();
                c_values.shrink_to_fit();
                ematrix.clear();

                break;
            }
            default:
                throw std::string("The variable '" + var_name + "' type has not been determined!");
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::fgetvar(const std::string &, float, std::vector<bool> &, effects_storage &, std::vector<std::string> &, std::string &)" << '\n';
            std::cerr << "Reason => " << e.what() << '\n';
            throw e;
        }
        catch (const std::string err)
        {
            std::cerr << "Exception in IOInterface::fgetvar(const std::string &, float, std::vector<bool> &, effects_storage &, std::vector<std::string> &, std::string &)" << '\n';
            std::cerr << "Reason => " << err << '\n';
            throw err;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::fgetvar(const std::string &, float, std::vector<bool> &, effects_storage &, std::vector<std::string> &, std::string &)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void IOInterface::fgetvar(const std::string &var_name, std::vector<int> &ref_vals_int, std::vector<std::string> &ref_vals_str, float miss_constant, std::vector<bool> &missing_vect, effects_storage &out_var, std::vector<std::string> &unique_levels, std::string &messege)
    {
        try
        {
            missing_constant = miss_constant;

            std::ifstream snpF;
            std::string line;
            std::string delimiter = " ";
            std::vector<std::string> vars_header;

            std::vector<std::string> data_str;
            std::vector<int> var_types;

            std::vector<std::string> data_str2;

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
                if (var_name == "1") // in the case of the intercept variable
                {
                    size_t num_rows = 0;
                    while (getline(snpF, line))
                        num_rows++;

                    compact_storage<int> fvalues(num_rows, 1);

                    for (size_t i = 0; i < num_rows; i++)
                        fvalues.append(1, i, 0);

                    fvalues.optimize();

                    out_var.set(fvalues);

                    fvalues.clear();

                    unique_levels.push_back("Intercept");

                    return;
                }

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
                        var_types.push_back(get_datatype(token,var_name));
                        data_str.push_back(token);
                    }

                    which_col++;
                }

                // for the very last column:
                if (which_col == (size_t)var_col) // do something if at right column
                {
                    var_types.push_back(get_datatype(line,var_name));
                    data_str.push_back(line);
                }

            }

            snpF.close();

            // ---------------- Detect types ------------------

            std::vector<bool> where_is_missing(var_types.size(), false);
            for (size_t i = 0; i < var_types.size(); i++)
                if (var_types[i] == 5)
                    where_is_missing[i] = true;

            int detected_type = define_vartype(var_types, var_name, messege);

            var_types.clear();
            var_types.shrink_to_fit();

            if ( detected_type == 2 )
                throw std::string("In attempt of extracting data for variable " + var_name + " using reference levels: there is only one level for a continious variable.");

            // ----------- Return specific data matrix --------

            if ( data_str.size() != where_is_missing.size() )
                throw std::string("data_str.size() != where_is_missing.size()");

            switch (detected_type)
            {
            case 2: // for continious variable, returns a vector
            {
                compact_storage<float> fvalues(data_str.size(), 1);

                std::vector<float> vect_fvalues;
                str_to_float(data_str, vect_fvalues);

                data_str.clear();
                data_str.shrink_to_fit();

                fvalues.append(vect_fvalues);

                out_var.set(fvalues);

                fvalues.clear();
                vect_fvalues.clear();

                unique_levels.push_back(var_name);

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
                    std::vector<int> ivalues;      // temporal container for integer-type values converted from string-type data
                    
                    str_to_int(data_str, where_is_missing, ivalues); // convert string-type data to integer (which will be processed further as a categorical-type data)

                    if ( ref_vals_int.empty() )
                        throw std::string("The expected integer-type vector of reference levels is empty.");

                    n_columns = int_to_cat(ivalues, ref_vals_int, c_values, unique_levels); // estimate the number of columns in an effect matrix, and convert integers to categorical data

                    ivalues.clear();
                    ivalues.shrink_to_fit();
                }
                else
                {
                    if ( ref_vals_str.empty() )
                        throw std::string("The expected string-type vector of reference levels is empty.");
                    
                    n_columns = int_to_cat(data_str, ref_vals_str, c_values, unique_levels); // estimate the number of columns in an effect matrix, and convert integers to categorical data
                }

                data_str.clear();
                data_str.shrink_to_fit();

                if (n_columns == 0)
                    throw std::string("The estimated number of columns of the effect matrix is zero!");

                if (c_values.size() == 0)
                    throw std::string("The estimated number of rows of the effect matrix is zero!");

                compact_storage<int> ematrix(c_values.size(), n_columns);

                if ( !missing_vect.empty() && ( c_values.size() != missing_vect.size() ) )
                    throw std::string("c_values.size() != missing_vect.size() while getting variable -> " + var_name);

                cat_to_effect(c_values, missing_vect, ematrix);

                out_var.set(ematrix);

                c_values.clear();
                c_values.shrink_to_fit();
                ematrix.clear();

                break;
            }
            default:
                throw std::string("The variable '" + var_name + "' type has not been determined!");
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::fgetvar(const std::string &, std::vector<int> &, std::vector<std::string> &, float, std::vector<bool> &, effects_storage &, std::vector<std::string> &, std::string &)" << '\n';
            std::cerr << "Reason => " << e.what() << '\n';
            throw e;
        }
        catch (const std::string err)
        {
            std::cerr << "Exception in IOInterface::fgetvar(const std::string &, std::vector<int> &, std::vector<std::string> &, float, std::vector<bool> &, effects_storage &, std::vector<std::string> &, std::string &)" << '\n';
            std::cerr << "Reason => " << err << '\n';
            throw err;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::fgetvar(const std::string &, std::vector<int> &, std::vector<std::string> &, float, std::vector<bool> &, effects_storage &, std::vector<std::string> &, std::string &)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void IOInterface::fget_var_levels(const std::string &var_name, float miss_constant, std::vector<int> &out_int, std::vector<std::string> &out_str, std::string &messege)
    {
        /*
                File format:

                [header]
                [list of data of different types with "space" delimiter]

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
            missing_constant = miss_constant;

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
                if (var_name == "1") // in the case of the intercept variable
                    throw std::string("Trying to obtain unique levels for the intercept.");

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
                        var_types.push_back(get_datatype(token,var_name));
                        data_str.push_back(token);
                    }

                    which_col++;
                }

                // for the very last column:
                if (which_col == (size_t)var_col) // do something if at right column
                {
                    var_types.push_back(get_datatype(line,var_name));
                    data_str.push_back(line);
                }
            }

            snpF.close();

            // ---------------- Detect types ------------------

            std::vector<bool> where_is_missing(var_types.size(), false);
            for (size_t i = 0; i < var_types.size(); i++)
                if (var_types[i] == 5)
                    where_is_missing[i] = true;

            int detected_type = define_vartype(var_types, var_name, messege);

            var_types.clear();
            var_types.shrink_to_fit();

            // ----------- Return specific data matrix --------

            switch (detected_type)
            {
            case 2: // for continious variable, returns a vector
            {
                // Do nothing because this type has only one level
                // check other output vectors, if they are empty - variable is continious
                break;
            }
            case 1:
            case 3:
            case 4: // for integer- and string-type variables (categorical), returns a matrix
            {
                if (detected_type == 3)
                {
                    std::vector<int> ivalues;      // temporal container for integer-type values converted from string-type data
                    str_to_int(data_str, where_is_missing, ivalues); // convert string-type data to integer (which will be processed further as a categorical-type data)
                    get_unique_levels(ivalues, out_int);
                }
                else
                    get_unique_levels(data_str, out_str); // estimate the number of columns in an effect matrix, and convert integers to categorical data

                data_str.clear();
                data_str.shrink_to_fit();

                break;
            }
            default:
                throw std::string("The variable '" + var_name + "' type has not been determined!");
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::fget_var_levels(const std::string &, float, std::vector<int> &, std::vector<std::string> &, std::string &)" << '\n';
            std::cerr << "Reason => " << e.what() << '\n';
            throw e;
        }
        catch (const std::string err)
        {
            std::cerr << "Exception in IOInterface::fget_var_levels(const std::string &, float, std::vector<int> &, std::vector<std::string> &, std::string &)" << '\n';
            std::cerr << "Reason => " << err << '\n';
            throw err;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::fget_var_levels(const std::string &, float, std::vector<int> &, std::vector<std::string> &, std::string &)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    bool IOInterface::is_var_in_header(const std::string &var_name)
    {
        bool present = false;

        try
        {
            std::ifstream snpF;
            std::string line;
            std::string delimiter = " ";
            std::vector<std::string> vars_header;

            // -------------- Opening file --------------------

            snpF.open(io_file.c_str(), std::fstream::in);
            if (!snpF.good())
                throw std::string("Cannot open file " + io_file + " for reading!");

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
                if (var_name == "1") // in the case of the intercept variable
                    throw std::string("Trying to obtain unique levels for intercept variable.");
            }
            else
                present = true;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::is_var_in_header(const std::string &)" << '\n';
            std::cerr << "Reason => " << e.what() << '\n';
            throw e;
        }
        catch (const std::string err)
        {
            std::cerr << "Exception in IOInterface::is_var_in_header(const std::string &)" << '\n';
            std::cerr << "Reason => " << err << '\n';
            throw err;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::is_var_in_header(const std::string &)" << '\n';
            throw;
        }

        return present;
    }
    //===============================================================================================================

    size_t IOInterface::int_to_cat(std::vector<int> &ivalues, std::vector<int> &cat_values, std::vector<std::string> &str_levels)
    {
        try
        {
            std::vector<int> in_values(ivalues);

            sort(in_values.begin(), in_values.end());

            if (adjacent_find(in_values.begin(), in_values.end()) != in_values.end())         // if not unique
                in_values.erase(unique(in_values.begin(), in_values.end()), in_values.end()); // make the vector unique

            if ( in_values.back() == missing_int )
                in_values.pop_back(); // remove missing constant to obtain correct number of levels in effect

            // recoding integer values in the 'ivalues' vector to the vcategorical data
            // for missing values: cat_values = -1
            for (size_t i = 0; i < ivalues.size(); i++)
            {
                int i_observation = find_value(in_values, ivalues[i]); // get categorical value (in the range [0,n]) of a specific observation 'ivalues[i]'
                cat_values.push_back(i_observation);                   // for each record write a specific category
            }

            for (auto const &v: in_values)
                str_levels.push_back(std::to_string(v));

            return in_values.size();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::int_to_cat(std::vector<int> &, std::vector<int> &, std::vector<std::string> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in IOInterface::int_to_cat(std::vector<int> &, std::vector<int> &, std::vector<std::string> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::int_to_cat(std::vector<int> &, std::vector<int> &, std::vector<std::string> &)" << '\n';
            throw;
        }
    }
    //===============================================================================================================

    size_t IOInterface::int_to_cat(std::vector<int> &ivalues, std::vector<int> &reference_ivalues, std::vector<int> &cat_values, std::vector<std::string> &str_levels)
    {
        try
        {
            std::vector<int> in_values(reference_ivalues);

            sort(in_values.begin(), in_values.end());

            if (adjacent_find(in_values.begin(), in_values.end()) != in_values.end())         // if not unique
                in_values.erase(unique(in_values.begin(), in_values.end()), in_values.end()); // make the vector unique

            if ( in_values.back() == missing_int )
                in_values.pop_back(); // remove missing constant to obtain correct number of levels in effect

            // recoding integer values in the 'ivalues' vector to the categorical data
            // for missing values: cat_values = -1
            for (size_t i = 0; i < ivalues.size(); i++)
            {
                int i_observation = find_value(in_values, ivalues[i]); // get categorical value (in the range [0,n]) of a specific observation 'ivalues[i]'
                cat_values.push_back(i_observation);                   // for each record write a specific category
            }

            for (auto const &v: in_values)
                str_levels.push_back(std::to_string(v));

            return in_values.size();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::int_to_cat(std::vector<int> &, std::vector<int> &, std::vector<int> &, std::vector<std::string> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in IOInterface::int_to_cat(std::vector<int> &, std::vector<int> &, std::vector<int> &, std::vector<std::string> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::int_to_cat(std::vector<int> &, std::vector<int> &, std::vector<int> &, std::vector<std::string> &)" << '\n';
            throw;
        }
    }
    //===============================================================================================================

    size_t IOInterface::int_to_cat(std::vector<std::string> &ivalues, std::vector<int> &cat_values, std::vector<std::string> &str_levels)
    {
        try
        {
            std::vector<std::string> in_values(ivalues);

            sort(in_values.begin(), in_values.end());

            if (adjacent_find(in_values.begin(), in_values.end()) != in_values.end())         // if not unique
                in_values.erase(unique(in_values.begin(), in_values.end()), in_values.end()); // make the vector unique

            std::string last = in_values.back();

            if ( is_number(last) )
            {
                int last_value = 0;
                last_value = stoi( last );

                if ( last_value == missing_int )
                    in_values.pop_back(); // remove missing constant to obtain correct number of levels in effect
            }

            // recoding integer values in the 'ivalues' vector to the vcategorical data
            // for missing values: cat_values = -1
            for (size_t i = 0; i < ivalues.size(); i++)
            {
                int i_observation = find_value(in_values, ivalues[i]); // get categorical value (in the range [0,n]) of a specific observation 'ivalues[i]'
                cat_values.push_back(i_observation);                   // for each record write a specific category
            }

            str_levels = in_values;

            return in_values.size();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::int_to_cat(std::vector<std::string> &, std::vector<int> &, std::vector<std::string> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in IOInterface::int_to_cat(std::vector<std::string> &, std::vector<int> &, std::vector<std::string> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::int_to_cat(std::vector<std::string> &, std::vector<int> &, std::vector<std::string> &)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    bool IOInterface::is_number(const std::string& s)
    {
        std::string::const_iterator it = s.begin();
        while (it != s.end() && std::isdigit(*it)) ++it;
        return !s.empty() && it == s.end();
    }
    //===============================================================================================================

    size_t IOInterface::int_to_cat(std::vector<std::string> &ivalues, std::vector<std::string> &reference_ivalues, std::vector<int> &cat_values, std::vector<std::string> &str_levels)
    {
        try
        {
            std::vector<std::string> in_values(reference_ivalues);

            sort(in_values.begin(), in_values.end());

            if (adjacent_find(in_values.begin(), in_values.end()) != in_values.end())         // if not unique
                in_values.erase(unique(in_values.begin(), in_values.end()), in_values.end()); // make the vector unique
            std::string last = in_values.back();

            if ( is_number(last) )
            {
                int last_value = 0;
                last_value = stoi( last );

                if ( last_value == missing_int )
                    in_values.pop_back(); // remove missing constant to obtain correct number of levels in effect
            }

            // recoding integer values in the 'ivalues' vector to the vcategorical data
            // for missing values: cat_values = -1
            for (size_t i = 0; i < ivalues.size(); i++)
            {
                int i_observation = find_value(in_values, ivalues[i]); // get categorical value (in the range [0,n]) of a specific observation 'ivalues[i]'
                cat_values.push_back(i_observation);                   // for each record write a specific category
            }

            str_levels = in_values;

            return in_values.size();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::int_to_cat(std::vector<std::string> &, std::vector<std::string> &, std::vector<int> &, std::vector<std::string> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in IOInterface::int_to_cat(std::vector<std::string> &, std::vector<std::string> &, std::vector<int> &, std::vector<std::string> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::int_to_cat(std::vector<std::string> &, std::vector<std::string> &, std::vector<int> &, std::vector<std::string> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    template <typename T>
    void IOInterface::get_unique_levels(std::vector<T> &ivalues, std::vector<T> &unique_ivalues)
    {
        try
        {
            std::vector<T> in_values(ivalues);

            sort(in_values.begin(), in_values.end());

            if (adjacent_find(in_values.begin(), in_values.end()) != in_values.end())         // if not unique
                in_values.erase(unique(in_values.begin(), in_values.end()), in_values.end()); // make the vector unique

            unique_ivalues = in_values;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::get_unique_levels(std::vector<T> &, std::vector<T> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in IOInterface::get_unique_levels(std::vector<T> &, std::vector<T> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::get_unique_levels(std::vector<T> &, std::vector<T> &)" << '\n';
            throw;
        }
    }

    template void IOInterface::get_unique_levels(std::vector<int> &ivalues, std::vector<int> &unique_ivalues);
    template void IOInterface::get_unique_levels(std::vector<std::string> &ivalues, std::vector<std::string> &unique_ivalues);

    //===============================================================================================================

    void IOInterface::cat_to_effect(std::vector<int> &cvalues, std::vector<bool> &missing_vect, compact_storage<int> &ematrix)
    {
        try
        {
            if ( missing_vect.empty() )
            {
                for (size_t i = 0; i < cvalues.size(); i++)
                    if ( cvalues[i] != -1 )
                        ematrix.append(1, i, cvalues[i]);
            }
            else
            {
                for (size_t i = 0; i < cvalues.size(); i++)
                {
                    if ( missing_vect[i] ) // keep the row empty if observation is missed
                        continue;
                    if ( cvalues[i] != -1 )
                        ematrix.append(1, i, cvalues[i]);
                }
            }

            ematrix.optimize();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::cat_to_effect(std::vector<int> &, compact_storage<int> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in IOInterface::cat_to_effect(std::vector<int> &, compact_storage<int> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::cat_to_effect(std::vector<int> &, compact_storage<int> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void IOInterface::str_to_int(std::vector<std::string> &data_str, std::vector<bool> &missing_pos, std::vector<int> &fvalues)
    {
        try
        {
            for (size_t i = 0; i < data_str.size(); i++)
                if ( missing_pos[i] )
                    fvalues.push_back( missing_int );
                else
                    fvalues.push_back(std::stoi(data_str[i]));

            std::vector<int>::iterator result = std::min_element(fvalues.begin(), fvalues.end());
            int min_value = *result;

            if (min_value < 0)
                throw std::string("There is negative value observed for the integer-type (categorical) data in a data file!");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::str_to_int(std::vector<std::string> &, std::vector<bool> &, std::vector<int> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in IOInterface::str_to_int(std::vector<std::string> &, std::vector<bool> &, std::vector<int> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::str_to_int(std::vector<std::string> &, std::vector<bool> &, std::vector<int> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void IOInterface::str_to_float(std::vector<std::string> &data_str, std::vector<float> &fvalues)
    {
        try
        {
            for (size_t i = 0; i < data_str.size(); i++)
                fvalues.push_back(std::stof(data_str[i]));
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::str_to_float(std::vector<std::string> &, std::vector<float> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in IOInterface::str_to_float(std::vector<std::string> &, std::vector<float> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::str_to_float(std::vector<std::string> &, std::vector<float> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    int IOInterface::define_vartype(std::vector<int> &types_vect, const std::string &var_name, std::string &messege)
    {
        try
        {
            sort(types_vect.begin(), types_vect.end());

            if (adjacent_find(types_vect.begin(), types_vect.end()) != types_vect.end())          // if not unique
                types_vect.erase(unique(types_vect.begin(), types_vect.end()), types_vect.end()); // make the vector unique

            if ( types_vect.back() == 5 )
                types_vect.pop_back();

            if (types_vect.size() > 2)
            {
                std::string detected_types = "";
                for (size_t i = 0; i < types_vect.size(); i++)
                {
                    if ( types_vect[i] == 1 )
                        detected_types = detected_types + "Logical; ";
                    if ( types_vect[i] == 2 )
                        detected_types = detected_types + "Real; ";
                    if ( types_vect[i] == 3 )
                        detected_types = detected_types + "Integer; ";
                    if ( types_vect[i] == 4 )
                        detected_types = detected_types + "String; ";
                }
                throw std::string("More than two different data types were detected: " + detected_types + "for the variable: " + var_name +"; in a data file: " + io_file + "!");
            }

            if (types_vect.size() == 0)
                throw std::string("Empty data types vector were detected for the variable: " + var_name +"; in a data file: " + io_file + "!");

            // only one clear type is detected
            if (types_vect.size() == 1)
                return types_vect[0];

            // if some of the integers have floating point representation, treat data vector as floating point
            if ( types_vect[0] == 2 && types_vect[1] == 3 )
            {
                messege = "WARNING: The mixed data types (Integer & Real) were detected for the variable: " + var_name +"; in a data file: " + io_file + ". The data vector will be considered as a Real type.";
                return types_vect[0];
            }

            // if a string types is mixed with floating point number
            if ((types_vect[0] == 4) || (types_vect[1] == 4))
            {
                if ((types_vect[0] == 2) || (types_vect[1] == 2))
                    throw std::string("The mixed types (String & Real) were detected for the variable: " + var_name +"; in a data file: " + io_file + "!");
            }

            // if a string types is mixed with floating integer number
            if ((types_vect[0] == 4) || (types_vect[1] == 4))
            {
                if ((types_vect[0] == 3) || (types_vect[1] == 3))
                    throw std::string("The mixed types (String & Integer) were detected for the variable: " + var_name +"; in a data file: " + io_file + "!");
            }

            // if a string types is mixed with floating point number
            if ((types_vect[0] == 1) || (types_vect[1] == 1))
            {
                if ((types_vect[0] == 2) || (types_vect[1] == 2))
                    throw std::string("The mixed types (Logical & Real) were detected for the variable: " + var_name +"; in a data file: " + io_file + "!");
            }

            // if a string types is mixed with floating integer number
            if ((types_vect[0] == 1) || (types_vect[1] == 1))
            {
                if ((types_vect[0] == 3) || (types_vect[1] == 3))
                    throw std::string("The mixed types (Logical & Integer) were detected for the variable: " + var_name +"; in a data file: " + io_file + "!");
            }

            // if a string types is mixed with floating integer number
            if ((types_vect[0] == 1) || (types_vect[1] == 1))
            {
                if ((types_vect[0] == 4) || (types_vect[1] == 4))
                    throw std::string("The mixed types (Logical & String) were detected for the variable: " + var_name +"; in a data file: " + io_file + "!");
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

    int IOInterface::get_datatype(std::string &str_token, const std::string &var_name)
    {
        try
        {
            std::regex boolean_expr = std::regex("^false|true$");                       // type 1, boolean
            //std::regex float_expr = std::regex("^[+-]?([0-9]+([.][0-9]*)|([.][0-9]+))$"); // type 2, float
            std::regex float_expr = std::regex("^([+-]?[0-9]+([.][0-9]*(?:e[+-]?[0-9]+)?))|([+-]?[0-9]+((?:e[+-]?[0-9]+){1}))$"); // type 2, float
            //std::regex integer_expr = std::regex("^-?[0-9]+$");                         // type 3, integer
            std::regex integer_expr = std::regex("^[+-]?[0-9]+$");                         // type 3, integer
            std::regex string_expr = std::regex("[a-zA-Z_#0-9]+");                      // type 4, string
            
            //std::string missing = std::to_string(missing_constant); This is unsreliable method of conversion !!!
            std::ostringstream ss; ss << missing_constant;
            std::string missing(ss.str()); // type 5, missing data

            int datatype = 0;

            if (std::regex_match(str_token, boolean_expr))
                datatype = 1;
            else if (std::regex_match(str_token, float_expr))
            {
                if ( missing.compare(str_token) == 0 ) // identical if strcmp == 0
                    datatype = 5;
                else
                    datatype = 2;
            }
            else if (std::regex_match(str_token, integer_expr))
            {
                if ( missing.compare(str_token) == 0 ) // identical if strcmp == 0
                    datatype = 5;
                else
                    datatype = 3;
            }
            else if (std::regex_match(str_token, string_expr))
                datatype = 4;
            else
                throw std::string("Cannot detect a data type for the token: " + str_token + "; of the variable: " + var_name + "; appeared in the file: " + io_file + "!");

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
    template int IOInterface::find_value(std::vector<std::int64_t> &where, std::int64_t what);
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
                throw std::string("Cannot open file " + io_file + " to read!");

            while (getline(snpF, line))
            {
                std::string delimiter = " ";
                size_t pos = 0;
                std::string token;

                // if (line.find(delimiter) == std::string::npos)
                //     delimiter = ",";
                
                while ((pos = line.find(delimiter)) != std::string::npos)
                {
                    if (pos == 0)
                        token = delimiter;
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
            std::cerr <<"Reason => "<< e.what() << '\n';
            std::cerr << " => operating on file "<< io_file << " which expected to consist strictly numerical data without header."<<"\n";
            throw e;
        }
        catch (const std::string &err)
        {
            std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<float>> &)." << '\n';
            std::cerr << "Error => " << err << '\n';
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

    template <typename T>
    void IOInterface::fgetdata(std::vector<T> &out, size_t start_row, size_t end_row, size_t col)
    {
        // Reads the specific part of general ASCII formated data

        // Reads file format:
        // [list of numbers with " " delimiter]
        // Example:
        //         12.2 20 51.1
        //         15.5 30 10
        //         21.0 45 562
        //         30.5 50 452
        //         40 61 231

        try
        {
            std::ifstream snpF;
            std::string line;
            std::vector<std::string> data_str;

            snpF.open(io_file.c_str(), std::fstream::in);

            if (!snpF.good())
                throw std::string("Cannot open file " + io_file + " to read!");

            long count_rows = -1;

            while (getline(snpF, line))
            {
                count_rows++;

                if ( count_rows < (long)start_row - 1 ) continue;
                if ( count_rows > (long)end_row - 1 ) break;

                std::string delimiter = " ";
                size_t pos = 0;
                std::string token;

                if (line.find(delimiter) == std::string::npos)
                    delimiter = ",";

                while ((pos = line.find(delimiter)) != std::string::npos)
                {                    
                    if (pos == 0)
                        token = delimiter;
                    else
                        token = line.substr(0, pos);

                    line.erase(0, pos + delimiter.length());

                    if (token.compare(delimiter) == 0) continue;

                    data_str.push_back(token);
                }

                data_str.push_back(line);

                if (data_str.size() > col)
                    throw std::string("The number of columns read from file " + io_file + " at row " + std::to_string(count_rows) + " is " + std::to_string(data_str.size()) + " which is lower than requested!");

                out.push_back(std::stof(data_str[col-1]));

                data_str.erase(data_str.begin(), data_str.end());
            }

            snpF.close();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<float>> &, size_t, size_t, size_t)." << '\n';
            std::cerr <<"Reason => "<< e.what() << '\n';
            std::cerr << " => operating on file "<< io_file << " which expected to consist strictly numerical data without header."<<"\n";
            throw e;
        }
        catch (const std::string &err)
        {
            std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<float>> &, size_t, size_t, size_t)." << '\n';
            std::cerr << "Error => " << err << '\n';
            throw err;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<float>> &, size_t, size_t, size_t)." << '\n';
            throw;
        }
    }

    template void IOInterface::fgetdata(std::vector<double> &out, size_t start_row, size_t end_row, size_t col);
    template void IOInterface::fgetdata(std::vector<float> &out, size_t start_row, size_t end_row, size_t col);

    //===============================================================================================================
    bool IOInterface::is_plink_file(const std::string &fname)
    {
        try
        {
            struct pio_file_t plink_file;

            if (pio_open(&plink_file, fname.c_str()) == PIO_OK)
            {
                pio_close(&plink_file);
                return true;
            }
            else
                return false;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::is_plink_file(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in IOInterface::is_plink_file(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::is_plink_file(const std::string &)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void IOInterface::scale_genotypes(std::vector<float> &values, size_t &nrows, size_t &ncols, std::vector<std::string> &snp_names)
    {
        try
        {
            evolm::matrix<int> M;
            evolm::matrix<float> Z;
            evolm::matrix<size_t> shape_of_Z;

            if (is_plink_file(io_file)) // the pipeline for binary (.bad) plink-formated data
            {
                get_m_matrix_plink(M, snp_names);
            }
            else
            {
                int file_format = detect_data_format_in_snp_txt_file(io_file); // Detect type of SNP text file                
                get_m_matrix_ascii(file_format,M, snp_names);
            }

            make_zmatrix(M,Z); // scalling SNPs
            M.clear();

            shape_of_Z = Z.shape();            
            nrows = shape_of_Z[0];
            ncols = shape_of_Z[1];
            shape_of_Z.clear();            
            Z.to_vector(values);
            Z.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::scale_genotypes(std::vector<T> &, size_t &, size_t &, std::vector<std::string> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in IOInterface::scale_genotypes(std::vector<T> &, size_t &, size_t &, std::vector<std::string> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::scale_genotypes(std::vector<T> &, size_t &, size_t &, std::vector<std::string> &)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void IOInterface::get_m_matrix_plink(evolm::matrix<int> &M, std::vector<std::string> &snp_names)
    {
        try
        {
            // use this elsewhere if needed !
            std::vector<std::int64_t> gmatID;                  // container for the list of G matrix IDs, initiated while reading pre-built G-matrix from file
            std::unordered_map<size_t, char *> samples_id_map; // key: assigned sample id, value: original id (as in .fam file)

            struct pio_file_t plink_file;
            snp_t *snp_buffer;
            size_t sample_id;

            if (pio_open(&plink_file, io_file.c_str()) != PIO_OK)
                throw std::string("Error: Could not open plink file!");

            if (!pio_one_locus_per_row(&plink_file))
                throw std::string("This script requires that snps are rows and samples columns.");

            // get number of samples and SNPs:
            size_t n_samples = pio_num_samples(&plink_file);
            size_t n_variants = pio_num_loci(&plink_file);

            M.resize(n_variants, n_samples);

            size_t locus_id = 0;

            snp_buffer = (snp_t *)malloc(pio_row_size(&plink_file));

            while (pio_next_row(&plink_file, snp_buffer) == PIO_OK)
            {
                for (sample_id = 0; sample_id < n_samples; sample_id++)
                {
                    //struct pio_sample_t *sample = pio_get_sample(&plink_file, sample_id);
                    M(locus_id, sample_id) = (int)snp_buffer[sample_id];
                }
                locus_id++;
            }

            // we need to extract and recode samples IDs, so we ge back to first row and read it onece
            pio_reset_row(&plink_file);
            if (pio_next_row(&plink_file, snp_buffer) == PIO_OK)
            {
                for (sample_id = 0; sample_id < n_samples; sample_id++)
                {
                    struct pio_sample_t *sample = pio_get_sample(&plink_file, sample_id);
                    samples_id_map[sample_id + 1] = sample->iid;
                    gmatID.push_back(std::stol(sample->iid));
                    std::string the_name = sample->iid;
                    snp_names.push_back(the_name);
                }
            }

            gmatID.shrink_to_fit();

            free(snp_buffer);
            pio_close(&plink_file);

            M.transpose();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::get_m_matrix(const std::string &, evolm::matrix<int> &, std::vector<std::string> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in IOInterface::get_m_matrix(const std::string &, evolm::matrix<int> &, std::vector<std::string> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::get_m_matrix(const std::string &, evolm::matrix<int> &, std::vector<std::string> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================
    void IOInterface::make_zmatrix( evolm::matrix<int> &M, evolm::matrix<float> &Z )
    {
        try
        {
            evolm::matrix<size_t> shp; // shape of snp_variants matrix
            shp = M.shape();

            if ( (shp[0] == 0) ||  (shp[1] == 0) )
                throw std::string("snp_variants matrix is empty!");

            size_t snpNum = shp[1];
            size_t sampleNum = shp[0];

            // vectors of SNPs frequences and missed values
            std::vector<float> P(snpNum, 0.0);
            std::vector<int> missed(snpNum, 0);
            std::vector<float> missed2pq(sampleNum, 0.0);

            // map of missed values locations
            std::vector<std::vector<int>> missedLocation;
            for (size_t i = 0; i < snpNum; i++)
                missedLocation.push_back(std::vector<int>());

            // count missing variants, and sums along columns (specific variants)
            for (size_t row = 0; row < sampleNum; row++)
            {
                for (size_t col = 0; col < snpNum; col++)
                {
                    if ( M(row, col) == 3 )
                    {
                        missed[col] += 1;
                        missedLocation[col].push_back(row);
                    }
                    else
                        P[col] += M(row, col);
                }
            }

            // finish to calculate allele frequences, additionally accounting missing values

            /*auto n_threads = std::thread::hardware_concurrency();
            auto block_size = static_cast<unsigned int>(P.size() / (n_threads));

            if (block_size < workload)
            {
                block_size = P.size();
                n_threads = 1;
            }*/

#pragma omp parallel for //schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < P.size(); i++)
            {
                P[i] = P[i] / (2.0 * (float)(sampleNum - missed[i]));
            }

            Z.resize(sampleNum, snpNum);

            for (size_t i = 0; i < sampleNum; i++)
            {
                for (size_t j = 0; j < snpNum; j++)
                {
                    Z(i, j) = (float)M(i, j) - 2.0 * P[j];
                }
            }

            // modify Z matrix, so instead of missing values we put population average (0.0)

#pragma omp parallel for //schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < missedLocation.size(); i++)
            {
                for (size_t j = 0; j < missedLocation[i].size(); j++)
                {
                    Z(missedLocation[i][j], i) = 0.0;
                    missed2pq[missedLocation[i][j]] = missed2pq[missedLocation[i][j]] + 2.0 * P[i] * (1.0 - P[i]);
                }
            }

            float freq = 0.0;
//#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t j = 0; j < P.size(); j++)
            {
                freq += P[j] * (1.0 - P[j]);
            }
            freq = 2.0 * freq;

#pragma omp parallel for //schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < sampleNum; i++)
            {
                missed2pq[i] = sqrt(freq / (freq - missed2pq[i]));
            }

            // After centering, adjust for missing markers for each animal;
            // adjust for sqrt[sum of 2pq over all loci /sum of 2pq over non-missing loci.

#pragma omp parallel for //schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < sampleNum; i++)
            {
                for (size_t j = 0; j < snpNum; j++)
                {
                    Z(i, j) = Z(i, j) * missed2pq[i];
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::make_zmatrix(evolm::matrix<int> &, evolm::matrix<float> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in IOInterface::make_zmatrix(evolm::matrix<int> &, evolm::matrix<float> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::make_zmatrix(evolm::matrix<int> &, evolm::matrix<float> &)" << '\n';
            throw;
        }
    }
    //===========================================================================================
    template <typename T>
    void IOInterface::fread_matrix(const std::string &fname, std::vector<T> &vals, std::vector<size_t> &keys, std::vector<std::int64_t> &ids)
    {
        std::fstream fA;
        //fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        std::string name_suffix(".corbin");
        size_t find = fname.find(name_suffix);
        std::string fname2(fname);
        if( find == std::string::npos ) // there is no suffix
            fname2 = fname + name_suffix;
        
        fA.open(fname2, fA.binary | fA.in);

        if (!fA.is_open())
            throw std::string("IOInterface::fread_matrix(const std::string &, std::vector<T> &, std::vector<size_t> &, std::vector<T2> &): Error while opening a binary file.");

        size_t B[6];

        fA.read( reinterpret_cast<char *>(&B), 6 * sizeof(size_t) ); // 0. reading a storage info

        // info about matrix values
        size_t var_inbytes = B[3];

        if ( var_inbytes != sizeof(T) )
            throw std::string("The matrix data type stored in file is not consistent with type of data trying to read in.");

        size_t vals_size;
        size_t keys_size;
        size_t ids_size;

        fA.read( reinterpret_cast<char *>( &vals_size ), sizeof(size_t) ); // 1. reading size of values
        vals.resize(vals_size);
        fA.read( reinterpret_cast<char *>( vals.data() ), vals_size * var_inbytes ); // 2. reading all values

        fA.read( reinterpret_cast<char *>( &keys_size ), sizeof(size_t) ); // 3. reading size of keys
        keys.resize(keys_size);
        fA.read( reinterpret_cast<char *>( keys.data() ), keys_size * sizeof(size_t) ); // 4. reading all keys

        fA.read( reinterpret_cast<char *>( &ids_size ), sizeof(size_t) ); // 3. reading size of keys
        ids.resize(ids_size);
        fA.read( reinterpret_cast<char *>( ids.data() ), ids_size * sizeof(std::int64_t) ); // 4. reading all keys

        fA.close();
    }
    template void IOInterface::fread_matrix(const std::string &fname, std::vector<float> &vals, std::vector<size_t> &keys, std::vector<std::int64_t> &ids);
    template void IOInterface::fread_matrix(const std::string &fname, std::vector<double> &vals, std::vector<size_t> &keys, std::vector<std::int64_t> &ids);
    //===========================================================================================
    template <typename T>
    void IOInterface::fread_matrix(const std::string &fname, std::vector<T> &vals, std::vector<size_t> &keys, std::vector<std::string> &ids)
    {
        std::fstream fA;
        //fA.exceptions(std::ifstream::failbit | std::ifstream::badbit);

        std::string name_suffix(".corbin");
        size_t find = fname.find(name_suffix);
        std::string fname2(fname);
        if( find == std::string::npos ) // there is no suffix
            fname2 = fname + name_suffix;

        fA.open(fname2, fA.binary | fA.in);

        if (!fA.is_open())
            throw std::string("IOInterface::fread_matrix(const std::string &, std::vector<T> &, std::vector<size_t> &, std::vector<T2> &): Error while opening a binary file.");

        size_t B[6];

        fA.read( reinterpret_cast<char *>(&B), 6 * sizeof(size_t) ); // 0. reading a storage info

        // info about matrix values
        size_t var_inbytes = B[3];

        if ( var_inbytes != sizeof(T) )
            throw std::string("The matrix data type stored in file is not consistent with type of data trying to read in.");

        size_t vals_size;
        size_t keys_size;
        size_t ids_size;

        fA.read( reinterpret_cast<char *>( &vals_size ), sizeof(size_t) ); // 1. reading size of values
        vals.resize(vals_size);
        fA.read( reinterpret_cast<char *>( vals.data() ), vals_size * var_inbytes ); // 2. reading all values

        fA.read( reinterpret_cast<char *>( &keys_size ), sizeof(size_t) ); // 3. reading size of keys
        keys.resize(keys_size);
        fA.read( reinterpret_cast<char *>( keys.data() ), keys_size * sizeof(size_t) ); // 4. reading all keys

        fA.read( reinterpret_cast<char *>( &ids_size ), sizeof(size_t) ); // 3. reading size of keys
        ids.resize(ids_size);
        for (size_t i = 0; i < ids_size; i++)
        {
            size_t str_size;
            std::string str;
            fA.read( reinterpret_cast<char *>( &str_size ), sizeof(size_t) ); // 5. writing size of ids
            str.resize(str_size);
            fA.read( reinterpret_cast<char *>( &str[0] ), str_size ); // 6. reading specific id
            ids[i] = str;
        }

        fA.close();
    }
    template void IOInterface::fread_matrix(const std::string &fname, std::vector<float> &vals, std::vector<size_t> &keys, std::vector<std::string> &ids);
    template void IOInterface::fread_matrix(const std::string &fname, std::vector<double> &vals, std::vector<size_t> &keys, std::vector<std::string> &ids);
    //===========================================================================================
    template <typename T>
    void IOInterface::fread_matrix(const std::string &fname, evolm::matrix<T> &matr, std::vector<std::int64_t> &ids)
    {
        std::string name_suffix(".dmbin");
        size_t find = fname.find(name_suffix);
        std::string fname2(fname);
        if( find == std::string::npos ) // there is no suffix
            fname2 = fname + name_suffix;

        matr.fread(fname2, ids);
    }
    template void IOInterface::fread_matrix(const std::string &fname, evolm::matrix<float> &matr, std::vector<std::int64_t> &ids);
    template void IOInterface::fread_matrix(const std::string &fname, evolm::matrix<double> &matr, std::vector<std::int64_t> &ids);
    //===========================================================================================
    template <typename T>
    void IOInterface::fread_matrix(const std::string &fname, evolm::matrix<T> &matr, std::vector<std::string> &ids)
    {
        std::string name_suffix(".dmbin");
        size_t find = fname.find(name_suffix);
        std::string fname2(fname);
        if( find == std::string::npos ) // there is no suffix
            fname2 = fname + name_suffix;
        
        matr.fread(fname2, ids);
    }
    template void IOInterface::fread_matrix(const std::string &fname, evolm::matrix<float> &matr, std::vector<std::string> &ids);
    template void IOInterface::fread_matrix(const std::string &fname, evolm::matrix<double> &matr, std::vector<std::string> &ids);
    //===========================================================================================
    void IOInterface::fread_matrix_info(const std::string &fname, size_t *info)
    {
        std::fstream fA;

        fA.open(fname, fA.binary | fA.in);

        if (!fA.is_open())
            throw std::string("IOInterface::fread_matrix_info(const std::string &, size_t &): Error while opening a binary file: " + fname + ".");

        fA.read( reinterpret_cast<char *>(info), 6 * sizeof(size_t) ); // reading a storage info

        fA.close();
    }

    //===============================================================================================================
    //===============================================================================================================
    //===============================================================================================================

}
