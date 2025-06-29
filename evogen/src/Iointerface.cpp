#include "Iointerface.hpp"

namespace evogen
{
    IOInterface::IOInterface()
    {
        use_genotypes = false; // by default we use haplotypes
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

    bool IOInterface::is_genotype_data()
    {
        try
        {
            return use_genotypes;
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

    void IOInterface::set_fname(std::string file, bool is_genotypes)
    {
        try
        {
            io_file = file;
            use_genotypes = is_genotypes;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::set_fname(std::string, bool)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::set_fname(std::string, bool)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    template <typename T>
    void IOInterface::str_parse(std::string &snpStr, std::vector<T> &markers)
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
                    int imarker = (int)cstr[i] - (int)48;
                    // std::cout<<"inside, imarker = "<<imarker<<", converted = "<<(T)imarker<<"\n";
                    // markers.push_back((T)cstr[i] - (T)48);
                    markers.push_back((T)imarker);
                }
            }

            delete[] cstr;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::str_parse(std::string &, std::vector<T> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::str_parse(std::string &, std::vector<T> &)." << '\n';
            throw;
        }
    }

    template void IOInterface::str_parse(std::string &snpStr, std::vector<short> &markers);
    template void IOInterface::str_parse(std::string &snpStr, std::vector<bool> &markers);

    //===============================================================================================================

    void IOInterface::fgetdata(std::vector<std::vector<unsigned long>> &id_map, std::vector<std::vector<bool>> &out)
    {
        /*
                Reads file format:
                [ID] [sex] [sere/dame] [list of SNPs with " " delimiter]

                Example:
                18 1 100 76 2 0 1 1 0 0 0 2 1 2
                19 0 100 90 5 0 0 0 0 2 0 2 1 0
                20 0 230 88 1 5 2 1 1 0 0 2 1 2
                21 1 120 90 0 0 2 1 0 1 0 2 2 1
        */

        try
        {
            std::ifstream snpF;
            std::string line;
            std::vector<std::string> data_list;

            snpF.open(io_file.c_str(), std::fstream::in);
            if (!snpF.good())
                throw std::string("Cannot open file: " + io_file);

            while (getline(snpF, line))
            {

                std::string delimiter = " ";
                size_t pos = 0;
                std::string token;
                
                if ( std::find(line.begin(), line.end(), ',') != line.end() )
                    delimiter = ",";

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

                    if (data_list.size() == 3)
                        break;
                }

                /* get the last element of the string */
                data_list.push_back(line);

                std::vector<bool> parsedMarkers;
                std::string strToParse = data_list[3];
                str_parse(strToParse, parsedMarkers);

                out.push_back(parsedMarkers);

                std::string s1, s2, s3;
                s1 = data_list[0];
                s2 = data_list[1];
                s2 = data_list[2];
                std::stringstream stream1(s1);
                std::stringstream stream2(s2);
                std::stringstream stream3(s3);
                unsigned long id1, id2, id3;
                stream1 >> id1;
                stream2 >> id2;
                stream3 >> id3;

                std::vector<unsigned long> ids;
                ids.push_back(stoi(data_list[0])); // id
                ids.push_back(stoi(data_list[1])); // sex
                ids.push_back(stoi(data_list[2])); // sire/dam

                id_map.push_back(ids);

                data_list.erase(data_list.begin(), data_list.end());
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

    void IOInterface::fgetdata(std::vector<std::vector<short>> &out)
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
                throw std::string("Cannot open file for reading!");

            while (getline(snpF, line))
            {
                std::vector<short> parsedMarkers;
                std::string strToParse = line;
                str_parse(strToParse, parsedMarkers);

                out.push_back(parsedMarkers);
            }

            snpF.close();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<short>> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (int err)
        {
            std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<short>> &)." << '\n';
            std::cerr << "Error code => " << err << '\n';
            throw err;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<short>> &)." << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void IOInterface::fgetdata(std::vector<std::vector<bool>> &out)
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
                throw std::string("Cannot open file " + io_file + " for reading!");

                while (getline(snpF, line))
            {
                std::vector<bool> parsedMarkers;
                std::string strToParse = line;
                str_parse(strToParse, parsedMarkers);

                out.push_back(parsedMarkers);
            }

            snpF.close();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<bool>> &)." << '\n';
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (int err)
        {
            std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<bool>> &)." << '\n';
            std::cerr << "Error code => " << err << '\n';
            throw err;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::fgetdata(std::vector<std::vector<bool>> &)." << '\n';
            throw;
        }
    }

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
                throw std::string("Cannot open file: " + io_file);

            while (getline(snpF, line))
            {

                std::string delimiter = " ";
                size_t pos = 0;
                std::string token;

                if ( std::find(line.begin(), line.end(), ',') != line.end() )
                    delimiter = ",";

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
    template void IOInterface::fgetdata(std::vector<std::vector<unsigned int>> &out);

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

            if ( std::find(line.begin(), line.end(), ',') != line.end() )
                delimiter = ",";

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
    void IOInterface::fgetvar(const std::string &var_name, std::vector<double> &vect_values)
    {
        /*
                Extract data for a specific variable accessed by the name 'var_name';

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
            std::ifstream snpF;
            std::string line;
            std::string delimiter = " ";
            std::vector<std::string> vars_header;

            std::vector<std::string> data_str;

            // -------------- Opening file --------------------

            snpF.open(io_file.c_str(), std::fstream::in);
            if (!snpF.good())
                throw std::string("Cannot open file " + io_file + " for reading!");

            // -------------- Read header ---------------------

            getline(snpF, line);

            if ( std::find(line.begin(), line.end(), ',') != line.end() )
                delimiter = ",";

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
                throw std::string("The following variable name is not in the data file header: " + var_name + "!");

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
                        data_str.push_back(token);
                    }

                    which_col++;
                }

                // for the very last column:
                if (which_col == (size_t)var_col) // do something if at right column
                {
                    data_str.push_back(line);
                }
            }

            snpF.close();

            // ----------- Return specific data matrix --------

            for (size_t i = 0; i < data_str.size(); i++)
            {
                if ( is_number(data_str[i]) )
                    vect_values.push_back(std::stod(data_str[i]));
                else
                    throw std::string("The following record: " + data_str[i] + ", is not a number. A number is required for the reading column in the data file!");
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::fgetvar(const std::string &, std::vector<double> &)" << '\n';
            std::cerr << "Reason => " << e.what() << '\n';
            throw e;
        }
        catch (const std::string err)
        {
            std::cerr << "Exception in IOInterface::fgetvar(const std::string &, std::vector<double> &)" << '\n';
            std::cerr << "Reason => " << err << '\n';
            throw err;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::fgetvar(const std::string &, std::vector<double> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================
    
    bool IOInterface::is_number(const std::string &s)
    {
        try
        {
            std::istringstream iss(s);
            double d;
            return iss >> std::noskipws >> d && iss.eof();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in IOInterface::is_number(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in IOInterface::is_number(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in IOInterface::is_number(const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    template <typename T>
    int IOInterface::find_value(std::vector<T> &where, T what)
    {
        try
        {
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
}
