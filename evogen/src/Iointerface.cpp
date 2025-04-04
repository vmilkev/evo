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
                                        int imarker = (int)cstr[i] - (int)48;
                                        //std::cout<<"inside, imarker = "<<imarker<<", converted = "<<(T)imarker<<"\n";
                                        //markers.push_back((T)cstr[i] - (T)48);
                                        markers.push_back( (T)imarker );
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

                        size_t nbytes = ceil(double(samples) / 4);

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
                catch (...)
                {
                        std::cerr << "Exception in IOInterface::read_bed(size_t, size_t, std::vector<std::vector<int>> &)." << '\n';
                        throw;
                }
        }

        //===============================================================================================================

        void IOInterface::fgetdata(std::vector< std::vector<unsigned long> > &id_map, std::vector<std::vector<bool>> &out)
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
                                throw 10;

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
                                throw 10;

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
        template void IOInterface::fgetdata(std::vector<std::vector<unsigned int>> &out);

        //===============================================================================================================
}
