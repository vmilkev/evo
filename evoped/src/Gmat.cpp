#include "Gmat.hpp"

namespace evoped
{
    //===============================================================================================================

    Gmat::Gmat()
    {
    }

    //===============================================================================================================

    void Gmat::get_matrix(evolm::matrix<double>& arr)
    {
        try
        {
            arr = G;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::get_matrix()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::get_matrix()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::get_matrix()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::bin_write()
    {
        try
        {
            G.fwrite();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::bin_write()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::bin_write()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::bin_write()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::bin_read()
    {
        try
        {
            G.fread();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::bin_read()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::bin_read()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::bin_read()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::make_matrix(const std::string &fname)
    {
        try
        {
            read_snp(fname);
            make_zmatrix();
            snp_map.clear();
            make_matrix();
            Z.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::make_matrix(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::make_matrix(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::make_matrix(const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::invert_matrix()
    {
        try
        {
            size_t sz = G.size();
            double *A;
            A = (double *)malloc(sz * sizeof(double));
            if (A == NULL)
            {
                free(A);
                throw std::string("Memory allocation error: invert_matrix()");
            }
            for (size_t i = 0; i < sz; i++)
            {
                A[i] = G[i];
            }

            lapack_int info = 0;
            lapack_int row = 9;
            lapack_int col = 9;
            int matrix_order = LAPACK_ROW_MAJOR;

            lapack_int *ipiv;

            ipiv = (lapack_int *)malloc(row * sizeof(lapack_int));

            if (ipiv == NULL)
            {
                free(ipiv);
                throw std::string("Memory allocation error: invert_matrix()");
            }

            for (lapack_int i = 0; i < (row); i++)
                ipiv[i] = 1;

            info = LAPACKE_dgetrf(matrix_order, row, col, A, col, ipiv);
        
            if (info != 0)
            {
                free(ipiv);
                throw std::string("Error during computation of the LU factorization of a general m-by-n matrix: invert_matrix()");
            }

            info = LAPACKE_dgetri(matrix_order, row, A, row, ipiv);

            if (info != 0)
            {
                free(ipiv);
                throw std::string("Error during computation the inverse of an LU-factored general matrix. matrix<T>::inv_rec(double *, MKL_INT, MKL_INT)");
            }
            free(ipiv);

            for (size_t i = 0; i < sz; i++)
            {
                G[i] = A[i];
            }

            //G.invert();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::invert_matrix()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::invert_matrix()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::invert_matrix()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::invert_matrix(std::vector<std::int64_t>& core_id)
    {
        // sparse inverse (APY)

        try
        {
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::invert_matrix(std::vector<std::int64_t>&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::invert_matrix(std::vector<std::int64_t>&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::invert_matrix(std::vector<std::int64_t>&)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    Gmat::~Gmat()
    {
        try
        {
            clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::~Gmat()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::~Gmat()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::~Gmat()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::clear()
    {
        try
        {
            G.clear();
            Z.clear();
            gmatID.clear();
            gmatID.shrink_to_fit();
            snp_map.clear();
            anim_id_map.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::clear()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::clear()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::clear()" << '\n';
            throw;
        }
    }

    //===============================================================================================================
    template <typename T>
    void Gmat::read_matrix(const std::string &gmat_file, std::vector<std::int64_t> &g_row, std::vector<std::int64_t> &g_col, std::vector<T> &g_val)
    {
        // reads G matrix from gmat.dat file into the linked list container

        std::ifstream ped;

        try
        {
            std::string where("Error during reading G matrix");
            std::string line;

            char *end;
            const char *p;
            std::vector<T> tmp_list;

            ped.exceptions(std::ifstream::failbit | std::ifstream::badbit);

            size_t diagonals = 0; /* for debugging */

            ped.open(gmat_file, std::fstream::in);

            if (ped)
            {
                while (getline(ped, line))
                {
                    p = line.c_str();
                    for (double f = std::strtod(p, &end); p != end; f = std::strtod(p, &end))
                    {
                        p = end;
                        if (errno == ERANGE)
                        {
                            throw std::string("Range error during reading G matrix");
                            errno = 0;
                        }
                        tmp_list.push_back(static_cast<T>(f));
                    }
                    g_row.push_back(static_cast<std::int64_t>(tmp_list[0]));
                    g_col.push_back(static_cast<std::int64_t>(tmp_list[1]));
                    g_val.push_back(tmp_list[2]);
                    gmatID.push_back(int(tmp_list[0]));
                    gmatID.push_back(int(tmp_list[1]));

                    if (static_cast<std::int64_t>(tmp_list[0]) == static_cast<std::int64_t>(tmp_list[1]))
                        diagonals++;

                    tmp_list.erase(tmp_list.begin(), tmp_list.end());

                    if (ped.eof())
                    {
                        ped.close();

                        if (!is_unique(gmatID))
                            gmatID.erase(unique(gmatID.begin(), gmatID.end()), gmatID.end()); // here the vector should be sorted and unique

                        if (diagonals != gmatID.size())
                            throw std::string("There are missing diagonals in G-matrix file.");
                        else
                            return;
                    }
                }
                ped.close();

                if (!is_unique(gmatID))
                    gmatID.erase(unique(gmatID.begin(), gmatID.end()), gmatID.end()); // here the vector should be sorted and unique
            }
            else
                throw std::string("Cannot open G-matrix file!");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::read_matrix(const std::string &, std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<T> &)" << '\n';
            std::cerr << e.what() << '\n';
            if (ped.eof())
                std::cerr << "Check for empty line(s) at the endd of the file!"
                          << "\n";
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::read_matrix(const std::string &, std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<T> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::read_matrix(const std::string &, std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<T> &)" << '\n';
            throw;
        }
    }

    template void Gmat::read_matrix(const std::string &fname, std::vector<std::int64_t> &g_row, std::vector<std::int64_t> &g_col, std::vector<float> &g_val);
    template void Gmat::read_matrix(const std::string &fname, std::vector<std::int64_t> &g_row, std::vector<std::int64_t> &g_col, std::vector<double> &g_val);

    //===============================================================================================================

    void Gmat::read_matrix(const std::string &gmat_file)
    {
        // NOTE: reads from file G mmatrix in compact format upper (lower) triangular part!
        try
        {
            std::vector<std::int64_t> g_row;
            std::vector<std::int64_t> g_col;
            std::vector<double> g_val;
            std::map<std::int64_t, std::int64_t> rid_map;

            read_matrix(gmat_file, g_row, g_col, g_val);

            if ( gmatID.empty() )
                throw std::string("Genotyped IDs vector is empty!");
                        
			get_RecodedIdMap(rid_map, gmatID);

            if ( rid_map.empty() )
                throw std::string("Recoded IDs map is empty!");

			//auto n_threads = std::thread::hardware_concurrency();
			//block_size = static_cast<unsigned int> (g_row.size()/(n_threads));

            //G.resize(gmatID.size(), gmatID.size());
            G.resize( gmatID.size() );

#pragma omp parallel for //schedule(static, block_size) num_threads(n_threads)
			for(size_t i = 0; i < g_row.size(); i++) {
				size_t r = rid_map[g_row[i]];
				size_t c = rid_map[g_col[i]];
				size_t ind = r*(r - 1)/2 + c - 1;
				if (c > r)
                    ind = c*(c - 1)/2 + r - 1;
                G[ind] = static_cast<double> (g_val[i]);
			}
G.symtorec();
			g_row.clear();
			g_row.shrink_to_fit();
			g_col.clear();
			g_col.shrink_to_fit();
			g_val.clear();
			g_val.shrink_to_fit();
			rid_map.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::read_matrix(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::read_matrix(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::read_matrix(const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::read_snp(const std::string &snp_file)
    {
        std::ifstream snpF;

        try
        {
            std::string line;
            std::vector<std::string> data_list;

            snpF.exceptions(std::ifstream::failbit | std::ifstream::badbit);

            snpF.open(snp_file, std::fstream::in);

            if (snpF)
            {
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
                            continue;

                        data_list.push_back(token);

                        if (data_list.size() == 2)
                            break;
                    }
                    /* get the last element of the string */
                    data_list.push_back(line);

                    /* now we have got the SNP data for one ID */
                    snp_map[stoi(data_list[0])] = data_list[2];
                    data_list.erase(data_list.begin(), data_list.end());

                    if (snpF.eof())
                    {
                        snpF.close();

                        /* build the map: <index, ID> */
                        size_t tmpInd = 0;
                        for (const auto &snp : snp_map)
                        {
                            anim_id_map[tmpInd] = snp.first;
                            tmpInd++;
                        }

                        return;
                    }
                }

                snpF.close();

                /* build the map: <index, ID> */
                size_t tmpInd = 0;
                for (const auto &snp : snp_map)
                {
                    anim_id_map[tmpInd] = snp.first;
                    tmpInd++;
                }
            }
            else
                throw std::string("Cannot open SNPs file!");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::read_snp(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            if (snpF.eof())
                std::cerr << "Check for empty line(s) at the endd of the file!"
                          << "\n";
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::read_snp(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::read_snp(const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    size_t Gmat::find_invect(std::vector<std::int64_t> &where, std::int64_t what)
    {
        try
        {
            std::vector<std::int64_t>::iterator it;
            it = find(where.begin(), where.end(), what);
            if (it != where.end())
                return it - where.begin();
            else
                return -1;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::find_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::find_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::find_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::find_RecodedIdMap(std::map<std::int64_t, std::int64_t> &id_map,
                                 std::vector<std::int64_t> &whereIidList,
                                 std::vector<std::int64_t> &whatIdList)
    {
        try
        {
            for (size_t i = 0; i < whatIdList.size(); i++)
                id_map[whatIdList[i]] = find_invect(whereIidList, whatIdList[i]) + 1;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::find_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::find_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::find_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::get_RecodedIdMap(std::map<std::int64_t, std::int64_t> &id_map,
                                std::vector<std::int64_t> &idVect)
    {
        try
        {
            size_t code_id = 1;
            for (auto const &elem : idVect)
            {
                id_map[elem] = code_id;
                code_id++;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::get_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::get_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::get_RecodedIdMap(std::map<std::int64_t, std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    bool Gmat::is_unique(std::vector<std::int64_t> &x)
    {
        try
        {
            bool out = false;
            sort(x.begin(), x.end());
            if (adjacent_find(x.begin(), x.end()) == x.end())
                out = true;
            return out;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::is_unique(std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::is_unique(std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::is_unique(std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::parse_string(std::string &snp_str, std::vector<int> &markers)
    {
        try
        {
            size_t sz = snp_str.length() + 1;

            char *cstr = new char[snp_str.length() + 1];
            std::strcpy(cstr, snp_str.c_str());

            for (size_t i = 0; i < sz; i++)
                if (isdigit(cstr[i]))
                    markers.push_back((int)cstr[i] - (int)48);
            
            delete[] cstr;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::parse_string(std::string& , std::vector<int>&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::parse_string(std::string& , std::vector<int>&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::parse_string(std::string& , std::vector<int>&)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::make_zmatrix()
    {
        try
        {
            /* get the number of SNPs */
            std::vector<int> tmpVect;
            
            auto it = snp_map.begin();
            std::string tmpStr = it->second;
            
            parse_string(tmpStr, tmpVect);
            
            size_t snpNum = tmpVect.size();
            
            tmpVect.clear();
            tmpVect.shrink_to_fit();
            tmpStr.clear();
            tmpStr.shrink_to_fit();
            // ------------------------

            // debug
            // std::cout<<"snpNum = "<<snpNum<<std::endl;
            // exit(1);
            // end debug

            /* declare matrix M*/
            evolm::matrix<double> M(snp_map.size(), snpNum);

            /* vector of SNPs frequences and missed values */
            std::vector<double> P(snpNum, 0.0);
            std::vector<int> missed(snpNum, 0);
            std::vector<double> missed2pq(snp_map.size(), 0.0);

            /* map of missed values locations */
            std::vector<std::vector<int>> missedLocation;
            for (size_t i = 0; i < snpNum; i++)
                missedLocation.push_back(std::vector<int>());

            /* parse SNPs and fill matrix M*/
            size_t rowI = 0;
            for (auto const &e : snp_map)
            {
                std::vector<int> parsedMarkers;
                std::string strToParse = e.second;

                parse_string(strToParse, parsedMarkers);

                for (size_t i = 0; i < parsedMarkers.size(); i++)
                {
                    M(rowI, i) = static_cast<double>(parsedMarkers[i]);
                    if (parsedMarkers[i] != 0 && parsedMarkers[i] != 1 && parsedMarkers[i] != 2)
                    {
                        missed[i] += 1;
                        missedLocation[i].push_back(rowI);
                    }
                    else
                        P[i] += static_cast<double>(parsedMarkers[i]);
                }
                rowI++;
            }

            /* finish to calculate allele frequences, additionally accounting missing values */

            auto n_threads = std::thread::hardware_concurrency();
            auto block_size = static_cast<unsigned int>(P.size() / (n_threads));

            if (block_size < workload)
            {
                block_size = P.size();
                n_threads = 1;
            }

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < P.size(); i++)
            {
                P[i] = P[i] / (2 * (snp_map.size() - missed[i]));
            }

            Z.resize(snp_map.size(), snpNum);

            for (size_t i = 0; i < snp_map.size(); i++)
            {
                for (size_t j = 0; j < snpNum; j++)
                {
                    Z(i, j) = M(i, j) - 2 * P[j];
                }
            }

            /* modify Z matrix, so instead of missing values we put population average (0.0) */

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < missedLocation.size(); i++)
            {
                for (size_t j = 0; j < missedLocation[i].size(); j++)
                {
                    Z(missedLocation[i][j], i) = 0.0;
                    missed2pq[missedLocation[i][j]] = missed2pq[missedLocation[i][j]] + 2 * P[i] * (1.0 - P[i]);
                }
            }

            freq = 0.0;
#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t j = 0; j < P.size(); j++)
            {
                freq += P[j] * (1.0 - P[j]);
            }
            freq = 2 * freq;

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < snp_map.size(); i++)
            {
                missed2pq[i] = sqrt(freq / (freq - missed2pq[i]));
            }

            /*
                After centering, adjust for missing markers for each animal;
                adjust for sqrt[sum of 2pq over all loci /sum of 2pq over non-missing loci.
            */
#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < snp_map.size(); i++)
            {
                for (size_t j = 0; j < snpNum; j++)
                {
                    Z(i, j) = Z(i, j) * missed2pq[i];
                }
            }

            M.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::make_zmatrix()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::make_zmatrix()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::make_zmatrix()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::make_matrix()
    {
        try
        {
            G = ( Z^2 ) * (1 / freq);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::make_matrix()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::make_matrix()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::make_matrix()" << '\n';
            throw;
        }
    }

    //===============================================================================================================
    
    template <typename T>
    void Gmat::get_gvalues(std::vector <std::int64_t>& row, std::vector <std::int64_t>& col, std::vector <T>& val, double diag_val)
    {
        try
        {
            for(size_t i = 0; i < anim_id_map.size(); i++)
            {
                for(size_t j = 0; j <= i; j++)
                {
                    col.push_back(anim_id_map[j]);
                    row.push_back(anim_id_map[i]);
                    
                    if ( i == j )
                        val.push_back( static_cast<T>( G(i,j) + diag_val ) );
                    else
                        val.push_back( static_cast<T>( G(i,j) ) );
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::get_gvalues(std::vector <std::int64_t>&, std::vector <std::int64_t>&, std::vector <T>&, double)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::get_gvalues(std::vector <std::int64_t>&, std::vector <std::int64_t>&, std::vector <T>&, double)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::get_gvalues(std::vector <std::int64_t>&, std::vector <std::int64_t>&, std::vector <T>&, double)" << '\n';
            throw;
        }
    }

    template void Gmat::get_gvalues(std::vector <std::int64_t>& row, std::vector <std::int64_t>& col, std::vector <float>& val, double diag_val);
    template void Gmat::get_gvalues(std::vector <std::int64_t>& row, std::vector <std::int64_t>& col, std::vector <double>& val, double diag_val);

    //===============================================================================================================

    void Gmat::get_gids(std::vector <std::int64_t>& ids)
    {
        try
        {
            if ( gmatID.empty() )
            {
            for(const auto& id: anim_id_map)
                ids.push_back(id.second);
            }
            else
                ids = gmatID;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::get_gids(std::vector <std::int64_t>&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::get_gids(std::vector <std::int64_t>&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::get_gids(std::vector <std::int64_t>&)" << '\n';
            throw;
        }
    }

    //===============================================================================================================
}