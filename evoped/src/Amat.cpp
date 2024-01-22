#include "Amat.hpp"

namespace evoped
{
    //===============================================================================================================

    Amat::Amat()
    {
        try
        {
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::Amat()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::Amat()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::Amat()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    Amat::~Amat()
    {
        try
        {
            clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::Amat()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::Amat()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::Amat()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::map_to_matr(std::map<PedPair, double> &amap, std::vector<std::int64_t> &ids, bool use_ainv)
    {
        try
        {
            Utilities2 u;

            std::map<std::int64_t, std::int64_t> rid_map;

            if ( ids.empty() )
                throw std::string("Empty traced pedigree IDs!");

            size_t limit = 0.5*(ids.size()-1)*ids.size()+ids.size();
            if ( !use_ainv )
                limit = ids.size() * ids.size();

            if ( limit < amap.size() )
                throw std::string("The number of elements in calculated A(-1) matrix is higher than the number of traced IDs!!");

            u.get_RecodedIdMap(rid_map, ids);

            if (rid_map.empty())
                throw std::string("Recoded IDs map is empty!");

            if ( !A.empty())
                A.clear();
            
            A.resize(ids.size());

//#pragma omp parallel for // this does not work ???
            for (auto const &elem : amap)
            {
                std::int64_t g_row = elem.first.val_1;
                std::int64_t g_col = elem.first.val_2;
                double g_val = elem.second;

                size_t r = rid_map[ g_row ];
                size_t c = rid_map[ g_col ];
                size_t ind = r * (r - 1) / 2 + c - 1;
                if (c > r)
                    ind = c * (c - 1) / 2 + r - 1;
                A[ind] = g_val;
            }

            rid_map.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::map_to_matr(std::map<PedPair, double> &, std::vector<std::int64_t> &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::map_to_matr(std::map<PedPair, double> &, std::vector<std::int64_t> &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::map_to_matr(std::map<PedPair, double> &, std::vector<std::int64_t> &, bool)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::make_matrix(const std::string &ped_file, bool use_ainv)
    {
        // Making A or A(-1) based on full pedigree
        try
        {
            std::vector<std::int64_t> pedID;
            std::map<PedPair, PedPair> pedigree_from_file;
            std::map<PedPair, PedPair> pedigree;

            fread_pedigree(ped_file, pedigree_from_file, pedID);

            if (pedID.empty())
                throw std::string("Empty pedigree IDs!");

            if (pedigree_from_file.empty())
                throw std::string("File provided pedigree is empty!");

            if (!traced_pedID.empty())
                traced_pedID.erase(traced_pedID.begin(), traced_pedID.end());

            trace_pedigree(pedigree_from_file, pedigree, pedID); // tracing full pedigree

            if ( !inbrF.empty() )
                inbrF.erase(inbrF.begin(), inbrF.end());

            std::map<PedPair, double> ainv;

            if ( use_ainv )
                get_ainv(pedigree, ainv, true);
            else
                get_a(pedigree, ainv);
 
            pedigree.clear();
            pedigree_from_file.clear();
            birth_id_map.clear();
            pedID.clear();
            pedID.shrink_to_fit();

            map_to_matr(ainv, traced_pedID, use_ainv); // converting ainv map to A or A(-1) matrix

            ainv.clear();

            id_iA = traced_pedID;
            A.fwrite();
            iA = A;            
            A.clear();
            traced_pedID.clear();
            traced_pedID.shrink_to_fit();

        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::make_matrix(std::string &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::make_matrix(std::string &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::make_matrix(std::string &, bool)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::make_matrix(const std::string &ped_file, const std::string &g_file, bool use_ainv)
    {
        // Making A(-1) based on reduced pedigree traced on IDs in the file 'g_file'
        try
        {
            std::vector<std::int64_t> pedID;
            std::vector<std::int64_t> genotypedID;
            std::map<PedPair, PedPair> pedigree_from_file;
            std::map<PedPair, PedPair> r_pedigree;
            
            fread_pedigree(ped_file, pedigree_from_file, pedID);

            if (pedigree_from_file.empty())
                throw std::string("File provided pedigree is empty!");

            fread_genotyped_id(g_file, genotypedID);

            if (genotypedID.empty())
                throw std::string("Cannot trace the reduced pedigree: ID's vector is empty!");
            
            if (!traced_pedID.empty())
                traced_pedID.erase(traced_pedID.begin(), traced_pedID.end());

            trace_pedigree(pedigree_from_file, r_pedigree, genotypedID); // tracing reduced pedigree for genotyped individuals

            if ( !inbrF.empty() )
                inbrF.erase(inbrF.begin(), inbrF.end());

            std::map<PedPair, double> r_ainv;

            if ( use_ainv )
                get_ainv(r_pedigree, r_ainv, true);
            else
                get_a(r_pedigree, r_ainv);

            r_pedigree.clear();
            pedigree_from_file.clear();
            birth_id_map.clear();
            genotypedID.clear();
            genotypedID.shrink_to_fit();
            pedID.clear();
            pedID.shrink_to_fit();

            map_to_matr(r_ainv, traced_pedID, use_ainv); // converting r_ainv map to A(-1) matrix

            r_ainv.clear();

            id_irA = traced_pedID;
            A.fwrite();
            irA = A;            
            A.clear();
            traced_pedID.clear();
            traced_pedID.shrink_to_fit();

        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::make_matrix(std::string &, std::string &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::make_matrix(std::string &, std::string &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::make_matrix(std::string &, std::string &, bool)" << '\n';
            throw;
        }
    }



    //===============================================================================================================

    void Amat::make_all(const std::string &ped_file, const std::string &g_file)
    {
        // Making all matrices required for ssBlup: A(-1), red_A(-1), A22, A22(-1)
        try
        {
            // -----------------------------------------------------
            // ---------- full A(-1) -------------------------------

            std::vector<std::int64_t> pedID;
            std::map<PedPair, PedPair> pedigree_from_file;
            std::map<PedPair, PedPair> pedigree;

            fread_pedigree(ped_file, pedigree_from_file, pedID);

            if (pedID.empty())
                throw std::string("Empty pedigree IDs!");

            if (pedigree_from_file.empty())
                throw std::string("File provided pedigree is empty!");

            if (!traced_pedID.empty())
                traced_pedID.erase(traced_pedID.begin(), traced_pedID.end());

            trace_pedigree(pedigree_from_file, pedigree, pedID); // tracing full pedigree

            if ( !inbrF.empty() )
                inbrF.erase(inbrF.begin(), inbrF.end());

            std::map<PedPair, double> ainv;

            get_ainv(pedigree, ainv, true);
 
            pedigree.clear();
            pedID.clear();
            pedID.shrink_to_fit();

            map_to_matr(ainv, traced_pedID, true); // converting ainv map to A or A(-1) matrix

            ainv.clear();

            id_iA = traced_pedID;
            
            A.fwrite(); // 1. Wrire to binary
            iA = A;     // 2. Copy matrix by just exchanging the internal binary file name
            //iA.fwrite();
            A.clear();  // 3. Clears the memory and gets new name for internal binary file

            // -----------------------------------------------------
            // ---------- reduced A(-1) ----------------------------

            std::vector<std::int64_t> genotypedID;
            std::map<PedPair, PedPair> r_pedigree;
          
            fread_genotyped_id(g_file, genotypedID);

            if (genotypedID.empty())
                throw std::string("Cannot trace the reduced pedigree: ID's vector is empty!");
            
            if (!traced_pedID.empty())
                traced_pedID.erase(traced_pedID.begin(), traced_pedID.end());

            trace_pedigree(pedigree_from_file, r_pedigree, genotypedID); // tracing reduced pedigree for genotyped individuals

            if ( !inbrF.empty() )
                inbrF.erase(inbrF.begin(), inbrF.end());

            std::map<PedPair, double> r_ainv;

            get_ainv(r_pedigree, r_ainv, true);

            pedigree_from_file.clear();

            map_to_matr(r_ainv, traced_pedID, true); // converting r_ainv map to A(-1) matrix

            r_ainv.clear();

            id_irA = traced_pedID;

            A.fwrite(); // 1. Wrire to binary
            irA = A;    // 2. Copy matrix by just exchanging the internal binary file name
            //irA.fwrite();

            A.clear();  // 3. Clears the memory and gets new name for internal binary file

            // -----------------------------------------------------
            // -------------------- A22 ----------------------------

            get_A22(r_pedigree, genotypedID);

            id_A22 = genotypedID;

            r_pedigree.clear();
            birth_id_map.clear();

            A.fwrite(); // 1. Wrire to binary
            A22 = A;    // 2. Copy matrix by just exchanging the internal binary file name
            //A22.fwrite();

            A.clear();  // 3. Clears the memory and gets new name for internal binary file

            // -----------------------------------------------------
            // -------------------- A22(-1) ------------------------

            irA.fread();

            get_iA22(irA, id_irA, genotypedID);

            A.fwrite(); // 1. Wrire to binary
            irA.fwrite();
            iA22 = A;   // 2. Copy matrix by just exchanging the internal binary file name
            //iA22.fwrite();
            A.clear();  // 3. Clears the memory and gets new name for internal binary file
            genotypedID.clear();
            genotypedID.shrink_to_fit();

        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::make_all(std::string &, std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::make_all(std::string &, std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::make_all(std::string &, std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::get_iA22(evolm::matrix<double>& full_matr, std::vector<std::int64_t>& matr_ids, std::vector<std::int64_t>& selected_ids)
    {
        try
        {
            Utilities2 u;

            evolm::matrix<size_t> shapeofh;
            shapeofh = full_matr.shape();

            size_t expected_size = 0.5 * (matr_ids.size() - 1) * matr_ids.size() + matr_ids.size();

            if ( shapeofh[0] != shapeofh[1] )
                throw std::string("The passed matrix has wrong dimension: number of raws is not the same as number of columns!");

            if ( shapeofh[0] != matr_ids.size() )
                throw std::string("The passed matrix has wrong dimension!");

            if ( full_matr.size() > expected_size )
                throw std::string("The number of elements in the passed matrix is greater then expected!");

            if ( selected_ids.size() > matr_ids.size() )
                throw std::string("The number of elements in the passed selected IDs array is greater then the number of IDs in the passed matrix!");

            if ( !u.is_value_in_vect(matr_ids, selected_ids) )
                throw std::string("There are IDs in the selected IDs array which are not part of the passed matrix!");

            // --------------------------------
            if ( !A.empty() )
            {
                A.fclear();
                A.clear();                
            }

            // Create the list of IDs which are in matr_ids but not in selected_ids vectors

            std::vector<std::int64_t> not_selected_ids;

            for ( size_t i = 0; i < matr_ids.size(); i++ )
            {
                int res = u.find_invect( selected_ids, matr_ids[i] );
                if ( res == -1 )
                    not_selected_ids.push_back( matr_ids[i] );
            }

            if ( ( selected_ids.size() + not_selected_ids.size() ) != matr_ids.size() )
                throw std::string("The sum of IDs from from two vectors is not equal to number of IDs in the matrix!");
            
            // ------------------------------------------------

            // Next steps: A22(-1) = A22 - A21 * A11(-1) * A12;

            evolm::matrix<double> a22; // in order to differrentiate from the 'global' A22
            evolm::matrix<double> A11;
            evolm::matrix<double> A21;
            evolm::matrix<double> A12;

            // -------------------- A11 -----------------------

            A11.resize( not_selected_ids.size() );

            std::vector<size_t> non_selected_pos;
            for (size_t i = 0; i < not_selected_ids.size(); i++)
                non_selected_pos.push_back( u.find_invect( matr_ids, not_selected_ids[i] ) );

            std::vector<size_t> selected_pos;
            for (size_t i = 0; i < selected_ids.size(); i++)
                selected_pos.push_back( u.find_invect( matr_ids, selected_ids[i] ) );

#pragma omp parallel for
            for (size_t i = 0; i < not_selected_ids.size(); i++)
            {
                //size_t pos_i = u.find_invect( matr_ids, not_selected_ids[i] );
                size_t pos_i = non_selected_pos[i];

                for (size_t j = 0; j <= i; j++)
                {
                    //size_t pos_j = u.find_invect( matr_ids, not_selected_ids[j] );
                    size_t pos_j = non_selected_pos[j];

                    if ( pos_i >= pos_j )
                        A11(i,j) = full_matr( pos_i, pos_j );
                    else
                        A11(i,j) = full_matr( pos_j, pos_i );
                }
            }

            A11.symtorec();

            A11.invert();

            A11.fwrite();
            // ------------------------------------------------
            //
            // -------------------- A21 -----------------------

            A21.resize(selected_ids.size(), not_selected_ids.size());

#pragma omp parallel for
            for (size_t i = 0; i < selected_ids.size(); i++)
            {
                //size_t pos_i = u.find_invect( matr_ids, selected_ids[i] );
                size_t pos_i = selected_pos[i];

                for (size_t j = 0; j < not_selected_ids.size(); j++)
                {
                    //size_t pos_j = u.find_invect( matr_ids, not_selected_ids[j] );
                    size_t pos_j = non_selected_pos[j];

                    if ( pos_i >= pos_j )
                        A21(i,j) = full_matr( pos_i, pos_j );
                    else
                        A21(i,j) = full_matr( pos_j, pos_i );
                }
            }
            // ------------------------------------------------
            //
            // -------------------- A12 -----------------------

            A12 = A21;

            A12.transpose();

            A12.fwrite();
            // ------------------------------------------------
            //
            // --------------- A21 * A11(-1) * A12 ------------

            A11.fread();
            evolm::matrix<double> res;

            res = A21 * A11;
            
            A11.fclear();
            A11.clear();

            A21.fclear();
            A21.clear();

            A12.fread();

            res = res * A12;

            A12.fclear();
            A12.clear();
            
            // ------------------------------------------------
            //
            // -------------------- A22 -----------------------

            a22.resize( selected_ids.size() );

#pragma omp parallel for
            for (size_t i = 0; i < selected_ids.size(); i++)
            {
                //size_t pos_i = u.find_invect( matr_ids, selected_ids[i] );
                size_t pos_i = selected_pos[i];

                for (size_t j = 0; j <= i; j++)
                {
                    //size_t pos_j = u.find_invect( matr_ids, selected_ids[j] );
                    size_t pos_j = selected_pos[j];

                    if ( pos_i >= pos_j )
                        a22(i,j) = full_matr( pos_i, pos_j );
                    else
                        a22(i,j) = full_matr( pos_j, pos_i );
                }
            }

            non_selected_pos.clear();
            non_selected_pos.shrink_to_fit();
            selected_pos.clear();
            selected_pos.shrink_to_fit();

            //a22.symtorec();

            // ------------------------------------------------
            //
            // ------------------ A22 - res -------------------

            evolm::matrix<size_t> shapeofa22;
            shapeofa22 = a22.shape();

            A.resize(shapeofa22[0]);

            res.rectosym();

#pragma omp parallel for
            for (size_t i = 0; i < a22.size(); i++)
                A[i] = a22[i] - res[i];
            
            //A = a22 - res;
            //A.rectosym();

            a22.fclear();
            a22.clear();            
            res.fclear();
            res.clear();            
            // ------------------------------------------------
            shapeofh.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_iA22(evolm::matrix<double>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_iA22(evolm::matrix<double>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_iA22(evolm::matrix<double>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::fread_pedigree(const std::string &ped_file, std::map<PedPair, PedPair> &out_ped, std::vector<std::int64_t> &out_ids)
    {
        try
        {
            Utilities2 u;

            std::string line;
            PedPair key;
            PedPair id_pair;

            char *end;
            const char *p;
            std::vector<double> tmp_list;

            std::ifstream ped;
            ped.open(ped_file, std::fstream::in);

            if ( !ped.good() )
                throw std::string("Cannot open pedigree file!");
            
            while (getline(ped, line))
            {
                p = line.c_str();
                for (double f = std::strtod(p, &end); p != end; f = std::strtod(p, &end))
                {
                    p = end;
                    if (errno == ERANGE)
                    {
                        errno = 0;
                        throw std::string("Range error during reading pedigree file");
                    }
                    tmp_list.push_back(f);
                }

                if (tmp_list.size() > 0)
                {
                    key.val_1 = static_cast<std::int64_t>(tmp_list[3]);
                    key.val_2 = static_cast<std::int64_t>(tmp_list[0]);
                    id_pair.val_1 = static_cast<std::int64_t>(tmp_list[1]);
                    id_pair.val_2 = static_cast<std::int64_t>(tmp_list[2]);
                    out_ped[key] = id_pair;
                    birth_id_map[static_cast<std::int64_t>(tmp_list[0])] = static_cast<std::int64_t>(tmp_list[3]);
                    out_ids.push_back(static_cast<std::int64_t>(tmp_list[0]));

                    tmp_list.erase(tmp_list.begin(), tmp_list.end());
                }
            }

            ped.close();

            if (!u.is_unique(out_ids))
                out_ids.erase(unique(out_ids.begin(), out_ids.end()), out_ids.end()); // here the vector should be sorted and unique
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::fread_pedigree(const std::string &, std::map<PedPair, PedPair> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::fread_pedigree(const std::string &, std::map<PedPair, PedPair> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::fread_pedigree(const std::string &, std::map<PedPair, PedPair> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::fread_genotyped_id(const std::string &g_file, std::vector<std::int64_t> &out_ids)
    {
        // reads genotyped and core IDs from typed file into the vectors: genotyped, core
        try
        {
            Utilities2 u;

            std::string line;
            std::int64_t t_gtyp, t_gcor;
            t_gtyp = t_gcor = 0;

            char *end;
            const char *p;
            std::vector<double> tmp_list;

            std::ifstream ped;
            ped.open(g_file, std::fstream::in);

            if ( !ped.good() )
                throw std::string("Cannot open genotyped ids file!");
            
            while (std::getline(ped, line))
            {
                p = line.c_str();

                for (double f = std::strtod(p, &end); p != end; f = std::strtod(p, &end))
                {
                    p = end;
                    if (errno == ERANGE)
                    {
                        errno = 0;
                        throw std::string("Range error during reading genotyped IDs file");
                    }
                    tmp_list.push_back(f);
                }

                if (tmp_list.size() > 0)
                {
                    if (static_cast<std::int64_t>(tmp_list[0]) != 0 && static_cast<std::int64_t>(tmp_list[0]) != t_gtyp)
                        out_ids.push_back(static_cast<std::int64_t>(tmp_list[0]));

                    //if (static_cast<std::int64_t>(tmp_list[1]) != 0 && static_cast<std::int64_t>(tmp_list[0]) != t_gtyp)
                    //    coreID.push_back(static_cast<std::int64_t>(tmp_list[0]));

                    t_gtyp = static_cast<std::int64_t>(tmp_list[0]);
                    //t_gcor = static_cast<std::int64_t>(tmp_list[1]);

                    tmp_list.erase(tmp_list.begin(), tmp_list.end());
                }
            }

            ped.close();

            if (!u.is_unique(out_ids))
                out_ids.erase(unique(out_ids.begin(), out_ids.end()), out_ids.end());

            //if (!u.is_unique(coreID))
            //    coreID.erase(unique(coreID.begin(), coreID.end()), coreID.end());
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::fread_genotyped_id(const std::string &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::fread_genotyped_id(const std::string &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::fread_genotyped_id(const std::string &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::trace_pedigree(std::map<PedPair, PedPair> &in_ped, std::map<PedPair, PedPair> &out_ped, std::vector<std::int64_t> &traced_id)
    {
        try
        {
            Utilities2 u;

            PedPair id_pair;
            PedPair key;
            std::vector<std::int64_t> gen_pedID(traced_id);
            std::int64_t elem_v;
            size_t iter_v = 0;
            bool exists;

            std::vector<std::int64_t> days;
            std::vector<std::int64_t> ids;
            std::vector<std::int64_t> sire;
            std::vector<std::int64_t> dame;

            for (auto const &elem_m : in_ped)
            {
                days.push_back(elem_m.first.val_1);
                ids.push_back(elem_m.first.val_2);
                sire.push_back(elem_m.second.val_1);
                dame.push_back(elem_m.second.val_2);
            }

            do
            {
                elem_v = gen_pedID[iter_v];

                exists = false;

                for (size_t i = 0; i < days.size(); i++)
                {

                    if (ids[i] == elem_v)
                    {

                        exists = true;
                        id_pair.val_1 = sire[i];
                        id_pair.val_2 = dame[i];

                        if ((birth_id_map[sire[i]] >= days[i]) || (birth_id_map[dame[i]] >= days[i]))
                            throw std::string("Pedigree is not correct: parents born before offspring!");

                        key.val_1 = days[i];
                        key.val_2 = elem_v;
                        out_ped[key] = id_pair;

                        if (id_pair.val_1 != 0)
                        {
                            if (!u.is_invect(gen_pedID, id_pair.val_1))
                                gen_pedID.push_back(id_pair.val_1);
                        }
                        if (id_pair.val_2 != 0)
                        {
                            if (!u.is_invect(gen_pedID, id_pair.val_2))
                                gen_pedID.push_back(id_pair.val_2);
                        }

                        break;
                    }
                }

                if (!exists)
                {
                    key.val_1 = -1;
                    key.val_2 = elem_v;
                    id_pair.val_1 = 0;
                    id_pair.val_2 = 0;
                    out_ped[key] = id_pair;
                    exists = false;
                }

                iter_v++;

            } while (iter_v < gen_pedID.size());
            
            if ( u.is_unique(gen_pedID) ) // Check if no repeated IDs appiar in the list of traced IDs
                traced_pedID = gen_pedID;
            else
                throw std::string("In the traced pedigree some IDs appiar more then one time!");

            gen_pedID.clear();
            days.clear();
            ids.clear();
            sire.clear();
            dame.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::trace_pedigree(std::map<PedPair, PedPair> &, std::map<PedPair, PedPair> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::trace_pedigree(std::map<PedPair, PedPair> &, std::map<PedPair, PedPair> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::trace_pedigree(std::map<PedPair, PedPair> &, std::map<PedPair, PedPair> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::int64_t Amat::pos_inped(std::map<std::int64_t, std::int64_t> &codemap, std::int64_t id)
    {
        try
        {
            auto pos = codemap.find(id);

            if (pos == codemap.end())
                return 0;
            else
                return pos->second;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::pos_inped(std::map<std::int64_t, std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::pos_inped(std::map<std::int64_t, std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::pos_inped(std::map<std::int64_t, std::int64_t> &, std::int64_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::get_dinv(std::map<PedPair, PedPair> &ped, std::vector<double> &dinv, bool inbreed)
    {
        // diagonal elements of A(-1)
        // Algorithm: M.Sargolzaei et al. 'A fast algorithm for computing inbreeding ...' (works correct!)
        try
        {
            if (inbreed)
            {
                size_t n = ped.size();
                size_t m = n + 1;

                std::vector<std::vector<std::int64_t>> Ped((m), std::vector<std::int64_t>(2, 0)); // main pedigree, defaults to zero initial value

                std::map<std::int64_t, std::int64_t> code_map;
                size_t code_id = 1;
                inbrF.push_back(0.0); // we shall be able to get the coeff. starting from index '1'

                for (auto const &elem : ped)
                {

                    code_map[elem.first.val_2] = code_id;

                    Ped[code_id][0] = pos_inped(code_map, elem.second.val_1);
                    Ped[code_id][1] = pos_inped(code_map, elem.second.val_2);

                    inbrF.push_back(0.0); // initialize vector of inbreeding coeff.

                    code_id++;
                }

                size_t i, j, k, rN, rS, S, D, MIP;

                std::vector<std::vector<std::int64_t>> rPed((m + 1), std::vector<std::int64_t>(2, 0)); // reduced pedigree
                std::vector<std::int64_t> SId(m, 0);                                                   // will contain the sorted animals ID based on the ID of their sires
                std::vector<std::int64_t> Link(m, 0);                                                  // will contain new ID of ancestors at position of their original ID
                std::vector<std::int64_t> MaxIdP(m + 1, 0);                                            // will contain maximum new ID of parents for each paternal group at position of the new ID of each sire

                std::vector<double> F((n + 1), 0.0); // inbreeding coefficients
                std::vector<double> B((m + 1), 0.0); // within family segregation variances
                std::vector<double> x((m + 1), 0.0); // x arrays

                // set values for the unknown parent
                F[0] = -1.0;
                x[0] = 0.0;
                Link[0] = 0;
                rN = 1;
                i = 1;

                for (; i <= n; ++i)
                { // extract and recode ancestors

                    SId[i] = i;
                    Link[i] = 0;
                    if (i <= m)
                        x[i] = 0.0;

                    S = Ped[i][0];
                    D = Ped[i][1];

                    if (S && !Link[S])
                    {
                        MaxIdP[rN] = Link[S] = rN;
                        rPed[rN][0] = Link[Ped[S][0]];
                        rPed[rN++][1] = Link[Ped[S][1]];
                    }
                    if (D && !Link[D])
                    {
                        /*MaxIdP[rN] =*/Link[D] = rN;
                        rPed[rN][0] = Link[Ped[D][0]];
                        rPed[rN++][1] = Link[Ped[D][1]];
                    }
                    if (MaxIdP[Link[S]] < Link[D])
                        MaxIdP[Link[S]] = Link[D]; // determine maximum ID of parents for each paternal group
                }

                // sort animals according to ID of their sires into SId
                // The pedigree is already sorted so the oldest animal IDs appears first
                for (size_t i = 1; i <= n; i++)
                {
                    SId[i] = i;
                }

                k = 1;
                i = 1;
                for (; i <= n;)
                { // do from the oldest sire to the youngest sire

                    int t = Ped[SId[i]][0];
                    if (!t)
                    {
                        F[SId[i++]] = 0.0; // sire is not known
                    }
                    else
                    {

                        S = Ped[SId[i]][0];
                        rS = Link[S];
                        MIP = MaxIdP[rS];
                        x[rS] = 1.0;

                        for (; k <= S; ++k)
                        { // compute within family segregation variances
                            if (Link[k])
                                B[Link[k]] = 0.5 - 0.25 * (F[Ped[k][0]] + F[Ped[k][1]]);
                        }

                        for (j = rS; j; --j)
                        { // trace back the reduced pedigree
                            if (x[j])
                            { // consider only ancestors of the sire

                                if (rPed[j][0])
                                    x[rPed[j][0]] += x[j] * 0.5;
                                if (rPed[j][1])
                                    x[rPed[j][1]] += x[j] * 0.5;
                                x[j] *= B[j];
                            }
                        }

                        for (j = 1; j <= MIP; ++j) // trace forth the reduced pedigree
                        {
                            x[j] += (x[rPed[j][0]] + x[rPed[j][1]]) * 0.5;
                        }

                        for (; i <= n; ++i) // obtain F for progeny of the current sire
                        {
                            if ((int)S != Ped[SId[i]][0])
                                break;
                            else
                                F[SId[i]] = x[Link[Ped[SId[i]][1]]] * 0.5;
                        }

                        for (j = 1; j <= MIP; ++j) // set to 0 for next evaluation of sire
                        {
                            x[j] = 0.0;
                        }
                    }
                }

                for (size_t i = 1; i <= n; i++)
                {
                    auto s = Ped[i][0];
                    auto d = Ped[i][1];

                    inbrF[i] = F[i]; // fill global vector of inbr. coeff.

                    if (s && d)
                    {
                        double g = 1.0 / (0.5 - 0.25 * (F[s] + F[d]));
                        dinv.push_back(g);
                    }
                    else if (!s && !d)
                    {
                        double g = 1.0;
                        dinv.push_back(g);
                    }
                    else
                    {
                        double g;
                        if (s)
                            g = 1.0 / (0.75 - 0.25 * (F[s]));
                        else
                            g = 1.0 / (0.75 - 0.25 * (F[d]));
                        dinv.push_back(g);
                    }
                }
            }
            else
            {
                for (size_t i = 0; i <= ped.size(); i++)
                    inbrF.push_back(0.0);

                std::int64_t s, d;
                for (auto const &elem : ped)
                {
                    s = elem.second.val_1;
                    d = elem.second.val_2;
                    if (s && d)
                        dinv.push_back(2.0);
                    else if (!s && !d)
                        dinv.push_back(1.0);
                    else
                        dinv.push_back(4.0 / 3.0);
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_dinv(std::map<PedPair, PedPair> &, std::vector<double> &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_dinv(std::map<PedPair, PedPair> &, std::vector<double> &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_dinv(std::map<PedPair, PedPair> &, std::vector<double> &, bool)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::get_ainv(std::map<PedPair, PedPair> &ped, std::map<PedPair, double> &ai, bool inbreed)
    {
        // Calculates A(-1) matrix (as map representation)
        try
        {
            PedPair akey;
            std::vector<double> di; // store D(-1)

            get_dinv(ped, di, inbreed);

            std::int64_t s, d, id;
            double dinv;

            for (auto const &elem : ped)
            {
                dinv = di.front();
                s = elem.second.val_1;
                d = elem.second.val_2;
                id = elem.first.val_2;

                if (s && d)
                {
                    akey.val_1 = id;
                    akey.val_2 = id;
                    ai[akey] = ai[akey] + dinv;

                    akey.val_1 = id;
                    akey.val_2 = s;
                    if (s > id)
                    {
                        akey.val_1 = s;
                        akey.val_2 = id;
                    }
                    ai[akey] = ai[akey] - dinv / 2.0;

                    akey.val_1 = id;
                    akey.val_2 = d;
                    if (d > id)
                    {
                        akey.val_1 = d;
                        akey.val_2 = id;
                    }
                    ai[akey] = ai[akey] - dinv / 2.0;

                    akey.val_1 = s;
                    akey.val_2 = s;
                    ai[akey] = ai[akey] + dinv / 4.0;

                    akey.val_1 = d;
                    akey.val_2 = d;
                    ai[akey] = ai[akey] + dinv / 4.0;

                    if (s >= d)
                    {
                        akey.val_1 = s;
                        akey.val_2 = d;
                    }
                    else
                    {
                        akey.val_1 = d;
                        akey.val_2 = s;
                    }
                    ai[akey] = ai[akey] + dinv / 4.0;
                }
                else if (!s && !d)
                {
                    akey.val_1 = id;
                    akey.val_2 = id;
                    ai[akey] = ai[akey] + dinv;
                }
                else
                {
                    akey.val_1 = id;
                    akey.val_2 = id;
                    ai[akey] = ai[akey] + dinv;

                    akey.val_1 = id;
                    if (s)
                    {
                        akey.val_2 = s;
                        if (s > id)
                        {
                            akey.val_1 = s;
                            akey.val_2 = id;
                        }
                    }
                    else
                    {
                        akey.val_2 = d;
                        if (d > id)
                        {
                            akey.val_1 = d;
                            akey.val_2 = id;
                        }
                    }
                    ai[akey] = ai[akey] - dinv / 2.0;

                    if (s)
                    {
                        akey.val_1 = s;
                        akey.val_2 = s;
                    }
                    else
                    {
                        akey.val_1 = d;
                        akey.val_2 = d;
                    }
                    ai[akey] = ai[akey] + dinv / 4.0;
                }
                di.erase(di.begin());
            }
            di.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_ainv(std::map<PedPair, PedPair> &, std::map<PedPair, double> &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_ainv(std::map<PedPair, PedPair> &, std::map<PedPair, double> &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_ainv(std::map<PedPair, PedPair> &, std::map<PedPair, double> &, bool)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::get_a(std::map<PedPair, PedPair> &in_ped, std::map<PedPair, double> &out_a)
    {
        // Calculates A matrix (as map representation)
        try
        {
            std::vector<std::int64_t> id;

            for (auto const &e : in_ped)
                id.push_back(e.first.val_2);

            PedPair akey;
            PedPair akey2;
            PedPair bkey;
            PedPair ckey;

            for (size_t i = 0; i < id.size(); i++)
            {
                std::int64_t id1 = id[i];
                
                for (size_t j = i; j < id.size(); j++)
                {
                    std::int64_t id2 = id[j];

                    akey2.val_1 = id2;
                    akey2.val_2 = id1;

                    akey.val_1 = id1;
                    akey.val_2 = id2;

                    auto begin = in_ped.begin();
                    std::advance(begin, j);

                    std::int64_t s2 = begin->second.val_1;
                    std::int64_t d2 = begin->second.val_2;

                    if ( s2 && d2 )
                    {
                        if ( id1 == id2 )
                        {
                            bkey.val_1 = s2;
                            bkey.val_2 = d2;
                            if (s2 > d2)
                            {
                                bkey.val_1 = d2;
                                bkey.val_2 = s2;
                            }
                            out_a[akey] = 1.0 + 0.5 * out_a[bkey];
                        }
                        else
                        {
                            bkey.val_1 = id1;
                            bkey.val_2 = s2;
                            if (id1 > s2)
                            {
                                bkey.val_1 = s2;
                                bkey.val_2 = id1;
                            }

                            ckey.val_1 = id1;
                            ckey.val_2 = d2;
                            if (id1 > d2)
                            {
                                ckey.val_1 = d2;
                                ckey.val_2 = id1;
                            }
                            out_a[akey] = 0.5 * ( out_a[bkey] + out_a[ckey] );
                            out_a[akey2] = out_a[akey];
                        }
                    }
                    if ( !s2 && !d2 )
                    {
                        if ( id1 == id2 )
                        {
                            out_a[akey] = 1.0;
                        }
                        else
                        {
                            out_a[akey] = 0.0;
                            out_a[akey2] = out_a[akey];
                        }
                    }
                    if ( s2 && !d2 )
                    {
                        if ( id1 == id2 )
                        {
                            out_a[akey] = 1.0;
                        }
                        else
                        {
                            bkey.val_1 = id1;
                            bkey.val_2 = s2;
                            if (id1 > s2)
                            {
                                bkey.val_1 = s2;
                                bkey.val_2 = id1;
                            }
                            out_a[akey] = 0.5 * out_a[bkey];
                            out_a[akey2] = out_a[akey];
                        }
                    }
                    if ( !s2 && d2 )
                    {
                        if ( id1 == id2 )
                        {
                            out_a[akey] = 1.0;
                        }
                        else
                        {
                            ckey.val_1 = id1;
                            ckey.val_2 = d2;
                            if (id1 > d2)
                            {
                                ckey.val_1 = d2;
                                ckey.val_2 = id1;
                            }
                            out_a[akey] = 0.5 * out_a[ckey];
                            out_a[akey2] = out_a[akey];
                        }
                    }

                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_a(std::map<PedPair, PedPair> &, std::map<PedPair, double> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_a(std::map<PedPair, PedPair> &, std::map<PedPair, double> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_a(std::map<PedPair, PedPair> &, std::map<PedPair, double> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::get_inbreeding(std::vector<double> &out)
    {
        try
        {
            out = inbrF;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_inbreeding(std::vector<double> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_inbreeding(std::vector<double> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_inbreeding(std::vector<double> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================
    
        void Amat::clear()
        {
            try
            {
                birth_id_map.clear();
                A.fclear();
                A.clear();                
                traced_pedID.clear();
                traced_pedID.shrink_to_fit();
                inbrF.clear();
                inbrF.shrink_to_fit();
                iA.fclear();
                iA.clear();
                irA.fclear();                
                irA.clear();
                iA22.fclear();
                iA22.clear();
                A22.fclear();
                A22.clear();                
                id_iA.clear();
                id_iA.shrink_to_fit();
                id_irA.clear();
                id_irA.shrink_to_fit();
                id_A22.clear();
                id_A22.shrink_to_fit();
            }
            catch (const std::exception &e)
            {
                std::cerr << "Exception in Amat::clear()" << '\n';
                std::cerr << e.what() << '\n';
                throw;
            }
            catch (const std::string &e)
            {
                std::cerr << "Exception in Amat::clear()" << '\n';
                std::cerr << "Reason: " << e << '\n';
                throw;
            }
            catch (...)
            {
                std::cerr << "Exception in Amat::clear()" << '\n';
                throw;
            }
        }

    //===============================================================================================================

    void Amat::get_matrix(const std::string &name, evolm::matrix<double> &arr, std::vector<std::int64_t> &out, bool keep_ondisk)
    {
        // Note, we operate with L-stored format, hence return lower triangular part
        try
        {
            if ( name == "A" )
            {
                // we use iA as a container for A as well
                iA.fread();
                arr = iA;
                out = id_iA;
                if (keep_ondisk)
                    iA.fwrite();
                else
                {
                    iA.fclear();
                    iA.clear();
                }
            }
            if ( name == "rA" )
            {
                // we use irA as a container for rA as well
                irA.fread();
                arr = irA;
                out = id_irA;
                if (keep_ondisk)
                    irA.fwrite();
                else
                {
                    irA.fclear();
                    irA.clear();
                }
            }
            if ( name == "iA" )
            {
                iA.fread();
                arr = iA;
                out = id_iA;
                if (keep_ondisk)
                    iA.fwrite();
                else
                {
                    iA.fclear();
                    iA.clear();
                }
            }
            if ( name == "irA" )
            {
                irA.fread();
                arr = irA;
                out = id_irA;
                if (keep_ondisk)
                    irA.fwrite();
                else
                {
                    irA.fclear();
                    irA.clear();
                }
            }
            if ( name == "iA22" )
            {
                iA22.fread();
                arr = iA22;
                out = id_A22;
                if (keep_ondisk)
                    iA22.fwrite();
                else
                {
                    iA22.fclear();
                    iA22.clear();
                }
            }
            if ( name == "A22" )
            {
                A22.fread();
                arr = A22;
                out = id_A22;
                if (keep_ondisk)
                    A22.fwrite();
                else
                {
                    A22.fclear();
                    A22.clear();
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_matrix(const std::string &, evolm::matrix<double> &, std::vector<std::int64_t> &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_matrix(const std::string &, evolm::matrix<double> &, std::vector<std::int64_t> &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_matrix(const std::string &, evolm::matrix<double> &, std::vector<std::int64_t> &, bool)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::get_matrix(const std::string &name, std::vector<double> &arr, std::vector<std::int64_t> &out, bool keep_ondisk)
    {
        // This method is for Python interfacing;
        // Note, we operate with L-stored format, hence return lower triangular part
        try
        {
            if ( name == "A" )
            {
                // we use iA as a container for A as well
                iA.fread();
                iA.to_vector(arr);
                out = id_iA;
                if (keep_ondisk)
                    iA.fwrite();
                else
                {
                    iA.fclear();
                    iA.clear();
                }
            }
            if ( name == "rA" )
            {
                // we use irA as a container for rA as well
                irA.fread();
                irA.to_vector(arr);
                out = id_irA;
                if (keep_ondisk)
                    irA.fwrite();
                else
                {
                    irA.fclear();
                    irA.clear();
                }
            }
            if ( name == "iA" )
            {
                iA.fread();
                iA.to_vector(arr);
                out = id_iA;
                if (keep_ondisk)
                    iA.fwrite();
                else
                {
                    iA.fclear();
                    iA.clear();
                }
            }
            if ( name == "irA" )
            {
                irA.fread();
                irA.to_vector(arr);
                out = id_irA;
                if (keep_ondisk)
                    irA.fwrite();
                else
                {
                    irA.fclear();
                    irA.clear();
                }
            }
            if ( name == "iA22" )
            {
                iA22.fread();
                iA22.to_vector(arr);
                out = id_A22;
                if (keep_ondisk)
                    iA22.fwrite();
                else
                {
                    iA22.fclear();
                    iA22.clear();
                }
            }
            if ( name == "A22" )
            {
                A22.fread();
                A22.to_vector(arr);
                out = id_A22;
                if (keep_ondisk)
                    A22.fwrite();
                else
                {
                    A22.fclear();
                    A22.clear();
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_matrix(const std::string &, std::vector<double> &, std::vector<std::int64_t> &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_matrix(const std::string &, std::vector<double> &, std::vector<std::int64_t> &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_matrix(const std::string &, std::vector<double> &, std::vector<std::int64_t> &, bool)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::get_A22(std::map <PedPair, PedPair> &ped, std::vector<std::int64_t> &genotypedID)
    {
        try
        {
            Utilities2 u;

            size_t n = ped.size(); // total amount of animals in pedigree
            std::vector<std::vector<std::int64_t> > Ped(n+1, std::vector<std::int64_t>(2, 0.0));
            std::vector<std::int64_t> GenID; // list of genotyped IDs

            std::map <std::int64_t, std::int64_t> code_map;
            std::map <std::int64_t, std::int64_t> gen_map;

            std::int64_t code_id = 1;
            for (auto const& elem : ped) {

                code_map[elem.first.val_2] = code_id;

                Ped[code_id][0] = pos_inped (code_map, elem.second.val_1);
                Ped[code_id][1] = pos_inped (code_map, elem.second.val_2);

                int pos = u.find_invect(genotypedID, elem.first.val_2);
                if (pos != -1) {
                    GenID.push_back(code_id);
                    gen_map[code_id] = pos;
                }

                code_id++;
            }

            size_t m = GenID.size(); // number of genotyped IDs

            if ( !A.empty() )
                A.clear();

            A.resize(m);

    #pragma omp parallel for
            for (size_t i = 0; i < m; i++)
            {
                std::vector<double> w(n+1, 0.0);
                std::vector<std::int64_t> v(n+1, 0);

                size_t ii = gen_map[GenID[i]]; // gives position of recoded (1...n) genotyped IDs

                v[GenID[i]] = 1;

                getA22vector(w, v, Ped);

                for (size_t j = 0; j < m; j++)
                {
                    if (ii >= (size_t)gen_map[GenID[j]])
                        A[ ii*(ii+1)/2+gen_map[ GenID[j] ] ] = w[ GenID[j] ];
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_A22(std::map <PedPai, PedPair> &, std::vector<std::int64_t> &, evolm::matrix<double> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_A22(std::map <PedPai, PedPai> &, std::vector<std::int64_t> &, evolm::matrix<double> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_A22(std::map <PedPai, PedPai> &, std::vector<std::int64_t> &, evolm::matrix<double> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::getA22vector(std::vector <double> &w, std::vector <std::int64_t> &v, std::vector<std::vector<std::int64_t> > &Ped)
    {
        try
        {
            size_t n = w.size()-1;
            std::vector<double> q(n+1, 0.0);

            for (size_t i = n; i >= 1; i--) {
                q[i] += v[i];
                auto s = Ped[i][0];
                auto d = Ped[i][1];
                if (s) q[s] += q[i]*0.5;
                if (d) q[d] += q[i]*0.5;
            }

            for (size_t i = 1; i <= n; i++)
            {
                auto s = Ped[i][0];
                auto d = Ped[i][1];
                auto di = (std::count (Ped[i].begin(), Ped[i].end(), 0)+2.0)/4.0 - 0.25 * (inbrF[s] + inbrF[d]);
                double temp = 0.0;
                if (s) temp += w[s];
                if (d) temp += w[d];
                w[i] = 0.5*temp;
                w[i] += di*q[i];
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::getA22vector(std::vector <double> &, std::vector <std::int64_t> &, std::vector<std::vector<std::int64_t> > &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::getA22vector(std::vector <double> &, std::vector <std::int64_t> &, std::vector<std::vector<std::int64_t> > &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::getA22vector(std::vector <double> &, std::vector <std::int64_t> &, std::vector<std::vector<std::int64_t> > &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

} // end of namespace evoped