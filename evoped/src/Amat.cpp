#include "Amat.hpp"

namespace evoped
{
    //===============================================================================================================

    template <typename T>
    Amat<T>::Amat()
    {
        try
        {
            data_sparsity = 0.0;
            use_sparse = false;
            sparsity_threshold = 90.0;

            IsEmpty.iA = true;
            IsEmpty.irA = true;
            IsEmpty.iA22 = true;
            IsEmpty.A22 = true;
            IsEmpty.iA_s = true;
            IsEmpty.irA_s = true;
            IsEmpty.iA22_s = true;
            IsEmpty.A22_s = true;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::Amat()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::Amat()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::Amat()" << '\n';
            throw;
        }
    }

    template Amat<float>::Amat();
    template Amat<double>::Amat();

    //===============================================================================================================

    template <typename T>
    Amat<T>::Amat(double threshold)
    {
        try
        {
            data_sparsity = 0.0;
            use_sparse = false;
            sparsity_threshold = threshold;

            IsEmpty.iA = true;
            IsEmpty.irA = true;
            IsEmpty.iA22 = true;
            IsEmpty.A22 = true;
            IsEmpty.iA_s = true;
            IsEmpty.irA_s = true;
            IsEmpty.iA22_s = true;
            IsEmpty.A22_s = true;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::Amat()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::Amat(double)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::Amat(double)" << '\n';
            throw;
        }
    }

    template Amat<float>::Amat(double threshold);
    template Amat<double>::Amat(double threshold);

    //===============================================================================================================

    template <typename T>
    Amat<T>::~Amat()
    {
        try
        {
            clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::Amat()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::Amat()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::Amat()" << '\n';
            throw;
        }
    }

    template Amat<float>::~Amat();
    template Amat<double>::~Amat();

    //===============================================================================================================

    template <typename T>
    void Amat<T>::set_sparsiity_threshold( double threshold )
    {
        try
        {
            /*
                Sparsity limit below which use_sparse == false;
                data will be processed through the dennse matrices operations.
                Higher the threshold - more zerros in a matrix should be in order
                to be considered sparse.
            */

            sparsity_threshold = threshold;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::Amat()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::set_sparsiity_threshold( double )" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::set_sparsiity_threshold( double )" << '\n';
            throw;
        }
    }

    template void Amat<float>::set_sparsiity_threshold( double threshold );
    template void Amat<double>::set_sparsiity_threshold( double threshold );

    //===============================================================================================================

    template <typename T>
    void Amat<T>::map_to_matr(std::map<PedPair, T> &in_amap, std::vector<std::int64_t> &in_ids, evolm::smatrix<T> &out_matr)
    {
        try
        {
            Utilities2 u;

            std::map<std::int64_t, std::int64_t> rid_map;

            if (in_ids.empty())
                throw std::string("Empty traced pedigree IDs!");

            u.get_RecodedIdMap(rid_map, in_ids); // list of real ids => std::map for the new consecutive list starting from 1

            if (rid_map.empty())
                throw std::string("Recoded IDs map is empty!");

            if (!out_matr.empty())
                out_matr.clear();

            out_matr.resize(in_ids.size());

            for (auto const &elem : in_amap)
            {
                std::int64_t g_row = elem.first.val_1; // id
                std::int64_t g_col = elem.first.val_2; // id
                T g_val = elem.second;

                if (std::abs(g_val) <= 0.000001) // because of the sparse container, we are not storing 0.0 values
                    continue;

                size_t r = rid_map[g_row];            // position in the list if ids, consecutive index of real id in the list of all ids
                size_t c = rid_map[g_col];            // position in the list if ids, consecutive index of real id in the list of all ids
                size_t ind = r * (r - 1) / 2 + c - 1; // r & c start from 1, but not from 0
                if (c > r)
                    ind = c * (c - 1) / 2 + r - 1;
                out_matr[ind] = g_val;
            }

            rid_map.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::map_to_matr(std::map<PedPair, T> &, std::vector<std::int64_t> &, evolm::smatrix<float> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::map_to_matr(std::map<PedPair, T> &, std::vector<std::int64_t> &, evolm::smatrix<float> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::map_to_matr(std::map<PedPair, T> &, std::vector<std::int64_t> &, evolm::smatrix<float> &)" << '\n';
            throw;
        }
    }

    template void Amat<float>::map_to_matr(std::map<PedPair, float> &in_amap, std::vector<std::int64_t> &in_ids, evolm::smatrix<float> &out_matr);
    template void Amat<double>::map_to_matr(std::map<PedPair, double> &in_amap, std::vector<std::int64_t> &in_ids, evolm::smatrix<double> &out_matr);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::map_to_matr(std::map<PedPair, T> &in_amap, std::vector<std::int64_t> &in_ids, evolm::matrix<T> &out_matr)
    {
        try
        {
            Utilities2 u;

            std::map<std::int64_t, std::int64_t> rid_map;

            if (in_ids.empty())
                throw std::string("Empty traced pedigree IDs!");

            u.get_RecodedIdMap(rid_map, in_ids); // list of real ids => std::map for the new consecutive list starting from 1

            if (rid_map.empty())
                throw std::string("Recoded IDs map is empty!");

            if (!out_matr.empty())
                out_matr.clear();

            out_matr.resize(in_ids.size());

            for (auto const &elem : in_amap)
            {
                std::int64_t g_row = elem.first.val_1; // id
                std::int64_t g_col = elem.first.val_2; // id
                T g_val = elem.second;

                size_t r = rid_map[g_row];            // position in the list if ids, consecutive index of real id in the list of all ids
                size_t c = rid_map[g_col];            // position in the list if ids, consecutive index of real id in the list of all ids
                size_t ind = r * (r - 1) / 2 + c - 1; // r & c start from 1, but not from 0
                if (c > r)
                    ind = c * (c - 1) / 2 + r - 1;
                out_matr[ind] = g_val;
            }

            rid_map.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::map_to_matr(std::map<PedPair, T> &, std::vector<std::int64_t> &, evolm::matrix<float> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::map_to_matr(std::map<PedPair, T> &, std::vector<std::int64_t> &, evolm::matrix<float> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::map_to_matr(std::map<PedPair, T> &, std::vector<std::int64_t> &, evolm::smatrix<float> &)" << '\n';
            throw;
        }
    }

    template void Amat<float>::map_to_matr(std::map<PedPair, float> &in_amap, std::vector<std::int64_t> &in_ids, evolm::matrix<float> &out_matr);
    template void Amat<double>::map_to_matr(std::map<PedPair, double> &in_amap, std::vector<std::int64_t> &in_ids, evolm::matrix<double> &out_matr);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::make_matrix(const std::string &ped_file, bool use_ainv)
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
        
            if (!inbrF.empty())
                inbrF.erase(inbrF.begin(), inbrF.end());

            std::map<PedPair, T> ainv;

            if (use_ainv)
                get_ainv(pedigree, ainv, true);
            else
                get_a(pedigree, ainv);

            pedigree.clear();

            pedigree_from_file.clear();

            birth_id_map.clear();

            pedID.clear();

            pedID.shrink_to_fit();

            size_t limit = 0.5 * (traced_pedID.size() - 1) * traced_pedID.size() + traced_pedID.size();
            if (!use_ainv)
                limit = traced_pedID.size() * traced_pedID.size();

            if (limit < ainv.size())
                throw std::string("The number of elements in calculated A(-1) matrix is higher than the number of traced IDs!!");

            data_sparsity = ( 1.0 - (double)ainv.size() / (double)limit ) * 100.0;
                        
            if ( data_sparsity >= sparsity_threshold )
                use_sparse = true; // the default is false

            if (use_sparse)
            {
                map_to_matr(ainv, traced_pedID, A_s);
                A_s.fwrite();
                iA_s = A_s;
                A_s.resize();

                IsEmpty.iA_s = false;
            }
            else
            {
                map_to_matr(ainv, traced_pedID, A);
                A.fwrite();
                iA = A;
                A.clear();

                IsEmpty.iA = false;
            }

            ainv.clear();

            id_iA = traced_pedID;

            traced_pedID.clear();
            traced_pedID.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::make_matrix(std::string &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::make_matrix(std::string &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::make_matrix(std::string &, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::make_matrix(const std::string &ped_file, bool use_ainv);
    template void Amat<double>::make_matrix(const std::string &ped_file, bool use_ainv);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::make_matrix(const std::string &ped_file, const std::string &g_file, bool use_ainv)
    {
        // Making A or A(-1) based on reduced pedigree traced on IDs in the file 'g_file'
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

            if (traced_pedID.empty())
                throw std::string("The list of traced IDs of the reduced pedigree is empty!");
            
            Utilities2 u;

            if ( !u.is_value_in_vect(traced_pedID, genotypedID) )
                throw std::string("There are genotyped IDs which are not part of the reduced A(-1) matrix!");

            if (!inbrF.empty())
                inbrF.erase(inbrF.begin(), inbrF.end());

            std::map<PedPair, T> r_ainv;

            if (use_ainv)
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

            size_t limit = 0.5 * (traced_pedID.size() - 1) * traced_pedID.size() + traced_pedID.size();
            if (!use_ainv)
                limit = traced_pedID.size() * traced_pedID.size();

            if (limit < r_ainv.size())
                throw std::string("The number of elements in calculated A(-1) matrix is higher than the number of traced IDs!!");

            data_sparsity = ( 1.0 - (double)r_ainv.size() / (double)limit ) * 100.0;

            if ( data_sparsity >= sparsity_threshold )
                use_sparse = true;

            if (use_sparse)
            {
                map_to_matr(r_ainv, traced_pedID, A_s);
                A_s.fwrite();
                irA_s = A_s;
                A_s.resize();

                IsEmpty.irA_s = false;
            }
            else
            {
                map_to_matr(r_ainv, traced_pedID, A);
                A.fwrite();
                irA = A;
                A.clear();

                IsEmpty.irA = false;
            }

            r_ainv.clear();

            id_irA = traced_pedID;

            traced_pedID.clear();
            traced_pedID.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::make_matrix(std::string &, std::string &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::make_matrix(std::string &, std::string &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::make_matrix(std::string &, std::string &, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::make_matrix(const std::string &ped_file, const std::string &g_file, bool use_ainv);
    template void Amat<double>::make_matrix(const std::string &ped_file, const std::string &g_file, bool use_ainv);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::make_all(const std::string &ped_file, const std::string &g_file)
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

            if (!inbrF.empty())
                inbrF.erase(inbrF.begin(), inbrF.end());

            std::map<PedPair, T> ainv;

            get_ainv(pedigree, ainv, true); // making A(-1)

            size_t limit = 0.5 * (traced_pedID.size() - 1) * traced_pedID.size() + traced_pedID.size();

            if (limit < ainv.size())
                throw std::string("The number of elements in calculated full A(-1) matrix is higher than the number of traced IDs!!");

            data_sparsity = ( 1.0 - (double)ainv.size() / (double)limit ) * 100.0;

            if ( data_sparsity >= sparsity_threshold )
                use_sparse = true;
            
            pedigree.clear();
            pedID.clear();
            pedID.shrink_to_fit();

            if (use_sparse)
            {
                map_to_matr(ainv, traced_pedID, A_s); // converting ainv map to sparse ( A or A(-1) ) matrix
                A_s.fwrite(); // 1. Wrire to binary
                iA_s = A_s;   // 2. Copy matrix by exchanging the internal binary file name
                A_s.resize();

                IsEmpty.iA_s = false;
            }
            else
            {
                map_to_matr(ainv, traced_pedID, A); // converting ainv map to ( A or A(-1) ) matrix
                A.fwrite(); // 1. Wrire to binary
                iA = A;     // 2. Copy matrix by exchanging the internal binary file name
                A.clear();  // 3. Clears the memory and gets new name for internal binary file

                IsEmpty.iA = false;
            }
            
            ainv.clear();

            id_iA = traced_pedID;

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

            if (traced_pedID.empty())
                throw std::string("The list of traced IDs of the reduced pedigree is empty!");
            
            Utilities2 u;

            if ( !u.is_value_in_vect(traced_pedID, genotypedID) )
                throw std::string("There are genotyped IDs which are not part of the reduced A(-1) matrix!");

            if (!inbrF.empty())
                inbrF.erase(inbrF.begin(), inbrF.end());

            std::map<PedPair, T> r_ainv;

            get_ainv(r_pedigree, r_ainv, true); // making reduced A(-1)

            if (r_ainv.size() == 0)
                throw std::string("The size of reduced pedigree for genotyped individuaals is empty!");

            pedigree_from_file.clear();

            limit = 0.5 * (traced_pedID.size() - 1) * traced_pedID.size() + traced_pedID.size();

            if (limit < r_ainv.size())
                throw std::string("The number of elements in calculated reduced A(-1) matrix is higher than the number of traced IDs!!");

            data_sparsity = ( 1.0 - (double)r_ainv.size() / (double)limit ) * 100.0;
            
            if ( data_sparsity >= sparsity_threshold )
                use_sparse = true;

            if (use_sparse)
            {
                map_to_matr(r_ainv, traced_pedID, A_s); // converting r_ainv map to A(-1) matrix
                A_s.fwrite(); // 1. Wrire to binary
                irA_s = A_s;  // 2. Copy matrix by just exchanging the internal binary file name
                A_s.resize();

                IsEmpty.irA_s = false;
            }
            else
            {
                map_to_matr(r_ainv, traced_pedID, A); // converting r_ainv map to A(-1) matrix
                A.fwrite(); // 1. Wrire to binary
                irA = A;    // 2. Copy matrix by just exchanging the internal binary file name
                A.clear();  // 3. Clears the memory and gets new name for internal binary file

                IsEmpty.irA = false;
            }

            r_ainv.clear();

            id_irA = traced_pedID;

            if ( traced_pedID.empty() )
                throw std::string("The list of traced IDs for genotyped individuals is empty!");

            // -----------------------------------------------------
            // This matrix assumed as always dense !
            // -------------------- A22 ----------------------------
            
            get_A22(r_pedigree, genotypedID);
            
            id_A22 = genotypedID;

            r_pedigree.clear();
            birth_id_map.clear();

            // Here we expect A22 is always dense, hense, the used matrix
            A.fwrite(); // 1. Wrire to binary
            A22 = A;    // 2. Copy matrix by just exchanging the internal binary file name
            A.clear();  // 3. Clears the memory and gets new name for internal binary file

            IsEmpty.A22 = false;

            // -----------------------------------------------------
            // The pathway depends on the sparsity of calculated irA_s or irA
            // -------------------- A22(-1) ------------------------

            if ( !IsEmpty.irA_s ) // irA_s is not empty (was processed through a sparse pathway)
            {
                irA_s.fread();

                get_iA22(irA_s, id_irA, genotypedID);

                irA_s.fwrite();

                A_s.fwrite(); // 1. Wrire to binary
                iA22_s = A_s; // 2. Copy matrix by just exchanging the internal binary file name
                A_s.resize(); // 3. Clears the memory and gets new name for internal binary file

                IsEmpty.iA22_s = false;
            }
            else
            {
                irA.fread();
                get_iA22(irA, id_irA, genotypedID);
                irA.fwrite();

                A.fwrite(); // 1. Wrire to binary
                iA22 = A;   // 2. Copy matrix by just exchanging the internal binary file name
                A.clear();  // 3. Clears the memory and gets new name for internal binary file

                IsEmpty.iA22 = false;
            }

            genotypedID.clear();
            genotypedID.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::make_all(std::string &, std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::make_all(std::string &, std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::make_all(std::string &, std::string &)" << '\n';
            throw;
        }
    }

    template void Amat<float>::make_all(const std::string &ped_file, const std::string &g_file);
    template void Amat<double>::make_all(const std::string &ped_file, const std::string &g_file);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::make_all(const std::string &ped_file, std::vector<std::int64_t> &g_ids)
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

            if (!inbrF.empty())
                inbrF.erase(inbrF.begin(), inbrF.end());

            std::map<PedPair, T> ainv;

            get_ainv(pedigree, ainv, true); // making A(-1)

            size_t limit = 0.5 * (traced_pedID.size() - 1) * traced_pedID.size() + traced_pedID.size();

            if (limit < ainv.size())
                throw std::string("The number of elements in calculated full A(-1) matrix is higher than the number of traced IDs!!");

            data_sparsity = ( 1.0 - (double)ainv.size() / (double)limit ) * 100.0;

            if ( data_sparsity >= sparsity_threshold )
                use_sparse = true;
            
            pedigree.clear();
            pedID.clear();
            pedID.shrink_to_fit();

            if (use_sparse)
            {
                map_to_matr(ainv, traced_pedID, A_s); // converting ainv map to sparse ( A or A(-1) ) matrix
                A_s.fwrite(); // 1. Wrire to binary
                iA_s = A_s;   // 2. Copy matrix by exchanging the internal binary file name
                A_s.resize();

                IsEmpty.iA_s = false;
            }
            else
            {
                map_to_matr(ainv, traced_pedID, A); // converting ainv map to ( A or A(-1) ) matrix
                A.fwrite(); // 1. Wrire to binary
                iA = A;     // 2. Copy matrix by exchanging the internal binary file name
                A.clear();  // 3. Clears the memory and gets new name for internal binary file

                IsEmpty.iA = false;
            }
            
            ainv.clear();

            id_iA = traced_pedID;

            // -----------------------------------------------------
            // ---------- reduced A(-1) ----------------------------

            std::vector<std::int64_t> genotypedID;
            std::map<PedPair, PedPair> r_pedigree;

            genotypedID = g_ids;

            if (genotypedID.empty())
                throw std::string("Cannot trace the reduced pedigree: ID's vector is empty!");

            if (!traced_pedID.empty())
                traced_pedID.erase(traced_pedID.begin(), traced_pedID.end());

            trace_pedigree(pedigree_from_file, r_pedigree, genotypedID); // tracing reduced pedigree for genotyped individuals

            if (traced_pedID.empty())
                throw std::string("The list of traced IDs of the reduced pedigree is empty!");
            
            Utilities2 u;

            if ( !u.is_value_in_vect(traced_pedID, genotypedID) )
                throw std::string("There are genotyped IDs which are not part of the reduced A(-1) matrix!");

            if (!inbrF.empty())
                inbrF.erase(inbrF.begin(), inbrF.end());

            std::map<PedPair, T> r_ainv;

            get_ainv(r_pedigree, r_ainv, true); // making reduced A(-1)

            if (r_ainv.size() == 0)
                throw std::string("The size of reduced pedigree for genotyped individuaals is empty!");

            pedigree_from_file.clear();

            limit = 0.5 * (traced_pedID.size() - 1) * traced_pedID.size() + traced_pedID.size();

            if (limit < r_ainv.size())
                throw std::string("The number of elements in calculated reduced A(-1) matrix is higher than the number of traced IDs!!");

            data_sparsity = ( 1.0 - (double)r_ainv.size() / (double)limit ) * 100.0;
            
            if ( data_sparsity >= sparsity_threshold )
                use_sparse = true;

            if (use_sparse)
            {
                map_to_matr(r_ainv, traced_pedID, A_s); // converting r_ainv map to A(-1) matrix
                A_s.fwrite(); // 1. Wrire to binary
                irA_s = A_s;  // 2. Copy matrix by just exchanging the internal binary file name
                A_s.resize();

                IsEmpty.irA_s = false;
            }
            else
            {
                map_to_matr(r_ainv, traced_pedID, A); // converting r_ainv map to A(-1) matrix
                A.fwrite(); // 1. Wrire to binary
                irA = A;    // 2. Copy matrix by just exchanging the internal binary file name
                A.clear();  // 3. Clears the memory and gets new name for internal binary file

                IsEmpty.irA = false;
            }

            r_ainv.clear();

            id_irA = traced_pedID;

            if ( traced_pedID.empty() )
                throw std::string("The list of traced IDs for genotyped individuals is empty!");

            // -----------------------------------------------------
            // This matrix assumed as always dense !
            // -------------------- A22 ----------------------------
            
            get_A22(r_pedigree, genotypedID);
            
            id_A22 = genotypedID;

            r_pedigree.clear();
            birth_id_map.clear();

            // Here we expect A22 is always dense, hense, the used matrix
            A.fwrite(); // 1. Wrire to binary
            A22 = A;    // 2. Copy matrix by just exchanging the internal binary file name
            A.clear();  // 3. Clears the memory and gets new name for internal binary file

            IsEmpty.A22 = false;

            // -----------------------------------------------------
            // The pathway depends on the sparsity of calculated irA_s or irA
            // -------------------- A22(-1) ------------------------

            if ( !IsEmpty.irA_s ) // irA_s is not empty (was processed through a sparse pathway)
            {
                irA_s.fread();

                get_iA22(irA_s, id_irA, genotypedID);

                irA_s.fwrite();

                A_s.fwrite(); // 1. Wrire to binary
                iA22_s = A_s; // 2. Copy matrix by just exchanging the internal binary file name
                A_s.resize(); // 3. Clears the memory and gets new name for internal binary file

                IsEmpty.iA22_s = false;
            }
            else
            {
                irA.fread();
                get_iA22(irA, id_irA, genotypedID);
                irA.fwrite();

                A.fwrite(); // 1. Wrire to binary
                iA22 = A;   // 2. Copy matrix by just exchanging the internal binary file name
                A.clear();  // 3. Clears the memory and gets new name for internal binary file

                IsEmpty.iA22 = false;
            }

            genotypedID.clear();
            genotypedID.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::make_all(std::string &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::make_all(std::string &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::make_all(std::string &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    template void Amat<float>::make_all(const std::string &ped_file, std::vector<std::int64_t> &g_ids);
    template void Amat<double>::make_all(const std::string &ped_file, std::vector<std::int64_t> &g_ids);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::make_matrix_forgenotyped(const std::string &ped_file, const std::string &g_file, bool use_ainv)
    {
        // Making matrices required for scalling G (GRM) matrix: A22 or A22(-1)
        try
        {
            // -----------------------------------------------------
            std::vector<std::int64_t> pedID;
            std::map<PedPair, PedPair> pedigree_from_file;

            fread_pedigree(ped_file, pedigree_from_file, pedID);

            if (pedID.empty())
                throw std::string("Empty pedigree IDs!");

            if (pedigree_from_file.empty())
                throw std::string("File provided pedigree is empty!");

            if (!traced_pedID.empty())
                traced_pedID.erase(traced_pedID.begin(), traced_pedID.end());

            // -----------------------------------------------------
            // ---------- reduced A(-1) ----------------------------

            std::vector<std::int64_t> genotypedID;
            std::map<PedPair, PedPair> r_pedigree;

            fread_genotyped_id(g_file, genotypedID);

            if (genotypedID.empty())
                throw std::string("Cannot trace the reduced pedigree: ID's vector is empty!");

            trace_pedigree(pedigree_from_file, r_pedigree, genotypedID); // tracing reduced pedigree for genotyped individuals

            if (traced_pedID.empty())
                throw std::string("The list of traced IDs of the reduced pedigree is empty!");
            
            Utilities2 u;

            if ( !u.is_value_in_vect(traced_pedID, genotypedID) )
                throw std::string("There are genotyped IDs which are not part of the reduced A(-1) matrix!");

            if (!inbrF.empty())
                inbrF.erase(inbrF.begin(), inbrF.end());

            std::map<PedPair, T> r_ainv;

            get_ainv(r_pedigree, r_ainv, true); // making reduced A(-1)

            pedigree_from_file.clear();

            size_t limit = 0.5 * (traced_pedID.size() - 1) * traced_pedID.size() + traced_pedID.size();

            if (limit < r_ainv.size())
                throw std::string("The number of elements in calculated reduced A(-1) matrix is higher than the number of traced IDs!!");

            data_sparsity = ( 1.0 - (double)r_ainv.size() / (double)limit ) * 100.0;

            if ( data_sparsity >= sparsity_threshold )
                use_sparse = true;

            if ( use_ainv )
            {
                if (use_sparse)
                {
                    map_to_matr(r_ainv, traced_pedID, A_s); // converting r_ainv map to A(-1) matrix
                    A_s.fwrite(); // 1. Wrire to binary
                    irA_s = A_s;  // 2. Copy matrix by just exchanging the internal binary file name
                    A_s.resize();

                    IsEmpty.irA_s = false;
                }
                else
                {
                    map_to_matr(r_ainv, traced_pedID, A); // converting r_ainv map to A(-1) matrix
                    A.fwrite(); // 1. Wrire to binary
                    irA = A;    // 2. Copy matrix by just exchanging the internal binary file name
                    A.clear();  // 3. Clears the memory and gets new name for internal binary file

                    IsEmpty.irA = false;
                }

                r_ainv.clear();

                id_irA = traced_pedID;

                // -----------------------------------------------------
                // The pathway depends on the sparsity of calculated irA_s or irA
                // -------------------- A22(-1) ------------------------

                if ( !IsEmpty.irA_s ) // irA_s is not empty (was processed through a sparse pathway)
                {
                    irA_s.fread();
                    
                    get_iA22(irA_s, id_irA, genotypedID);

                    irA_s.fwrite();

                    A_s.fwrite(); // 1. Wrire to binary
                    iA22_s = A_s; // 2. Copy matrix by just exchanging the internal binary file name
                    
                    A_s.resize(); // 3. Clears the memory and gets new name for internal binary file

                    IsEmpty.iA22_s = false;
                }
                else
                {
                    irA.fread();

                    get_iA22(irA, id_irA, genotypedID);
                    irA.fwrite();

                    A.fwrite(); // 1. Wrire to binary
                    iA22 = A;   // 2. Copy matrix by just exchanging the internal binary file name
                    A.clear();  // 3. Clears the memory and gets new name for internal binary file

                    IsEmpty.iA22 = false;
                }

                id_A22 = genotypedID;
            }
            else
            {
                // -----------------------------------------------------
                // This matrix assumed as always dense !
                // -------------------- A22 ----------------------------
                
                get_A22(r_pedigree, genotypedID);
                
                id_A22 = genotypedID;

                r_pedigree.clear();
                birth_id_map.clear();

                // Here we expect A22 is always dense, hense, the used matrix
                A.fwrite(); // 1. Wrire to binary
                A22 = A;    // 2. Copy matrix by just exchanging the internal binary file name
                A.clear();  // 3. Clears the memory and gets new name for internal binary file

                IsEmpty.A22 = false;
            }

            genotypedID.clear();
            genotypedID.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::make_matrix_forgenotyped(std::string &, std::string &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::make_matrix_forgenotyped(std::string &, std::string &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::make_matrix_forgenotyped(std::string &, std::string &, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::make_matrix_forgenotyped(const std::string &ped_file, const std::string &g_file, bool use_ainv);
    template void Amat<double>::make_matrix_forgenotyped(const std::string &ped_file, const std::string &g_file, bool use_ainv);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::make_matrix_forgenotyped(const std::string &ped_file, std::vector<std::int64_t> &genotyped_ids, bool use_ainv)
    {
        // Making matrices required for scalling G (GRM) matrix: A22 or A22(-1)
        try
        {
            // -----------------------------------------------------
            std::vector<std::int64_t> pedID;
            std::map<PedPair, PedPair> pedigree_from_file;

            fread_pedigree(ped_file, pedigree_from_file, pedID);

            if (pedID.empty())
                throw std::string("Empty pedigree IDs!");

            if (pedigree_from_file.empty())
                throw std::string("File provided pedigree is empty!");

            if (!traced_pedID.empty())
                traced_pedID.erase(traced_pedID.begin(), traced_pedID.end());

            // -----------------------------------------------------
            // ---------- reduced A(-1) ----------------------------

            std::vector<std::int64_t> genotypedID;
            std::map<PedPair, PedPair> r_pedigree;

            genotypedID = genotyped_ids;

            if (genotypedID.empty())
                throw std::string("Cannot trace the reduced pedigree: ID's vector is empty!");

            trace_pedigree(pedigree_from_file, r_pedigree, genotypedID); // tracing reduced pedigree for genotyped individuals

            if (traced_pedID.empty())
                throw std::string("The list of traced IDs of the reduced pedigree is empty!");
            
            Utilities2 u;

            if ( !u.is_value_in_vect(traced_pedID, genotypedID) )
                throw std::string("There are genotyped IDs which are not part of the reduced A(-1) matrix!");

            if (!inbrF.empty())
                inbrF.erase(inbrF.begin(), inbrF.end());

            std::map<PedPair, T> r_ainv;

            get_ainv(r_pedigree, r_ainv, true); // making reduced A(-1)

            pedigree_from_file.clear();

            size_t limit = 0.5 * (traced_pedID.size() - 1) * traced_pedID.size() + traced_pedID.size();

            if (limit < r_ainv.size())
                throw std::string("The number of elements in calculated reduced A(-1) matrix is higher than the number of traced IDs!!");

            data_sparsity = ( 1.0 - (double)r_ainv.size() / (double)limit ) * 100.0;

            if ( data_sparsity >= sparsity_threshold )
                use_sparse = true;

            if ( use_ainv )
            {
                if (use_sparse)
                {
                    map_to_matr(r_ainv, traced_pedID, A_s); // converting r_ainv map to A(-1) matrix
                    A_s.fwrite(); // 1. Wrire to binary
                    irA_s = A_s;  // 2. Copy matrix by just exchanging the internal binary file name
                    A_s.resize();

                    IsEmpty.irA_s = false;
                }
                else
                {
                    map_to_matr(r_ainv, traced_pedID, A); // converting r_ainv map to A(-1) matrix
                    A.fwrite(); // 1. Wrire to binary
                    irA = A;    // 2. Copy matrix by just exchanging the internal binary file name
                    A.clear();  // 3. Clears the memory and gets new name for internal binary file

                    IsEmpty.irA = false;
                }

                r_ainv.clear();

                id_irA = traced_pedID;

                // -----------------------------------------------------
                // The pathway depends on the sparsity of calculated irA_s or irA
                // -------------------- A22(-1) ------------------------

                if ( !IsEmpty.irA_s ) // irA_s is not empty (was processed through a sparse pathway)
                {
                    irA_s.fread();

                    get_iA22(irA_s, id_irA, genotypedID);

                    irA_s.fwrite();

                    A_s.fwrite(); // 1. Wrire to binary
                    iA22_s = A_s; // 2. Copy matrix by just exchanging the internal binary file name
                    A_s.resize(); // 3. Clears the memory and gets new name for internal binary file

                    IsEmpty.iA22_s = false;
                }
                else
                {
                    irA.fread();

                    get_iA22(irA, id_irA, genotypedID);
                    irA.fwrite();

                    A.fwrite(); // 1. Wrire to binary
                    iA22 = A;   // 2. Copy matrix by just exchanging the internal binary file name
                    A.clear();  // 3. Clears the memory and gets new name for internal binary file

                    IsEmpty.iA22 = false;
                }

                id_A22 = genotypedID;
            }
            else
            {
                // -----------------------------------------------------
                // This matrix assumed as always dense !
                // -------------------- A22 ----------------------------
                
                get_A22(r_pedigree, genotypedID);
                
                id_A22 = genotypedID;

                r_pedigree.clear();
                birth_id_map.clear();

                // Here we expect A22 is always dense, hense, the used matrix
                A.fwrite(); // 1. Wrire to binary
                A22 = A;    // 2. Copy matrix by just exchanging the internal binary file name
                A.clear();  // 3. Clears the memory and gets new name for internal binary file

                IsEmpty.A22 = false;
            }

            genotypedID.clear();
            genotypedID.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::make_matrix_forgenotyped(std::string &, std::vector<std::int64_t> &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::make_matrix_forgenotyped(std::string &, std::vector<std::int64_t> &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::make_matrix_forgenotyped(std::string &, std::vector<std::int64_t> &, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::make_matrix_forgenotyped(const std::string &ped_file, std::vector<std::int64_t> &genotyped_ids, bool use_ainv);
    template void Amat<double>::make_matrix_forgenotyped(const std::string &ped_file, std::vector<std::int64_t> &genotyped_ids, bool use_ainv);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::get_iA22(evolm::matrix<T> &full_matr, std::vector<std::int64_t> &matr_ids, std::vector<std::int64_t> &selected_ids)
    {
        try
        {
            Utilities2 u;

            evolm::matrix<size_t> shapeofh;
            shapeofh = full_matr.shape();

            size_t expected_size = 0.5 * (matr_ids.size() - 1) * matr_ids.size() + matr_ids.size();

            if (shapeofh[0] != shapeofh[1])
                throw std::string("The passed matrix has wrong dimension: number of raws is not the same as number of columns!");

            if (shapeofh[0] != matr_ids.size())
                throw std::string("The passed matrix has wrong dimension!");

            if (full_matr.size() > expected_size)
                throw std::string("The number of elements in the passed matrix is greater then expected!");

            if (selected_ids.size() > matr_ids.size())
                throw std::string("The number of elements in the passed selected IDs array is greater then the number of IDs in the passed matrix!");

            // --------------------------------
            if (!A.empty())
            {
                A.fclear();
                A.clear();
            }

            // Create the list of IDs which are in matr_ids but not in selected_ids vectors

            std::vector<std::int64_t> not_selected_ids;

            for (size_t i = 0; i < matr_ids.size(); i++)
            {
                int res = u.find_invect(selected_ids, matr_ids[i]);
                if (res == -1)
                    not_selected_ids.push_back(matr_ids[i]);
            }

            if ( not_selected_ids.empty() )
                throw std::string("The vector of individuals which are in pedigree but not in genotyped IDs list is empty!");

            if ((selected_ids.size() + not_selected_ids.size()) != matr_ids.size())
                throw std::string("The sum of IDs from from two vectors is not equal to number of IDs in the matrix!");

            // ------------------------------------------------

            // Next steps: A22(-1) = A22 - A21 * A11(-1) * A12;

            evolm::matrix<T> a22; // in order to differrentiate from the 'global' A22
            evolm::matrix<T> A11;
            evolm::matrix<T> A21;
            evolm::matrix<T> A12;

            // -------------------- A11 -----------------------

            A11.resize(not_selected_ids.size());

            std::vector<size_t> non_selected_pos;
            for (size_t i = 0; i < not_selected_ids.size(); i++)
                non_selected_pos.push_back(u.find_invect(matr_ids, not_selected_ids[i]));

            std::vector<size_t> selected_pos;
            for (size_t i = 0; i < selected_ids.size(); i++)
                selected_pos.push_back(u.find_invect(matr_ids, selected_ids[i]));

#pragma omp parallel for
            for (size_t i = 0; i < not_selected_ids.size(); i++)
            {
                size_t pos_i = non_selected_pos[i];

                for (size_t j = 0; j <= i; j++)
                {
                    size_t pos_j = non_selected_pos[j];

                    if (pos_i >= pos_j)
                        A11(i, j) = full_matr(pos_i, pos_j);
                    else
                        A11(i, j) = full_matr(pos_j, pos_i);
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
                size_t pos_i = selected_pos[i];

                for (size_t j = 0; j < not_selected_ids.size(); j++)
                {
                    size_t pos_j = non_selected_pos[j];

                    if (pos_i >= pos_j)
                        A21(i, j) = full_matr(pos_i, pos_j);
                    else
                        A21(i, j) = full_matr(pos_j, pos_i);
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
            evolm::matrix<T> res;

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

            a22.resize(selected_ids.size());

#pragma omp parallel for
            for (size_t i = 0; i < selected_ids.size(); i++)
            {
                size_t pos_i = selected_pos[i];

                for (size_t j = 0; j <= i; j++)
                {
                    size_t pos_j = selected_pos[j];

                    if (pos_i >= pos_j)
                        a22(i, j) = full_matr(pos_i, pos_j);
                    else
                        a22(i, j) = full_matr(pos_j, pos_i);
                }
            }

            non_selected_pos.clear();
            non_selected_pos.shrink_to_fit();
            selected_pos.clear();
            selected_pos.shrink_to_fit();

            // ------------------------------------------------
            //
            // ------------------ A22 - res -------------------

            evolm::matrix<size_t> shapeofa22;
            shapeofa22 = a22.shape();

            // A.resize(shapeofa22[0]); 1

            res.rectosym();

#pragma omp parallel for
            for (size_t i = 0; i < a22.size(); i++)
            {
                // A[i] = a22[i] - res[i]; 2
                res[i] = a22[i] - res[i]; // 2
            }

            a22.fclear();
            a22.clear();

            A = res; // 1

            res.fclear();
            res.clear();
            // ------------------------------------------------
            shapeofh.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::get_iA22(evolm::matrix<T>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::get_iA22(evolm::matrix<T>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::get_iA22(evolm::matrix<T>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::get_iA22(evolm::matrix<float> &full_matr, std::vector<std::int64_t> &matr_ids, std::vector<std::int64_t> &selected_ids);
    template void Amat<double>::get_iA22(evolm::matrix<double> &full_matr, std::vector<std::int64_t> &matr_ids, std::vector<std::int64_t> &selected_ids);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::get_iA22(evolm::smatrix<T> &full_matr, std::vector<std::int64_t> &matr_ids, std::vector<std::int64_t> &selected_ids)
    {
        try
        {
            // high memory use, though, no more than required for making G matrix; allows fast computations
            making_iA22_d(full_matr, matr_ids, selected_ids); // initial data and final output in SPARSE format, all internal operations in DENSE format

            // low memory use, though, takes much more computation time than making_iA22_d(...) method
            //making_iA22_s(full_matr, matr_ids, selected_ids); // initial data and final output in SPARSE format, as well as all internal operations
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::get_iA22(evolm::smatrix<T>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::get_iA22(evolm::smatrix<T>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::get_iA22(evolm::smatrix<T>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::get_iA22(evolm::smatrix<float> &full_matr, std::vector<std::int64_t> &matr_ids, std::vector<std::int64_t> &selected_ids);
    template void Amat<double>::get_iA22(evolm::smatrix<double> &full_matr, std::vector<std::int64_t> &matr_ids, std::vector<std::int64_t> &selected_ids);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::making_iA22_s(evolm::smatrix<T> &full_matr, std::vector<std::int64_t> &matr_ids, std::vector<std::int64_t> &selected_ids)
    {
        try
        {
            Utilities2 u;

            size_t n_rows = full_matr.nrows();
            size_t n_cols = full_matr.ncols();

            if ( matr_ids.size() == selected_ids.size() )
                throw std::string("The genotyped individuals have no ancestors in pedigree => the list of population IDs is equal to the list of genotyped IDs!");

            if (n_rows != n_cols)
                throw std::string("The passed matrix has wrong dimension: number of raws is not the same as number of columns!");

            if (n_rows != matr_ids.size())
                throw std::string("The passed matrix has wrong dimension!");

            if (full_matr.size() > full_matr.max_key() + 1)
                throw std::string("The number of elements in the passed matrix is greater then expected!");

            if (selected_ids.size() > matr_ids.size())
                throw std::string("The number of elements in the passed selected IDs array is greater then the number of IDs in the passed matrix!");

            if (!u.is_value_in_vect(matr_ids, selected_ids))
                throw std::string("There are IDs in the selected IDs array which are not part of the passed matrix!");

            // --------------------------------
            // Create the list of IDs which are in matr_ids but not in selected_ids vectors

            std::vector<std::int64_t> not_selected_ids;

            for (size_t i = 0; i < matr_ids.size(); i++)
            {
                int res = u.find_invect(selected_ids, matr_ids[i]);
                if (res == -1)
                    not_selected_ids.push_back(matr_ids[i]);
            }

            if ( not_selected_ids.empty() )
                throw std::string("The list of individuals which are in pedigree but not in genotyped IDs list is empty!");

            if ((selected_ids.size() + not_selected_ids.size()) != matr_ids.size())
                throw std::string("The sum of IDs from from two vectors is not equal to number of IDs in the matrix!");

            // ------------------------------------------------
            
            // Next steps: A22(-1) = A22 - A21 * A11(-1) * A12;
            
            // -------------------- A11 -----------------------
            
            evolm::smatrix<T> a22; // in order to differrentiate from the 'global' A22
            evolm::smatrix<T> A11;
            evolm::smatrix<T> A21;
            evolm::smatrix<T> A12;
            evolm::matrix<T> dense_A11; // this is the temporal storage for making inverse
            
            dense_A11.resize(not_selected_ids.size());
            
            std::vector<size_t> non_selected_pos;
            for (size_t i = 0; i < not_selected_ids.size(); i++)
                non_selected_pos.push_back(u.find_invect(matr_ids, not_selected_ids[i]));

            std::vector<size_t> selected_pos;
            for (size_t i = 0; i < selected_ids.size(); i++)
                selected_pos.push_back(u.find_invect(matr_ids, selected_ids[i]));

            T zerro_value = (T)0;

#pragma omp parallel for
            for (size_t i = 0; i < not_selected_ids.size(); i++)
            {
                size_t pos_i = non_selected_pos[i];
                
                T value = (T)0;
                
                for (size_t j = 0; j <= i; j++)
                {
                    size_t pos_j = non_selected_pos[j];

                    if (pos_i >= pos_j)
                    {
                        value = full_matr.get_nonzero(pos_i, pos_j);

                        if (value != zerro_value)
                            dense_A11(i, j) = value;
                        else
                            dense_A11(i, j) = zerro_value;                        
                    }
                    else
                    {
                        value = full_matr.get_nonzero(pos_j, pos_i);

                        if (value != zerro_value)
                            dense_A11(i, j) = value;
                        else
                            dense_A11(i, j) = zerro_value;                        
                    }
                }
            }
            
            dense_A11.symtorec();

            dense_A11.invert();

            // converting A11 to sparse matrix

            A11.resize(not_selected_ids.size(), not_selected_ids.size());

            for (size_t i = 0; i < not_selected_ids.size(); i++)
            {
                for (size_t j = 0; j < not_selected_ids.size(); j++)
                {
                    if (dense_A11(i, j) != 0.0)
                        A11(i, j) = dense_A11(i, j);
                }
            }

            dense_A11.fclear(); // we don't need the dense dense_A11 any more
            dense_A11.clear();

            A11.fwrite();

            // ------------------------------------------------
            // -------------------- A21 -----------------------

            A21.resize(selected_ids.size(), not_selected_ids.size());

            for (size_t i = 0; i < selected_ids.size(); i++)
            {
                size_t pos_i = selected_pos[i];
                
                T value = (T)0;

                for (size_t j = 0; j < not_selected_ids.size(); j++)
                {
                    size_t pos_j = non_selected_pos[j];

                    if (pos_i >= pos_j)
                    {
                        value = full_matr.get_nonzero(pos_i, pos_j);

                        if (value != zerro_value)
                            A21(i, j) = value;
                    }
                    else
                    {
                        value = full_matr.get_nonzero(pos_j, pos_i);

                        if (value != zerro_value)
                            A21(i, j) = value;
                    }
                }
            }

            // ------------------------------------------------
            // -------------------- A12 -----------------------

            A12 = A21;

            A12.transpose();

            A12.fwrite();

            // ------------------------------------------------
            // --------------- A21 * A11(-1) * A12 ------------

            A11.fread();
            
            evolm::smatrix<T> res;

            res = A21 * A11;

            A11.fclear();
            A11.clear();

            A21.fclear();
            A21.clear();

            A12.fread();

            res = res * A12;

            A12.fclear();
            A12.clear();

            res.rectosym();

            // ------------------------------------------------
            // -------------------- A22 -----------------------

            a22.resize(selected_ids.size());

            for (size_t i = 0; i < selected_ids.size(); i++)
            {
                size_t pos_i = selected_pos[i];
                T value = (T)0;
                T zerro_value = (T)0;

                for (size_t j = 0; j <= i; j++)
                {
                    size_t pos_j = selected_pos[j];

                    if (pos_i >= pos_j)
                    {
                        value = full_matr.get_nonzero(pos_i, pos_j);

                        if (value != zerro_value)
                            a22(i, j) = value;
                    }
                    else
                    {
                        value = full_matr.get_nonzero(pos_j, pos_i);

                        if (value != zerro_value)
                            a22(i, j) = value;
                    }
                }
            }

            non_selected_pos.clear();
            non_selected_pos.shrink_to_fit();
            selected_pos.clear();
            selected_pos.shrink_to_fit();

            // ------------------------------------------------
            // ------------------ A22 - res -------------------

            if (!A_s.empty())
                A_s.resize();

            A_s = a22 - res;

            a22.fclear();
            a22.clear();

            res.fclear();
            res.clear();

        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::making_iA22_s(evolm::smatrix<T>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::making_iA22_s(evolm::smatrix<T>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::making_iA22_s(evolm::smatrix<T>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::making_iA22_s(evolm::smatrix<float> &full_matr, std::vector<std::int64_t> &matr_ids, std::vector<std::int64_t> &selected_ids);
    template void Amat<double>::making_iA22_s(evolm::smatrix<double> &full_matr, std::vector<std::int64_t> &matr_ids, std::vector<std::int64_t> &selected_ids);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::making_iA22_d(evolm::smatrix<T> &full_matr, std::vector<std::int64_t> &matr_ids, std::vector<std::int64_t> &selected_ids)
    {
        try
        {
            Utilities2 u;

            size_t n_rows = full_matr.nrows();
            size_t n_cols = full_matr.ncols();

            if ( matr_ids.size() == selected_ids.size() )
                throw std::string("The genotyped individuals have no ancestors in pedigree => the list of population IDs is equal to the list of genotyped IDs!");

            if (n_rows != n_cols)
                throw std::string("The passed matrix has wrong dimension: number of raws is not the same as number of columns!");

            if (n_rows != matr_ids.size())
                throw std::string("The passed matrix has wrong dimension!");

            if (full_matr.size() > full_matr.max_key() + 1)
                throw std::string("The number of elements in the passed matrix is greater then expected!");

            if (selected_ids.size() > matr_ids.size())
                throw std::string("The number of elements in the passed selected IDs array is greater then the number of IDs in the passed matrix!");

            if (!u.is_value_in_vect(matr_ids, selected_ids))
                throw std::string("There are IDs in the selected IDs array which are not part of the passed matrix!");

            // --------------------------------
            
            // Create the list of IDs which are in matr_ids but not in selected_ids vectors

            std::vector<std::int64_t> not_selected_ids;

            for (size_t i = 0; i < matr_ids.size(); i++)
            {
                int res = u.find_invect(selected_ids, matr_ids[i]);
                if (res == -1)
                    not_selected_ids.push_back(matr_ids[i]);
            }

            if ( not_selected_ids.empty() )
                throw std::string("The list of individuals which are in pedigree but not in genotyped IDs list is empty!");

            if ((selected_ids.size() + not_selected_ids.size()) != matr_ids.size())
                throw std::string("The sum of IDs from from two vectors is not equal to number of IDs in the matrix!");

            // ------------------------------------------------
            
            // Next steps: A22(-1) = A22 - A21 * A11(-1) * A12;
            
            // -------------------- A11 -----------------------            
            
            evolm::matrix<T> a22; // in order to differrentiate from the 'global' A22
            evolm::matrix<T> A11;
            evolm::matrix<T> A21;
            evolm::matrix<T> A12;

            A11.resize(not_selected_ids.size());

            std::vector<size_t> non_selected_pos;
            for (size_t i = 0; i < not_selected_ids.size(); i++)
                non_selected_pos.push_back(u.find_invect(matr_ids, not_selected_ids[i]));

            std::vector<size_t> selected_pos;
            for (size_t i = 0; i < selected_ids.size(); i++)
                selected_pos.push_back(u.find_invect(matr_ids, selected_ids[i]));

#pragma omp parallel for
            for (size_t i = 0; i < not_selected_ids.size(); i++)
            {
                size_t pos_i = non_selected_pos[i];
                T value = (T)0;
                //T zerro_value = (T)0;

                for (size_t j = 0; j <= i; j++)
                {
                    size_t pos_j = non_selected_pos[j];

                    if (pos_i >= pos_j)
                        value = full_matr.get_nonzero(pos_i, pos_j);
                    else
                        value = full_matr.get_nonzero(pos_j, pos_i);

                    A11(i, j) = value;
                }
            }

            A11.symtorec();

            A11.invert();

            A11.fwrite();

            // ------------------------------------------------
            // -------------------- A21 -----------------------

            A21.resize(selected_ids.size(), not_selected_ids.size());

#pragma omp parallel for
            for (size_t i = 0; i < selected_ids.size(); i++)
            {
                size_t pos_i = selected_pos[i];
                T value = (T)0;
                //T zerro_value = (T)0;

                for (size_t j = 0; j < not_selected_ids.size(); j++)
                {
                    size_t pos_j = non_selected_pos[j];

                    if (pos_i >= pos_j)
                        value = full_matr.get_nonzero(pos_i, pos_j);
                    else
                        value = full_matr.get_nonzero(pos_j, pos_i);
                    
                    A21(i, j) = value;
                }
            }

            // ------------------------------------------------
            // -------------------- A12 -----------------------

            A12 = A21;

            A12.transpose();

            A12.fwrite();

            // ------------------------------------------------
            // --------------- A21 * A11(-1) * A12 ------------

            A11.fread();
            evolm::matrix<T> res;

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
            // -------------------- A22 -----------------------

            a22.resize(selected_ids.size());

#pragma omp parallel for
            for (size_t i = 0; i < selected_ids.size(); i++)
            {
                size_t pos_i = selected_pos[i];
                T value = (T)0;

                for (size_t j = 0; j <= i; j++)
                {
                    size_t pos_j = selected_pos[j];

                    if (pos_i >= pos_j)
                        value = full_matr.get_nonzero(pos_i, pos_j);
                    else
                        value = full_matr.get_nonzero(pos_j, pos_i);
                    
                    a22(i, j) = value;
                }
            }

            non_selected_pos.clear();
            non_selected_pos.shrink_to_fit();
            selected_pos.clear();
            selected_pos.shrink_to_fit();

            // ------------------------------------------------
            // ------------------ A22 - res -------------------

            evolm::matrix<size_t> shapeofa22;
            shapeofa22 = a22.shape();

            res.rectosym();

#pragma omp parallel for
            for (size_t i = 0; i < a22.size(); i++)
                res[i] = a22[i] - res[i];

            a22.fclear();
            a22.clear();

            A_s.resize(shapeofa22[0]);
            for(size_t i = 0; i < res.size(); i++)
                if( res[i] != (T)0 )
                    A_s[i] = res[i];

            res.fclear();
            res.clear();

        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::making_iA22_d(evolm::smatrix<T>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::making_iA22_d(evolm::smatrix<T>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::making_iA22_d(evolm::smatrix<T>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::making_iA22_d(evolm::smatrix<float> &full_matr, std::vector<std::int64_t> &matr_ids, std::vector<std::int64_t> &selected_ids);
    template void Amat<double>::making_iA22_d(evolm::smatrix<double> &full_matr, std::vector<std::int64_t> &matr_ids, std::vector<std::int64_t> &selected_ids);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::fread_pedigree(const std::string &ped_file, std::map<PedPair, PedPair> &out_ped, std::vector<std::int64_t> &out_ids)
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

            if (!ped.good())
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
            std::cerr << "Exception in Amat<T>::fread_pedigree(const std::string &, std::map<PedPair, PedPair> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::fread_pedigree(const std::string &, std::map<PedPair, PedPair> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::fread_pedigree(const std::string &, std::map<PedPair, PedPair> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    template void Amat<float>::fread_pedigree(const std::string &ped_file, std::map<PedPair, PedPair> &out_ped, std::vector<std::int64_t> &out_ids);
    template void Amat<double>::fread_pedigree(const std::string &ped_file, std::map<PedPair, PedPair> &out_ped, std::vector<std::int64_t> &out_ids);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::fread_genotyped_id(const std::string &g_file, std::vector<std::int64_t> &out_ids)
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

            if (!ped.good())
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

                    // if (static_cast<std::int64_t>(tmp_list[1]) != 0 && static_cast<std::int64_t>(tmp_list[0]) != t_gtyp)
                    //     coreID.push_back(static_cast<std::int64_t>(tmp_list[0]));

                    t_gtyp = static_cast<std::int64_t>(tmp_list[0]);
                    // t_gcor = static_cast<std::int64_t>(tmp_list[1]);

                    tmp_list.erase(tmp_list.begin(), tmp_list.end());
                }
            }

            ped.close();

            if (!u.is_unique(out_ids))
                out_ids.erase(unique(out_ids.begin(), out_ids.end()), out_ids.end());

            // if (!u.is_unique(coreID))
            //     coreID.erase(unique(coreID.begin(), coreID.end()), coreID.end());
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::fread_genotyped_id(const std::string &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::fread_genotyped_id(const std::string &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::fread_genotyped_id(const std::string &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    template void Amat<float>::fread_genotyped_id(const std::string &g_file, std::vector<std::int64_t> &out_ids);
    template void Amat<double>::fread_genotyped_id(const std::string &g_file, std::vector<std::int64_t> &out_ids);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::trace_pedigree(std::map<PedPair, PedPair> &in_ped, std::map<PedPair, PedPair> &out_ped, std::vector<std::int64_t> &traced_id)
    {
        try
        {
            Utilities2 u;

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

            std::vector<size_t> loads;
            thread_loads(traced_id, loads);

            std::vector<std::vector<std::int64_t>> res_id;
            std::vector<std::map<PedPair, PedPair>> res_ped;
            std::vector<std::thread> vect_thr;

            for(size_t i = 0; i < loads.size(); i++)
            {
                std::vector<std::int64_t> id;
                std::map<PedPair, PedPair> ped;
                res_id.push_back(id);
                res_ped.push_back(ped);
            }

            for(size_t i = 0; i < loads.size(); i++)
            {
                vect_thr.emplace_back( &Amat::trace_operation,
                                        this,
                                        std::ref(days),
                                        std::ref(ids),
                                        std::ref(sire),
                                        std::ref(dame),
                                        std::ref(traced_id),
                                        std::ref(res_id[i]),
                                        std::ref(res_ped[i]),
                                        std::ref(loads),
                                        i);
            }

            for(size_t i = 0; i < loads.size(); i++)
                vect_thr[i].join();

            for(size_t i = 0; i < loads.size(); i++)
            {
                traced_pedID.insert( traced_pedID.end(), res_id[i].begin(), res_id[i].end() );
                res_id[i].clear();
            }

            for(size_t i = 0; i < loads.size(); i++)
            {
                out_ped.insert( res_ped[i].begin(), res_ped[i].end() );
                res_ped[i].clear();
            }

            if (!u.is_unique(traced_pedID)) // due to multiple threads traced_pedID is not unique
                traced_pedID.erase(unique(traced_pedID.begin(), traced_pedID.end()), traced_pedID.end()); // here the vector should be sorted and unique

            days.clear();
            ids.clear();
            sire.clear();
            dame.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::trace_pedigree(std::map<PedPair, PedPair> &, std::map<PedPair, PedPair> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::trace_pedigree(std::map<PedPair, PedPair> &, std::map<PedPair, PedPair> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::trace_pedigree(std::map<PedPair, PedPair> &, std::map<PedPair, PedPair> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    template void Amat<float>::trace_pedigree(std::map<PedPair, PedPair> &in_ped, std::map<PedPair, PedPair> &out_ped, std::vector<std::int64_t> &traced_id);
    template void Amat<double>::trace_pedigree(std::map<PedPair, PedPair> &in_ped, std::map<PedPair, PedPair> &out_ped, std::vector<std::int64_t> &traced_id);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::trace_operation(std::vector<std::int64_t> &days,
                                  std::vector<std::int64_t> &ids,
                                  std::vector<std::int64_t> &sire,
                                  std::vector<std::int64_t> &dame,
                                  std::vector<std::int64_t> &in_traced_id,
                                  std::vector<std::int64_t> &out_traced_id,
                                  std::map<PedPair, PedPair> &out_ped,
                                  std::vector<size_t> &loads_vect,
                                  size_t thr_id)
    {

        try
        {
            size_t i_first = 0; // very first element in a threads range

            if (thr_id != 0)
                i_first = loads_vect[thr_id - 1] + 1;

            size_t i_last = loads_vect[thr_id]; // very last element in a threads range

            std::vector<std::int64_t> gen_pedID; // working vector of IDs to be traced

            for (size_t i = i_first; i <= i_last; i++)
                gen_pedID.push_back(in_traced_id[i]); // get only those IDs assigned to the current thread

            std::int64_t elem_v;
            size_t iter_v = 0;
            bool exists;

            Utilities2 u;
            PedPair id_pair;
            PedPair key;

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

                        if ((birth_id_map[sire[i]] > days[i]) || (birth_id_map[dame[i]] > days[i]))
                        {
                            std::string s1("id: "+std::to_string(ids[i])+", birth: "+std::to_string(days[i])+"; ");
                            std::string s2("sire: "+std::to_string(sire[i])+", birth: "+std::to_string(birth_id_map[sire[i]])+"; ");
                            std::string s3("dame: "+std::to_string(dame[i])+", birth: "+std::to_string(birth_id_map[dame[i]])+".");
                            throw std::string("Pedigree is not correct, parents born before offspring! "+s1+s2+s3);
                        }

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

            if (u.is_unique(gen_pedID)) // Check if no repeated IDs appiar in the list of traced IDs
                out_traced_id = gen_pedID;
            else
                throw std::string("In the traced pedigree some IDs appiar more then one time!");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::trace_operation( ... ): ";
            std::cerr << e.what() << '\n';
            throw e;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::trace_operation( ... ): ";
            std::cerr << e << '\n';
            throw e;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::trace_operation( ... )" << '\n';
            throw;
        }
    }

    template void Amat<double>::trace_operation(std::vector<std::int64_t> &days,
                                                std::vector<std::int64_t> &ids,
                                                std::vector<std::int64_t> &sire,
                                                std::vector<std::int64_t> &dame,
                                                std::vector<std::int64_t> &in_traced_id,
                                                std::vector<std::int64_t> &out_traced_id,
                                                std::map<PedPair, PedPair> &out_ped,
                                                std::vector<size_t> &loads_vect,
                                                size_t thr_id);
    template void Amat<float>::trace_operation(std::vector<std::int64_t> &days,
                                               std::vector<std::int64_t> &ids,
                                               std::vector<std::int64_t> &sire,
                                               std::vector<std::int64_t> &dame,
                                               std::vector<std::int64_t> &in_traced_id,
                                               std::vector<std::int64_t> &out_traced_id,
                                               std::map<PedPair, PedPair> &out_ped,
                                               std::vector<size_t> &loads_vect,
                                               size_t thr_id);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::thread_loads(std::vector<std::int64_t> &in, std::vector<size_t> &out)
    {
        try
        {
            const auto processor_count = std::thread::hardware_concurrency(); // may return 0 when not able to detect

            size_t n_threads = processor_count;
            size_t vect_size = in.size();
            size_t work_load = 0;
            work_load = (size_t)(vect_size / n_threads); // expected work load per thread
            size_t max_load = 1;                         // max load (elements in range) per thread (should be probably changed!)

            if (work_load <= max_load) // correct the number of assigned threads (is more the case for a small data or large max_load)
            {
                n_threads = (size_t)n_threads / 2.0;
                work_load = (size_t)(vect_size / n_threads);

                if (work_load <= max_load)
                {
                    n_threads = 1;
                    work_load = vect_size;
                }
            }

            size_t last_index = 0; // approx. load

            for (size_t i = 0; i < n_threads - 1; i++)
            {
                last_index += work_load;
                out.push_back(last_index);
            }
            out.push_back(vect_size - 1);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::thread_loads(std::vector<std::int64_t> &, std::vector<size_t> &)" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::thread_loads(std::vector<std::int64_t> &, std::vector<size_t> &)" << '\n';
        }
    }

    template void Amat<float>::thread_loads(std::vector<std::int64_t> &in, std::vector<size_t> &out);
    template void Amat<double>::thread_loads(std::vector<std::int64_t> &in, std::vector<size_t> &out);

    //===============================================================================================================

    template <typename T>
    std::int64_t Amat<T>::pos_inped(std::map<std::int64_t, std::int64_t> &codemap, std::int64_t id)
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
            std::cerr << "Exception in Amat<T>::pos_inped(std::map<std::int64_t, std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::pos_inped(std::map<std::int64_t, std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::pos_inped(std::map<std::int64_t, std::int64_t> &, std::int64_t)" << '\n';
            throw;
        }
    }

    template std::int64_t Amat<float>::pos_inped(std::map<std::int64_t, std::int64_t> &codemap, std::int64_t id);
    template std::int64_t Amat<double>::pos_inped(std::map<std::int64_t, std::int64_t> &codemap, std::int64_t id);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::get_dinv(std::map<PedPair, PedPair> &ped, std::vector<T> &dinv, bool inbreed)
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

                std::vector<T> F((n + 1), 0.0); // inbreeding coefficients
                std::vector<T> B((m + 1), 0.0); // within family segregation variances
                std::vector<T> x((m + 1), 0.0); // x arrays

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
                        T g = 1.0 / (0.5 - 0.25 * (F[s] + F[d]));
                        dinv.push_back(g);
                    }
                    else if (!s && !d)
                    {
                        T g = 1.0;
                        dinv.push_back(g);
                    }
                    else
                    {
                        T g;
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
            std::cerr << "Exception in Amat<T>::get_dinv(std::map<PedPair, PedPair> &, std::vector<T> &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::get_dinv(std::map<PedPair, PedPair> &, std::vector<T> &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::get_dinv(std::map<PedPair, PedPair> &, std::vector<T> &, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::get_dinv(std::map<PedPair, PedPair> &ped, std::vector<float> &dinv, bool inbreed);
    template void Amat<double>::get_dinv(std::map<PedPair, PedPair> &ped, std::vector<double> &dinv, bool inbreed);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::get_ainv(std::map<PedPair, PedPair> &ped, std::map<PedPair, T> &ai, bool inbreed)
    {
        // Calculates A(-1) matrix (as map representation)
        try
        {
            PedPair akey;
            std::vector<T> di; // store D(-1)

            get_dinv(ped, di, inbreed);

            std::int64_t s, d, id;
            T dinv;

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
            std::cerr << "Exception in Amat<T>::get_ainv(std::map<PedPair, PedPair> &, std::map<PedPair, T> &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::get_ainv(std::map<PedPair, PedPair> &, std::map<PedPair, T> &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::get_ainv(std::map<PedPair, PedPair> &, std::map<PedPair, T> &, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::get_ainv(std::map<PedPair, PedPair> &ped, std::map<PedPair, float> &ai, bool inbreed);
    template void Amat<double>::get_ainv(std::map<PedPair, PedPair> &ped, std::map<PedPair, double> &ai, bool inbreed);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::get_a(std::map<PedPair, PedPair> &in_ped, std::map<PedPair, T> &out_a)
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

                    if (s2 && d2)
                    {
                        if (id1 == id2)
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
                            out_a[akey] = 0.5 * (out_a[bkey] + out_a[ckey]);
                            out_a[akey2] = out_a[akey];
                        }
                    }
                    if (!s2 && !d2)
                    {
                        if (id1 == id2)
                        {
                            out_a[akey] = 1.0;
                        }
                        else
                        {
                            out_a[akey] = 0.0;
                            out_a[akey2] = out_a[akey];
                        }
                    }
                    if (s2 && !d2)
                    {
                        if (id1 == id2)
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
                    if (!s2 && d2)
                    {
                        if (id1 == id2)
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
            std::cerr << "Exception in Amat<T>::get_a(std::map<PedPair, PedPair> &, std::map<PedPair, T> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::get_a(std::map<PedPair, PedPair> &, std::map<PedPair, T> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::get_a(std::map<PedPair, PedPair> &, std::map<PedPair, T> &)" << '\n';
            throw;
        }
    }

    template void Amat<float>::get_a(std::map<PedPair, PedPair> &in_ped, std::map<PedPair, float> &out_a);
    template void Amat<double>::get_a(std::map<PedPair, PedPair> &in_ped, std::map<PedPair, double> &out_a);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::get_inbreeding(std::vector<T> &out)
    {
        try
        {
            out = inbrF;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::get_inbreeding(std::vector<T> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::get_inbreeding(std::vector<T> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::get_inbreeding(std::vector<T> &)" << '\n';
            throw;
        }
    }

    template void Amat<float>::get_inbreeding(std::vector<float> &out);
    template void Amat<double>::get_inbreeding(std::vector<double> &out);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::get_inbreeding(const std::string &fname)
    {
        try
        {
            if ( traced_pedID.size() != inbrF.size() )
                throw std::string("traced_pedID.size() != inbrF.size()");
            
            std::ofstream result(fname);

            if (result.is_open())
            {
                for (size_t i = 0; i < inbrF.size(); i++)
                    result << traced_pedID[i] << ", " << inbrF[i] << "\n";

                result.close();
            }
            else
                throw std::string("Unable to open file for output!");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::get_inbreeding(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::get_inbreeding(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::get_inbreeding(const std::string &)" << '\n';
            throw;
        }
    }

    template void Amat<float>::get_inbreeding(const std::string &fname);
    template void Amat<double>::get_inbreeding(const std::string &fname);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::clear()
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
            
            A_s.clean();
            iA_s.clean();
            irA_s.clean();
            iA22_s.clean();
            A22_s.clean();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::clear()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::clear()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::clear()" << '\n';
            throw;
        }
    }

    template void Amat<float>::clear();
    template void Amat<double>::clear();

    //===============================================================================================================

    template <typename T>
    void Amat<T>::get_matrix(const std::string &name, evolm::matrix<T> &arr, std::vector<std::int64_t> &out)
    {
        /*
            Note, we operate with L-stored format, hence return lower triangular part.
            Note: in this current implementation copied-to-outside matrices remain linked with a class instance
            because the same file names which become common for internal containers and copied matrices. When
            class instance calls clear() or destructor - the binary files are cleaned from disk.
        */

        try
        {
            Utilities2 u;

            if (name == "A") // we use iA as a container for A as well
            {
                if ( IsEmpty.iA ) // was used sparse pipeline
                {
                    iA.resize( id_iA.size() );
                    iA_s.fread();                    
                    u.sparse_to_dense( iA_s, iA );
                    iA_s.fwrite();
                    iA.fwrite();
                }
               
                arr = iA;
                out = id_iA;
            }

            if (name == "rA") // we use irA as a container for rA as well
            {
                if ( IsEmpty.irA ) // was used sparse pipeline
                {
                    irA.resize( id_irA.size() );
                    irA_s.fread();                    
                    u.sparse_to_dense( irA_s, irA );
                    irA_s.fwrite();
                    irA.fwrite();
                }

                arr = irA;
                out = id_irA;
            }

            if (name == "iA")
            {
                if ( IsEmpty.iA ) // was used sparse pipeline
                {
                    iA.resize( id_iA.size() );
                    iA_s.fread();                    
                    u.sparse_to_dense( iA_s, iA );
                    iA_s.fwrite();
                    iA.fwrite();
                }

                arr = iA;
                out = id_iA;
            }

            if (name == "irA")
            {
                if ( IsEmpty.irA ) // was used sparse pipeline
                {
                    irA.resize( id_irA.size() );
                    irA_s.fread();                    
                    u.sparse_to_dense( irA_s, irA );
                    irA_s.fwrite();
                    irA.fwrite();
                }

                arr = irA;
                out = id_irA;
            }

            if (name == "iA22")
            {
                if ( IsEmpty.iA22 ) // was used sparse pipeline
                {
                    iA22.resize( id_A22.size() );
                    iA22_s.fread();                    
                    u.sparse_to_dense( iA22_s, iA22 );
                    iA22_s.fwrite();
                    iA22.fwrite();
                }

                arr = iA22;
                out = id_A22;
            }

            if (name == "A22")
            {
                if ( IsEmpty.A22 ) // was used sparse pipeline
                {
                    A22.resize( id_A22.size() );
                    A22_s.fread();                    
                    u.sparse_to_dense( A22_s, A22 );
                    A22_s.fwrite();
                    A22.fwrite();
                }

                arr = A22;
                out = id_A22;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::get_matrix(const std::string &, evolm::matrix<T> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::get_matrix(const std::string &, evolm::matrix<T> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::get_matrix(const std::string &, evolm::matrix<T> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    template void Amat<float>::get_matrix(const std::string &name, evolm::matrix<float> &arr, std::vector<std::int64_t> &out);
    template void Amat<double>::get_matrix(const std::string &name, evolm::matrix<double> &arr, std::vector<std::int64_t> &out);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::get_matrix(const std::string &name, evolm::smatrix<T> &arr, std::vector<std::int64_t> &out)
    {
        /*
            Note, we operate with L-stored format, hence return lower triangular part.
            Here we copy not by value but by binary file: arr matrix inherits all propertie of copying matrix,
            if it is on disc, only file name is copied.

            Note: in this current implementation copied-to-outside matrices remain linked with a class instance
            because the same file names which become common for internal containers and copied matrices. When
            class instance calls clear() or destructor - the binary files are cleaned from disk. 
        */
        try
        {
            Utilities2 u;

            if (name == "A") // we use iA as a container for A as well
            {
                if ( IsEmpty.iA_s ) // was used dense pipeline
                {
                    iA_s.resize( id_iA.size() );
                    iA.fread();                    
                    u.dense_to_sparse( iA, iA_s );
                    iA_s.fwrite();
                    iA.fwrite();
                }

                arr = iA_s;
                out = id_iA;
            }

            if (name == "rA") // we use irA as a container for rA as well
            {
                if ( IsEmpty.irA_s ) // was used dense pipeline
                {
                    irA_s.resize( id_irA.size() );
                    irA.fread();                    
                    u.dense_to_sparse( irA, irA_s );
                    irA_s.fwrite();
                    irA.fwrite();
                }

                arr = irA_s;
                out = id_irA;
            }

            if (name == "iA")
            {
                if ( IsEmpty.iA_s ) // was used dense pipeline
                {
                    iA_s.resize( id_iA.size() );
                    iA.fread();                    
                    u.dense_to_sparse( iA, iA_s );
                    iA_s.fwrite();
                    iA.fwrite();
                }

                arr = iA_s;
                out = id_iA;
            }

            if (name == "irA")
            {
                if ( IsEmpty.irA_s ) // was used dense pipeline
                {
                    irA_s.resize( id_irA.size() );
                    irA.fread();                    
                    u.dense_to_sparse( irA, irA_s );
                    irA_s.fwrite();
                    irA.fwrite();
                }

                arr = irA_s;
                out = id_irA;
            }

            if (name == "iA22")
            {
                if ( IsEmpty.iA22_s ) // was used dense pipeline
                {
                    iA22_s.resize( id_A22.size() );
                    iA22.fread();                    
                    u.dense_to_sparse( iA22, iA22_s );
                    iA22_s.fwrite();
                    iA22.fwrite();
                }

                arr = iA22_s;
                out = id_A22;
            }

            if (name == "A22") // was used dense pipeline, this is assumed to be always dense
            {
                if ( IsEmpty.A22_s )
                {
                    A22_s.resize( id_A22.size() );
                    A22.fread();
                    u.dense_to_sparse( A22, A22_s );
                    A22_s.fwrite();
                    A22.fwrite();
                }
 
                arr = A22_s;
                out = id_A22;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::get_matrix(const std::string &, evolm::smatrix<T> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::get_matrix(const std::string &, evolm::smatrix<T> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::get_matrix(const std::string &, evolm::smatrix<T> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    template void Amat<float>::get_matrix(const std::string &name, evolm::smatrix<float> &arr, std::vector<std::int64_t> &out);
    template void Amat<double>::get_matrix(const std::string &name, evolm::smatrix<double> &arr, std::vector<std::int64_t> &out);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::save_matrix(const std::string &name, const std::string &out_fname)
    {
        /* 
            This method is for Python interfacing;
            Note, we operate with L-stored format, hence return lower triangular part.
            
            Note: in this current implementation copied-to-outside matrices remain linked with a class instance
            because the same file names which become common for internal containers and copied matrices. When
            class instance calls clear() or destructor - the binary files are cleaned from disk. 
        */
        try
        {
            Utilities2 u;

            std::vector<T> values;
            std::vector<size_t> keys;

            if (name == "A") // we use iA as a container for A as well
            {
                if ( IsEmpty.iA_s ) // was used dense pipeline
                {
                    // convert iA to vect, use id_iA for ids
                    iA.fread();
                    iA.to_vector(values);
                    for (size_t i = 0; i < id_iA.size(); i++)
                        for (size_t j = 0; j <= i; j++)
                            keys.push_back(i*(i+1)/2 + j);
                    iA.fwrite();
                }
                else
                {
                    // convert iA_s to vect, use id_iA for ids
                    iA_s.fread();
                    iA_s.to_vect(values, keys);
                    iA_s.clear();
                }
                u.fwrite_matrix(out_fname, values, keys, id_iA);
            }

            if (name == "rA") // we use irA as a container for rA as well
            {
                if ( IsEmpty.irA_s ) // was used dense pipeline
                {
                    // convert irA to vect, use id_irA for ids
                    irA.fread();
                    irA.to_vector(values);
                    for (size_t i = 0; i < id_irA.size(); i++)
                        for (size_t j = 0; j <= i; j++)
                            keys.push_back(i*(i+1)/2 + j);
                    irA.fwrite();
                }
                else
                {
                    // convert irA_s to vect, use id_irA for ids
                    irA_s.fread();
                    irA_s.to_vect(values, keys);
                    irA_s.clear();
                }
                u.fwrite_matrix(out_fname, values, keys, id_irA);
            }

            if (name == "iA")
            {
                if ( IsEmpty.iA_s ) // was used dense pipeline
                {
                    // convert iA to vect, use id_iA for ids
                    iA.fread();
                    iA.to_vector(values);
                    for (size_t i = 0; i < id_iA.size(); i++)
                        for (size_t j = 0; j <= i; j++)
                            keys.push_back(i*(i+1)/2 + j);
                    iA.fwrite();
                }
                else
                {
                    // convert iA_s to vect, use id_iA for ids
                    iA_s.fread();
                    iA_s.to_vect(values, keys);
                    iA_s.clear();
                }
                u.fwrite_matrix(out_fname, values, keys, id_iA);
            }

            if (name == "irA")
            {
                if ( IsEmpty.irA_s ) // was used dense pipeline
                {
                    // convert irA to vect, use id_irA for ids
                    irA.fread();
                    irA.to_vector(values);
                    for (size_t i = 0; i < id_irA.size(); i++)
                        for (size_t j = 0; j <= i; j++)
                            keys.push_back(i*(i+1)/2 + j);
                    irA.fwrite();
                }
                else
                {
                    // convert irA_s to vect, use id_irA for ids
                    irA_s.fread();
                    irA_s.to_vect(values, keys);
                    irA_s.clear();
                }
                u.fwrite_matrix(out_fname, values, keys, id_irA);
            }

            if (name == "iA22")
            {
                if ( IsEmpty.iA22_s ) // was used dense pipeline
                {
                    // convert iA22 to vect, use id_A22 for ids
                    iA22.fread();
                    iA22.to_vector(values);
                    for (size_t i = 0; i < id_A22.size(); i++)
                        for (size_t j = 0; j <= i; j++)
                            keys.push_back(i*(i+1)/2 + j);
                    iA22.fwrite();
                }
                else
                {
                    // convert iA22_s to vect, use id_A22 for ids
                    iA22_s.fread();
                    iA22_s.to_vect(values, keys);
                    iA22_s.clear();
                }
                u.fwrite_matrix(out_fname, values, keys, id_A22);
            }

            if (name == "A22") // was used dense pipeline, this is assumed to be always dense
            {
                if ( IsEmpty.A22_s )
                {
                    // convert A22 to vect, use id_A22 for ids
                    A22.fread();
                    A22.to_vector(values);
                    for (size_t i = 0; i < id_A22.size(); i++)
                        for (size_t j = 0; j <= i; j++)
                            keys.push_back(i*(i+1)/2 + j);
                    A22.fwrite();
                }
                else
                {
                    // convert A22_s to vect, use id_A22 for ids
                    A22_s.fread();
                    A22_s.to_vect(values, keys);
                    A22_s.clear();
                }
                u.fwrite_matrix(out_fname, values, keys, id_A22);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::save_matrix(const std::string &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::save_matrix(const std::string &, const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::save_matrix(const std::string &, const std::string &)" << '\n';
            throw;
        }
    }

    template void Amat<float>::save_matrix(const std::string &name, const std::string &out_fname);
    template void Amat<double>::save_matrix(const std::string &name, const std::string &out_fname);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::save_ids(const std::string &name, const std::string &out_fname)
    {
        try
        {
            std::ofstream out;
            out.open (out_fname, std::ofstream::out | std::ofstream::trunc);
            
            if (!out.is_open()) throw std::string("Cannot open file "+ out_fname + "for writing!");

            out << "ref_ids"<<'\n';
            
            if ( name == "A" || name == "iA" )
            {
                if (id_iA.empty()) throw std::string("The requested vector of ids is empty!");
                for (auto const &v: id_iA) out << v <<'\n';
            }

            if ( name == "rA" || name == "irA" )
            {
                if (id_irA.empty()) throw std::string("The requested vector of ids is empty!");
                for (auto const &v: id_irA) out << v <<'\n';
            }

            if ( name == "iA22" || name == "A22" )
            {
                if (id_A22.empty()) throw std::string("The requested vector of ids is empty!");
                for (auto const &v: id_A22) out << v <<'\n';
            }

            out.close();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::save_ids(const std::string &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::save_ids(const std::string &, const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::save_ids(const std::string &, const std::string &)" << '\n';
            throw;
        }
    }

    template void Amat<float>::save_ids(const std::string &name, const std::string &out_fname);
    template void Amat<double>::save_ids(const std::string &name, const std::string &out_fname);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::get_A22(std::map<PedPair, PedPair> &ped, std::vector<std::int64_t> &genotypedID)
    {
        try
        {
            Utilities2 u;

            size_t n = ped.size(); // total amount of animals in pedigree
            std::vector<std::vector<std::int64_t>> Ped(n + 1, std::vector<std::int64_t>(2, 0.0));
            std::vector<std::int64_t> GenID; // list of genotyped IDs

            std::map<std::int64_t, std::int64_t> code_map;
            std::map<std::int64_t, std::int64_t> gen_map;

            std::int64_t code_id = 1;
            for (auto const &elem : ped)
            {

                code_map[elem.first.val_2] = code_id;

                Ped[code_id][0] = pos_inped(code_map, elem.second.val_1);
                Ped[code_id][1] = pos_inped(code_map, elem.second.val_2);

                int pos = u.find_invect(genotypedID, elem.first.val_2);
                if (pos != -1)
                {
                    GenID.push_back(code_id);
                    gen_map[code_id] = pos;
                }

                code_id++;
            }

            size_t m = GenID.size(); // number of genotyped IDs

            if (!A.empty())
                A.clear();

            A.resize(m);

#pragma omp parallel for
            for (size_t i = 0; i < m; i++)
            {
                std::vector<T> w(n + 1, 0.0);
                std::vector<std::int64_t> v(n + 1, 0);

                size_t ii = gen_map[GenID[i]]; // gives position of recoded (1...n) genotyped IDs

                v[GenID[i]] = 1;

                getA22vector(w, v, Ped);

                for (size_t j = 0; j < m; j++)
                {
                    if (ii >= (size_t)gen_map[GenID[j]])
                        A[ii * (ii + 1) / 2 + gen_map[GenID[j]]] = w[GenID[j]];
                }
            }

        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::get_A22(std::map <PedPai, PedPair> &, std::vector<std::int64_t> &, evolm::matrix<T> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::get_A22(std::map <PedPai, PedPai> &, std::vector<std::int64_t> &, evolm::matrix<T> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::get_A22(std::map <PedPai, PedPai> &, std::vector<std::int64_t> &, evolm::matrix<T> &)" << '\n';
            throw;
        }
    }

    template void Amat<float>::get_A22(std::map<PedPair, PedPair> &ped, std::vector<std::int64_t> &genotypedID);
    template void Amat<double>::get_A22(std::map<PedPair, PedPair> &ped, std::vector<std::int64_t> &genotypedID);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::getA22vector(std::vector<T> &w, std::vector<std::int64_t> &v, std::vector<std::vector<std::int64_t>> &Ped)
    {
        try
        {
            size_t n = w.size() - 1;
            std::vector<T> q(n + 1, 0.0);

            for (size_t i = n; i >= 1; i--)
            {
                q[i] += v[i];
                auto s = Ped[i][0];
                auto d = Ped[i][1];
                if (s)
                    q[s] += q[i] * 0.5;
                if (d)
                    q[d] += q[i] * 0.5;
            }

            for (size_t i = 1; i <= n; i++)
            {
                auto s = Ped[i][0];
                auto d = Ped[i][1];
                auto di = (std::count(Ped[i].begin(), Ped[i].end(), 0) + 2.0) / 4.0 - 0.25 * (inbrF[s] + inbrF[d]);
                T temp = (T)0;
                if (s)
                    temp += w[s];
                if (d)
                    temp += w[d];
                w[i] = 0.5 * temp;
                w[i] += di * q[i];
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::getA22vector(std::vector <T> &, std::vector <std::int64_t> &, std::vector<std::vector<std::int64_t> > &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::getA22vector(std::vector <T> &, std::vector <std::int64_t> &, std::vector<std::vector<std::int64_t> > &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::getA22vector(std::vector <T> &, std::vector <std::int64_t> &, std::vector<std::vector<std::int64_t> > &)" << '\n';
            throw;
        }
    }

    template void Amat<float>::getA22vector(std::vector<float> &w, std::vector<std::int64_t> &v, std::vector<std::vector<std::int64_t>> &Ped);
    template void Amat<double>::getA22vector(std::vector<double> &w, std::vector<std::int64_t> &v, std::vector<std::vector<std::int64_t>> &Ped);

    //===============================================================================================================

} // end of namespace evoped