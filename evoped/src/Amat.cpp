#include "Amat.hpp"

namespace evoped
{
    //===============================================================================================================

    template <typename T>
    Amat<T>::Amat()
    {
        try
        {
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
    void Amat<T>::map_to_matr(std::map<PedPair, T> &amap, std::vector<std::int64_t> &ids, bool use_ainv, bool use_large)
    {
        try
        {
            Utilities2 u;

            std::map<std::int64_t, std::int64_t> rid_map;

            if (ids.empty())
                throw std::string("Empty traced pedigree IDs!");

            size_t limit = 0.5 * (ids.size() - 1) * ids.size() + ids.size();
            if (!use_ainv)
                limit = ids.size() * ids.size();

            if (limit < amap.size())
                throw std::string("The number of elements in calculated A(-1) matrix is higher than the number of traced IDs!!");

            u.get_RecodedIdMap(rid_map, ids); // list of real ids => std::map for the new consecutive list starting from 1

            if (rid_map.empty())
                throw std::string("Recoded IDs map is empty!");

            if (!A.empty())
                A.clear();

            if (!A_s.empty())
                A_s.clear();

            if (use_large)
                A_s.resize(ids.size());
            else
                A.resize(ids.size());

            if (use_large)
            {
                for (auto const &elem : amap)
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
                    A_s[ind] = g_val;
                }
            }
            else
            {
                for (auto const &elem : amap)
                {
                    std::int64_t g_row = elem.first.val_1; // id
                    std::int64_t g_col = elem.first.val_2; // id
                    T g_val = elem.second;

                    size_t r = rid_map[g_row];            // position in the list if ids, consecutive index of real id in the list of all ids
                    size_t c = rid_map[g_col];            // position in the list if ids, consecutive index of real id in the list of all ids
                    size_t ind = r * (r - 1) / 2 + c - 1; // r & c start from 1, but not from 0
                    if (c > r)
                        ind = c * (c - 1) / 2 + r - 1;
                    A[ind] = g_val;
                }
            }

            rid_map.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::map_to_matr(std::map<PedPair, T> &, std::vector<std::int64_t> &, bool, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::map_to_matr(std::map<PedPair, T> &, std::vector<std::int64_t> &, bool, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::map_to_matr(std::map<PedPair, T> &, std::vector<std::int64_t> &, bool, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::map_to_matr(std::map<PedPair, float> &amap, std::vector<std::int64_t> &ids, bool use_ainv, bool use_large);
    template void Amat<double>::map_to_matr(std::map<PedPair, double> &amap, std::vector<std::int64_t> &ids, bool use_ainv, bool use_large);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::make_matrix(const std::string &ped_file, bool use_ainv, bool use_large)
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

            map_to_matr(ainv, traced_pedID, use_ainv, use_large); // converting ainv map to A or A(-1) matrix

            ainv.clear();

            id_iA = traced_pedID;

            if (use_large)
            {
                A_s.fwrite();
                iA_s = A_s;
                A_s.resize();
            }
            else
            {
                A.fwrite();
                iA = A;
                A.clear();
            }

            traced_pedID.clear();
            traced_pedID.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::make_matrix(std::string &, bool, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::make_matrix(std::string &, bool, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::make_matrix(std::string &, bool, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::make_matrix(const std::string &ped_file, bool use_ainv, bool use_large);
    template void Amat<double>::make_matrix(const std::string &ped_file, bool use_ainv, bool use_large);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::make_matrix(const std::string &ped_file, const std::string &g_file, bool use_ainv, bool use_large)
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

            map_to_matr(r_ainv, traced_pedID, use_ainv, use_large); // converting r_ainv map to A(-1) matrix

            r_ainv.clear();

            id_irA = traced_pedID;

            if (use_large)
            {
                A_s.fwrite();
                irA_s = A_s;
                A_s.resize();
            }
            else
            {
                A.fwrite();
                irA = A;
                A.clear();
            }

            traced_pedID.clear();
            traced_pedID.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::make_matrix(std::string &, std::string &, bool, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::make_matrix(std::string &, std::string &, bool, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::make_matrix(std::string &, std::string &, bool, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::make_matrix(const std::string &ped_file, const std::string &g_file, bool use_ainv, bool use_large);
    template void Amat<double>::make_matrix(const std::string &ped_file, const std::string &g_file, bool use_ainv, bool use_large);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::make_all(const std::string &ped_file, const std::string &g_file, bool use_large)
    {
        // Making all matrices required for ssBlup: A(-1), red_A(-1), A22, A22(-1)
        try
        {
            // -----------------------------------------------------
            // ---------- full A(-1) -------------------------------
            std::cout << "full A(-1) ..."
                      << "\n";
            std::vector<std::int64_t> pedID;
            std::map<PedPair, PedPair> pedigree_from_file;
            std::map<PedPair, PedPair> pedigree;
            std::cout << "reading pedigree ..."
                      << "\n";
            auto start = std::chrono::high_resolution_clock::now();

            fread_pedigree(ped_file, pedigree_from_file, pedID);

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout << "fread_pedigree() duration (milliseconds): " << duration.count() << std::endl;

            if (pedID.empty())
                throw std::string("Empty pedigree IDs!");

            if (pedigree_from_file.empty())
                throw std::string("File provided pedigree is empty!");

            if (!traced_pedID.empty())
                traced_pedID.erase(traced_pedID.begin(), traced_pedID.end());
            std::cout << "tracing pedigree ..."
                      << "\n";

            stop = std::chrono::high_resolution_clock::now();

            trace_pedigree(pedigree_from_file, pedigree, pedID); // tracing full pedigree

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout << "trace_pedigree() duration (milliseconds): " << duration.count() << std::endl;

            if (!inbrF.empty())
                inbrF.erase(inbrF.begin(), inbrF.end());

            std::map<PedPair, T> ainv;
            std::cout << "making A(-1) ..."
                      << "\n";

            stop = std::chrono::high_resolution_clock::now();

            get_ainv(pedigree, ainv, true); // making A(-1)
            
            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout << "get_ainv() duration (milliseconds): " << duration.count() << std::endl;

            pedigree.clear();
            pedID.clear();
            pedID.shrink_to_fit();

            std::cout << "map_to_matr ..."
                      << "\n";

            map_to_matr(ainv, traced_pedID, true, use_large); // converting ainv map to A or A(-1) matrix
            
            std::cout << "Done."
                      << "\n";
            ainv.clear();

            id_iA = traced_pedID;

            if (use_large)
            {
                A_s.fwrite(); // 1. Wrire to binary
                iA_s = A_s;   // 2. Copy matrix by exchanging the internal binary file name
                A_s.resize();
            }
            else
            {
                A.fwrite(); // 1. Wrire to binary
                iA = A;     // 2. Copy matrix by exchanging the internal binary file name
                A.clear();  // 3. Clears the memory and gets new name for internal binary file
            }
            std::cout << "Done."
                      << "\n";
            // -----------------------------------------------------
            // ---------- reduced A(-1) ----------------------------
            std::cout << "reduced A(-1)"
                      << "\n";
            std::vector<std::int64_t> genotypedID;
            std::map<PedPair, PedPair> r_pedigree;

            fread_genotyped_id(g_file, genotypedID);

            if (genotypedID.empty())
                throw std::string("Cannot trace the reduced pedigree: ID's vector is empty!");

            if (!traced_pedID.empty())
                traced_pedID.erase(traced_pedID.begin(), traced_pedID.end());

            trace_pedigree(pedigree_from_file, r_pedigree, genotypedID); // tracing reduced pedigree for genotyped individuals

            if (!inbrF.empty())
                inbrF.erase(inbrF.begin(), inbrF.end());

            std::map<PedPair, T> r_ainv;

            get_ainv(r_pedigree, r_ainv, true); // making reduced A(-1)

            pedigree_from_file.clear();

            map_to_matr(r_ainv, traced_pedID, true, use_large); // converting r_ainv map to A(-1) matrix

            r_ainv.clear();

            id_irA = traced_pedID;

            if (use_large)
            {
                A_s.fwrite(); // 1. Wrire to binary
                irA_s = A_s;  // 2. Copy matrix by just exchanging the internal binary file name
                A_s.resize();
            }
            else
            {
                A.fwrite(); // 1. Wrire to binary
                irA = A;    // 2. Copy matrix by just exchanging the internal binary file name
                A.clear();  // 3. Clears the memory and gets new name for internal binary file
            }
            std::cout << "Done."
                      << "\n";
            // -----------------------------------------------------
            // These two matrices are considered always dense !
            // -------------------- A22 ----------------------------
            std::cout << "A22"
                      << "\n";
            
            stop = std::chrono::high_resolution_clock::now();

            get_A22(r_pedigree, genotypedID);
            
            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout << "get_A22() duration (milliseconds): " << duration.count() << std::endl;

            id_A22 = genotypedID;

            r_pedigree.clear();
            birth_id_map.clear();

            size_t non_serro = 0;
            for (size_t i = 0; i < A.size(); i++)
            {
                if (A[i] != 0.0)
                    non_serro++;
            }
            std::cout << "Sparsity of A22: " << (1.0 - (double)non_serro / A.size()) * 100.0 << "\n";

            // Here we are using always dense matrix
            A.fwrite(); // 1. Wrire to binary
            A22 = A;    // 2. Copy matrix by just exchanging the internal binary file name
            A.clear();  // 3. Clears the memory and gets new name for internal binary file
            std::cout << "Done."
                      << "\n";
            // -----------------------------------------------------
            // -------------------- A22(-1) ------------------------
            std::cout << "A22(-1)"
                      << "\n";
            if (use_large)
            {
                irA_s.fread();

                stop = std::chrono::high_resolution_clock::now();

                get_iA22(irA_s, id_irA, genotypedID);

                stop = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
                std::cout << "get_A22() duration (milliseconds): " << duration.count() << std::endl;

                irA_s.fwrite();

                A_s.fwrite(); // 1. Wrire to binary
                iA22_s = A_s; // 2. Copy matrix by just exchanging the internal binary file name
                A_s.resize(); // 3. Clears the memory and gets new name for internal binary file

                // A.fwrite(); // 1. Wrire to binary
                // iA22 = A;   // 2. Copy matrix by just exchanging the internal binary file name
                // A.clear();  // 3. Clears the memory and gets new name for internal binary file
            }
            else
            {
                irA.fread();
                get_iA22(irA, id_irA, genotypedID);
                irA.fwrite();

                A.fwrite(); // 1. Wrire to binary
                iA22 = A;   // 2. Copy matrix by just exchanging the internal binary file name
                A.clear();  // 3. Clears the memory and gets new name for internal binary file
            }

            genotypedID.clear();
            genotypedID.shrink_to_fit();
            std::cout << "Done."
                      << "\n";
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::make_all(std::string &, std::string &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::make_all(std::string &, std::string &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::make_all(std::string &, std::string &, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::make_all(const std::string &ped_file, const std::string &g_file, bool use_large);
    template void Amat<double>::make_all(const std::string &ped_file, const std::string &g_file, bool use_large);

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

            if (!u.is_value_in_vect(matr_ids, selected_ids))
                throw std::string("There are IDs in the selected IDs array which are not part of the passed matrix!");

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
            std::cout << "Entering for making iA22"
                      << "\n";
            Utilities2 u;

            size_t n_rows = full_matr.nrows();
            size_t n_cols = full_matr.ncols();

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
            if (!A_s.empty())
            {
                A_s.resize();
                // A.fclear();
                // A.clear();
            }

            // Create the list of IDs which are in matr_ids but not in selected_ids vectors

            std::vector<std::int64_t> not_selected_ids;

            for (size_t i = 0; i < matr_ids.size(); i++)
            {
                int res = u.find_invect(selected_ids, matr_ids[i]);
                if (res == -1)
                    not_selected_ids.push_back(matr_ids[i]);
            }

            if ((selected_ids.size() + not_selected_ids.size()) != matr_ids.size())
                throw std::string("The sum of IDs from from two vectors is not equal to number of IDs in the matrix!");

            // ------------------------------------------------

            // Next steps: A22(-1) = A22 - A21 * A11(-1) * A12;

            evolm::smatrix<T> a22; // in order to differrentiate from the 'global' A22
            evolm::smatrix<T> A11;
            evolm::smatrix<T> A21;
            evolm::smatrix<T> A12;
            evolm::matrix<T> A11d; // this is the temporal storage for making inverse

            // -------------------- A11 -----------------------
            std::cout << "   Making A11 ..."
                      << "\n";
            A11d.resize(not_selected_ids.size());

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
                T zerro_value = (T)0;

                for (size_t j = 0; j <= i; j++)
                {
                    size_t pos_j = non_selected_pos[j];

                    if (pos_i >= pos_j)
                    {
                        value = full_matr.get_nonzero(pos_i, pos_j);

                        if (value != zerro_value)
                            A11d(i, j) = value;
                        else
                            A11d(i, j) = zerro_value;
                    }
                    else
                    {
                        value = full_matr.get_nonzero(pos_j, pos_i);

                        if (value != zerro_value)
                            A11d(i, j) = value;
                        else
                            A11d(i, j) = zerro_value;
                    }
                }
            }
            std::cout << "   Get sparsity of A11 ..."
                      << "\n";
            size_t non_zeros = 0;
            for (size_t i = 0; i < A11d.size(); i++)
            {
                if (A11d[i] != 0.0)
                    non_zeros++;
            }
            std::cout << "        A11d: not_selected_ids.size() is " << not_selected_ids.size() << "\n";
            std::cout << "        non zeros = " << non_zeros << ", expected all in rect: " << not_selected_ids.size() * not_selected_ids.size() << "\n";
            std::cout << "        sparsity of symmetric = " << (1.0 - (T)non_zeros / (T)A11d.size()) * 100.0 << "\n";
            std::cout << "        sparsity of rect = " << (1.0 - (T)non_zeros / (T)(not_selected_ids.size() * not_selected_ids.size())) * 100.0 << "\n";

            auto start = std::chrono::high_resolution_clock::now();

            A11d.symtorec();

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout << "A11d.symtorec() duration (milliseconds): " << duration.count()  << std::endl;

            std::cout << "   Get sparsity of rectangular A11 ..."
                      << "\n";
            non_zeros = 0;
            for (size_t i = 0; i < A11d.size(); i++)
            {
                if (A11d[i] != 0.0)
                    non_zeros++;
            }
            std::cout << "        non zeros = " << non_zeros << ", expected all in rect: " << A11d.size() << "\n";
            std::cout << "        sparsity = " << (1.0 - (T)non_zeros / (T)A11d.size()) * 100.0 << "\n";

            std::cout << "      inverting A11 ..."
                      << "\n";
            A11d.invert();

            std::cout << "   Get sparsity of inverted A11 ..."
                      << "\n";
            non_zeros = 0;
            for (size_t i = 0; i < A11d.size(); i++)
            {
                if (A11d[i] != 0.0)
                    non_zeros++;
            }
            std::cout << "        non zeros = " << non_zeros << ", expected all in rect: " << A11d.size() << "\n";
            std::cout << "        sparsity = " << (1.0 - (T)non_zeros / (T)A11d.size()) * 100.0 << "\n";

            std::cout << "A11d => A11 ..."
                      << "\n";

            start = std::chrono::high_resolution_clock::now();

            A11.resize(not_selected_ids.size(), not_selected_ids.size());
            for (size_t i = 0; i < not_selected_ids.size(); i++)
            {
                for (size_t j = 0; j < not_selected_ids.size(); j++)
                {
                    if (A11d(i, j) != 0.0)
                        A11(i, j) = A11d(i, j);
                }
            }

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout << "A11d => A11 duration (milliseconds): " << duration.count()  << std::endl;

            // A11d.print("A11d");
            // A11.print("A11");
            A11d.fclear();
            A11d.clear();

            A11.fwrite();
            std::cout << "   Done."
                      << "\n";
            // ------------------------------------------------
            //
            // -------------------- A21 -----------------------
            std::cout << "   Making A21 ..."
                      << "\n";

            start = std::chrono::high_resolution_clock::now();

            A21.resize(selected_ids.size(), not_selected_ids.size());

            // #pragma omp parallel for
            for (size_t i = 0; i < selected_ids.size(); i++)
            {
                size_t pos_i = selected_pos[i];
                T value = (T)0;
                T zerro_value = (T)0;

                for (size_t j = 0; j < not_selected_ids.size(); j++)
                {
                    size_t pos_j = non_selected_pos[j];

                    if (pos_i >= pos_j)
                    {
                        value = full_matr.get_nonzero(pos_i, pos_j);

                        if (value != zerro_value)
                            A21(i, j) = value;
                        // else
                        // A21(i,j) = 0.0;
                    }
                    else
                    {
                        value = full_matr.get_nonzero(pos_j, pos_i);

                        if (value != zerro_value)
                            A21(i, j) = value;
                        // else
                        // A21(i,j) = 0.0;
                    }
                }
            }

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout << "assigning A21 duration (milliseconds): " << duration.count()  << std::endl;

            // A21.print("A21");
            std::cout << "   Done."
                      << "\n";

            std::cout << "        sparsity of A21 = " << (1.0 - (T)A21.size() / (T)A21.max_key()) * 100.0 << "\n";
            // ------------------------------------------------
            //
            // -------------------- A12 -----------------------
            std::cout << "   Making A12 ..."
                      << "\n";
            std::cout << "   A12 = A21"
                      << "\n";
            A12 = A21;
            std::cout << "   A12.transpose()"
                      << "\n";

            start = std::chrono::high_resolution_clock::now();

            A12.transpose();
            // A12.print("A12");

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout << "A12.transpose() duration (milliseconds): " << duration.count()  << std::endl;

            std::cout << "   Done."
                      << "\n";

            std::cout << "   Writing."
                      << "\n";

            A12.fwrite();
            std::cout << "   Done."
                      << "\n";
            // ------------------------------------------------
            //
            // --------------- A21 * A11(-1) * A12 ------------
            start = std::chrono::high_resolution_clock::now();

            std::cout << "   Some reading and multiplication: reading ..."
                      << "\n";
            A11.fread();
            evolm::smatrix<T> res;
            std::cout << "      A21 * A11 ... "
                      << "elements in A21: " << A21.size() << ", A11: " << A11.size() << "\n";
            res = A21 * A11;
            // res.print("A21 * A11");
            std::cout << "fclearing ..."
                      << "\n";
            A11.fclear();
            A11.clear();
            // A11.resize();

            A21.fclear();
            A21.clear();
            std::cout << "resizing of A21."
                      << "\n";
            // A21.resize();

            A12.fread();
            std::cout << "      res * A12 ..."
                      << "\n";
            // res = res * A12; <= this does not work!!! Clean res afterwards!
            // evolm::smatrix<T> res2;
            res = res * A12;
            // res.print("res * A12");
            std::cout << "fclear of A12."
                      << "\n";
            A12.fclear();
            A12.clear();
            // A12.resize();

            std::cout << "rectosym()."
                      << "\n";
            res.rectosym();

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout << "A21 * A11(-1) * A12 duration (milliseconds): " << duration.count()  << std::endl;

            std::cout << "   Done."
                      << "\n";
            // ------------------------------------------------
            //
            // -------------------- A22 -----------------------
            std::cout << "   Making A22 ..."
                      << "\n";
            start = std::chrono::high_resolution_clock::now();

            a22.resize(selected_ids.size());

            // #pragma omp parallel for
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
                        // else
                        // a22(i,j) = 0.0;
                    }
                    else
                    {
                        value = full_matr.get_nonzero(pos_j, pos_i);

                        if (value != zerro_value)
                            a22(i, j) = value;
                        // else
                        // a22(i,j) = 0.0;
                    }
                }
            }

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout << "assigning A22 duration (milliseconds): " << duration.count()  << std::endl;

            // a22.print("a22");
            non_selected_pos.clear();
            non_selected_pos.shrink_to_fit();
            selected_pos.clear();
            selected_pos.shrink_to_fit();
            std::cout << "   Done."
                      << "\n";
            // ------------------------------------------------
            //
            // ------------------ A22 - res -------------------
            std::cout << "   Finalising ..."
                      << "\n";
            // evolm::matrix<size_t> shapeofa22;
            // shapeofa22 = a22.shape();

            std::cout << "        sparsity of a22 = " << (1.0 - (T)a22.size() / (T)a22.max_key()) * 100.0 << "\n";

            // #pragma omp parallel for
            // for (size_t i = 0; i < a22.size(); i++)
            // res[i] = a22[i] - res[i];
            // res.resize();
            std::cout << "   a22 - res ..."
                      << "\n";

            start = std::chrono::high_resolution_clock::now();

            res = a22 - res;

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
            std::cout << "a22 - res duration (milliseconds): " << duration.count()  << std::endl;

            std::cout << "   Done."
                      << "\n";
            // res.print("a22[i] - res[i] d");

            // a22.fclear();
            // a22.clear();
            a22.resize();

            // A = res;
            A_s = res;

            // res.fclear();
            // res.clear();
            res.resize();

            std::cout << "   Done."
                      << "\n";

            std::cout << "        sparsity of A22(-1) = " << (1.0 - (T)A_s.size() / (T)A_s.max_key()) * 100.0 << ", rows & cols: " << A_s.nrows() << " " << A_s.ncols() << "\n";

            // ------------------------------------------------
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
            std::vector<T> tmp_list;

            std::ifstream ped;
            ped.open(ped_file, std::fstream::in);

            if (!ped.good())
                throw std::string("Cannot open pedigree file!");

            while (getline(ped, line))
            {
                p = line.c_str();
                for (T f = std::strtod(p, &end); p != end; f = std::strtod(p, &end))
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
            std::vector<T> tmp_list;

            std::ifstream ped;
            ped.open(g_file, std::fstream::in);

            if (!ped.good())
                throw std::string("Cannot open genotyped ids file!");

            while (std::getline(ped, line))
            {
                p = line.c_str();

                for (T f = std::strtod(p, &end); p != end; f = std::strtod(p, &end))
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

            if (u.is_unique(gen_pedID)) // Check if no repeated IDs appiar in the list of traced IDs
                out_traced_id = gen_pedID;
            else
                throw std::string("In the traced pedigree some IDs appiar more then one time!");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in smatrix<T>::trace_operation( ... )" << '\n';
            std::cerr << e.what() << '\n';
        }
        catch (...)
        {
            std::cerr << "Exception in smatrix<T>::trace_operation( ... )" << '\n';
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
    void Amat<T>::get_matrix(const std::string &name, evolm::matrix<T> &arr, std::vector<std::int64_t> &out, bool keep_ondisk)
    {
        // Note, we operate with L-stored format, hence return lower triangular part
        try
        {
            if (name == "A")
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
            if (name == "rA")
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
            if (name == "iA")
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
            if (name == "irA")
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
            if (name == "iA22")
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
            if (name == "A22")
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
            std::cerr << "Exception in Amat<T>::get_matrix(const std::string &, evolm::matrix<T> &, std::vector<std::int64_t> &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::get_matrix(const std::string &, evolm::matrix<T> &, std::vector<std::int64_t> &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::get_matrix(const std::string &, evolm::matrix<T> &, std::vector<std::int64_t> &, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::get_matrix(const std::string &name, evolm::matrix<float> &arr, std::vector<std::int64_t> &out, bool keep_ondisk);
    template void Amat<double>::get_matrix(const std::string &name, evolm::matrix<double> &arr, std::vector<std::int64_t> &out, bool keep_ondisk);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::get_matrix(const std::string &name, evolm::smatrix<T> &arr, std::vector<std::int64_t> &out, bool keep_ondisk)
    {
        // Note, we operate with L-stored format, hence return lower triangular part
        try
        {
            if (name == "A")
            {
                // we use iA as a container for A as well
                iA_s.fread();
                arr = iA_s;
                out = id_iA;
                if (keep_ondisk)
                    iA_s.fwrite();
                else
                {
                    iA_s.fclear();
                    iA_s.clear();
                }
            }
            if (name == "rA")
            {
                // we use irA as a container for rA as well
                irA_s.fread();
                arr = irA_s;
                out = id_irA;
                if (keep_ondisk)
                    irA_s.fwrite();
                else
                {
                    irA_s.fclear();
                    irA_s.clear();
                }
            }
            if (name == "iA")
            {
                iA_s.fread();
                arr = iA_s;
                out = id_iA;
                if (keep_ondisk)
                    iA_s.fwrite();
                else
                {
                    iA_s.fclear();
                    iA_s.clear();
                }
            }
            if (name == "irA")
            {
                irA_s.fread();
                arr = irA_s;
                out = id_irA;
                if (keep_ondisk)
                    irA_s.fwrite();
                else
                {
                    irA_s.fclear();
                    irA_s.clear();
                }
            }
            if (name == "iA22")
            {
                iA22_s.fread();
                arr = iA22_s;
                out = id_A22;
                if (keep_ondisk)
                    iA22_s.fwrite();
                else
                {
                    iA22_s.fclear();
                    iA22_s.clear();
                }
            }
            if (name == "A22") // this is always dense
            {
                throw std::string("The A22 matrix is always dense. Use the same but overloaded method which allows dense matrices!");
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat<T>::get_matrix(const std::string &, evolm::smatrix<T> &, std::vector<std::int64_t> &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::get_matrix(const std::string &, evolm::smatrix<T> &, std::vector<std::int64_t> &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::get_matrix(const std::string &, evolm::smatrix<T> &, std::vector<std::int64_t> &, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::get_matrix(const std::string &name, evolm::smatrix<float> &arr, std::vector<std::int64_t> &out, bool keep_ondisk);
    template void Amat<double>::get_matrix(const std::string &name, evolm::smatrix<double> &arr, std::vector<std::int64_t> &out, bool keep_ondisk);

    //===============================================================================================================

    template <typename T>
    void Amat<T>::get_matrix(const std::string &name, std::vector<T> &arr, std::vector<std::int64_t> &out, bool keep_ondisk)
    {
        // This method is for Python interfacing;
        // Note, we operate with L-stored format, hence return lower triangular part
        try
        {
            if (name == "A")
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
            if (name == "rA")
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
            if (name == "iA")
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
            if (name == "irA")
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
            if (name == "iA22")
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
            if (name == "A22")
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
            std::cerr << "Exception in Amat<T>::get_matrix(const std::string &, std::vector<T> &, std::vector<std::int64_t> &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat<T>::get_matrix(const std::string &, std::vector<T> &, std::vector<std::int64_t> &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat<T>::get_matrix(const std::string &, std::vector<T> &, std::vector<std::int64_t> &, bool)" << '\n';
            throw;
        }
    }

    template void Amat<float>::get_matrix(const std::string &name, std::vector<float> &arr, std::vector<std::int64_t> &out, bool keep_ondisk);
    template void Amat<double>::get_matrix(const std::string &name, std::vector<double> &arr, std::vector<std::int64_t> &out, bool keep_ondisk);

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

auto start = std::chrono::high_resolution_clock::now();
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
auto stop = std::chrono::high_resolution_clock::now();
auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
std::cout << "in get_A22()->1 duration (milliseconds): " << duration.count() << std::endl;

            size_t m = GenID.size(); // number of genotyped IDs

            if (!A.empty())
                A.clear();

            A.resize(m);
start = std::chrono::high_resolution_clock::now();

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

stop = std::chrono::high_resolution_clock::now();
duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
std::cout << "in get_A22()->2 duration (milliseconds): " << duration.count() << std::endl;

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