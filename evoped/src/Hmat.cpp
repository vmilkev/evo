#include "Hmat.hpp"
#include "Gmat.hpp"
#include "Amat.hpp"

namespace evoped
{
    //===============================================================================================================
    template <typename T>
    Hmat<T>::Hmat()
    {
        //IsEmpty.H = true;
        //IsEmpty.H_s = true;
    }
    template Hmat<float>::Hmat();
    template Hmat<double>::Hmat();

    //===============================================================================================================
    template <typename T>
    Hmat<T>::~Hmat()
    {
        try
        {
            clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat<T>::~Hmat()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat<T>::~Hmat()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat<T>::~Hmat()" << '\n';
            throw;
        }
    }
    template Hmat<float>::~Hmat();
    template Hmat<double>::~Hmat();
    //===============================================================================================================
    template <typename T>
    void Hmat<T>::make_matrix(const std::string &g_file, const std::string &ped_file, const std::string &out_file)
    {
        try
        {
            // check extension of g_file to determine if it is ready-to-use matrix or SNP file
            std::string f_suffix(".gmatrix");
            size_t found_suffix = g_file.find(f_suffix);
            bool is_gmatrix = false;
            if (found_suffix != std::string::npos)
                is_gmatrix = true;
            
            Gmat<T> g;
            evolm::matrix<T> iG;
            std::vector<std::int64_t> g_id;
            
            Amat<T> a;
            evolm::matrix<T> A22;
            std::vector<int64_t> a22_id;
            evolm::smatrix<T> iA;
            std::vector<int64_t> iA_id;
            evolm::matrix<T> iA22;
            std::vector<int64_t> iA22_id;

            if ( is_gmatrix )          // read G matrix from ASCII file with IDs
                g.read_matrix(g_file);
            else                       // build G matrix from ASCII SNP file with IDs
                g.make_matrix(g_file);

            g.scale_diag(0.01);
            g.get_ids(g_id); // here we need just IDs

            a.make_matrix_forgenotyped( ped_file, g_id, false );
            a.get_matrix("A22", A22, a22_id);

            A22.fread();
            
            g.scale_matrix(A22,0.25);

            A22.clear();
            a.clear();

            g.invert_matrix(true); // inverting as a full-store            
            g.save_matrix2("iG.dmat");
            g.clear();

            a.make_matrix( ped_file, true ); // A(-1)
            a.make_matrix_forgenotyped( ped_file, g_id, true ); // A22(-1)            
            a.get_matrix("iA", iA, iA_id);            
            a.get_matrix("iA22", iA22, iA22_id);

            iA22.fread();
            iA.fread();

            iG.fread("iG.dmat");

            make_matrix(iA, iA_id, iA22, iA22_id, iG, g_id);            
            save_matrix(out_file);

            iA.fclear(); iA.clear();
            iA_id.clear(); iA_id.shrink_to_fit();
            iA22.fclear(); iA22.clear();
            iA22_id.clear(); iA22_id.shrink_to_fit();
            iG.fclear("iG.dmat"); iG.fclear(); iG.clear();
            g_id.clear(); g_id.shrink_to_fit();

            clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(const std::string &, const std::string &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(const std::string &, const std::string &, const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(const std::string &, const std::string &, const std::string &)" << '\n';
            throw;
        }
    }
    template void Hmat<float>::make_matrix(const std::string &g_file, const std::string &ped_file, const std::string &out_file);
    template void Hmat<double>::make_matrix(const std::string &g_file, const std::string &ped_file, const std::string &out_file);

    //===============================================================================================================
    template <typename T>
    void Hmat<T>::make_matrix(const std::string &g_file, const std::string &g_id_file, const std::string &ped_file, const std::string &out_file)
    {
        try
        {
            Gmat<T> g;
            evolm::matrix<T> iG;
            std::vector<std::int64_t> g_id;
            
            Amat<T> a;
            evolm::matrix<T> A22;
            std::vector<int64_t> a22_id;
            evolm::smatrix<T> iA;
            std::vector<int64_t> iA_id;
            evolm::matrix<T> iA22;
            std::vector<int64_t> iA22_id;

            g.make_matrix(g_file, g_id_file); // build G matrix from ASCII SNP file with IDs in sepaarate file
            g.scale_diag(0.01);
            g.get_ids(g_id); // here we need just IDs

            a.make_matrix_forgenotyped( ped_file, g_id, false );
            a.get_matrix("A22", A22, a22_id);
            A22.fread();
            
            g.scale_matrix(A22,0.25);

            A22.clear();
            a.clear();

            g.invert_matrix(true); // inverting as a full-store            
            g.save_matrix("iG.dmat");            
            g.clear();

            a.make_matrix( ped_file, true ); // A(-1)
            a.make_matrix_forgenotyped( ped_file, g_id, true ); // A22(-1)            
            a.get_matrix("iA", iA, iA_id);            
            a.get_matrix("iA22", iA22, iA22_id);

            iA22.fread();
            iA.fread();

            iG.fread("iG.dmat");

            make_matrix(iA, iA_id, iA22, iA22_id, iG, g_id);            
            save_matrix(out_file);

            iA.fclear(); iA.clear();
            iA_id.clear(); iA_id.shrink_to_fit();
            iA22.fclear(); iA22.clear();
            iA22_id.clear(); iA22_id.shrink_to_fit();
            iG.fclear("iG.dmat"); iG.fclear(); iG.clear();
            g_id.clear(); g_id.shrink_to_fit();

            clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(const std::string &, const std::string &, const std::string &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(const std::string &, const std::string &, const std::string &, const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(const std::string &, const std::string &, const std::string &, const std::string &)" << '\n';
            throw;
        }
    }
    template void Hmat<float>::make_matrix(const std::string &g_file, const std::string &g_id_file, const std::string &ped_file, const std::string &out_file);
    template void Hmat<double>::make_matrix(const std::string &g_file, const std::string &g_id_file, const std::string &ped_file, const std::string &out_file);

    //===============================================================================================================
    template <typename T>
    void Hmat<T>::clear()
    {
        try
        {
            h_values.clear(); h_values.shrink_to_fit();
            h_keys.clear(); h_keys.shrink_to_fit();
            hmat_id.clear(); hmat_id.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat<T>::clear()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat<T>::clear()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat<T>::clear()" << '\n';
            throw;
        }
    }
    template void Hmat<float>::clear();
    template void Hmat<double>::clear();
    //===============================================================================================================
    template <typename T>
    void Hmat<T>::make_matrix(const std::string &a_matr, // A(-1)
                              const std::string &a_ids,
                              const std::string &a_red_matr, // A22(-1)
                              const std::string &a_red_ids,
                              const std::string &g_matr, // G(-1)
                              const std::string &g_ids)
    {
        // Making single-step matriix, H(-1);
        // Assume all input matrices are symmetric and in L-store format
        try
        {
            Utilities2 u;

            const size_t dense_matr = 1001;
            const size_t sparse_matr = 2002;

            size_t a_matr_kind = u.fget_matrix_kind(a_matr);
            size_t a_red_matr_kind = u.fget_matrix_kind(a_red_matr);
            size_t g_matr_kind = u.fget_matrix_kind(g_matr);

            if ( g_matr_kind != 0 )
                throw std::string("Cannot correctly identify the matrix type!");

            if ( a_red_matr_kind != 0 )
                throw std::string("Cannot correctly identify the matrix type!");

            if ( a_matr_kind != 0 )
                throw std::string("Cannot correctly identify the matrix type!");

            if ( g_matr_kind != dense_matr )
                throw std::string("G matrix expected to be dense!");

            evolm::matrix<T> g_mat;
            g_mat.fread(g_matr);

            std::vector<std::int64_t> g_id;
            u.vect_from_binary(g_id, g_ids);

            std::vector<std::int64_t> a_id;
            u.vect_from_binary(a_id, a_ids);

            std::vector<std::int64_t> ar_id;
            u.vect_from_binary(ar_id, a_red_ids);
            
            if ( (a_matr_kind == dense_matr) && (a_red_matr_kind == dense_matr) )
            {
                evolm::matrix<T> a_mat;
                evolm::matrix<T> ar_mat;
                a_mat.fread(a_matr);
                ar_mat.fread(a_red_matr);
                make_matrix(a_mat, a_id, ar_mat, ar_id, g_mat, g_id);
            }
            
            if ( (a_matr_kind == sparse_matr) && (a_red_matr_kind == dense_matr) )
            {
                evolm::smatrix<T> a_mat;
                evolm::matrix<T> ar_mat;
                a_mat.fread(a_matr);
                ar_mat.fread(a_red_matr);
                make_matrix(a_mat, a_id, ar_mat, ar_id, g_mat, g_id);
            }

            if ( (a_matr_kind == dense_matr) && (a_red_matr_kind == sparse_matr) )
            {
                evolm::matrix<T> a_mat;
                evolm::smatrix<T> ar_mat;
                a_mat.fread(a_matr);
                ar_mat.fread(a_red_matr);
                make_matrix(a_mat, a_id, ar_mat, ar_id, g_mat, g_id);
            }

            if ( (a_matr_kind == sparse_matr) && (a_red_matr_kind == sparse_matr) )
            {
                evolm::smatrix<T> a_mat;
                evolm::smatrix<T> ar_mat;
                a_mat.fread(a_matr);
                ar_mat.fread(a_red_matr);
                make_matrix(a_mat, a_id, ar_mat, ar_id, g_mat, g_id);
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(const std::string&, const std::string&, const std::string&, const std::string&, const std::string&, const std::string&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(const std::string&, const std::string&, const std::string&, const std::string&, const std::string&, const std::string&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(const std::string&, const std::string&, const std::string&, const std::string&, const std::string&, const std::string&)" << '\n';
            throw;
        }
    }
    template void Hmat<float>::make_matrix(const std::string &a_matr,
                                           const std::string &a_ids,
                                           const std::string &a_red_matr,
                                           const std::string &a_red_ids,
                                           const std::string &g_matr,
                                           const std::string &g_ids);
    template void Hmat<double>::make_matrix(const std::string &a_matr,
                                            const std::string &a_ids,
                                            const std::string &a_red_matr,
                                            const std::string &a_red_ids,
                                            const std::string &g_matr,
                                            const std::string &g_ids);

    //===============================================================================================================
    template <typename T>
    void Hmat<T>::make_matrix(evolm::matrix<T> &a_matr, // A(-1)
                              std::vector<std::int64_t> &a_ids,
                              evolm::matrix<T> &a_red_matr, // A22(-1)
                              std::vector<std::int64_t> &a_red_ids,
                              evolm::matrix<T> &g_matr, // G(-1)
                              std::vector<std::int64_t> &g_ids)
    {
        // Making single-step matriix, H(-1);
        // Assume all input matrices are symmetric and in L-store format
        try
        {
            size_t expected_size = 0.5 * (a_ids.size() - 1) * a_ids.size() + a_ids.size();

            if (a_matr.size() > expected_size)
                throw std::string("The number of elements in the passed A(-1) matrix is greater then expected!");

            expected_size = 0.5 * (a_red_ids.size() - 1) * a_red_ids.size() + a_red_ids.size();

            if (a_red_matr.size() > expected_size)
                throw std::string("The number of elements in the passed redused A(-1) matrix is greater then expected!");

            expected_size = 0.5 * (g_ids.size() - 1) * g_ids.size() + g_ids.size();

            if (g_matr.size() > expected_size)
                throw std::string("The number of elements in the passed G(-1) matrix is greater then expected!");

            if (g_ids.size() > a_ids.size())
                throw std::string("The number of elements in the passed G(-1) maatrix is greater then the number of IDs in full A(-1) matrix!");

            if (g_ids.size() != a_red_ids.size())
                throw std::string("The number of elements in the passed G(-1) maatrix is not equal to the number of IDs in A22(-1) matrix!");

            Utilities2 u;

            if (!u.is_value_in_vect(a_red_ids, g_ids))
                throw std::string("There are IDs in G(-1) matrix which are not part of the reduced A(-1) matrix!");

            if (!u.is_value_in_vect(a_ids, g_ids))
                throw std::string("There are IDs in G(-1) matrix which are not part of the A(-1) matrix!");

            evolm::matrix<size_t> shapeofa;
            shapeofa = a_matr.shape();

            evolm::matrix<size_t> shapeofra;
            shapeofra = a_red_matr.shape();

            evolm::matrix<size_t> shapeofg;
            shapeofg = g_matr.shape();

            if ((shapeofa[0] != a_ids.size()) || (shapeofa[0] != shapeofa[1]))
                throw std::string("The dimension of A(-1) is not correct!");

            if ((shapeofra[0] != a_red_ids.size()) || (shapeofra[0] != shapeofra[1]))
                throw std::string("The dimension of reduced A(-1) is not correct!");

            if ((shapeofg[0] != g_ids.size()) || (shapeofg[0] != shapeofg[1]))
                throw std::string("The dimension of G(-1) is not correct!");

            if ((shapeofg[0] != shapeofra[0]) || (shapeofg[1] != shapeofra[1]))
                throw std::string("The dimension of G(-1) is not the same as A22(-1)!");

            shapeofa.clear();
            shapeofra.clear();
            shapeofg.clear();

            // ----------------- G22(-1) = G(-1) - A22(-1) -----------------

            std::vector<size_t> pos_map;
            for (size_t i = 0; i < g_ids.size(); i++)
                pos_map.push_back(u.find_invect(a_red_ids, g_ids[i]));

#pragma omp parallel for
            for (size_t i = 0; i < g_ids.size(); i++)
            {
                size_t row = pos_map[i];

                for (size_t j = 0; j <= i; j++)
                {
                    size_t col = pos_map[j];

                    if (row >= col)
                        g_matr(i, j) = g_matr(i, j) - a_red_matr(row, col);
                    else
                        g_matr(i, j) = g_matr(i, j) - a_red_matr(col, row);
                }
            }

            a_red_matr.clear();

            // --------------------- A(-1) + G22(-1) -----------------------

            pos_map.clear();
            pos_map.shrink_to_fit();

            for (size_t i = 0; i < g_ids.size(); i++)
                pos_map.push_back(u.find_invect(a_ids, g_ids[i]));

            hmat_id = a_ids;

            // define size of h-vectors to resize
            size_t h_size = a_ids.size() * ( a_ids.size() + 1 ) / 2;

            h_values.resize(h_size);
            h_keys.resize(h_size);

#pragma omp parallel for
            for (size_t i = 0; i < g_ids.size(); i++)
            {
                size_t row = pos_map[i];

                for (size_t j = 0; j <= i; j++)
                {
                    size_t col = pos_map[j];

                    if (row >= col)
                    {
                        a_matr(row, col) = a_matr(row, col) + g_matr(i, j);
                    }
                    else
                    {
                        a_matr(col, row) = a_matr(col, row) + g_matr(i, j);
                    }
                }
            }

#pragma omp parallel for
            for (size_t i = 0; i < a_matr.size(); i++)
            {
                h_values[i] = a_matr[i];
                h_keys[i] = i;
            }

            a_matr.clear();
            g_matr.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(evolm::matrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(evolm::matrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(evolm::matrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&)" << '\n';
            throw;
        }
    }
    template void Hmat<float>::make_matrix(evolm::matrix<float> &a_matr,
                                           std::vector<std::int64_t> &a_ids,
                                           evolm::matrix<float> &a_red_matr,
                                           std::vector<std::int64_t> &a_red_ids,
                                           evolm::matrix<float> &g_matr,
                                           std::vector<std::int64_t> &g_ids);
    template void Hmat<double>::make_matrix(evolm::matrix<double> &a_matr,
                                            std::vector<std::int64_t> &a_ids,
                                            evolm::matrix<double> &a_red_matr,
                                            std::vector<std::int64_t> &a_red_ids,
                                            evolm::matrix<double> &g_matr,
                                            std::vector<std::int64_t> &g_ids);

    //===============================================================================================================
    template <typename T>
    void Hmat<T>::make_matrix(evolm::smatrix<T> &a_matr, // A(-1)
                              std::vector<std::int64_t> &a_ids,
                              evolm::matrix<T> &a_red_matr, // A22(-1)
                              std::vector<std::int64_t> &a_red_ids,
                              evolm::matrix<T> &g_matr, // G(-1)
                              std::vector<std::int64_t> &g_ids)
    {
        // Making single-step matriix, H(-1);
        // Assume all input matrices are symmetric and in L-store format
        try
        {
            size_t expected_size = 0.5 * (a_ids.size() - 1) * a_ids.size() + a_ids.size();

            if (a_matr.size() > expected_size)
                throw std::string("The number of elements in the passed A(-1) matrix is greater then expected!");

            expected_size = 0.5 * (a_red_ids.size() - 1) * a_red_ids.size() + a_red_ids.size();

            if (a_red_matr.size() > expected_size)
                throw std::string("The number of elements in the passed redused A(-1) matrix is greater then expected!");

            expected_size = 0.5 * (g_ids.size() - 1) * g_ids.size() + g_ids.size();

            if (g_matr.size() > expected_size)
                throw std::string("The number of elements in the passed G(-1) matrix is greater then expected!");

            if (g_ids.size() > a_ids.size())
                throw std::string("The number of elements in the passed G(-1) maatrix is greater then the number of IDs in full A(-1) matrix!");

            if (g_ids.size() != a_red_ids.size())
                throw std::string("The number of elements in the passed G(-1) maatrix is not equal to the number of IDs in A22(-1) matrix!");

            Utilities2 u;

            if (!u.is_value_in_vect(a_red_ids, g_ids))
                throw std::string("There are IDs in G(-1) matrix which are not part of the reduced A(-1) matrix!");

            if (!u.is_value_in_vect(a_ids, g_ids))
                throw std::string("There are IDs in G(-1) matrix which are not part of the A(-1) matrix!");

            evolm::matrix<size_t> shapeofra;
            shapeofra = a_red_matr.shape();

            evolm::matrix<size_t> shapeofg;
            shapeofg = g_matr.shape();

            if ((a_matr.nrows() != a_ids.size()) || (a_matr.nrows() != a_matr.ncols()))
                throw std::string("The dimension of A(-1) is not correct!");

            if ((shapeofra[0] != a_red_ids.size()) || (shapeofra[0] != shapeofra[1]))
                throw std::string("The dimension of reduced A(-1) is not correct!");

            if ((shapeofg[0] != g_ids.size()) || (shapeofg[0] != shapeofg[1]))
                throw std::string("The dimension of G(-1) is not correct!");

            if ((shapeofg[0] != shapeofra[0]) || (shapeofg[1] != shapeofra[1]))
                throw std::string("The dimension of G(-1) is not the same as A22(-1)!");

            shapeofra.clear();
            shapeofg.clear();

            // ----------------- G22(-1) = G(-1) - A22(-1) -----------------

            std::vector<size_t> pos_map;
            for (size_t i = 0; i < g_ids.size(); i++)
                pos_map.push_back(u.find_invect(a_red_ids, g_ids[i]));

#pragma omp parallel for
            for (size_t i = 0; i < g_ids.size(); i++)
            {
                size_t row = pos_map[i];

                for (size_t j = 0; j <= i; j++)
                {
                    size_t col = pos_map[j];

                    if (row >= col)
                        g_matr(i, j) = g_matr(i, j) - a_red_matr(row, col);
                    else
                        g_matr(i, j) = g_matr(i, j) - a_red_matr(col, row);
                }
            }

            a_red_matr.clear();

            // ----------------- H = A(-1) not in G22(-1) -----------------

            std::vector<std::int64_t> a_id_notin_g;
            u.check_id( g_ids, a_ids, a_id_notin_g ); // find ids of A(-1) which are not in G(-1)

            std::vector<size_t> matr_keys; // keys which are not part of genotyped
            a_matr.get_keyslist(matr_keys);

            std::vector<size_t> matr_keys_in_g; // keys which are part of genotyped

            std::vector<T> t_values(matr_keys.size());
            std::vector<std::int64_t> t_keys(matr_keys.size(), -2);
            std::vector<std::int64_t> t_matr_keys_in_g(matr_keys.size(), -2);

#pragma omp parallel for
            for (size_t i = 0; i < matr_keys.size(); i++) // write keys and values for non-genotyped
            {
                size_t irow = a_matr.row_insym(matr_keys[i]);
                size_t icol = a_matr.col_insym(matr_keys[i], irow);
                size_t real_row = a_ids[irow];
                size_t real_col = a_ids[icol];
                int found_row_pos = u.find_invect( a_id_notin_g, (std::int64_t)real_row );
                int found_col_pos = u.find_invect( a_id_notin_g, (std::int64_t)real_col );
                if ( found_col_pos != -1 ) // all ids which are not genotyped
                {
                    t_values[i] = a_matr[ matr_keys[i] ];
                    t_keys[i] = matr_keys[i];
                }
                else if ( found_row_pos == -1 && found_col_pos == -1) // store all genotyped ids
                {
                    t_matr_keys_in_g[i] = matr_keys[i];
                }
            }

            for (size_t i = 0; i < t_values.size(); i++)
            {
                if ( t_keys[i] != -2 )
                {
                    h_values.push_back(t_values[i]);
                    h_keys.push_back(t_keys[i]);
                }

                if ( t_matr_keys_in_g[i] != -2 )
                    matr_keys_in_g.push_back(t_matr_keys_in_g[i]);
            }

            t_values.clear(); t_values.shrink_to_fit();
            t_keys.clear(); t_keys.shrink_to_fit();
            t_matr_keys_in_g.clear(); t_matr_keys_in_g.shrink_to_fit();

            // ------------------- H = A(-1) + G22(-1) ---------------------
           
            for (size_t i = 0; i < matr_keys_in_g.size(); i++) // correct genotyped values by related A-matrix values
            {
                size_t ikey = matr_keys_in_g[i];
                size_t irow = a_matr.row_insym(ikey);
                size_t icol = a_matr.col_insym(ikey, irow);
                size_t real_row = a_ids[irow];
                size_t real_col = a_ids[icol];

                size_t g_row = u.find_invect( g_ids, (std::int64_t)real_row );
                size_t g_col = u.find_invect( g_ids, (std::int64_t)real_col );
                
                if (g_row >= g_col)
                    g_matr(g_row, g_col) = a_matr[ikey] + g_matr(g_row, g_col);
                else
                    g_matr(g_col, g_row) = a_matr[ikey] + g_matr(g_col, g_row);
            }

            pos_map.clear();
            pos_map.shrink_to_fit();

            for (size_t i = 0; i < g_ids.size(); i++)
                pos_map.push_back(u.find_invect(a_ids, g_ids[i]));

            hmat_id = a_ids;

            size_t key = 0;
            size_t row = 0;
            size_t col = 0;

            for (size_t i = 0; i < g_ids.size(); i++) // write keys and values for all genotyped including corrected
            {
                row = pos_map[i];

                for (size_t j = 0; j <= i; j++)
                {
                    col = pos_map[j];

                    if (row >= col)
                        key = row * (row + 1) / 2 + col;
                    else
                        key = col * (col + 1) / 2 + row;

                    h_values.push_back( g_matr(i, j) );
                    h_keys.push_back( key );
                }
            }

            a_matr.clear();
            g_matr.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(evolm::smatrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(evolm::smatrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(evolm::smatrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&)" << '\n';
            throw;
        }
    }
    template void Hmat<float>::make_matrix(evolm::smatrix<float> &a_matr,
                                           std::vector<std::int64_t> &a_ids,
                                           evolm::matrix<float> &a_red_matr,
                                           std::vector<std::int64_t> &a_red_ids,
                                           evolm::matrix<float> &g_matr,
                                           std::vector<std::int64_t> &g_ids);
    template void Hmat<double>::make_matrix(evolm::smatrix<double> &a_matr,
                                            std::vector<std::int64_t> &a_ids,
                                            evolm::matrix<double> &a_red_matr,
                                            std::vector<std::int64_t> &a_red_ids,
                                            evolm::matrix<double> &g_matr,
                                            std::vector<std::int64_t> &g_ids);

    //===============================================================================================================
    template <typename T>
    void Hmat<T>::make_matrix(evolm::smatrix<T> &a_matr, // A(-1)
                              std::vector<std::int64_t> &a_ids,
                              evolm::smatrix<T> &a_red_matr, // A22(-1)
                              std::vector<std::int64_t> &a_red_ids,
                              evolm::matrix<T> &g_matr, // G(-1)
                              std::vector<std::int64_t> &g_ids)
    {
        // Making single-step matriix, H(-1);
        // Assume all input matrices are symmetric and in L-store format
        try
        {
            size_t expected_size = 0.5 * (a_ids.size() - 1) * a_ids.size() + a_ids.size();

            if (a_matr.size() > expected_size)
                throw std::string("The number of elements in the passed A(-1) matrix is greater then expected!");

            expected_size = 0.5 * (a_red_ids.size() - 1) * a_red_ids.size() + a_red_ids.size();

            if (a_red_matr.size() > expected_size)
                throw std::string("The number of elements in the passed redused A(-1) matrix is greater then expected!");

            expected_size = 0.5 * (g_ids.size() - 1) * g_ids.size() + g_ids.size();

            if (g_matr.size() > expected_size)
                throw std::string("The number of elements in the passed G(-1) matrix is greater then expected!");

            if (g_ids.size() > a_ids.size())
                throw std::string("The number of elements in the passed G(-1) maatrix is greater then the number of IDs in full A(-1) matrix!");

            if (g_ids.size() != a_red_ids.size())
                throw std::string("The number of elements in the passed G(-1) maatrix is not equal to the number of IDs in A22(-1) matrix!");

            Utilities2 u;

            if (!u.is_value_in_vect(a_red_ids, g_ids))
                throw std::string("There are IDs in G(-1) matrix which are not part of the reduced A(-1) matrix!");

            if (!u.is_value_in_vect(a_ids, g_ids))
                throw std::string("There are IDs in G(-1) matrix which are not part of the A(-1) matrix!");

            evolm::matrix<size_t> shapeofg;
            shapeofg = g_matr.shape();

            if ((a_matr.nrows() != a_ids.size()) || (a_matr.nrows() != a_matr.ncols()))
                throw std::string("The dimension of A(-1) is not correct!");

            if ((a_red_matr.nrows() != a_red_ids.size()) || (a_red_matr.nrows() != a_red_matr.ncols()))
                throw std::string("The dimension of reduced A(-1) is not correct!");

            if ((shapeofg[0] != g_ids.size()) || (shapeofg[0] != shapeofg[1]))
                throw std::string("The dimension of G(-1) is not correct!");

            if ((shapeofg[0] != a_red_matr.nrows()) || (shapeofg[1] != a_red_matr.ncols()))
                throw std::string("The dimension of G(-1) is not the same as A22(-1)!");

            shapeofg.clear();

            // ----------------- G22(-1) = G(-1) - A22(-1) -----------------

            std::vector<size_t> pos_map;
            for (size_t i = 0; i < g_ids.size(); i++)
                pos_map.push_back(u.find_invect(a_red_ids, g_ids[i]));

            T zerro_val = (T)0;

#pragma omp parallel for
            for (size_t i = 0; i < g_ids.size(); i++)
            {
                T a_val = (T)0;

                size_t row = pos_map[i];

                for (size_t j = 0; j <= i; j++)
                {
                    size_t col = pos_map[j];

                    if (row >= col)
                        a_val = a_red_matr.get_nonzero(row, col);
                    else
                        a_val = a_red_matr.get_nonzero(col, row);
                    
                    if ( a_val != zerro_val )
                        g_matr(i, j) = g_matr(i, j) - a_val;
                }
            }

            a_red_matr.clear();

            // ----------------- H = A(-1) not in G22(-1) -----------------

            std::vector<std::int64_t> a_id_notin_g;
            u.check_id( g_ids, a_ids, a_id_notin_g ); // find ids of A(-1) which are not in G(-1)

            size_t key = 0;
            T a_value = (T)0;

            for (size_t i = 0; i < a_ids.size(); i++)
            {
                for (size_t j = 0; j < a_id_notin_g.size(); j++)
                {
                    if ( j > i )
                        continue;

                    key = i * (i + 1) / 2 + j;

                    a_value = a_matr.get_nonzero( key );

                    if ( a_value != (T)0 )
                    {
                        h_values.push_back( a_matr[key] );
                        h_keys.push_back( key );
                    }
                }
            }

            // --------------------- A(-1) + G22(-1) -----------------------

            pos_map.clear();
            pos_map.shrink_to_fit();

            for (size_t i = 0; i < g_ids.size(); i++)
                pos_map.push_back(u.find_invect(a_ids, g_ids[i]));

            hmat_id = a_ids;

            key = 0;
            a_value = (T)0;
            size_t row = 0;
            size_t col = 0;

            for (size_t i = 0; i < g_ids.size(); i++)
            {
                row = pos_map[i];

                for (size_t j = 0; j <= i; j++)
                {
                    col = pos_map[j];

                    if (row >= col)
                        key = row * (row + 1) / 2 + col;
                    else
                        key = col * (col + 1) / 2 + row;

                    a_value = a_matr.get_nonzero( key );
                    h_values.push_back( a_value + g_matr(i, j) );
                    h_keys.push_back( key );
                }
            }

            a_matr.clear();
            g_matr.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(evolm::smatrix<T>&, std::vector<std::int64_t>&, evolm::smatrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(evolm::smatrix<T>&, std::vector<std::int64_t>&, evolm::smatrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(evolm::smatrix<T>&, std::vector<std::int64_t>&, evolm::smatrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&)" << '\n';
            throw;
        }
    }
    template void Hmat<float>::make_matrix(evolm::smatrix<float> &a_matr,
                                           std::vector<std::int64_t> &a_ids,
                                           evolm::smatrix<float> &a_red_matr,
                                           std::vector<std::int64_t> &a_red_ids,
                                           evolm::matrix<float> &g_matr,
                                           std::vector<std::int64_t> &g_ids);
    template void Hmat<double>::make_matrix(evolm::smatrix<double> &a_matr,
                                            std::vector<std::int64_t> &a_ids,
                                            evolm::smatrix<double> &a_red_matr,
                                            std::vector<std::int64_t> &a_red_ids,
                                            evolm::matrix<double> &g_matr,
                                            std::vector<std::int64_t> &g_ids);

    //===============================================================================================================
    template <typename T>
    void Hmat<T>::make_matrix(evolm::matrix<T> &a_matr, // A(-1)
                              std::vector<std::int64_t> &a_ids,
                              evolm::smatrix<T> &a_red_matr, // A22(-1)
                              std::vector<std::int64_t> &a_red_ids,
                              evolm::matrix<T> &g_matr, // G(-1)
                              std::vector<std::int64_t> &g_ids)
    {
        // Making single-step matriix, H(-1);
        // Assume all input matrices are symmetric and in L-store format
        try
        {
            size_t expected_size = 0.5 * (a_ids.size() - 1) * a_ids.size() + a_ids.size();

            if (a_matr.size() > expected_size)
                throw std::string("The number of elements in the passed A(-1) matrix is greater then expected!");

            expected_size = 0.5 * (a_red_ids.size() - 1) * a_red_ids.size() + a_red_ids.size();

            if (a_red_matr.size() > expected_size)
                throw std::string("The number of elements in the passed redused A(-1) matrix is greater then expected!");

            expected_size = 0.5 * (g_ids.size() - 1) * g_ids.size() + g_ids.size();

            if (g_matr.size() > expected_size)
                throw std::string("The number of elements in the passed G(-1) matrix is greater then expected!");

            if (g_ids.size() > a_ids.size())
                throw std::string("The number of elements in the passed G(-1) maatrix is greater then the number of IDs in full A(-1) matrix!");

            if (g_ids.size() != a_red_ids.size())
                throw std::string("The number of elements in the passed G(-1) maatrix is not equal to the number of IDs in A22(-1) matrix!");

            Utilities2 u;

            if (!u.is_value_in_vect(a_red_ids, g_ids))
                throw std::string("There are IDs in G(-1) matrix which are not part of the reduced A(-1) matrix!");

            if (!u.is_value_in_vect(a_ids, g_ids))
                throw std::string("There are IDs in G(-1) matrix which are not part of the A(-1) matrix!");

            evolm::matrix<size_t> shapeofa;
            shapeofa = a_matr.shape();

            evolm::matrix<size_t> shapeofg;
            shapeofg = g_matr.shape();

            if ((shapeofa[0] != a_ids.size()) || (shapeofa[0] != shapeofa[1]))
                throw std::string("The dimension of A(-1) is not correct!");

            if ((a_red_matr.nrows() != a_red_ids.size()) || (a_red_matr.nrows() != a_red_matr.ncols()))
                throw std::string("The dimension of reduced A(-1) is not correct!");

            if ((shapeofg[0] != g_ids.size()) || (shapeofg[0] != shapeofg[1]))
                throw std::string("The dimension of G(-1) is not correct!");

            if ((shapeofg[0] != a_red_matr.nrows()) || (shapeofg[1] != a_red_matr.ncols()))
                throw std::string("The dimension of G(-1) is not the same as A22(-1)!");

            shapeofg.clear();

            // ----------------- G22(-1) = G(-1) - A22(-1) -----------------

            std::vector<size_t> pos_map;
            for (size_t i = 0; i < g_ids.size(); i++)
                pos_map.push_back(u.find_invect(a_red_ids, g_ids[i]));

            T zerro_val = (T)0;

#pragma omp parallel for
            for (size_t i = 0; i < g_ids.size(); i++)
            {
                T a_val = (T)0;

                size_t row = pos_map[i];

                for (size_t j = 0; j <= i; j++)
                {
                    size_t col = pos_map[j];

                    if (row >= col)
                        a_val = a_red_matr.get_nonzero(row, col);
                    else
                        a_val = a_red_matr.get_nonzero(col, row);
                    
                    if ( a_val != zerro_val )
                        g_matr(i, j) = g_matr(i, j) - a_val;
                }
            }

            a_red_matr.clear();

            // --------------------- A(-1) + G22(-1) -----------------------

            pos_map.clear();
            pos_map.shrink_to_fit();

            for (size_t i = 0; i < g_ids.size(); i++)
                pos_map.push_back(u.find_invect(a_ids, g_ids[i]));

            hmat_id = a_ids;

            size_t h_size = a_matr.size() * ( a_matr.size() + 1 )/2;

            h_values.resize(h_size);
            h_keys.resize(h_size);

#pragma omp parallel for
            for (size_t i = 0; i < g_ids.size(); i++)
            {
                size_t row = pos_map[i];

                for (size_t j = 0; j <= i; j++)
                {
                    size_t col = pos_map[j];

                    if (row >= col)
                        a_matr(row, col) = a_matr(row, col) + g_matr(i, j);
                    else
                        a_matr(col, row) = a_matr(col, row) + g_matr(i, j);
                }
            }

#pragma omp parallel for
            for (size_t i = 0; i < a_matr.size(); i++)
            {
                h_values[i] = a_matr[i];
                h_keys[i] = i;
            }

            a_matr.clear();
            g_matr.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(evolm::matrix<T>&, std::vector<std::int64_t>&, evolm::smatrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(evolm::matrix<T>&, std::vector<std::int64_t>&, evolm::smatrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat<T>::make_matrix(evolm::matrix<T>&, std::vector<std::int64_t>&, evolm::smatrix<T>&, std::vector<std::int64_t>&, evolm::matrix<T>&, std::vector<std::int64_t>&)" << '\n';
            throw;
        }
    }
    template void Hmat<float>::make_matrix(evolm::matrix<float> &a_matr,
                                           std::vector<std::int64_t> &a_ids,
                                           evolm::smatrix<float> &a_red_matr,
                                           std::vector<std::int64_t> &a_red_ids,
                                           evolm::matrix<float> &g_matr,
                                           std::vector<std::int64_t> &g_ids);
    template void Hmat<double>::make_matrix(evolm::matrix<double> &a_matr,
                                            std::vector<std::int64_t> &a_ids,
                                            evolm::smatrix<double> &a_red_matr,
                                            std::vector<std::int64_t> &a_red_ids,
                                            evolm::matrix<double> &g_matr,
                                            std::vector<std::int64_t> &g_ids);

    //===============================================================================================================
    template <typename T>
    void Hmat<T>::save_matrix(const std::string &out_fname)
    {
        // Note, we operate with L-stored format, hence return lower triangular part
        try
        {
            Utilities2 u;
            u.fwrite_matrix(out_fname, h_values, h_keys, hmat_id);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat<T>::get_matrix(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat<T>::get_matrix(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat<T>::get_matrix(const std::string &)" << '\n';
            throw;
        }
    }
    template void Hmat<float>::save_matrix(const std::string &out_fname);
    template void Hmat<double>::save_matrix(const std::string &out_fname);

    //===============================================================================================================
    template <typename T>
    void Hmat<T>::get_matrix(std::vector<T> &values, std::vector<size_t> &keys, std::vector<std::int64_t> &id)
    {
        // Note, we operate with L-stored format, hence return lower triangular part
        try
        {
            values = h_values;
            keys = h_keys;
            id = hmat_id;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat<T>::get_matrix(std::vector<T> &, std::vector<size_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat<T>::get_matrix(std::vector<T> &, std::vector<size_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat<T>::get_matrix(std::vector<T> &, std::vector<size_t> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }
    template void Hmat<float>::get_matrix(std::vector<float> &values, std::vector<size_t> &keys, std::vector<std::int64_t> &id);
    template void Hmat<double>::get_matrix(std::vector<double> &values, std::vector<size_t> &keys, std::vector<std::int64_t> &id);

    //===============================================================================================================
    /*template <typename T>
    void Hmat<T>::get_matrix(evolm::matrix<T> &arr, std::vector<std::int64_t> &ids)
    {
        // Note, we operate with L-stored format, hence return lower triangular part
        try
        {
            Utilities2 u;

            if (IsEmpty.H) // was used sparse pipeline
            {
                H.resize(hmat_id.size());
                if ( H_s.is_ondisk() )
                    H_s.fread();
                u.sparse_to_dense(H_s, H);
                H_s.fwrite();
            }

            arr = H;
            ids = hmat_id;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat<T>::get_matrix(evolm::matrix<T> &, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat<T>::get_matrix(evolm::matrix<T> &, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat<T>::get_matrix(evolm::matrix<T> &, std::vector<std::int64_t>&, bool)" << '\n';
            throw;
        }
    }
    template void Hmat<float>::get_matrix(evolm::matrix<float> &arr, std::vector<std::int64_t> &ids);
    template void Hmat<double>::get_matrix(evolm::matrix<double> &arr, std::vector<std::int64_t> &ids);

    //===============================================================================================================
    template <typename T>
    void Hmat<T>::get_matrix(evolm::smatrix<T> &arr, std::vector<std::int64_t> &ids)
    {
        // Note, we operate with L-stored format, hence return lower triangular part
        try
        {
            Utilities2 u;

            if (IsEmpty.H_s) // was used dense pipeline
            {
                H_s.resize(hmat_id.size());
                if ( H.is_ondisk() )
                    H.fread();
                u.dense_to_sparse(H, H_s);
                H.fwrite();
            }

            arr = H_s;
            ids = hmat_id;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat<T>::get_matrix(evolm::smatrix<T> &, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat<T>::get_matrix(evolm::smatrix<T> &, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat<T>::get_matrix(evolm::smatrix<T> &, std::vector<std::int64_t>&, bool)" << '\n';
            throw;
        }
    }
    template void Hmat<float>::get_matrix(evolm::smatrix<float> &arr, std::vector<std::int64_t> &ids);
    template void Hmat<double>::get_matrix(evolm::smatrix<double> &arr, std::vector<std::int64_t> &ids);

    //===============================================================================================================
    template <typename T>
    void Hmat<T>::save_matrix(const std::string &arr, const std::string &ids)
    {
        // This method is for Python interfacing;
        // Note, we operate with L-stored format, hence return lower triangular part
        try
        {
            Utilities2 u;

            if (IsEmpty.H_s) // was used dense pipeline
            {
                if ( H.is_ondisk() )
                    H.fread();
                H.fwrite(arr);
            }
            else if (IsEmpty.H) // was used sparse pipeline
            {
                if ( H_s.is_ondisk() )
                    H_s.fread();
                H_s.fwrite(arr);
            }
            else
                throw std::string("H matrix is empty!");

            u.vect_to_binary(hmat_id, ids);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat<T>::save_matrix(std::vector<T> &, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat<T>::save_matrix(std::vector<T> &, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat<T>::save_matrix(std::vector<T> &, std::vector<std::int64_t>&, bool)" << '\n';
            throw;
        }
    }
    template void Hmat<float>::save_matrix(const std::string &arr, const std::string &ids);
    template void Hmat<double>::save_matrix(const std::string &arr, const std::string &ids);
    */
    //===============================================================================================================


    //===============================================================================================================


    //===============================================================================================================

}