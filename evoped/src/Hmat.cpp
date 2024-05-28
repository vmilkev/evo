#include "Hmat.hpp"

namespace evoped
{
    //===============================================================================================================
    template <typename T>
    Hmat<T>::Hmat()
    {
        IsEmpty.H = true;
        IsEmpty.H_s = true;
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
    void Hmat<T>::clear()
    {
        try
        {
            H.fclear();
            H.clear();
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
            {
                pos_map.push_back(u.find_invect(a_red_ids, g_ids[i]));
            }

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

            H = a_matr;

            a_matr.clear();

            hmat_id = a_ids;

#pragma omp parallel for
            for (size_t i = 0; i < g_ids.size(); i++)
            {
                size_t row = pos_map[i];

                for (size_t j = 0; j <= i; j++)
                {
                    size_t col = pos_map[j];

                    if (row >= col)
                        H(row, col) = H(row, col) + g_matr(i, j);
                    else
                        H(col, row) = H(col, row) + g_matr(i, j);
                }
            }

            //H.fwrite();

            IsEmpty.H = false;

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
            {
                pos_map.push_back(u.find_invect(a_red_ids, g_ids[i]));
            }

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

            //--------------------------------------
            /*
            H_s.resize( a_ids.size() );

            std::vector<std::int64_t> gmat_keys;
            std::vector<std::int64_t> amat_in_gmat_keys; // whose keys which coincide with gmat_keys
            std::vector<std::int64_t> amat_not_gmat_keys;
            std::vector<std::int64_t> amat_keys;
            
            size_t key = 0;

            for (size_t i = 0; i < g_ids.size(); i++)
            {
                size_t row = pos_map[i];

                for (size_t j = 0; j <= i; j++)
                {
                    size_t col = pos_map[j];

                    if (row >= col)
                    {
                        key = row * (row + 1) / 2 + col;
                        H_s[key] = g_matr(i, j);
                    }
                    else
                    {
                        key = col * (col + 1) / 2 + row;
                        H_s[key] = g_matr(i, j);
                    }
                    gmat_keys.push_back(key);
                }
            }

            a_matr.get_keyslist( amat_keys );
            u.get_what_in_vect( gmat_keys, amat_keys, amat_in_gmat_keys);
            u.get_missing_in_vect( amat_in_gmat_keys, amat_keys, amat_not_gmat_keys);

            std::cout<<"amat_keys.size(): "<<amat_keys.size()<<"\n";
            std::cout<<"gmat_keys.size(): "<<gmat_keys.size()<<"\n";
            std::cout<<"amat_in_gmat_keys.size(): "<<amat_in_gmat_keys.size()<<"\n";
            std::cout<<"amat_not_gmat_keys.size(): "<<amat_not_gmat_keys.size()<<"\n";

            for (size_t i = 0; i < amat_in_gmat_keys.size(); i++)
                H_s[ amat_in_gmat_keys[i] ] = a_matr[ amat_in_gmat_keys[i] ] + H_s[ amat_in_gmat_keys[i] ];

            for (size_t i = 0; i < amat_not_gmat_keys.size(); i++)
                H_s[ amat_not_gmat_keys[i] ] = a_matr[ amat_not_gmat_keys[i] ];

            hmat_id = a_ids;
            */
            //--------------------------------------

            /**/
            H_s = a_matr;

            a_matr.clear();

            hmat_id = a_ids;

            size_t key = 0;

            for (size_t i = 0; i < g_ids.size(); i++)
            {
                size_t row = pos_map[i];

                for (size_t j = 0; j <= i; j++)
                {
                    size_t col = pos_map[j];

                    if (row >= col)
                    {
                        key = row * (row + 1) / 2 + col;
                        //H_s(row, col) = H_s(row, col) + g_matr(i, j);
                    }
                    else
                    {
                        key = col * (col + 1) / 2 + row;
                        //H_s(col, row) = H_s(col, row) + g_matr(i, j);
                    }
                    H_s[key] = H_s[key] + g_matr(i, j);
                }
            }

            //H_s.fwrite();

            IsEmpty.H_s = false;

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
            {
                pos_map.push_back(u.find_invect(a_red_ids, g_ids[i]));
            }

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

            H_s = a_matr;

            a_matr.clear();

            hmat_id = a_ids;

            for (size_t i = 0; i < g_ids.size(); i++)
            {
                size_t row = pos_map[i];

                for (size_t j = 0; j <= i; j++)
                {
                    size_t col = pos_map[j];

                    if (row >= col)
                        H_s(row, col) = H_s(row, col) + g_matr(i, j);
                    else
                        H_s(col, row) = H_s(col, row) + g_matr(i, j);
                }
            }

            //H_s.fwrite();

            IsEmpty.H_s = false;

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
            {
                pos_map.push_back(u.find_invect(a_red_ids, g_ids[i]));
            }

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

            H = a_matr;

            a_matr.clear();

            hmat_id = a_ids;

#pragma omp parallel for
            for (size_t i = 0; i < g_ids.size(); i++)
            {
                size_t row = pos_map[i];

                for (size_t j = 0; j <= i; j++)
                {
                    size_t col = pos_map[j];

                    if (row >= col)
                        H(row, col) = H(row, col) + g_matr(i, j);
                    else
                        H(col, row) = H(col, row) + g_matr(i, j);
                }
            }

            //H.fwrite();

            IsEmpty.H = false;

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

    //===============================================================================================================


    //===============================================================================================================


    //===============================================================================================================

}