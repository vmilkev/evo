#include "Hmat.hpp"

namespace evoped
{
    //===============================================================================================================

    Hmat::Hmat()
    {
    }

    //===============================================================================================================

    Hmat::~Hmat()
    {
        try
        {
            clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat::~Hmat()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::~Hmat()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::~Hmat()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Hmat::clear()
    {
        try
        {
            H.clear();
            H.fclear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat::clear()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::clear()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::clear()" << '\n';
            throw;
        }
    }




    //===============================================================================================================

    void Hmat::make_matrix(std::vector<double>& a_matr,
                           std::vector<std::int64_t>& a_ids,
                           std::vector<double>& a_red_matr,
                           std::vector<std::int64_t>& a_red_ids,
                           std::vector<double>& g_matr,
                           std::vector<std::int64_t>& g_ids)
    {
        // Making single-step matriix, H(-1);
        // Assume all input matrices are symmetric and in L-store format
        try
        {
            size_t expected_size = 0.5 * (a_ids.size() - 1) * a_ids.size() + a_ids.size();

            if ( a_matr.size() > expected_size )
                throw std::string("The number of elements in the passed A(-1) matrix is greater then expected!");

            expected_size = 0.5 * (a_red_ids.size() - 1) * a_red_ids.size() + a_red_ids.size();

            if ( a_red_matr.size() > expected_size )
                throw std::string("The number of elements in the passed redused A(-1) matrix is greater then expected!");

            expected_size = 0.5 * (g_ids.size() - 1) * g_ids.size() + g_ids.size();

            if ( g_matr.size() > expected_size )
                throw std::string("The number of elements in the passed G(-1) matrix is greater then expected!");

            if ( g_ids.size() > a_ids.size() )
                throw std::string("The number of elements in the passed G(-1) maatrix is greater then the number of IDs in full A(-1) matrix!");

            if ( !is_value_in_vect(a_red_ids, g_ids) )
                throw std::string("There are IDs in G(-1) matrix which are not part of the reduced A(-1) matrix!");

            if ( !is_value_in_vect(a_ids, g_ids) )
                throw std::string("There are IDs in G(-1) matrix which are not part of the A(-1) matrix!");

            evolm::matrix<double> amatr( a_ids.size() );
            amatr.from_vector(a_matr);
            //amatr.fwrite(); ???

            evolm::matrix<double> ramatr( a_red_ids.size() );
            ramatr.from_vector(a_red_matr);
            //ramatr.fwrite(); ???

            evolm::matrix<double> gmatr( g_ids.size() );
            gmatr.from_vector(g_matr);
            //gmatr.fwrite(); ???

            make_matrix(amatr, a_ids, ramatr, a_red_ids, gmatr, g_ids);

            amatr.clear();
            amatr.fclear();
            ramatr.clear();
            ramatr.fclear();
            gmatr.clear();
            gmatr.fclear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat::make_matrix(std::vector<double>&, std::vector<std::int64_t>&, std::vector<double>&, std::vector<std::int64_t>&, std::vector<double>&, std::vector<std::int64_t>&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::make_matrix(std::vector<double>&, std::vector<std::int64_t>&, std::vector<double>&, std::vector<std::int64_t>&, std::vector<double>&, std::vector<std::int64_t>&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::make_matrix(std::vector<double>&, std::vector<std::int64_t>&, std::vector<double>&, std::vector<std::int64_t>&, std::vector<double>&, std::vector<std::int64_t>&)" << '\n';
            throw;
        }
    }




    //===============================================================================================================

    void Hmat::make_matrix(evolm::matrix<double>& a_matr,
                           std::vector<std::int64_t>& a_ids,
                           evolm::matrix<double>& a_red_matr,
                           std::vector<std::int64_t>& a_red_ids,
                           evolm::matrix<double>& g_matr,
                           std::vector<std::int64_t>& g_ids)
    {
        // Making single-step matriix, H(-1);
        // Assume all input matrices are symmetric and in L-store format
        try
        {
            size_t expected_size = 0.5 * (a_ids.size() - 1) * a_ids.size() + a_ids.size();

            if ( a_matr.size() > expected_size )
                throw std::string("The number of elements in the passed A(-1) matrix is greater then expected!");

            expected_size = 0.5 * (a_red_ids.size() - 1) * a_red_ids.size() + a_red_ids.size();

            if ( a_red_matr.size() > expected_size )
                throw std::string("The number of elements in the passed redused A(-1) matrix is greater then expected!");

            expected_size = 0.5 * (g_ids.size() - 1) * g_ids.size() + g_ids.size();

            if ( g_matr.size() > expected_size )
                throw std::string("The number of elements in the passed G(-1) matrix is greater then expected!");

            if ( g_ids.size() > a_ids.size() )
                throw std::string("The number of elements in the passed G(-1) maatrix is greater then the number of IDs in full A(-1) matrix!");

            if ( !is_value_in_vect(a_red_ids, g_ids) )
                throw std::string("There are IDs in G(-1) matrix which are not part of the reduced A(-1) matrix!");

            if ( !is_value_in_vect(a_ids, g_ids) )
                throw std::string("There are IDs in G(-1) matrix which are not part of the A(-1) matrix!");

            evolm::matrix<size_t> shapeofa;
            shapeofa = a_matr.shape();

            evolm::matrix<size_t> shapeofra;
            shapeofra = a_red_matr.shape();

            evolm::matrix<size_t> shapeofg;
            shapeofg = g_matr.shape();

            if ( ( shapeofa[0] != a_ids.size() ) || ( shapeofa[0] != shapeofa[1] ) )
                throw std::string("The dimension of A(-1) is not correct!");

            if ( ( shapeofra[0] != a_red_ids.size() ) || ( shapeofra[0] != shapeofra[1] ) )
                throw std::string("The dimension of reduced A(-1) is not correct!");

            if ( ( shapeofg[0] != g_ids.size() ) || ( shapeofg[0] != shapeofg[1] ) )
                throw std::string("The dimension of G(-1) is not correct!");

            shapeofa.clear();
            shapeofra.clear();
            shapeofg.clear();

            // --------------------- A22(-1) -------------------------
            evolm::matrix<double> iAgen( g_ids.size() ); // A22(-1)

            get_submatrix( a_red_matr, a_red_ids, g_ids, true);

            evolm::matrix<size_t> shapeofh;
            shapeofh = H.shape();

            if ( ( shapeofh[0] != g_ids.size() ) || ( shapeofh[0] != shapeofh[1] ) )
                throw std::string("The dimension of H == A22(-1) is not correct!");
            
            iAgen = H;

            H.clear();
            H.fclear();

            // ----------------- G22(-1) = G(-1) - A22(-1) -----------------

            //std::cout<<"reading G"<<"\n";
            //g_matr.fread();

            if ( g_matr.size() != iAgen.size() )
                throw  std::string("The number of elements in the G(-1) matrix is not the same as the numbr of elements in the A22(-1) matrix!");

            for (size_t i = 0; i < g_matr.size(); i++)
                g_matr[i] = g_matr[i] - iAgen[i];

            iAgen.clear();
            iAgen.fclear();

            // --------------------- A(-1) + G22(-1) -----------------------
            //std::cout<<"reading A"<<"\n";
            //a_matr.fread();

            H.resize( a_ids.size() );

            H = a_matr;

            a_matr.clear();
            a_matr.fclear();

            for (size_t i = 0; i < g_ids.size(); i++)
            {
                size_t pos_i = find_invect( a_ids, g_ids[i] );

                for (size_t j = 0; j <= i; j++)
                {
                    size_t pos_j = find_invect( a_ids, g_ids[j] );
                    H(pos_i,pos_j) = H(pos_i,pos_j) + g_matr(i,j);
                }
            }

            hmat_id = a_ids;

            g_matr.clear();
            g_matr.fclear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat::make_matrix(evolm::matrix<double>&, std::vector<std::int64_t>&, evolm::matrix<double>&, std::vector<std::int64_t>&, evolm::matrix<double>&, std::vector<std::int64_t>&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::make_matrix(evolm::matrix<double>&, std::vector<std::int64_t>&, evolm::matrix<double>&, std::vector<std::int64_t>&, evolm::matrix<double>&, std::vector<std::int64_t>&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::make_matrix(evolm::matrix<double>&, std::vector<std::int64_t>&, evolm::matrix<double>&, std::vector<std::int64_t>&, evolm::matrix<double>&, std::vector<std::int64_t>&)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Hmat::get_submatrix(std::vector<double>& full_matr, std::vector<std::int64_t>& matr_ids, std::vector<std::int64_t>& selected_ids, bool use_invverse)
    {
        // Extract a sub-matriix corresponding to the 'selected_ids' list from the matrix 'full_matr' with IDs list 'matr_ids';
        // It is expected the 'full_matr' is symmetrixand in L-store fromat
        try
        {
            size_t expected_size = 0.5 * (matr_ids.size() - 1) * matr_ids.size() + matr_ids.size();

            if ( full_matr.size() > expected_size )
                throw std::string("The number of elements in the passed matrix is greater then expected!");

            if ( selected_ids.size() > matr_ids.size() )
                throw std::string("The number of elements in the passed selected IDs array is greater then the number of IDs in the passed matrix!");

            if ( !is_value_in_vect(matr_ids, selected_ids) )
                throw std::string("There are IDs in the selected IDs array which are not part of the passed matrix!");
            
            evolm::matrix<double> hmatr( selected_ids.size() );

            hmatr.from_vector(full_matr);

            get_submatrix(hmatr, matr_ids, selected_ids, use_invverse);

            hmatr.clear();
            hmatr.fclear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat::get_submatrix(std::vector<double>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::get_submatrix(std::vector<double>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::get_submatrix(std::vector<double>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Hmat::get_submatrix(evolm::matrix<double>& full_matr, std::vector<std::int64_t>& matr_ids, std::vector<std::int64_t>& selected_ids, bool use_invverse)
    {
        // Extract a sub-matriix corresponding to the 'selected_ids' list from the matrix 'full_matr' with IDs list 'matr_ids';
        // It is expected the 'full_matr' is symmetrixand in L-store fromat
        try
        {
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

            if ( !is_value_in_vect(matr_ids, selected_ids) )
                throw std::string("There are IDs in the selected IDs array which are not part of the passed matrix!");
            
            // --------------------------------
            if ( !H.empty() )
            {
                H.clear();
                H.fclear();
            }

            if ( !use_invverse )
            {
                H.resize( selected_ids.size() );

                for (size_t i = 0; i < selected_ids.size(); i++)
                {
                    size_t pos_i = find_invect( matr_ids, selected_ids[i] );

                    for (size_t j = 0; j <= i; j++)
                    {
                        size_t pos_j = find_invect( matr_ids, selected_ids[j] );

                        H(i,j) = full_matr( pos_i, pos_j );
                    }
                } 
            }
            else
            {
                // Create the list of IDs which are in matr_ids but not in selected_ids vectors

                std::vector<std::int64_t> not_selected_ids;

                for ( size_t i = 0; i < matr_ids.size(); i++ )
                {
                    int res = find_invect( selected_ids, matr_ids[i] );
                    if ( res == -1 )
                        not_selected_ids.push_back( matr_ids[i] );
                }

                if ( ( selected_ids.size() + not_selected_ids.size() ) != matr_ids.size() )
                    throw std::string("The sum of IDs from from two vectors is not equal to number of IDs in the matrix!");
                
                // ------------------------------------------------

                // Next steps: A22(-1) = A22 - A21 * A11(-1) * A12;

                evolm::matrix<double> A22;
                evolm::matrix<double> A11;
                evolm::matrix<double> A21;
                evolm::matrix<double> A12;

                // -------------------- A11 -----------------------

                A11.resize( not_selected_ids.size() );

                for (size_t i = 0; i < not_selected_ids.size(); i++)
                {
                    size_t pos_i = find_invect( matr_ids, not_selected_ids[i] );

                    for (size_t j = 0; j <= i; j++)
                    {
                        size_t pos_j = find_invect( matr_ids, not_selected_ids[j] );

                        A11(i,j) = full_matr( pos_i, pos_j );
                    }
                }

                A11.symtorec();

                A11.invert();

                A11.fwrite();
                // ------------------------------------------------
                //
                // -------------------- A21 -----------------------
                A21.resize(selected_ids.size(), not_selected_ids.size());

                for (size_t i = 0; i < selected_ids.size(); i++)
                {
                    size_t pos_i = find_invect( matr_ids, selected_ids[i] );

                    for (size_t j = 0; j < not_selected_ids.size(); j++)
                    {
                        size_t pos_j = find_invect( matr_ids, not_selected_ids[j] );

                        A21(i,j) = full_matr( pos_i, pos_j );
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

                A11.clear();
                A11.fclear();
                A21.clear();
                A21.fclear();
                A12.fread();

                res = res * A12;

                A12.clear();
                A12.fclear();
                // ------------------------------------------------
                //
                // -------------------- A22 -----------------------
                A22.resize(selected_ids.size(), selected_ids.size());

                for (size_t i = 0; i < selected_ids.size(); i++)
                {
                    size_t pos_i = find_invect( matr_ids, selected_ids[i] );

                    for (size_t j = 0; j < selected_ids.size(); j++)
                    {
                        size_t pos_j = find_invect( matr_ids, selected_ids[j] );

                        A22(i,j) = full_matr( pos_i, pos_j );
                    }
                }

                // ------------------------------------------------
                //
                // ------------------ A22 - res -------------------
                H = A22 - res;

                H.rectosym();

                A22.clear();
                A22.fclear();
                res.clear();
                res.fclear();
                // ------------------------------------------------
                shapeofh.clear();
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat::get_submatrix(evolm::matrix<double>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::get_submatrix(evolm::matrix<double>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::get_submatrix(evolm::matrix<double>&, std::vector<std::int64_t>&, std::vector<std::int64_t>&, bool)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    int Hmat::find_invect(std::vector<std::int64_t> &where, std::int64_t what)
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
            std::cerr << "Exception in Hmat::find_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::find_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::find_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Hmat::get_matrix(evolm::matrix<double> &arr)
    {
        // Note, we operate with L-stored format, hence return lower triangular part
        try
        {
            arr = H;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat::get_matrix(evolm::matrix<double> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::get_matrix(evolm::matrix<double> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::get_matrix(evolm::matrix<double> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Hmat::get_matrix(std::vector<double> &arr)
    {
        // This method is for Python interfacing;
        // Note, we operate with L-stored format, hence return lower triangular part
        try
        {
            H.to_vector(arr);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat::get_matrix(std::vector<double> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::get_matrix(std::vector<double> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::get_matrix(std::vector<double> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Hmat::get_ids(std::vector<std::int64_t> &ids)
    {
        try
        {
            ids = hmat_id;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat::get_ids(std::vector<std::int64_t>&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::get_ids(std::vector<std::int64_t>&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::get_ids(std::vector<std::int64_t>&)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Hmat::bin_write()
    {
        try
        {
            H.fwrite();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat::bin_write()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::bin_write()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::bin_write()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Hmat::bin_read()
    {
        try
        {
            H.fread();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat::bin_read()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::bin_read()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::bin_read()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Hmat::check_id(std::vector<std::int64_t> &id_list, std::vector<std::int64_t> &checked_id, std::vector<std::int64_t> &missing_id)
    {
        try
        {
            // id_list - vector of ids from where to check
            // checked_id - vector of ids to check
            // missing_id - vector of missing ids (whose not found in id vector)

            for (size_t i = 0; i < checked_id.size(); i++)
                if (!std::binary_search(id_list.begin(), id_list.end(), checked_id[i]))
                    missing_id.push_back(checked_id[i]);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat::check_id(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::check_id(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::check_id(std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    bool Hmat::is_value_in_vect(std::vector<std::int64_t> &where_tocheck, std::vector<std::int64_t> &what_tocheck)
    {
        // are values from what_tocheck in where_tocheck?
        try
        {
            bool out = true;

            std::vector<std::int64_t> missing;
            check_id(where_tocheck, what_tocheck, missing);

            if (!missing.empty())
            {
                out = false;
                missing.clear();
            }

            return out;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat::is_value_in_vect(std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::is_value_in_vect(std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::is_value_in_vect(std::vector<std::int64_t> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

}