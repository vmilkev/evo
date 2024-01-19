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
            H.fclear();
            H.clear();            
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

    void Hmat::make_matrix(std::vector<double>& a_matr,          // A(-1)
                           std::vector<std::int64_t>& a_ids,
                           std::vector<double>& a_red_matr,      // A22(-1)
                           std::vector<std::int64_t>& a_red_ids,
                           std::vector<double>& g_matr,          // G(-1)
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

            amatr.fclear();
            amatr.clear();
            
            ramatr.fclear();
            ramatr.clear();
            
            gmatr.fclear();
            gmatr.clear();            
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

    void Hmat::make_matrix(evolm::matrix<double>& a_matr,     // A(-1)
                           std::vector<std::int64_t>& a_ids,
                           evolm::matrix<double>& a_red_matr, // A22(-1)
                           std::vector<std::int64_t>& a_red_ids,
                           evolm::matrix<double>& g_matr,     // G(-1)
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

            if ( ( shapeofg[0] != shapeofra[0] ) || ( shapeofg[1] != shapeofra[1] ) )
                throw std::string("The dimension of G(-1) is not the same as A22(-1)!");

            shapeofa.clear();
            shapeofra.clear();
            shapeofg.clear();

            // ----------------- G22(-1) = G(-1) - A22(-1) -----------------

            std::vector<size_t> pos_map;
            for (size_t i = 0; i < g_ids.size(); i++)
            {
                pos_map.push_back( find_invect( a_red_ids, g_ids[i] ) );
            }

#pragma omp parallel for
			for (size_t i = 0; i < g_ids.size(); i++)
            {
                size_t row = pos_map[ i ];

				for (size_t j = 0; j <= i; j++)
                {					
					size_t col = pos_map[ j ];

                    if ( row >= col )
                        g_matr(i,j) = g_matr(i,j) - a_red_matr(row,col);
                    else
                        g_matr(i,j) = g_matr(i,j) - a_red_matr(col,row);
				}
			}

            a_red_matr.fclear();
            a_red_matr.clear();

            // --------------------- A(-1) + G22(-1) -----------------------

            pos_map.clear();
            pos_map.shrink_to_fit();

            for (size_t i = 0; i < g_ids.size(); i++)
                pos_map.push_back( find_invect( a_ids, g_ids[i] ) );

#pragma omp parallel for
            for (size_t i = 0; i < g_ids.size(); i++)
            {
                size_t row = pos_map[ i ];

                for (size_t j = 0; j <= i; j++)
                {
                    size_t col = pos_map[ j ];

                    if (row >= col)
                        a_matr(row,col) = a_matr(row,col) + g_matr(i,j);
                    else
                        a_matr(col,row) = a_matr(col,row) + g_matr(i,j);
                }
            }

            a_matr.fwrite();

            H = a_matr;
            hmat_id = a_ids;

            a_matr.clear();

            g_matr.fclear();
            g_matr.clear();
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

    void Hmat::get_matrix(evolm::matrix<double> &arr, std::vector<std::int64_t>& ids, bool keep_ondisk)
    {
        // Note, we operate with L-stored format, hence return lower triangular part
        try
        {
            H.fread();

            arr = H;
            ids = hmat_id;

            if (keep_ondisk)
                H.fwrite();
            else
            {
                H.fclear();
                H.clear();
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat::get_matrix(evolm::matrix<double> &, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::get_matrix(evolm::matrix<double> &, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::get_matrix(evolm::matrix<double> &, std::vector<std::int64_t>&, bool)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Hmat::get_matrix(std::vector<double> &arr, std::vector<std::int64_t>& ids, bool keep_ondisk)
    {
        // This method is for Python interfacing;
        // Note, we operate with L-stored format, hence return lower triangular part
        try
        {
            H.fread();
            
            H.to_vector(arr);
            ids = hmat_id;

            if (keep_ondisk)
                H.fwrite();
            else
            {
                H.fclear();
                H.clear();
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Hmat::get_matrix(std::vector<double> &, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Hmat::get_matrix(std::vector<double> &, std::vector<std::int64_t>&, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Hmat::get_matrix(std::vector<double> &, std::vector<std::int64_t>&, bool)" << '\n';
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