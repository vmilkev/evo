#include "Gmat.hpp"

namespace evoped
{
    //===============================================================================================================

    Gmat::Gmat()
    {
    }

    //===============================================================================================================

    void Gmat::get_matrix(evolm::matrix<double> &arr)
    {
        // Note, we operate with L-stored format, hence return lower triangular part
        try
        {
            arr = G;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::get_matrix(evolm::matrix<double> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::get_matrix(evolm::matrix<double> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::get_matrix(evolm::matrix<double> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::get_matrix(std::vector<double> &arr)
    {
        // This method is for Python interfacing;
        // Note, we operate with L-stored format, hence return lower triangular part
        try
        {
            G.to_vector(arr);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::get_matrix(std::vector<double> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::get_matrix(std::vector<double> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::get_matrix(std::vector<double> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::get_ids(std::vector<std::int64_t> &ids)
    {
        try
        {
            ids = gmatID;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::get_ids(std::vector<std::int64_t>&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::get_ids(std::vector<std::int64_t>&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::get_ids(std::vector<std::int64_t>&)" << '\n';
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

    void Gmat::make_matrix(const std::string &fname, const std::string &fname_ids)
    {
        try
        {
            read_snp(fname, fname_ids);
            make_zmatrix();
            snp_map.clear();
            make_matrix();
            Z.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::make_matrix(const std::string &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::make_matrix(const std::string &, const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::make_matrix(const std::string &, const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::scale_genotypes(const std::string &fname)
    {
        try
        {
            read_snp(fname);
            make_zmatrix();
            G = Z;
            snp_map.clear();
            Z.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::scale_genotypes(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::scale_genotypes(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::scale_genotypes(const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::scale_genotypes(const std::string &fname, const std::string &fname_ids)
    {
        try
        {
            read_snp(fname, fname_ids);
            make_zmatrix();
            G = Z;
            snp_map.clear();
            Z.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::scale_genotypes(const std::string &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::scale_genotypes(const std::string &, const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::scale_genotypes(const std::string &, const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::invert_matrix()
    {
        // Invert matrix as it is despite the conditions of PD or PSD may not be sutisfied;
        // therefore, it is recommended to call the scale_matrix(...) methods first
        try
        {
            G.invert();
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

    void Gmat::invert_matrix(bool full_store)
    {
        // Invert matrix as it is despite the conditions of PD or PSD may not be sutisfied;
        // therefore, it is recommended to call the scale_matrix(...) methods first
        try
        {
            if (full_store)
            {
                G.symtorec();
                G.invert();

                // To be consistent with G-matrix format,
                // transform it back to L-stored fromat

                G.rectosym();
            }
            else
                G.invert();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::invert_matrix(bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::invert_matrix(bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::invert_matrix(bool)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::invert_matrix(std::vector<std::int64_t> &core_id)
    {
        // sparse inverse (APY)

        try
        {
            // Check if all IDs in the core_id vector are in the gmatID vector

            int is_invector;
            for (size_t i = 0; i < core_id.size(); i++)
            {
                is_invector = find_invect(gmatID, core_id[i]);
                if (is_invector == -1)
                {
                    std::cerr << "Core ID: " << core_id[i] << "\n";
                    throw std::string("The core ID is not part of the G-matrix IDs list!");
                }
            }

            // ----------------------------------------------
            //                 Gcc
            // ----------------------------------------------

            size_t r_gcc, c_gcc;
            r_gcc = c_gcc = core_id.size();

            evolm::matrix<double> Gcc(r_gcc, c_gcc);

            // make the list of positions of coreIDs in genotypedIDs
            std::map<std::int64_t, std::int64_t> corePositions;

            find_RecodedIdMap(corePositions, gmatID, core_id);

            // n_threads = std::thread::hardware_concurrency();
            // block_size = static_cast<unsigned int> (coreID.size()/(1*n_threads));

#pragma omp parallel for // schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < core_id.size(); i++)
            {
                size_t r = corePositions[core_id[i]];
                for (size_t j = 0; j <= i; j++)
                {
                    size_t c = corePositions[core_id[j]];
                    Gcc(i, j) = Gcc(j, i) = G(r, c); // for half-store Gf
                }
            }

            corePositions.clear();

            // Inverting Gcc
            Gcc.invert();

            // Write out Gcc to a file
            Gcc.fwrite();

            // ----------------------------------------------
            //                 Gnc
            // ----------------------------------------------

            size_t r_gnc, c_gnc;

            r_gnc = gmatID.size() - core_id.size();
            ;
            c_gnc = c_gcc;

            evolm::matrix<double> Gnc(r_gnc, c_gnc);

            // Make a vector of non-core IDs
            std::vector<int> noncoreID;

            for (size_t i = 0; i < gmatID.size(); i++)
            {
                std::int64_t id = gmatID[i];
                if (!find_invect(core_id, id))
                    noncoreID.push_back(id);
            }

            std::sort(noncoreID.begin(), noncoreID.end());

            // Build Gnc from G matrix

            // make combined vector coreNoncoreIDs
            std::vector<std::int64_t> coreNonCoreIDs;

            coreNonCoreIDs.insert(coreNonCoreIDs.end(), core_id.begin(), core_id.end());

            coreNonCoreIDs.insert(coreNonCoreIDs.end(), noncoreID.begin(), noncoreID.end());

            find_RecodedIdMap(corePositions, gmatID, coreNonCoreIDs);

            // n_threads = std::thread::hardware_concurrency();
            // block_size = static_cast<unsigned int> (noncoreID.size()/(1*n_threads));

#pragma omp parallel for // schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < noncoreID.size(); i++)
            {
                size_t r = corePositions[noncoreID[i]];
                for (size_t j = 0; j < core_id.size(); j++)
                {
                    size_t c = corePositions[core_id[j]];
                    size_t ind;
                    if (noncoreID[i] >= core_id[j])
                        ind = r * (r - 1) / 2 + c - 1;
                    else
                        ind = c * (c - 1) / 2 + r - 1;
                    Gnc[i * core_id.size() + j] = G[ind]; // for half-stored Gf
                }
            }

            coreNonCoreIDs.clear();

            // Write out Gnc to a file
            // Gnc.fwrite();

            // ----------------------------------------------
            //                 Gcn
            // ----------------------------------------------

            // Make Gcn (transpose Gnc to get Gcn) and write it to a file

            evolm::matrix<double> Gcn(c_gnc, r_gnc);

            // n_threads = std::thread::hardware_concurrency();
            // block_size = static_cast<unsigned int> (r_gnc/(1*n_threads));

#pragma omp parallel for // schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < r_gnc; i++)
            {
                for (size_t j = 0; j < c_gnc; j++)
                {
                    Gcn[j * r_gnc + i] = Gnc[i * c_gnc + j];
                }
            }

            Gcn.fwrite();

            // ----------------------------------------------
            //                 Gnn
            // ----------------------------------------------

            // Calculate Gnn and then - invert it

            evolm::matrix<double> GncGcc(r_gnc, c_gnc);

            // restore Gcc from file:
            Gcc.fread();

            GncGcc = Gnc * Gcc;

            Gnc.fwrite();
            Gcc.fwrite();

            // restore Gcn from file:
            Gcn.fread();

            evolm::matrix<double> Gnnfull(r_gnc, r_gnc);

            Gnnfull = GncGcc * Gcn;

            Gcn.fwrite();
            GncGcc.fwrite();

            size_t r_gnn, c_gnn;
            r_gnn = c_gnn = gmatID.size() - core_id.size();

            evolm::matrix<double> Gnn(r_gnn, 1);

            // n_threads = std::thread::hardware_concurrency();
            // block_size = static_cast<unsigned int> (noncoreID.size()/(1*n_threads));

            // Here is  the inversion of Gnn;
            // According to the applied method, Gnn is diagonal

#pragma omp parallel for // schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 1; i <= noncoreID.size(); i++)
            {
                size_t r = corePositions[noncoreID[i - 1]];
                Gnn[i - 1] = 1.0 / (G[r * (r - 1) / 2 + r - 1] - Gnnfull[(i - 1) * r_gnc + (i - 1)]); // half-store case
            }

            Gnnfull.fwrite();
            G.fwrite();
            Gnn.fwrite();

            // ----------------------------------------------
            //                 G12
            //    G12 = - Gcc_i * Gcn * Gnn_i;
            //
            // 1) G12 = Gcc_i * Gcn;
            // 2) G12 = G12 * Gnn_i;
            // 3) G12 = -1 * G12;
            // ----------------------------------------------

            // Restore Gcn and Gcc from the files

            Gcn.fread();
            Gcc.fread();

            // Matrix product, produce Gcci*Gcn
            evolm::matrix<double> G12(r_gcc, r_gnc);

            G12 = Gcc * Gcn;

            Gcc.fwrite();
            Gcn.fwrite();

            size_t r_g12 = r_gcc;
            size_t c_g12 = c_gnn;

            // Matrix product, produce Gcc_cn*Gnn (G12 matrix)

            // n_threads = std::thread::hardware_concurrency();
            // block_size = static_cast<unsigned int> (r_gcc/(1*n_threads));

            Gnn.fread();

#pragma omp parallel for // schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < r_gcc; i++)
            {
                for (size_t j = 0; j < r_gnc; j++)
                {
                    G12[j + i * r_gnc] = G12[j + i * r_gnc] * Gnn[j];
                }
            }

            Gnn.fwrite();

            // ----------------------------------------------
            //                 G11
            // G11 = (I + G12 * Gnc) * Gcc_i.
            // 1) G12_nc = G12 * Gnc.
            // 2) G12_nc = G12_nc + I.
            //							not exists anymore 2) G12_nc_cc = G12_nc * Gcc_i.
            // 3) G11 = G12_nc * Gcc_i.
            // 							not exists anymore 3) G11 = Gcc_i + G12_nc_cc.
            // ----------------------------------------------

            Gnc.fread();

			evolm::matrix<double> G12_nc(r_g12, c_gnc);

            G12_nc = G12 * Gnc;

            Gnc.fwrite();

            // complete making G12 by mult. by -1.
#pragma omp parallel for
			for (size_t i = 0; i < (r_g12 * c_g12); i++)
				G12[i] = G12[i] * (-1.0);

            G12.fwrite();

			// adding identity matrix to G12_nc

			//n_threads = std::thread::hardware_concurrency();
			//block_size = static_cast<unsigned int> (r_g12/(1*n_threads));

#pragma omp parallel for //schedule(static, block_size) num_threads(n_threads)
			for (size_t i = 0; i < r_g12; i++)
				G12_nc[ i + i * c_gnc ] = G12_nc[ i + i * c_gnc ] + 1.0;

			// restore Gcc
            Gcc.fread();

			evolm::matrix<double> G11(r_gcc, c_gcc);

            G11 = G12_nc * Gcc;

            Gcc.fwrite();
            G12_nc.fwrite();

            // ----------------------------------------------
            //                 G(-1)
            // ----------------------------------------------
            evolm::matrix<double> G_11_12;
            evolm::matrix<double> G_21_22;

            G12.fread();
            
            G_11_12 = G11 << G12;
            
            G_11_12.fwrite();
            
            G11.clear();
            G11.fclear();

            G12.transpose();
            Gnn.fread();

            G_21_22 = G12 << Gnn;

            G12.clear();
            G12.fclear();
            Gnn.clear();
            Gnn.fclear();
            G.clear();
            G.fclear();
            Gcc.clear();
            Gcc.fclear();
            Gcn.clear();
            Gcn.fclear();
            Gnc.clear();
            Gnc.fclear();

            G_11_12.fread();
            G.resize(gmatID.size(), gmatID.size());

            G = G_11_12 >> G_21_22;

            G_11_12.clear();
            G_11_12.fclear();
            G_21_22.clear();
            G_21_22.fclear();

            G.rectosym();

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

    void Gmat::scale_matrix(double scale_coef)
    {
        // Scale diagonal elements of G matrix by the scale_coef, then invert it
        try
        {
            evolm::matrix<size_t> shapeofg;
            shapeofg = G.shape();

            for (size_t i = 0; i < shapeofg[0]; i++)
                G(i, i) = G(i, i) * scale_coef;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::scale_matrix(double)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::scale_matrix(double)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::scale_matrix(double)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::scale_matrix(std::vector<double> &scale_matr, double scaling_weight)
    {
        // Use the matrix scale_matr to scale G matrix, then invert it;
        // This method is just to fit to the Python interface;
        // Note, scale_matr is symetric and consists of only L-part (L-store format)
        try
        {
            // It is required the scaling matrix is the same dimension as the inverting G

            evolm::matrix<size_t> shapeofg;
            shapeofg = G.shape();
            if (shapeofg[0] != shapeofg[1])
                throw std::string("G matrix has wrong dimension: number of raws is not the same as number of columns!");

            evolm::matrix<double> amat; // assumed it is in L-stored format
            amat.resize(shapeofg[0]);

            if (amat.size() != scale_matr.size())
                throw std::string("The passed scaleing matrix has wrong dimension: it is not the same as the dimension of inverting G matrix!");

            // build the scalling matrix from the passed vector and call inversion method:
            amat.from_vector(scale_matr);

            scale_matrix(amat, scaling_weight);

            amat.clear();
            shapeofg.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::scale_matrix(std::vector<double>&, double)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::scale_matrix(std::vector<double>&, double)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::scale_matrix(std::vector<double>&, double)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Gmat::scale_matrix(evolm::matrix<double> &scale_matr, double scaling_weight)
    {
        // Direct inverse of G with prior scaling by the matrix 'scale_matr':
        // G = (1-w)*Ga + wA; Ga = beta*G + alpha;
        // mean( diag(G) )*beta + alpha = mean( diag(A) );
        // mean( G )*beta + alpha = mean( A )

        try
        {
            evolm::matrix<size_t> shapeofa;
            shapeofa = scale_matr.shape();

            evolm::matrix<size_t> shapeofg;
            shapeofg = G.shape();

            if ( shapeofa[0] != shapeofg[0] )
                throw std::string("The passed scaleing matrix has wrong dimension: the number of rows is not the same as in the inverting G matrix!");

            if ( shapeofa[1] != shapeofg[1] )
                throw std::string("The passed scaleing matrix has wrong dimension: the number of columns is not the same as in the inverting G matrix!");

            // Get mean values of the scaler matrix:

            double a_ofd_mean = 0.0;
            double a_all_mean = 0.0;
            double a_diag_mean = 0.0;

#pragma omp parallel for reduction(+ : a_diag_mean)
            for (size_t i = 0; i < shapeofa[0]; i++)
            {
                a_diag_mean += scale_matr(i, i);
            }

            a_all_mean = a_diag_mean;
            a_diag_mean = a_diag_mean / ((double)shapeofa[0]);

#pragma omp parallel for reduction(+ : a_ofd_mean)
            for (size_t i = 0; i < shapeofa[0]; i++)
            {
                for (size_t j = 0; j < i; j++)
                {
                    a_ofd_mean += scale_matr(i, j);
                }
            }

            a_all_mean = (a_all_mean + 2 * a_ofd_mean) / ((double)scale_matr.size());
            a_ofd_mean = 2 * a_ofd_mean / ((double)scale_matr.size() - (double)shapeofa[0]);

            // Get mean values of G matrix

            double g_ofd_mean = 0.0;
            double g_all_mean = 0.0;
            double g_diag_mean = 0.0;

            double alpha = 0.0;
            double betha = 1.0;

#pragma omp parallel for reduction(+ : g_diag_mean)
            for (size_t i = 0; i < shapeofg[0]; i++)
            {
                g_diag_mean += G(i, i);
            }

#pragma omp parallel for reduction(+ : g_ofd_mean)
            for (size_t i = 0; i < shapeofg[0]; i++)
            {
                for (size_t j = 0; j < i; j++)
                {
                    g_ofd_mean += G(i, j);
                }
            }

            g_all_mean = (g_diag_mean + 2 * g_ofd_mean) / ((double)G.size());
            g_diag_mean = g_diag_mean / ((double)shapeofg[0]);

            betha = (a_all_mean - a_diag_mean) / (g_all_mean - g_diag_mean);
            alpha = a_diag_mean - g_diag_mean * betha;

            g_ofd_mean = 2 * g_ofd_mean / ((double)G.size() - (double)shapeofg[0]);

            shapeofa.clear();
            shapeofg.clear();

#pragma omp parallel for
            for (size_t i = 0; i < G.size(); i++)
            {
                G[i] = (G[i] * betha + alpha) * (1 - scaling_weight) + scale_matr[i] * scaling_weight;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::scale_matrix(evolm::matrix<double>&, double)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::scale_matrix(evolm::matrix<double>&, double)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::scale_matrix(evolm::matrix<double>&, double)" << '\n';
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
            G.fclear();
            Z.clear();
            Z.fclear();
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
        try
        {
            std::string line;

            char *end;
            const char *p;
            std::vector<T> tmp_list;
            size_t diagonals = 0;

            std::ifstream ped;
            ped.open(gmat_file, std::fstream::in);

            if (!ped.good())
                throw std::string("Cannot open G-matrix file!");

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
            }

            ped.close();

            if (!is_unique(gmatID))
                gmatID.erase(unique(gmatID.begin(), gmatID.end()), gmatID.end()); // here the vector should be sorted and unique

            if (diagonals != gmatID.size())
                throw std::string("There are missing diagonals in G-matrix file.");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::read_matrix(const std::string &, std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<T> &)" << '\n';
            std::cerr << e.what() << '\n';
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
        //       The function should accept a full store format as well, though will write
        //       only to the lower triangular part assuming the G matrix is always symmetric.
        try
        {
            std::vector<std::int64_t> g_row;
            std::vector<std::int64_t> g_col;
            std::vector<double> g_val;
            std::map<std::int64_t, std::int64_t> rid_map;

            read_matrix(gmat_file, g_row, g_col, g_val);

            if (gmatID.empty())
                throw std::string("Genotyped IDs vector is empty!");

            get_RecodedIdMap(rid_map, gmatID); // here indexing starts from 1 (but not from 0) !

            if (rid_map.empty())
                throw std::string("Recoded IDs map is empty!");

            // auto n_threads = std::thread::hardware_concurrency();
            // block_size = static_cast<unsigned int> (g_row.size()/(n_threads));

            G.resize(gmatID.size());

#pragma omp parallel for // schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < g_row.size(); i++) // here indexing fro r & c starts from 1 (but not from 0), because of rid_map coding !
            {
                size_t r = rid_map[g_row[i]];
                size_t c = rid_map[g_col[i]];
                size_t ind = r * (r - 1) / 2 + c - 1;
                if (c > r)
                    ind = c * (c - 1) / 2 + r - 1;
                G[ind] = g_val[i];
            }

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
        /*
                Reads file format:
                [observation ID] [devise code] [list of SNPs with " " delimiter]

                Example:
                18 1000 2 0 1 1 0 0 0 2 1 2
                19 1000 5 0 0 0 0 2 0 2 1 0
                20 1000 1 5 2 1 1 0 0 2 1 2
                21 1000 0 0 2 1 0 1 0 2 2 1
        */
        try
        {
            std::string line;
            std::vector<std::string> data_list;

            std::ifstream snpF;
            snpF.open(snp_file, std::fstream::in);

            if (!snpF.good())
                throw std::string("Cannot open SNPs file!");

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
                // get the last element of the string
                data_list.push_back(line);

                // now we have got the SNP data for one ID
                snp_map[stoi(data_list[0])] = data_list[2];

                gmatID.push_back(stoi(data_list[0]));

                data_list.erase(data_list.begin(), data_list.end());
            }

            snpF.close();

            // build the map: <index, ID>
            size_t tmpInd = 0;
            for (const auto &snp : snp_map)
            {
                anim_id_map[tmpInd] = snp.first;
                tmpInd++;
            }

            // Check if gmat IDs are unique
            std::vector<std::int64_t> t_gmatID(gmatID);

            if (!is_unique(t_gmatID))
                throw std::string("Thhere are repeated IDs in the processed SNPs file!");

            t_gmatID.clear();
            t_gmatID.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::read_snp(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
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

    void Gmat::read_snp(const std::string &snp_file, const std::string &ids_file)
    {
        /*
                Reads SNPs file of the format:
                [list of SNPs with " " delimiter]

                Example:
                2 0 1 1 0 0 0 2 1 2
                5 0 0 0 0 2 0 2 1 0
                1 5 2 1 1 0 0 2 1 2
                0 0 2 1 0 1 0 2 2 1

                Reads IDs file of the format:
                [vector of IDs of type integer]

                Example:
                101
                543
                20987
                345
        */
        try
        {
            std::string line;
            std::vector<std::string> data_list;
            std::ifstream snpF;

            // (1) reading the IDs first:

            snpF.open(ids_file, std::fstream::in);

            if (!snpF.good())
                throw std::string("Cannot open G-matrix IDs file!");

            if (!gmatID.empty())
            {
                gmatID.clear();
                gmatID.shrink_to_fit();
            }

            std::int64_t id;
            while (snpF >> id)
                gmatID.push_back(id);

            snpF.close();

            // Check if gmat IDs are unique
            std::vector<std::int64_t> t_gmatID(gmatID);

            if (!is_unique(t_gmatID))
                throw std::string("Thhere are repeated IDs in the processed SNPs file!");

            t_gmatID.clear();
            t_gmatID.shrink_to_fit();

            // (2) reading the SNPs:

            snpF.open(snp_file, std::fstream::in);

            if (!snpF.good())
                throw std::string("Cannot open SNPs file!");

            while (getline(snpF, line))
                data_list.push_back(line);

            snpF.close();

            if (gmatID.size() != data_list.size())
                throw std::string("The number of IDs in the G-matrix IDs file is not the same as the number of genotypes in the SNPs file!");

            if (!snp_map.empty())
                snp_map.clear();

            for (size_t i = 0; i < gmatID.size(); i++)
                snp_map[gmatID[i]] = data_list[i];

            data_list.clear();
            data_list.shrink_to_fit();

            if (!anim_id_map.empty())
                anim_id_map.clear();

            // build the map: <index, ID>
            size_t tmpInd = 0;
            for (const auto &snp : snp_map)
            {
                anim_id_map[tmpInd] = snp.first;
                tmpInd++;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat::read_snp(const std::string &, const std::string&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat::read_snp(const std::string &, const std::string&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat::read_snp(const std::string &, const std::string&)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    int Gmat::find_invect(std::vector<std::int64_t> &where, std::int64_t what)
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
            // get the number of SNPs
            std::vector<int> tmpVect;

            auto it = snp_map.begin();
            std::string tmpStr = it->second;

            parse_string(tmpStr, tmpVect);

            size_t snpNum = tmpVect.size();

            tmpVect.clear();
            tmpVect.shrink_to_fit();
            tmpStr.clear();
            tmpStr.shrink_to_fit();

            // declare the matrix M
            evolm::matrix<double> M(snp_map.size(), snpNum);

            /* vector of SNPs frequences and missed values */
            std::vector<double> P(snpNum, 0.0);
            std::vector<int> missed(snpNum, 0);
            std::vector<double> missed2pq(snp_map.size(), 0.0);

            // map of missed values locations
            std::vector<std::vector<int>> missedLocation;
            for (size_t i = 0; i < snpNum; i++)
                missedLocation.push_back(std::vector<int>());

            // parse SNPs and fill matrix M
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

            // finish to calculate allele frequences, additionally accounting missing values

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

            // modify Z matrix, so instead of missing values we put population average (0.0)

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

            // After centering, adjust for missing markers for each animal;
            // adjust for sqrt[sum of 2pq over all loci /sum of 2pq over non-missing loci.

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
            G = (Z ^ 2) * (1 / freq);

            // Because G is symmetric,
            // make it L-stored to save some memory

            G.rectosym();
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
    // ??? not used so far !
    template <typename T>
    void Gmat::get_gvalues(std::vector<std::int64_t> &row, std::vector<std::int64_t> &col, std::vector<T> &val, double diag_val)
    {
        try
        {
            for (size_t i = 0; i < anim_id_map.size(); i++)
            {
                for (size_t j = 0; j <= i; j++)
                {
                    col.push_back(anim_id_map[j]);
                    row.push_back(anim_id_map[i]);

                    if (i == j)
                        val.push_back(static_cast<T>(G(i, j) + diag_val));
                    else
                        val.push_back(static_cast<T>(G(i, j)));
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

    template void Gmat::get_gvalues(std::vector<std::int64_t> &row, std::vector<std::int64_t> &col, std::vector<float> &val, double diag_val);
    template void Gmat::get_gvalues(std::vector<std::int64_t> &row, std::vector<std::int64_t> &col, std::vector<double> &val, double diag_val);

    //===============================================================================================================

    /*void Gmat::get_gids(std::vector<std::int64_t> &ids)
    {
        try
        {
            if (gmatID.empty())
            {
                for (const auto &id : anim_id_map)
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
    }*/

    //===============================================================================================================
}