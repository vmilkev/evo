/**
 * @file Gmat.cpp
 * @author Viktor Milkevych
 * @brief 
 * @version 0.1
 * @date 2024-05-14
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#include "Gmat.hpp"

namespace evoped
{
    //===============================================================================================================
    /**
     * @brief Construct a new Gmat< T>:: Gmat object
     * 
     * @tparam T 
     */
    template <typename T> Gmat<T>::
    Gmat()
    { }
    template Gmat<float>::Gmat();
    template Gmat<double>::Gmat();
    //===============================================================================================================
    /**
     * @brief I/O interface for accessing the internal container storing results.
     *          Note, we operate with L-stored format, hence return lower triangular part
     * 
     * @param arr empty dense matrix class object where the internal
     *            storage conttainer (holding result) will be copied
     * @param ids empty std vector where samples (individuals) IDs
     *            corresponding to the data in arr will be copied
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    get_matrix(evolm::matrix<T> &arr, std::vector<std::int64_t> &ids)
    {        
        try
        {
            //if ( !G.is_ondisk() )
            //    G.fwrite();
            //if ( G.is_ondisk() )
            //    G.fread();

            arr = G;

            if (!ids.empty())
                ids.clear();
            
            ids = gmatID;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::get_matrix(evolm::matrix<T> &, std::vector<std::int64_t>&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::get_matrix(evolm::matrix<T> &, std::vector<std::int64_t>&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::get_matrix(evolm::matrix<T> &, std::vector<std::int64_t>&)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::get_matrix(evolm::matrix<float> &arr, std::vector<std::int64_t> &ids);
    template void Gmat<double>::get_matrix(evolm::matrix<double> &arr, std::vector<std::int64_t> &ids);
    //===============================================================================================================
    /**
     * @brief I/O interface for saving on disk the internal matrix container holding results.
     *          Note, we operate with L-stored format, hence return lower triangular part
     * 
     * @param arr file name for saving a dense matrix class object where the internal
     *            storage conttainer (holding result) will be copied
     * @param ids file name for saving std vector where samples (individuals) IDs
     *            corresponding to the data in arr will be copied
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    save_matrix(const std::string &arr, const std::string &ids)
    {        
        try
        {
            Utilities2 u;

            //if ( G.is_ondisk() )
            //    G.fread();

            G.fwrite(arr);
            //G.fwrite();

            u.vect_to_binary(gmatID, ids);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::save_matrix(const std::string &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::save_matrix(const std::string &, const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::save_matrix(const std::string &, const std::string &)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::save_matrix(const std::string &arr, const std::string &ids);
    template void Gmat<double>::save_matrix(const std::string &arr, const std::string &ids);
    //===============================================================================================================
    /**
     * @brief I/O interface for saving on disk the internal matrix container holding results.
     *          Note, we operate with L-stored format, hence return lower triangular part
     * 
     * @param arr file name for saving a dense matrix class object where the internal
     *            storage conttainer (holding result) will be copied
     *
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    save_matrix(const std::string &arr)
    {        
        try
        {
            //if ( G.is_ondisk() )
            //    G.fread();

            G.fwrite(arr);
            //G.fwrite();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::save_matrix(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::save_matrix(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::save_matrix(const std::string &)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::save_matrix(const std::string &arr);
    template void Gmat<double>::save_matrix(const std::string &arr);
    //===============================================================================================================
    /**
     * @brief I/O interface for saving on disk the list of IDs of the internal matrix container.
     * 
     * @param ids file name for saving std vector where samples (individuals) IDs
     *            corresponding to the data in arr will be copied
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    save_ids(const std::string &ids)
    {        
        try
        {
            Utilities2 u;

            u.vect_to_binary(gmatID, ids);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::save_ids(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::save_ids(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::save_ids(const std::string &)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::save_ids(const std::string &ids);
    template void Gmat<double>::save_ids(const std::string &ids);

    //===============================================================================================================
    /**
     * @brief I/O interface for accessing the internal container storing results.
     *          Note, we operate with L-stored format, hence return lower triangular part
     * 
     * @param arr empty dense matrix class object where the internal
     *            storage conttainer (holding result) will be copied
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    get_matrix(evolm::matrix<T> &arr)
    {        
        try
        {
            //if ( !G.is_ondisk() )
            //    G.fwrite();

            arr = G;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::get_matrix(evolm::matrix<T> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::get_matrix(evolm::matrix<T> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::get_matrix(evolm::matrix<T> &)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::get_matrix(evolm::matrix<float> &arr);
    template void Gmat<double>::get_matrix(evolm::matrix<double> &arr);
    //===============================================================================================================
    /**
     * @brief I/O interface for accessing the internal container storing results.
     *          Note, we operate with L-stored format, hence return lower triangular part
     * 
     * @param ids empty std vector where samples (individuals) IDs
     *            corresponding to the data in arr will be copied
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    get_ids(std::vector<std::int64_t> &ids)
    {        
        try
        {
             if (!ids.empty())
                ids.clear();
            
            ids = gmatID;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::get_ids(std::vector<std::int64_t>&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::get_ids(std::vector<std::int64_t>&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::get_matrix(std::vector<std::int64_t>&)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::get_ids(std::vector<std::int64_t> &ids);
    template void Gmat<double>::get_ids(std::vector<std::int64_t> &ids);
    //===============================================================================================================
    /**
     * @brief I/O interface for accessing the internal container storing the results.
     *          Note, we operate with L-stored format, hence return lower triangular part
     * 
     * @param arr empty std vector object where the internal storage conttainer (holding result) will be copied
     * @param ids empty std vector where samples (individuals) IDs corresponding to the data in arr will be copied
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    get_matrix(const std::string &out_fname)
    {
        try
        {
            std::vector<T> values;
            std::vector<size_t> keys;

            G.to_vector(values);

            for (size_t i = 0; i < gmatID.size(); i++)
                for (size_t j = 0; j <= i; j++)
                    keys.push_back(i*(i+1)/2 + j);

            Utilities2 u;
            u.fwrite_matrix(out_fname, values, keys, gmatID);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::get_matrix(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::get_matrix(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::get_matrix(const std::string &)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::get_matrix(const std::string &out_fname);
    template void Gmat<double>::get_matrix(const std::string &out_fname);
    //===============================================================================================================
    /**
     * @brief Constructs G matrix by reading a text file consisting of samples and snp variants.
     * 
     * @tparam T defines type, float or double
     * @param fname text file name where the first two cols are variant ID and chip ID, the rest is SNP variants
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    make_matrix(const std::string &fname)
    {
        try
        {
            if ( is_plink_file(fname) ) // the pipeline for binary (.bad) plink-formated data
            {
std::cout<<"Pass 1"<<"\n";
                evolm::matrix<int> M;
                get_m_matrix(fname, M);
std::cout<<"Pass 2"<<"\n";
                make_zmatrix(M); // scalling SNPs
std::cout<<"Pass 3"<<"\n";
                M.clear();
            }
            else // the pipeline for text (.ped) plink-formated data !!! not implemented yet
            {
                read_snp(fname); // reads SNPs with variant IDs from ASCII fiele and output to the snp_map
                make_zmatrix(); // scalling SNPs
                snp_map.clear();
            }
std::cout<<"Pass 4"<<"\n";
            make_matrix(); // making G matrix
std::cout<<"Pass 5"<<"\n";
            Z.fclear();
            Z.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::make_matrix(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::make_matrix(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::make_matrix(const std::string &)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::make_matrix(const std::string &fname);
    template void Gmat<double>::make_matrix(const std::string &fname);
    //===============================================================================================================
    /**
     * @brief Constructs G matrix by reading two text files,
     *        one consisting of snp variants, and other consisting of samples IDs.
     * 
     * @tparam T defines type, float or double
     * @param fname the text file with snp variants data
     * @param fname_ids the text file with samples IDs
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    make_matrix(const std::string &fname, const std::string &fname_ids)
    {
        try
        {
            read_snp(fname, fname_ids); // reads SNPs and its IDs from ASCII fieles (one for SNPs, another for IDs)
            make_zmatrix(); // scalling SNPs
            snp_map.clear();
            make_matrix(); // making G matrix
            Z.fclear();
            Z.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::make_matrix(const std::string &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::make_matrix(const std::string &, const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::make_matrix(const std::string &, const std::string &)" << '\n';
            throw;
        }
    }

    template void Gmat<float>::make_matrix(const std::string &fname, const std::string &fname_ids);
    template void Gmat<double>::make_matrix(const std::string &fname, const std::string &fname_ids);
    //===============================================================================================================
    /**
     * @brief Check if provided file is a family of plink files (.bed, .bim, .fam).
     * 
     * @tparam T defines type, float or double
     * @param fname the data file name
     * 
     * @returns false if the data file is not a plink family, true otherwise
     * 
     */
    template <typename T> bool Gmat<T>::
    is_plink_file(const std::string &fname)
    {
        try
        {
            struct pio_file_t plink_file;

            if( pio_open( &plink_file, fname.c_str() ) == PIO_OK )
            {
                pio_close( &plink_file );
                return true;
            }
            else
                return false;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::is_plink_file(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::is_plink_file(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::is_plink_file(const std::string &)" << '\n';
            throw;
        }
    }
    template bool Gmat<float>::is_plink_file(const std::string &fname);
    template bool Gmat<double>::is_plink_file(const std::string &fname);
    //===============================================================================================================
    /**
     * @brief Extract variants and samples from plink files (.bed, .bim, .fam).
     * 
     * @tparam T defines type, float or double
     * @param fname the data file name
     * @param M empty object of dense matrix class where the resulting vriants will be copied
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    get_m_matrix(const std::string &fname, evolm::matrix<int> &M)
    {
        try
        {
            struct pio_file_t plink_file;
            snp_t *snp_buffer;
            size_t sample_id;

            if( pio_open( &plink_file, fname.c_str() ) != PIO_OK )
                throw std::string( "Error: Could not open plink file!" );

            if( !pio_one_locus_per_row( &plink_file ) )
                throw std::string( "This script requires that snps are rows and samples columns." );

            // get number of samples and SNPs:
            size_t n_samples = pio_num_samples(&plink_file);
            size_t n_variants = pio_num_loci(&plink_file);

            M.resize(n_variants, n_samples);

            size_t locus_id = 0;

            snp_buffer = (snp_t *) malloc( pio_row_size( &plink_file ) );

            while( pio_next_row( &plink_file, snp_buffer ) == PIO_OK )
            {
                for( sample_id = 0; sample_id < n_samples; sample_id++)
                {
                    struct pio_sample_t *sample = pio_get_sample( &plink_file, sample_id );
                    M(locus_id, sample_id) = (int)snp_buffer[ sample_id ];
                }
                locus_id++;
            }

            // we need to extract and recode samples IDs, so we ge back to first row and read it onece
            pio_reset_row(&plink_file);
            if( pio_next_row( &plink_file, snp_buffer ) == PIO_OK )
            {
                for( sample_id = 0; sample_id < n_samples; sample_id++)
                {
                    struct pio_sample_t *sample = pio_get_sample( &plink_file, sample_id );
                    samples_id_map[sample_id+1] = sample->iid;
                    //gmatID.push_back(sample_id+1);
                    gmatID.push_back( std::stol(sample->iid) );
                }
            }

            gmatID.shrink_to_fit();

            free( snp_buffer );
            pio_close( &plink_file );

            M.transpose();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::get_m_matrix(const std::string &, evolm::matrix<int> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::get_m_matrix(const std::string &, evolm::matrix<int> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::get_m_matrix(const std::string &, evolm::matrix<int> &)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::get_m_matrix(const std::string &fname, evolm::matrix<int> &M);
    template void Gmat<double>::get_m_matrix(const std::string &fname, evolm::matrix<int> &M);
    //===============================================================================================================
    /**
     * @brief Scales SNPs from the text file.
     * 
     * @tparam T defines type, float or double
     * @param fname the text file name consisting the snp variants with samples information (IDs)
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    scale_genotypes(const std::string &fname)
    {
        try
        {
            if ( is_plink_file(fname) ) // the pipeline for binary (.bad) plink-formated data
            {
                evolm::matrix<int> M;
                get_m_matrix(fname, M);
                make_zmatrix(M); // scalling SNPs
                M.clear();
                G = Z; // copy to the main container
                Z.fclear();
                Z.clear();
            }
            else
            {
                read_snp(fname); // reads SNPs with IDs from ASCII fiele
                make_zmatrix(); // scalling
                //Z.fwrite(); // move data to a binary file
                G = Z; // copy to the main container
                snp_map.clear();
                Z.fclear();
                Z.clear();
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::scale_genotypes(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::scale_genotypes(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::scale_genotypes(const std::string &)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::scale_genotypes(const std::string &fname);
    template void Gmat<double>::scale_genotypes(const std::string &fname);
    //===============================================================================================================
    /**
     * @brief Scales SNPs from the text file.
     * 
     * @tparam T defines type, float or double
     * @param fname the text file name consisting the snp variants without samples information (IDs)
     * @param fname_ids the text file name with samples IDs
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    scale_genotypes(const std::string &fname, const std::string &fname_ids)
    {
        try
        {
            read_snp(fname, fname_ids); // reads SNPs and its IDs from ASCII fieles (one for SNPs, another for IDs)
            make_zmatrix(); // scalling
            //Z.fwrite(); // move data to a binary file
            G = Z; // copy to the main container
            snp_map.clear();
            Z.fclear();
            Z.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::scale_genotypes(const std::string &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::scale_genotypes(const std::string &, const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::scale_genotypes(const std::string &, const std::string &)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::scale_genotypes(const std::string &fname, const std::string &fname_ids);
    template void Gmat<double>::scale_genotypes(const std::string &fname, const std::string &fname_ids);
    //===============================================================================================================
    /**
     * @brief Inverting G matrix regardless of whether it is in compact or full storage format.
     *          Invert matrix as it is despite the conditions of PD or PSD may not be satisfied;
     *          therefore, it is recommended to call the scale_matrix(...) methods first.
     * 
     * @tparam T defines type, float or double
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    invert_matrix()
    {
        try
        {
            G.invert();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::invert_matrix()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::invert_matrix()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::invert_matrix()" << '\n';
            throw;
        }
    }
    template void Gmat<float>::invert_matrix();
    template void Gmat<double>::invert_matrix();
    //===============================================================================================================
    /**
     * @brief Invert G matrix as in full storage format.
     *          Invert matrix as it is despite the conditions of
     *          PD or PSD may not be sutisfied;therefore,
     *          it is recommended to call the scale_matrix(...) methods first.
     * 
     * @tparam T defines type, float or double
     * @param full_store boolean parameter indicating a lower triangular or
     *                   full store format inversion should be provided
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    invert_matrix(bool full_store)
    {
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
            std::cerr << "Exception in Gmat<T>::invert_matrix(bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::invert_matrix(bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::invert_matrix(bool)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::invert_matrix(bool full_store);
    template void Gmat<double>::invert_matrix(bool full_store);
    //===============================================================================================================
    /**
     * @brief Invert G matrix through its sparse approximation, the APY method
     * 
     * @tparam T defines type, float or double
     * @param core_id std vector of samples IDs which are core-id according to the APY method
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    invert_matrix(std::vector<std::int64_t> &core_id)
    {
        try
        {
            Utilities2 u;

            if (G.empty())
                throw std::string("There is no G matrix which needs to be inverted!");

            if (gmatID.empty())
                throw std::string("The vector of G matrix IDs is empty!");

            // Check if all IDs in the core_id vector are in the gmatID vector

            int is_invector;
            for (size_t i = 0; i < core_id.size(); i++)
            {
                is_invector = u.find_invect(gmatID, core_id[i]);
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

            evolm::matrix<T> Gcc(r_gcc, c_gcc);

            // make the list of positions of coreIDs in genotypedIDs
            std::map<std::int64_t, std::int64_t> corePositions;

            u.find_RecodedIdMap(corePositions, gmatID, core_id);

            // n_threads = std::thread::hardware_concurrency();
            // block_size = static_cast<unsigned int> (coreID.size()/(1*n_threads));

#pragma omp parallel for // schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < core_id.size(); i++)
            {
                size_t r = corePositions[core_id[i]];
                for (size_t j = 0; j <= i; j++)
                {
                    size_t c = corePositions[core_id[j]];
                    Gcc(i, j) = Gcc(j, i) = G(r - 1, c - 1); // for half-store Gf
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
            c_gnc = c_gcc;

            evolm::matrix<T> Gnc(r_gnc, c_gnc);

            // Make a vector of non-core IDs
            std::vector<int> noncoreID;

            for (size_t i = 0; i < gmatID.size(); i++)
            {
                std::int64_t id = gmatID[i];
                if (u.find_invect(core_id, id) == -1)
                    noncoreID.push_back(id);
            }

            std::sort(noncoreID.begin(), noncoreID.end());

            // Build Gnc from G matrix

            // make combined vector coreNoncoreIDs
            std::vector<std::int64_t> coreNonCoreIDs;

            coreNonCoreIDs.insert(coreNonCoreIDs.end(), core_id.begin(), core_id.end());

            coreNonCoreIDs.insert(coreNonCoreIDs.end(), noncoreID.begin(), noncoreID.end());

            u.find_RecodedIdMap(corePositions, gmatID, coreNonCoreIDs);

            // n_threads = std::thread::hardware_concurrency();
            // block_size = static_cast<unsigned int> (noncoreID.size()/(1*n_threads));

#pragma omp parallel for // schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < noncoreID.size(); i++)
            {
                size_t r = corePositions[noncoreID[i]];

                for (size_t j = 0; j < core_id.size(); j++)
                {
                    size_t c = corePositions[core_id[j]];

                    if (r >= c)
                        Gnc(i, j) = G(r - 1, c - 1);
                    else
                        Gnc(i, j) = G(c - 1, r - 1);
                    /*size_t ind;
                     if (noncoreID[i] >= core_id[j])
                         ind = r * (r - 1) / 2 + c - 1;
                     else
                         ind = c * (c - 1) / 2 + r - 1;

                     Gnc[i * core_id.size() + j] = G[ind]; // for half-stored Gf*/
                }
            }

            // Write out Gnc to a file
            // Gnc.fwrite();

            // ----------------------------------------------
            //                 Gcn
            // ----------------------------------------------

            // Make Gcn (transpose Gnc to get Gcn) and write it to a file

            evolm::matrix<T> Gcn(c_gnc, r_gnc);

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

            evolm::matrix<T> GncGcc(r_gnc, c_gnc);

            // restore Gcc from file:
            Gcc.fread();

            GncGcc = Gnc * Gcc;

            Gnc.fwrite();
            Gcc.fwrite();

            // restore Gcn from file:
            Gcn.fread();

            evolm::matrix<T> Gnnfull(r_gnc, r_gnc);

            Gnnfull = GncGcc * Gcn;

            Gcn.fwrite();

            // GncGcc.fwrite();
            GncGcc.fclear();
            GncGcc.clear();

            size_t r_gnn, c_gnn;
            r_gnn = c_gnn = gmatID.size() - core_id.size();

            evolm::matrix<T> Gnn(r_gnn, 1);

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

            // Gnnfull.fwrite();
            Gnnfull.fclear();
            Gnnfull.clear();

            // G.fwrite();
            G.fclear();
            G.clear();

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
            evolm::matrix<T> G12(r_gcc, r_gnc);

            G12 = Gcc * Gcn;

            Gcc.fwrite();
            // Gcn.fwrite();
            Gcn.fclear();
            Gcn.clear();

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

            evolm::matrix<T> G12_nc(r_g12, c_gnc);

            G12_nc = G12 * Gnc;

            // Gnc.fwrite();
            Gnc.fclear();
            Gnc.clear();

            // complete making G12 by mult. by -1.
#pragma omp parallel for
            for (size_t i = 0; i < (r_g12 * c_g12); i++)
                G12[i] = G12[i] * (-1.0);

            G12.fwrite();

            // adding identity matrix to G12_nc

            // n_threads = std::thread::hardware_concurrency();
            // block_size = static_cast<unsigned int> (r_g12/(1*n_threads));

#pragma omp parallel for // schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < r_g12; i++)
                G12_nc[i + i * c_gnc] = G12_nc[i + i * c_gnc] + 1.0;

            // restore Gcc
            Gcc.fread();

            evolm::matrix<T> G11(r_gcc, c_gcc);

            G11 = G12_nc * Gcc;

            // Gcc.fwrite();
            // G12_nc.fwrite();
            Gcc.fclear();
            Gcc.clear();
            G12_nc.fclear();
            G12_nc.clear();

            // ----------------------------------------------
            //                 G(-1)
            // ----------------------------------------------
            evolm::matrix<T> G_11_12;
            evolm::matrix<T> G_21_22;
            evolm::matrix<T> Gnn_full;

            G12.fread();

            G_11_12 = G11 << G12;

            G_11_12.fwrite();

            G11.fclear();
            G11.clear();

            G12.transpose();
            Gnn.fread();

            Gnn_full.resize(Gnn.size(), Gnn.size());

            for (size_t i = 0; i < Gnn.size(); i++)
                Gnn_full(i, i) = Gnn[i];

            G_21_22 = G12 << Gnn_full;

            G12.fclear();
            G12.clear();
            Gnn.fclear();
            Gnn.clear();
            Gnn_full.fclear();
            Gnn_full.clear();

            G_11_12.fread();
            G.resize(gmatID.size(), gmatID.size());

            G = G_11_12 >> G_21_22;

            G_11_12.fclear();
            G_11_12.clear();
            G_21_22.fclear();
            G_21_22.clear();

            G.rectosym();

            //G.fwrite();

            gmatID.clear();
            gmatID = coreNonCoreIDs;

            coreNonCoreIDs.clear();
            coreNonCoreIDs.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::invert_matrix(std::vector<std::int64_t>&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::invert_matrix(std::vector<std::int64_t>&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::invert_matrix(std::vector<std::int64_t>&)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::invert_matrix(std::vector<std::int64_t> &core_id);
    template void Gmat<double>::invert_matrix(std::vector<std::int64_t> &core_id);
    //===============================================================================================================
    /**
     * @brief Add scalar value to the diagonal elements of matrix
     * 
     * @tparam T defines type, float or double
     * @param scale_coef floating point value
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    scale_diag(T scale_coef)
    {
        try
        {
            evolm::matrix<size_t> shapeofg;
            shapeofg = G.shape();

            for (size_t i = 0; i < shapeofg[0]; i++)
                G(i, i) = G(i, i) + scale_coef;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::scale_diag(T)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::scale_diag(T)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::scale_diag(T)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::scale_diag(float scale_coef);
    template void Gmat<double>::scale_diag(double scale_coef);
    //===============================================================================================================
    /**
     * @brief Scaling G matrix by the A matrix:
     *          G = (1-w)*Ga + wA;
     *          Ga = beta*G + alpha;
     *          mean( diag(G) )*beta + alpha = mean( diag(A) );
     *          mean( G )*beta + alpha = mean( A ).
     *          Note, scale_matr is symetric and consists of only L-part (L-store format).
     * 
     * @tparam T defines type, float or double
     * @param scale_matr std vector of the scaler A matrix
     * @param scaling_weight floatin point scalar, w
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    scale_matrix(std::vector<T> &scale_matr, T scaling_weight)
    {
        try
        {
            // It is required the scaling matrix is the same dimension as the inverting G

            evolm::matrix<size_t> shapeofg;
            shapeofg = G.shape();
            if (shapeofg[0] != shapeofg[1])
                throw std::string("G matrix has wrong dimension: number of raws is not the same as number of columns!");

            evolm::matrix<T> amat; // assumed it is in L-stored format
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
            std::cerr << "Exception in Gmat<T>::scale_matrix(std::vector<T>&, T)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::scale_matrix(std::vector<T>&, T)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::scale_matrix(std::vector<T>&, T)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::scale_matrix(std::vector<float> &scale_matr, float scaling_weight);
    template void Gmat<double>::scale_matrix(std::vector<double> &scale_matr, double scaling_weight);
    //===============================================================================================================
    /**
     * @brief Scaling G matrix by the A matrix:
     *          G = (1-w)*Ga + wA;
     *          Ga = beta*G + alpha;
     *          mean( diag(G) )*beta + alpha = mean( diag(A) );
     *          mean( G )*beta + alpha = mean( A ).
     *          Note, scale_matr is symetric and consists of only L-part (L-store format).
     * 
     * @tparam T defines type, float or double
     * @param scale_matr dense matrix class object representing the scaler A matrix
     * @param scaling_weight floatin point scalar, w
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    scale_matrix(evolm::matrix<T> &scale_matr, T scaling_weight)
    {
        try
        {
            evolm::matrix<size_t> shapeofa;
            shapeofa = scale_matr.shape();

            evolm::matrix<size_t> shapeofg;
            shapeofg = G.shape();

            if (shapeofa[0] != shapeofg[0])
                throw std::string("The passed scaleing matrix has wrong dimension: the number of rows is not the same as in the inverting G matrix!");

            if (shapeofa[1] != shapeofg[1])
                throw std::string("The passed scaleing matrix has wrong dimension: the number of columns is not the same as in the inverting G matrix!");

            // Get mean values of the scaler matrix:

            T a_ofd_mean = 0.0;
            T a_all_mean = 0.0;
            T a_diag_mean = 0.0;

#pragma omp parallel for reduction(+ : a_diag_mean)
            for (size_t i = 0; i < shapeofa[0]; i++)
            {
                a_diag_mean += scale_matr(i, i);
            }

            a_all_mean = a_diag_mean;
            a_diag_mean = a_diag_mean / ((T)shapeofa[0]);

#pragma omp parallel for reduction(+ : a_ofd_mean)
            for (size_t i = 1; i < shapeofa[0]; i++)
            {
                for (size_t j = 0; j < i; j++)
                {
                    a_ofd_mean += scale_matr(i, j);
                }
            }

            a_all_mean = (a_all_mean + 2 * a_ofd_mean) / ((T)shapeofa[0] * (T)shapeofa[0]);
            a_ofd_mean = 2 * a_ofd_mean / ((T)shapeofa[0] * (T)shapeofa[0] - (T)shapeofa[0]);

            // Get mean values of G matrix

            T g_ofd_mean = 0.0;
            T g_all_mean = 0.0;
            T g_diag_mean = 0.0;

            T alpha = 0.0;
            T betha = 1.0;

#pragma omp parallel for reduction(+ : g_diag_mean)
            for (size_t i = 0; i < shapeofg[0]; i++)
            {
                g_diag_mean += G(i, i);
            }

#pragma omp parallel for reduction(+ : g_ofd_mean)
            for (size_t i = 1; i < shapeofg[0]; i++)
            {
                for (size_t j = 0; j < i; j++)
                {
                    g_ofd_mean += G(i, j);
                }
            }

            g_all_mean = (g_diag_mean + 2 * g_ofd_mean) / ((T)shapeofg[0] * (T)shapeofg[0]);
            g_diag_mean = g_diag_mean / ((T)shapeofg[0]);

            betha = (a_all_mean - a_diag_mean) / (g_all_mean - g_diag_mean);
            alpha = a_diag_mean - g_diag_mean * betha;

            g_ofd_mean = 2 * g_ofd_mean / ((T)shapeofg[0] * (T)shapeofg[0] - (T)shapeofg[0]);

            shapeofa.clear();
            shapeofg.clear();

#ifdef UTEST
            scaling_a = alpha;
            scaling_b = betha;
            scaling_a_diag = a_diag_mean;
            scaling_a_ofd = a_ofd_mean;
            scaling_g_diag = g_diag_mean;
            scaling_g_ofd = g_ofd_mean;
#endif

            if (scaling_weight == 0.0)
            {
#pragma omp parallel for
                for (size_t i = 0; i < G.size(); i++)
                    G[i] = G[i] * betha + alpha;
            }
            else
            {
#pragma omp parallel for
                for (size_t i = 0; i < G.size(); i++)
                    G[i] = (G[i] * betha + alpha) * (1 - scaling_weight) + scale_matr[i] * scaling_weight;
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::scale_matrix(evolm::matrix<T>&, T)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::scale_matrix(evolm::matrix<T>&, T)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::scale_matrix(evolm::matrix<T>&, T)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::scale_matrix(evolm::matrix<float> &scale_matr, float scaling_weight);
    template void Gmat<double>::scale_matrix(evolm::matrix<double> &scale_matr, double scaling_weight);
    //===============================================================================================================
    /**
     * @brief Destroy the Gmat< T>:: Gmat object
     * 
     * @tparam T defines type, float or double
     * 
     */
    template <typename T> Gmat<T>::
    ~Gmat()
    {
        try
        {
            clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::~Gmat()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::~Gmat()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::~Gmat()" << '\n';
            throw;
        }
    }
    template Gmat<float>::~Gmat();
    template Gmat<double>::~Gmat();
    //===============================================================================================================
    /**
     * @brief Clears the internal storage containers
     * 
     * @tparam T defines type, float or double
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    clear()
    {
        try
        {
            G.fclear();
            G.clear();
            Z.fclear();
            Z.clear();

            gmatID.clear();
            gmatID.shrink_to_fit();
            snp_map.clear();
            anim_id_map.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::clear()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::clear()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::clear()" << '\n';
            throw;
        }
    }
    template void Gmat<float>::clear();
    template void Gmat<double>::clear();
    //===============================================================================================================
    /**
     * @brief Reads G matrix from a text file into std::vector containers
     * 
     * @tparam T defines type, float or double
     * @param gmat_file the text file in [row col value] format consisting G matrix
     * @param g_row empty std vector where G-matrix row identifiers will be storeed
     * @param g_col empty std vector where G-matrix col identifiers will be storeed
     * @param g_val empty std vector where G-matrix values will be storeed
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    read_matrix(const std::string &gmat_file, std::vector<std::int64_t> &g_row, std::vector<std::int64_t> &g_col, std::vector<T> &g_val)
    {
        try
        {
            Utilities2 u;

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
                // gmatID.push_back(int(tmp_list[1]));

                if (static_cast<std::int64_t>(tmp_list[0]) == static_cast<std::int64_t>(tmp_list[1]))
                    diagonals++;

                tmp_list.erase(tmp_list.begin(), tmp_list.end());
            }

            ped.close();

            if (!u.is_unique(gmatID))
                gmatID.erase(unique(gmatID.begin(), gmatID.end()), gmatID.end()); // here the vector should be sorted and unique

            if (diagonals != gmatID.size())
                throw std::string("There are missing diagonals in G-matrix file.");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::read_matrix(const std::string &, std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<T> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::read_matrix(const std::string &, std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<T> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::read_matrix(const std::string &, std::vector<std::int64_t> &, std::vector<std::int64_t> &, std::vector<T> &)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::read_matrix(const std::string &fname, std::vector<std::int64_t> &g_row, std::vector<std::int64_t> &g_col, std::vector<float> &g_val);
    template void Gmat<double>::read_matrix(const std::string &fname, std::vector<std::int64_t> &g_row, std::vector<std::int64_t> &g_col, std::vector<double> &g_val);
    //===============================================================================================================
    /**
     * @brief Reads G mmatrix in compact format (upper/lower triangular part) from a text file.
     *          The method should accept a full store format as well, though will write
     *          only to the lower triangular part considering the matrix is always symmetric.
     * 
     * @tparam T defines type, float or double
     * @param gmat_file the text file in [row col value] format consisting G matrix
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    read_matrix(const std::string &gmat_file)
    {
        try
        {
            Utilities2 u;

            std::vector<std::int64_t> g_row;
            std::vector<std::int64_t> g_col;
            std::vector<T> g_val;
            std::map<std::int64_t, std::int64_t> rid_map;

            read_matrix(gmat_file, g_row, g_col, g_val);

            if (gmatID.empty())
                throw std::string("Genotyped IDs vector is empty!");

            u.get_RecodedIdMap(rid_map, gmatID); // here indexing starts from 1 (but not from 0) !

            if (rid_map.empty())
                throw std::string("Recoded IDs map is empty!");

            // auto n_threads = std::thread::hardware_concurrency();
            // block_size = static_cast<unsigned int> (g_row.size()/(n_threads));

            G.resize(gmatID.size());

#pragma omp parallel for                              // schedule(static, block_size) num_threads(n_threads)
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
            std::cerr << "Exception in Gmat<T>::read_matrix(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::read_matrix(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::read_matrix(const std::string &)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::read_matrix(const std::string &gmat_file);
    template void Gmat<double>::read_matrix(const std::string &gmat_file);
    //===============================================================================================================
    /**
     * @brief Reads variants and samples data from a text file
     *        Reads file format:
     *        [observation ID] [devise code] [list of SNPs with " " delimiter]
     * 
     *          Example:
     *          18 1000 2 0 1 1 0 0 0 2 1 2
     *          19 1000 5 0 0 0 0 2 0 2 1 0
     *          20 1000 1 5 2 1 1 0 0 2 1 2
     *          21 1000 0 0 2 1 0 1 0 2 2 1
     * 
     * @tparam T defines type, float or double
     * @param snp_file the name of SNPs text file
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    read_snp(const std::string &snp_file)
    {
        try
        {
            Utilities2 u;

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

            if (!u.is_unique(t_gmatID))
                throw std::string("Thhere are repeated IDs in the processed SNPs file!");

            t_gmatID.clear();
            t_gmatID.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::read_snp(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::read_snp(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::read_snp(const std::string &)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::read_snp(const std::string &snp_file);
    template void Gmat<double>::read_snp(const std::string &snp_file);
    //===============================================================================================================
    /**
     * @brief Reads variants and samples data from a text file
     *        Reads file format:
     *        [list of SNPs with " " delimiter]
     * 
     *          Example:
     *          2 0 1 1 0 0 0 2 1 2
     *          5 0 0 0 0 2 0 2 1 0
     *          1 5 2 1 1 0 0 2 1 2
     *          0 0 2 1 0 1 0 2 2 1
     * 
     *          Reads IDs file of the format:
     *          [vector of IDs of type integer]
     * 
     * @tparam T defines type, float or double
     * @param snp_file the name of SNPs text file
     * @param ids+file text file name for samples identifiers (IDs)
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    read_snp(const std::string &snp_file, const std::string &ids_file)
    {
        try
        {
            Utilities2 u;

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

            if (!u.is_unique(t_gmatID))
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
            std::cerr << "Exception in Gmat<T>::read_snp(const std::string &, const std::string&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::read_snp(const std::string &, const std::string&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::read_snp(const std::string &, const std::string&)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::read_snp(const std::string &snp_file, const std::string &ids_file);
    template void Gmat<double>::read_snp(const std::string &snp_file, const std::string &ids_file);
    //===============================================================================================================
    /**
     * @brief Parse std string of snp variants to integers
     * 
     * @tparam T defines type, float or double
     * @param snp_str std string of snp variants
     * @param markers empty std vector where parsed variants will be copied
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    parse_string(std::string &snp_str, std::vector<int> &markers)
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
            std::cerr << "Exception in Gmat<T>::parse_string(std::string& , std::vector<int>&)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::parse_string(std::string& , std::vector<int>&)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::parse_string(std::string& , std::vector<int>&)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::parse_string(std::string &snp_str, std::vector<int> &markers);
    template void Gmat<double>::parse_string(std::string &snp_str, std::vector<int> &markers);
    //===============================================================================================================
    /**
     * @brief Construct scaled variants (snps) matrix
     * 
     * @tparam T defines type, float or double
     * 
     * @returns none
     */
    template <typename T> void Gmat<T>::
    make_zmatrix()
    {
        try
        {
            // get the number of SNPs
            std::vector<int> tmpVect;

            if ( snp_map.empty() )
                throw std::string("snp_map is empty!");

            auto it = snp_map.begin();
            std::string tmpStr = it->second;

            parse_string(tmpStr, tmpVect);

            size_t snpNum = tmpVect.size();

            tmpVect.clear();
            tmpVect.shrink_to_fit();
            tmpStr.clear();
            tmpStr.shrink_to_fit();

            // declare the matrix M
            evolm::matrix<int> M(snp_map.size(), snpNum);

            /* vector of SNPs frequences and missed values */
            std::vector<T> P(snpNum, 0.0);
            std::vector<int> missed(snpNum, 0);
            std::vector<T> missed2pq(snp_map.size(), 0.0);

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
                    M(rowI, i) = parsedMarkers[i];
                    if (parsedMarkers[i] != 0 && parsedMarkers[i] != 1 && parsedMarkers[i] != 2)
                    {
                        missed[i] += 1;
                        missedLocation[i].push_back(rowI);
                    }
                    else
                        P[i] += static_cast<T>(parsedMarkers[i]);
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
                P[i] = P[i] / (2.0 * (T)(snp_map.size() - missed[i]));
            }

            Z.resize(snp_map.size(), snpNum);

            for (size_t i = 0; i < snp_map.size(); i++)
            {
                for (size_t j = 0; j < snpNum; j++)
                {
                    Z(i, j) = (T)M(i, j) - 2.0 * P[j];
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
//#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
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
            std::cerr << "Exception in Gmat<T>::make_zmatrix()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::make_zmatrix()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::make_zmatrix()" << '\n';
            throw;
        }
    }
    template void Gmat<float>::make_zmatrix();
    template void Gmat<double>::make_zmatrix();
    //===============================================================================================================
    /**
     * @brief Construct scaled variants (snps) matrix
     * 
     * @tparam T defines type, float or double
     * @param M object of the dense matrix class with snp variants;
     *          samples are in rows, variants are in columns;
     *          M should be of dim:(sampleNum, snpNum)
     * 
     * @returns none
     * 
     */
    template <typename T> void Gmat<T>::
    make_zmatrix( evolm::matrix<int> &M )
    {
        try
        {
            evolm::matrix<size_t> shp; // shape of snp_variants matrix
            shp = M.shape();

            if ( (shp[0] == 0) ||  (shp[1] == 0) )
                throw std::string("snp_variants matrix is empty!");

            size_t snpNum = shp[1];
            size_t sampleNum = shp[0];

            // vectors of SNPs frequences and missed values
            std::vector<T> P(snpNum, 0.0);
            std::vector<int> missed(snpNum, 0);
            std::vector<T> missed2pq(sampleNum, 0.0);

            // map of missed values locations
            std::vector<std::vector<int>> missedLocation;
            for (size_t i = 0; i < snpNum; i++)
                missedLocation.push_back(std::vector<int>());

            // count missing variants, and sums along columns (specific variants)
            for (size_t row = 0; row < sampleNum; row++)
            {
                for (size_t col = 0; col < snpNum; col++)
                {
                    if ( M(row, col) == 3 )
                    {
                        missed[col] += 1;
                        missedLocation[col].push_back(row);
                    }
                    else
                        P[col] += M(row, col);
                }
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
                P[i] = P[i] / (2.0 * (T)(sampleNum - missed[i]));
            }

            Z.resize(sampleNum, snpNum);

            for (size_t i = 0; i < sampleNum; i++)
            {
                for (size_t j = 0; j < snpNum; j++)
                {
                    Z(i, j) = (T)M(i, j) - 2.0 * P[j];
                }
            }

            // modify Z matrix, so instead of missing values we put population average (0.0)

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < missedLocation.size(); i++)
            {
                for (size_t j = 0; j < missedLocation[i].size(); j++)
                {
                    Z(missedLocation[i][j], i) = 0.0;
                    missed2pq[missedLocation[i][j]] = missed2pq[missedLocation[i][j]] + 2.0 * P[i] * (1.0 - P[i]);
                }
            }

            freq = 0.0;
//#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t j = 0; j < P.size(); j++)
            {
                freq += P[j] * (1.0 - P[j]);
            }
            freq = 2.0 * freq;

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < sampleNum; i++)
            {
                missed2pq[i] = sqrt(freq / (freq - missed2pq[i]));
            }

            // After centering, adjust for missing markers for each animal;
            // adjust for sqrt[sum of 2pq over all loci /sum of 2pq over non-missing loci.

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < sampleNum; i++)
            {
                for (size_t j = 0; j < snpNum; j++)
                {
                    Z(i, j) = Z(i, j) * missed2pq[i];
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::make_zmatrix(evolm::matrix<int> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::make_zmatrix(evolm::matrix<int> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::make_zmatrix(evolm::matrix<int> &)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::make_zmatrix(evolm::matrix<int> &M);
    template void Gmat<double>::make_zmatrix(evolm::matrix<int> &M);
    //===============================================================================================================
    /**
     * @brief Constructs G matrix: G = (Z ^ 2) * (1 / freq)
     * 
     * @tparam T defines type, float or double
     * 
     * @returns none
     */
    template <typename T> void Gmat<T>::
    make_matrix()
    {
        try
        {
std::cout<<"G = Z"<<"\n";
            G = Z;
std::cout<<"Z.transpose()"<<"\n";
            Z.transpose();
std::cout<<"G = G * Z"<<"\n";
            G = G * Z;
std::cout<<"G.rectosym()"<<"\n";
            G.rectosym();
std::cout<<"G.scale(1.0 / freq)"<<"\n";
            G.scale(1.0 / freq);

            //G = (Z ^ 2) * (1 / freq);

            // Because G is symmetric,
            // make it L-stored to save some memory

            //G.rectosym();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::make_matrix()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::make_matrix()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::make_matrix()" << '\n';
            throw;
        }
    }
    template void Gmat<float>::make_matrix();
    template void Gmat<double>::make_matrix();
    //===============================================================================================================
#ifdef UTEST
    /**
     * @brief I/O interface for utest, returns some scaling properties of G matrix
     * 
     * @tparam T 
     * @param alpha 
     * @param beta 
     * @param a_diag 
     * @param a_ofd 
     * @param g_diag 
     * @param g_ofd 
     */
    template <typename T> void Gmat<T>::
    get_alpha_beta(T &alpha, T &beta, T &a_diag, T &a_ofd, T &g_diag, T &g_ofd)
    {
        try
        {
            alpha = scaling_a;
            beta = scaling_b;
            a_diag = scaling_a_diag;
            a_ofd = scaling_a_ofd;
            g_diag = scaling_g_diag;
            g_ofd = scaling_g_ofd;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Gmat<T>::get_alpha_beta(T &, T &, T &, T &, T &, T &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Gmat<T>::get_alpha_beta(T &, T &, T &, T &, T &, T &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Gmat<T>::get_alpha_beta(T &, T &, T &, T &, T &, T &)" << '\n';
            throw;
        }
    }
    template void Gmat<float>::get_alpha_beta(float &alpha, float &beta, float &a_diag, float &a_ofd, float &g_diag, float &g_ofd);
    template void Gmat<double>::get_alpha_beta(double &alpha, double &beta, double &a_diag, double &a_ofd, double &g_diag, double &g_ofd);
#endif
    //===============================================================================================================
} // end of namespace evoped