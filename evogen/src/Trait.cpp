#include "Trait.hpp"

namespace evogen
{
    //===============================================================================================================
    Trait::Trait()
    {
        try
        {
            cleared = true;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::Trait()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::Trait()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::Trait()" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    Trait::~Trait()
    {
        try
        {
            clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::~Trait()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::~Trait()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::~Trait()" << '\n';
            throw;
        }
    }
    //===============================================================================================================
#ifdef PYBIND

    std::vector<float> Trait::py_tovect(pybind11::array_t<float> py_vect)
    {
        try
        {
            std::vector<float> out_vect;

            pybind11::buffer_info buf1 = py_vect.request();

            if (buf1.ndim != 1)
                throw std::runtime_error("Number of dimensions must be one");

            float *ptr1 = static_cast<float *>(buf1.ptr);

            for (pybind11::ssize_t i = 0; i < buf1.shape[0]; i++)
                out_vect.push_back(ptr1[i]);

            return out_vect;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::py_tovect(pybind11::array_t<float>)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::py_tovect(pybind11::array_t<float>)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::py_tovect(pybind11::array_t<float>)" << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    std::vector<std::vector<float>> Trait::py_tovect2d(pybind11::array_t<float> py_vect)
    {
        try
        {
            std::vector<std::vector<float>> out_vect;

            pybind11::buffer_info buf1 = py_vect.request();

            if (buf1.ndim != 2)
                throw std::runtime_error("Number of dimensions must be two");

            float *ptr1 = static_cast<float *>(buf1.ptr);

            for (pybind11::ssize_t i = 0; i < buf1.shape[0]; i++)
            {
                std::vector<float> vect;

                for (pybind11::ssize_t j = 0; j < buf1.shape[1]; j++)
                {
                    vect.push_back( ptr1[ i*buf1.shape[1]+j ] );
                }
                out_vect.push_back(vect);
            }

            return out_vect;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::py_tovect2d(pybind11::array_t<float>)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::py_tovect2d(pybind11::array_t<float>)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::py_tovect2d(pybind11::array_t<float>)" << '\n';
            throw;
        }        
    }
    //===============================================================================================================
    Trait::Trait(Population &pop,
                 pybind11::array_t<float> py_trmean,
                 pybind11::array_t<float> py_qtl_prop_chrom,
                 pybind11::array_t<float> py_corr_g, // 2d
                 pybind11::array_t<float> py_varr_g,
                 pybind11::array_t<float> py_corr_e, // 2d
                 pybind11::array_t<float> py_varr_e,
                 pybind11::array_t<float> py_envr,
                 size_t dist_model,
                 pybind11::array_t<float> py_dist_par)
    {
        try
        {
            std::vector<float> trmean = py_tovect(py_trmean);
            std::vector<float> qtl_prop_chrom = py_tovect(py_qtl_prop_chrom);
            std::vector<std::vector<float>> corr_g = py_tovect2d(py_corr_g);
            std::vector<float> varr_g = py_tovect(py_varr_g);
            std::vector<std::vector<float>> corr_e = py_tovect2d(py_corr_e);
            std::vector<float> varr_e = py_tovect(py_varr_e);
            std::vector<float> envr = py_tovect(py_envr);
            std::vector<float> dist_par = py_tovect(py_dist_par);

            // call the native c++ method

            cleared = true;
            set_trait(pop, trmean, qtl_prop_chrom, corr_g, varr_g, corr_e, varr_e, envr, dist_model, dist_par);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::Trait(...overloaded...)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::Trait(...overloaded...)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::Trait(...overloaded...)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Trait::set_trait(Population &pop,
                        pybind11::array_t<float> py_trmean,
                        pybind11::array_t<float> py_qtl_prop_chrom,
                        pybind11::array_t<float> py_corr_g, // 2d
                        pybind11::array_t<float> py_varr_g,
                        pybind11::array_t<float> py_corr_e, // 2d
                        pybind11::array_t<float> py_varr_e,
                        pybind11::array_t<float> py_envr,
                        size_t dist_model,
                        pybind11::array_t<float> py_dist_par)
    {
        try
        {
            std::vector<float> trmean = py_tovect(py_trmean);
            std::vector<float> qtl_prop_chrom = py_tovect(py_qtl_prop_chrom);
            std::vector<std::vector<float>> corr_g = py_tovect2d(py_corr_g);
            std::vector<float> varr_g = py_tovect(py_varr_g);
            std::vector<std::vector<float>> corr_e = py_tovect2d(py_corr_e);
            std::vector<float> varr_e = py_tovect(py_varr_e);
            std::vector<float> envr = py_tovect(py_envr);
            std::vector<float> dist_par = py_tovect(py_dist_par);
            
            // call the native c++ method
            set_trait(pop, trmean, qtl_prop_chrom, corr_g, varr_g, corr_e, varr_e, envr, dist_model, dist_par);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::set_trait(...python overloaded...)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::set_trait(...python overloaded...)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::set_trait(...python overloaded...)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Trait::reset_trait(Population &pop,
                            pybind11::array_t<float> py_trmean,
                            pybind11::array_t<float> py_qtl_prop_chrom,
                            pybind11::array_t<float> py_corr_g, // 2d
                            pybind11::array_t<float> py_varr_g,
                            pybind11::array_t<float> py_corr_e, // 2d
                            pybind11::array_t<float> py_varr_e,
                            pybind11::array_t<float> py_envr,
                            size_t dist_model,
                            pybind11::array_t<float> py_dist_par)
    {
        try
        {
            std::vector<float> trmean = py_tovect(py_trmean);
            std::vector<float> qtl_prop_chrom = py_tovect(py_qtl_prop_chrom);
            std::vector<std::vector<float>> corr_g = py_tovect2d(py_corr_g);
            std::vector<float> varr_g = py_tovect(py_varr_g);
            std::vector<std::vector<float>> corr_e = py_tovect2d(py_corr_e);
            std::vector<float> varr_e = py_tovect(py_varr_e);
            std::vector<float> envr = py_tovect(py_envr);
            std::vector<float> dist_par = py_tovect(py_dist_par);

            // call the native c++ method

            clear();
            set_trait(pop, trmean, qtl_prop_chrom, corr_g, varr_g, corr_e, varr_e, envr, dist_model, dist_par);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::reset_trait(...python overloaded...)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::reset_trait(...python overloaded...)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::reset_trait(...python overloaded...)" << '\n';
            throw;
        }
    }

#endif
    //===============================================================================================================
    Trait::Trait(Population &pop,
                 std::vector<float> &trmean,
                 std::vector<float> &qtl_prop_chrom,
                 std::vector<std::vector<float>> &corr_g,
                 std::vector<float> &varr_g,
                 std::vector<std::vector<float>> &corr_e,
                 std::vector<float> &varr_e,
                 std::vector<float> &envr,
                 size_t dist_model,
                 std::vector<float> &dist_par)
    {
        try
        {
            cleared = true;
            set_trait(pop, trmean, qtl_prop_chrom, corr_g, varr_g, corr_e, varr_e, envr, dist_model, dist_par);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::Trait(...overloaded...)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::Trait(...overloaded...)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::Trait(...overloaded...)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Trait::set_trait(Population &pop,
                          std::vector<float> &trmean,
                          std::vector<float> &qtl_prop_chrom,
                          std::vector<std::vector<float>> &corr_g,
                          std::vector<float> &varr_g,
                          std::vector<std::vector<float>> &corr_e,
                          std::vector<float> &varr_e,
                          std::vector<float> &envr,
                          size_t dist_model,
                          std::vector<float> &dist_par)
    {
        try
        {
            // --------- Initial set-up and input check -------------
            size_t n_trait = trmean.size();

            if (pop.size() < 1)
                throw std::string("Cannot set the trait for the empty base population!");

            std::vector<std::vector<unsigned long>> snp_table = pop.get_genome_table();
            size_t n_chr = snp_table.size();

            if (n_chr != qtl_prop_chrom.size())
                throw std::string("Number of rows in qtl_prop_chrom array does not corresponds to the number of chromosomes in the populations's genome!");

            std::vector<size_t> n_qtl_in_chr; // number of QTLs for each chromosome
            size_t n_all_qtls = 0;
            size_t n_all_snps = 0;
            for (size_t i = 0; i < n_chr; i++)
            {
                if ((qtl_prop_chrom[i] >= 0) && (qtl_prop_chrom[i] <= 1))
                {
                    size_t n_snps = snp_table[i][1] - snp_table[i][0] + 1;
                    size_t iqtls = (size_t)std::roundf(n_snps * qtl_prop_chrom[i]);

                    // std::cout<<"n_snps: "<<n_snps<<", iqtls: "<<iqtls<<", not rounded: "<<n_snps * qtl_prop_chrom[i]<<", qtl_prop_chrom[i]: "<<qtl_prop_chrom[i]<<"\n";
                    n_qtl_in_chr.push_back(iqtls);
                    n_all_qtls = n_all_qtls + iqtls;
                    n_all_snps = n_all_snps + n_snps;
                }
                else
                    throw std::string("The values of proportion of snps specified for one of the chromosomes is not eligible!");
            }

            if (corr_g.size() != n_trait)
                throw std::string("The number of rows in the genetic correlation matrix is not equal the number of correlated traits!");

            if (corr_g[0].size() != n_trait)
                throw std::string("The number of cols in the genetic correlation matrix is not equal the number of correlated traits!");

            if (corr_e.size() != n_trait)
                throw std::string("The number of rows in the residual correlation matrix is not equal the number of correlated traits!");

            if (corr_e[0].size() != n_trait)
                throw std::string("The number of cols in the residual correlation matrix is not equal the number of correlated traits!");

            if (varr_g.size() != n_trait)
                throw std::string("The number of genetic variances is not equal the number of correlated traits!");

            if (varr_e.size() != n_trait)
                throw std::string("The number of residual variances is not equal the number of correlated traits!");

            if (envr.size() != n_trait)
                throw std::string("The number of values in the ENVR matrix is not equal the number of correlated traits!");

            // --------- Sample genes (qtls) responsible for the trait ------------------------

            qtls.resize(n_all_qtls, 1);
            sample_genes(n_qtl_in_chr, snp_table, qtl_prop_chrom);

            // qtls.print("sampled genes"); // debugging
            // std::cout << "n_snps in pop: " << n_all_snps << ", n_qtls in pop: " << n_all_qtls << "\n";

            // --------- Calculate upper Cholesky decomposition of correlation matrices -------

            evolm::matrix<float> Ug;
            evolm::matrix<float> Ue;

            Ug.from_vector2d(corr_g);
            Ue.from_vector2d(corr_e);
            Ug.lchol();
            Ug.transpose();
            Ue.lchol();
            Ue.transpose();

            // Ug.printf("Ug.dat", false); // debugging
            // Ue.print("Ue"); // debugging

            // --------- Sample effects (a, e, k) ---------------------------------------------

            sample_effects(n_trait);
            sample_dom(dist_model, dist_par);

            // a.print("a1"); // debugging
            // e.print("e"); // debugging
            // k.print("k"); // debugging
            // a.printf("a0.dat", false); // debugging

            // --------- Adjust effects acording to correlation matrices ---------------------

            a = a * Ug;
            e = e * Ue;

            // a.print("a2"); // debugging
            // e.print("e"); // debugging
            //a.printf("a.dat", false); // debugging
            //e.printf("e.dat", false); // debugging
            //k.printf("k.dat", false); // debugging

            // --------- Allocating memory for traits values containers ----------------------

            realloc_traits(pop, n_trait);

            // --------- Calculate un-adjusted (to required variances) traits ----------------

            calculate_trait(pop, envr, n_trait);

            // ta.printf("ta1.dat", false); // debugging
            // te.printf("te1.dat", false); // debugging

            // --------- Calculate scaling (diagonal) matrices -------------------------------
            // square root of diag_matr of required variances * square root of inverse of diag_matr of current variances

            evolm::matrix<float> scal_a = get_scaler(varr_g, ta);
            evolm::matrix<float> scal_e = get_scaler(varr_e, te);

            // scal_a.printf("scal_a.dat", false); // debugging
            // scal_a.printf("scal_e.dat", false); // debugging

            // --------- Adjust QTL effects to match required the variance -------------------

            a = a * scal_a;
            e = e * scal_e;

            // --------- Calculate re-adjusted traits ----------------------------------------

            realloc_traits(pop, n_trait);

            calculate_trait(pop, envr, n_trait);

            //ta.printf("ta2.dat", false); // debugging
            //te.printf("te2.dat", false); // debugging

            // --------- Calculate current mean of traits -------------------------------------

            calculate_correction_mean(trmean);

            //t_mean.printf("tmean.dat", false); // debugging

            // --------- Clean some arrays ----------------------------------------------------

            ta.clear();
            te.clear();

            // --------- Store genome structure of base population -----------------------------

            base_genome_structure = pop.get_genome_structure();

            cleared = false;
        }
        catch (const std::exception &e)
        {
            std::string msg = "Exception in Trait::set_trait(Population &,\n"
                              "std::vector<float> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<std::vector<float>> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<std::vector<float>> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<float> &,\n"
                              "size_t\n"
                              "std::vector<float> &,\n)";
            std::cerr << msg << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::string msg = "Exception in Trait::set_trait(Population &,\n"
                              "std::vector<float> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<std::vector<float>> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<std::vector<float>> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<float> &,\n"
                              "size_t\n"
                              "std::vector<float> &,\n)";
            std::cerr << msg << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::string msg = "Exception in Trait::set_trait(Population &,\n"
                              "std::vector<float> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<std::vector<float>> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<std::vector<float>> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<float> &,\n"
                              "size_t\n"
                              "std::vector<float> &,\n)";
            std::cerr << msg << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Trait::reset_trait(Population &pop,
                            std::vector<float> &trmean,
                            std::vector<float> &qtl_prop_chrom,
                            std::vector<std::vector<float>> &corr_g,
                            std::vector<float> &varr_g,
                            std::vector<std::vector<float>> &corr_e,
                            std::vector<float> &varr_e,
                            std::vector<float> &envr,
                            size_t dist_model,
                            std::vector<float> &dist_par)
    {
        try
        {
            clear();
            set_trait(pop, trmean, qtl_prop_chrom, corr_g, varr_g, corr_e, varr_e, envr, dist_model, dist_par);
        }
        catch (const std::exception &e)
        {
            std::string msg = "Exception in Trait::reset_trait(Population &,\n"
                              "std::vector<float> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<std::vector<float>> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<std::vector<float>> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<float> &,\n"
                              "size_t\n"
                              "std::vector<float> &,\n)";
            std::cerr << msg << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::string msg = "Exception in Trait::reset_trait(Population &,\n"
                              "std::vector<float> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<std::vector<float>> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<std::vector<float>> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<float> &,\n"
                              "size_t\n"
                              "std::vector<float> &,\n)";
            std::cerr << msg << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::string msg = "Exception in Trait::reset_trait(Population &,\n"
                              "std::vector<float> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<std::vector<float>> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<std::vector<float>> &\n"
                              "std::vector<float> &,\n"
                              "std::vector<float> &,\n"
                              "size_t\n"
                              "std::vector<float> &,\n)";
            std::cerr << msg << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Trait::clear()
    {
        try
        {
            a.clear();
            e.clear();
            k.clear();
            qtls.clear();
            ta.clear();
            te.clear();

            base_genome_structure.clear();
            base_genome_structure.shrink_to_fit();

            cleared = true;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::clear()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::clear()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::clear()" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    bool Trait::is_cleared()
    {
        try
        {
            return cleared;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::is_cleared()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::is_cleared()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::is_cleared()" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Trait::sample_dom(size_t which_dist, std::vector<float> &dist_param)
    {
        try
        {
            Utilites u;

            std::vector<float> v;

            switch (which_dist)
            {
            case 1 /* Uniform distribution */:
                if (dist_param.size() != 2)
                    throw std::string("The number of provided uniform distribution parameters is not equal 2!");

                v = u.get_runi_rand(qtls.size(), dist_param[0], dist_param[1], false);
                break;

            case 2 /* Normal distribution */:
                if (dist_param.size() != 2)
                    throw std::string("The number of provided normal distribution parameters is not equal 2!");

                v = u.get_norm_rand(qtls.size(), dist_param[0], dist_param[1], false);
                break;

            case 3 /* Gamma distribution */:
                if (dist_param.size() != 1)
                    throw std::string("The number of provided gamma distribution parameters is not equal 1!");

                v = u.get_gamma_rand(qtls.size(), 2.0f, dist_param[0], false);
                break;

            default:
                break;
            }

            k.from_vector(v);

            v.clear();
            v.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::sample_genes(size_t, std::vector<float> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::sample_genes(size_t, std::vector<float> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::sample_genes(size_t, std::vector<float> &)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Trait::sample_genes(std::vector<size_t> &n_qtl, std::vector<std::vector<unsigned long>> &stable, std::vector<float> &qtl_prop)
    {
        /* We are sampling not physical locations but indexes pointing to the std::vector of maarkers in the genome. */
        try
        {
            Utilites u;
            size_t first_qtl = 0;

            for (size_t i = 0; i < n_qtl.size(); i++)
            {
                if (qtl_prop[i] < 0.5) // this condition avoids many repetitive genes
                {
                    std::vector<size_t> igenes = u.get_uni_rand(n_qtl[i], stable[i][0], stable[i][1], false);
                    std::sort(igenes.begin(), igenes.end());

                    for (size_t j = 0; j < igenes.size(); j++)
                    {
                        qtls[first_qtl + j] = igenes[j];
                    }
                    first_qtl = first_qtl + igenes.size();
                }
                else
                {
                    for (size_t j = 0; j < n_qtl[i]; j++)
                    {
                        qtls[first_qtl + j] = stable[i][0] + j;
                    }
                    first_qtl = first_qtl + n_qtl[i];
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::sample_genes(std::vector<size_t> &, std::vector<std::vector<unsigned long>> &, std::vector<float> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::sample_genes(std::vector<size_t> &, std::vector<std::vector<unsigned long>> &, std::vector<float> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::sample_genes(std::vector<size_t> &, std::vector<std::vector<unsigned long>> &, std::vector<float> &)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Trait::sample_effects(size_t n_trate)
    {
        try
        {
            Utilites u;

            std::vector<std::vector<float>> v;

            for (size_t i = 0; i < n_trate; i++)
            {
                std::vector<float> v1 = u.get_norm_rand(qtls.size(), 0.0f, 1.0f, true);
                v.push_back(v1);
            }

            a.from_vector2d(v);
            a.transpose();

            v.clear();
            v.shrink_to_fit();

            for (size_t i = 0; i < n_trate; i++)
            {
                std::vector<float> v1 = u.get_norm_rand(qtls.size(), 0.0f, 1.0f, true);
                v.push_back(v1);
            }

            e.from_vector2d(v);
            e.transpose();

            v.clear();
            v.shrink_to_fit();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::sample_effects(size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::sample_effects(size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::sample_effects(size_t)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Trait::calculate_trait(Population &in_pop, std::vector<float> &envr, size_t n_trate)
    {
        try
        {
            size_t n_ploidy = in_pop.get_ploidy();
            size_t n_indiv = in_pop.size();

            float p_eff = ploidy_effect(n_ploidy, 2.0f);

            for (size_t trait = 0; trait < n_trate; trait++)
            {
                for (size_t individ = 0; individ < n_indiv; individ++)
                {
                    for (size_t iqtl = 0; iqtl < qtls.size(); iqtl++)
                    {
                        std::vector<int> locus_state = get_locus_state(in_pop, individ, qtls(iqtl, 0));

                        float qtl_val = (float)locus_state[0];   // qtl value of genotype, number of ref alleles
                        float dom_cases = (float)locus_state[1]; // cases of dominance

                        ta(individ, trait) = ta(individ, trait) + p_eff * (1.0f + dom_cases * k(iqtl, 0)) * a(iqtl, trait) * qtl_val;
                        te(individ, trait) = te(individ, trait) + p_eff * (e(iqtl, trait) + 0.5f * envr[trait] * std::abs(e(iqtl, trait))) * qtl_val;
                    }
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::calculate_trait(Population &, std::vector<float> &, size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::calculate_trait(Population &, std::vector<float> &, size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::calculate_trait(Population &, std::vector<float> &, size_t)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    void Trait::calculate_trait(Population &in_pop, std::vector<size_t> &ind_list, std::vector<float> &envr, size_t n_trate)
    {
        try
        {
            size_t n_ploidy = in_pop.get_ploidy();

            size_t n_indiv = ind_list.size();

            float p_eff = ploidy_effect(n_ploidy, 2.0f);

            for (size_t trait = 0; trait < n_trate; trait++)
            {
                for (size_t individ = 0; individ < n_indiv; individ++)
                {
                    for (size_t iqtl = 0; iqtl < qtls.size(); iqtl++)
                    {
                        std::vector<int> locus_state = get_locus_state(in_pop, ind_list[individ], qtls(iqtl, 0)); // access only active individuals

                        float qtl_val = (float)locus_state[0];   // qtl value of genotype, number of ref alleles
                        float dom_cases = (float)locus_state[1]; // cases of dominance

                        ta(individ, trait) = ta(individ, trait) + p_eff * (1.0f + dom_cases * k(iqtl, 0)) * a(iqtl, trait) * qtl_val;
                        te(individ, trait) = te(individ, trait) + p_eff * (e(iqtl, trait) + 0.5f * envr[trait] * std::abs(e(iqtl, trait))) * qtl_val;
                    }
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::calculate_trait(Population &, std::vector<size_t> &, std::vector<float> &, size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::calculate_trait(Population &, std::vector<size_t> &, std::vector<float> &, size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::calculate_trait(Population &, std::vector<size_t> &, std::vector<float> &, size_t)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    float Trait::ploidy_effect(size_t n_ploidy, float degree)
    {
        try
        {
            return (1.0f / (1.0f + std::pow((1.0f - 2.0f / n_ploidy), degree)));
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::ploidy_effect(size_t, float)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::ploidy_effect(size_t, float)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::ploidy_effect(size_t, float)" << '\n';
            throw;
        }
    }
    //===============================================================================================================
    std::vector<int> Trait::get_locus_state(Population &in_pop, size_t individ, size_t qtl_position)
    {
        try
        {
            // std::cout<<"Animal: "<<individ<<", pos: "<<qtl_position<<"\n";

            std::vector<int> lstate;

            int ref_alleles_count = 0; // number of reference alleles for all pairs of haplotypes
            size_t alleles_count = 0;  // number of reference alleles for each pair of haplotypes
            int dom_cases = 0;         // how many dominance cases in a locus; 1 < dom_cases < inf

            std::vector<short> val = in_pop.get_genome_at(individ, qtl_position);

            for (size_t i = 0; i < val.size(); i++)
            {
                ref_alleles_count = ref_alleles_count + (int)val[i]; // sum for entire DNA strands
                alleles_count = alleles_count + (size_t)val[i];      // sum for pair of DNA strands

                // std::cout<<"strand: "<<i+1<<", ref allele: "<<(int)val[i]<<"\n";

                if (((i + 1) % 2) == 0) // consider complete pair of haplotypes
                {
                    if (alleles_count == 1) // evaluate dominance condition on pair of haplotypes
                        dom_cases = dom_cases + 1;
                    alleles_count = 0;
                }
            }

            lstate.push_back(ref_alleles_count);
            lstate.push_back(dom_cases);

            // std::cout<<"ref_alleles_count: "<<ref_alleles_count<<", dom_cases: "<<dom_cases<<"\n";

            return lstate;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::get_locus_state(Population &, size_t, size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::get_locus_state(Population &, size_t, size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::get_locus_state(Population &, size_t, size_t)" << '\n';
            throw;
        }
    }
    //===============================================================================================================

    evolm::matrix<float> Trait::var_diag(evolm::matrix<float> &arr)
    {
        try
        {
            evolm::matrix<size_t> shape(2, 1);
            shape = arr.shape();

            size_t row = shape[0];
            size_t col = shape[1];

            evolm::matrix<float> out_var(col, col);

            // calculate mean
            evolm::matrix<float> mean(col, col);

            for (size_t i = 0; i < col; i++)
            {
                for (size_t j = 0; j < row; j++)
                {
                    mean(i, i) += arr(j, i);
                }
            }
            mean = mean * (1.0f / (float)row);

            // calculate variannce
            for (size_t i = 0; i < col; i++)
            {
                for (size_t j = 0; j < row; j++)
                {
                    out_var(i, i) += (arr(j, i) - mean(i, i)) * (arr(j, i) - mean(i, i));
                }
            }

            return out_var * (1.0f / (float)row);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::var_diag( evolm::matrix<float> & )" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::var_diag( evolm::matrix<float> & )" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::var_diag( evolm::matrix<float> & )" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    evolm::matrix<float> Trait::get_scaler(std::vector<float> &in_var, evolm::matrix<float> &in_arr)
    {
        try
        {
            // evolm::matrix<float> out;
            evolm::matrix<float> requested_std(in_var.size(), in_var.size());
            evolm::matrix<float> current_std(in_var.size(), in_var.size());

            evolm::matrix<float> var = var_diag(in_arr);

            for (size_t i = 0; i < in_var.size(); i++)
            {
                requested_std(i, i) = std::sqrt(in_var[i]);
                current_std(i, i) = 1.0f / std::sqrt(var(i, i));
            }

            // requested_std.printf("req_var.dat", false);
            // current_std.printf("cur_var.dat", false);

            return requested_std * current_std;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::get_scaler( std::vector<float> &, evolm::matrix<float> & )" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::get_scaler( std::vector<float> &, evolm::matrix<float> & )" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::get_scaler( std::vector<float> &, evolm::matrix<float> & )" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Trait::realloc_traits(Population &in_pop, size_t n_trait)
    {
        try
        {
            ta.clear();
            te.clear();

            ta.resize(in_pop.size(), n_trait);

            if (ta.size() != in_pop.size() * n_trait)
                throw std::string("The memory for the genotypic trait values (ta array) was not allocated properly!");

            te.resize(in_pop.size(), n_trait);

            if (te.size() != in_pop.size() * n_trait)
                throw std::string("The memory for the environmental trait values (te array) was not allocated properly!");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::realloc_traits(Population &,size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::realloc_traits(Population &,size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::realloc_traits(Population &,size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Trait::realloc_traits(size_t pop_size, size_t n_trait)
    {
        try
        {
            ta.clear();
            te.clear();

            ta.resize(pop_size, n_trait);

            if (ta.size() != pop_size * n_trait)
                throw std::string("The memory for the genotypic trait values (ta array) was not allocated properly!");

            te.resize(pop_size, n_trait);

            if (te.size() != pop_size * n_trait)
                throw std::string("The memory for the environmental trait values (te array) was not allocated properly!");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::realloc_traits(size_t,size_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::realloc_traits(size_t,size_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::realloc_traits(size_t,size_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Trait::calculate_correction_mean(std::vector<float> &in_mean)
    {
        try
        {
            evolm::matrix<size_t> shape(2, 1);
            shape = ta.shape();

            size_t row = shape[0];
            size_t col = shape[1];

            t_mean.resize(col, 1);

            for (size_t i = 0; i < col; i++)
            {
                for (size_t j = 0; j < row; j++)
                {
                    t_mean(i, 0) += ta(j, i) + te(j, i);
                }
                // std::cout<<"in_mean[i] => "<<in_mean[i]<<", t_mean(i, 0) / (float)row => "<<t_mean(i, 0) / (float)row<<"\n";
                t_mean(i, 0) = in_mean[i] - t_mean(i, 0) / (float)row;
                // std::cout<<"t_mean(i, 0) => "<<t_mean(i, 0)<<"\n";
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::calculate_correction_mean(std::vector<float> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::calculate_correction_mean(std::vector<float> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::calculate_correction_mean(std::vector<float> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

#ifdef PYBIND

    void Trait::get_observations(Population &in_pop, pybind11::array_t<float> py_env)
    {
        try
        {
            std::vector<float> env = py_tovect(py_env);

            if (cleared)
                throw std::string("The trait is configured, hence the call of set_trait() is required!");

            evolm::matrix<size_t> shape(2, 1);
            shape = a.shape();

            size_t n_individuals = in_pop.size();
            size_t n_traits = shape[1];

            if (n_individuals < 1)
                throw std::string("Cannot provide observation on empty population!");

            if (env.size() != n_traits)
                throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

            if (base_genome_structure != in_pop.get_genome_structure())
                throw std::string("The genome structure of the base population is not the same as in the population being observed!");

            realloc_traits(in_pop, n_traits);
            calculate_trait(in_pop, env, n_traits);

            evolm::matrix<float> t(n_individuals, n_traits);

            for (size_t i = 0; i < n_traits; i++)
            {
                for (size_t j = 0; j < n_individuals; j++)
                {
                    t(j, i) = ta(j, i) + te(j, i) + t_mean(i, 0);
                }
            }

            // register the observed phenotypes for corresponding individual
            for (size_t i = 0; i < n_individuals; i++)
            {
                std::vector<float> obs;
                for (size_t j = 0; j < n_traits; j++)
                {
                    obs.push_back(t(i, j));
                }
                in_pop.phenotype_at_cpp(i, obs);
            }

            ta.clear();
            te.clear();
            t.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, pybind11::array_t<float>)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, pybind11::array_t<float>)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, pybind11::array_t<float>)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Trait::get_observations(Population &in_pop, pybind11::array_t<float> py_env, pybind11::array_t<float> py_t)
    {
        try
        {
            std::vector<float> env = py_tovect(py_env);

            if (cleared)
                throw std::string("The trait is configured, hence the call of set_trait() is required!");

            evolm::matrix<size_t> shape(2, 1);
            shape = a.shape();

            size_t n_individuals = in_pop.size();
            size_t n_traits = shape[1];

            if (n_individuals < 1)
                throw std::string("Cannot provide observation on empty population!");

            if (env.size() != n_traits)
                throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

            if (base_genome_structure != in_pop.get_genome_structure())
                throw std::string("The genome structure of the base population is not the same as in the population being observed!");

            realloc_traits(in_pop, n_traits);
            calculate_trait(in_pop, env, n_traits);

            evolm::matrix<float> t(n_individuals, n_traits);

            for (size_t i = 0; i < n_traits; i++)
            {
                for (size_t j = 0; j < n_individuals; j++)
                {
                    t(j, i) = ta(j, i) + te(j, i) + t_mean(i, 0);
                }
            }

            // register the observed phenotypes for corresponding individual
            for (size_t i = 0; i < n_individuals; i++)
            {
                std::vector<float> obs;
                for (size_t j = 0; j < n_traits; j++)
                {
                    obs.push_back(t(i, j));
                }
                in_pop.phenotype_at_cpp(i, obs);
            }

            //------------------------------
            std::vector<std::vector<float>> out_t;

            t.to_vector2d(out_t);

            size_t N = out_t.size();
            size_t M = out_t[0].size();

            py_t.resize( {N,M}, false );

            //  allocate the buffer
            pybind11::buffer_info buf2 = py_t.request();

            if (buf2.ndim != 2)
                throw std::runtime_error("Number of dimensions must be two");

            float *ptr2 = (float *) buf2.ptr;

            for (size_t i = 0; i < N; i++) {
                for (size_t j = 0; j < M; j++) {
                    ptr2[ i * M + j ] = out_t[i][j];
                }
            }
            //------------------------------

            ta.clear();
            te.clear();
            t.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, pybind11::array_t<float>, pybind11::array_t<float>)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, pybind11::array_t<float>, pybind11::array_t<float>)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, pybind11::array_t<float>, pybind11::array_t<float>)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Trait::get_observations(Population &in_pop, pybind11::array_t<float> py_env, pybind11::array_t<float> py_t, pybind11::array_t<int> py_g)
    {
        try
        {
            std::vector<float> env = py_tovect(py_env);

            if (cleared)
                throw std::string("The trait is configured, hence the call of set_trait() is required!");

            evolm::matrix<size_t> shape(2, 1);
            shape = a.shape();

            size_t n_individuals = in_pop.size();
            size_t n_traits = shape[1];

            if (n_individuals < 1)
                throw std::string("Cannot provide observation on empty population!");

            if (env.size() != n_traits)
                throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

            if (base_genome_structure != in_pop.get_genome_structure())
                throw std::string("The genome structure of the base population is not the same as in the population being observed!");

            realloc_traits(in_pop, n_traits);
            calculate_trait(in_pop, env, n_traits);

            evolm::matrix<float> t(n_individuals, n_traits);

            for (size_t i = 0; i < n_traits; i++)
            {
                for (size_t j = 0; j < n_individuals; j++)
                {
                    t(j, i) = ta(j, i) + te(j, i) + t_mean(i, 0);
                }
            }

            // register the observed phenotypes for corresponding individual
            for (size_t i = 0; i < n_individuals; i++)
            {
                std::vector<float> obs;
                for (size_t j = 0; j < n_traits; j++)
                {
                    obs.push_back(t(i, j));
                }
                in_pop.phenotype_at_cpp(i, obs);
            }

            //------------------------------
            std::vector<std::vector<float>> out_t;

            t.to_vector2d(out_t);

            size_t N = out_t.size();
            size_t M = out_t[0].size();

            py_t.resize( {N,M}, false );

            //  allocate the buffer
            pybind11::buffer_info buf2 = py_t.request();

            if (buf2.ndim != 2)
                throw std::runtime_error("Number of dimensions must be two");

            float *ptr2 = (float *) buf2.ptr;

            for (size_t i = 0; i < N; i++) {
                for (size_t j = 0; j < M; j++) {
                    ptr2[ i * M + j ] = out_t[i][j];
                }
            }
            //------------------------------

            ta.clear();
            te.clear();
            t.clear();

            //--------------------------------------
            std::vector<std::vector<short>> out_g;

            in_pop.get_all_genotypes(out_g);

            N = M = 0;

            N = out_g.size();
            M = out_g[0].size();

            py_g.resize( {N,M}, false );

            //  allocate the buffer
            pybind11::buffer_info buf3 = py_g.request();

            if (buf3.ndim != 2)
                throw std::runtime_error("Number of dimensions must be two");

            int *ptr3 = (int *) buf3.ptr;

            for (size_t i = 0; i < N; i++) {
                for (size_t j = 0; j < M; j++) {
                    ptr3[ i * M + j ] = (int)out_g[i][j];
                }
            }
            //--------------------------------------
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, pybind11::array_t<float>, pybind11::array_t<float>, pybind11::array_t<float>)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Trait::get_observations(Population &in_pop, pybind11::array_t<float> py_env, const std::string &out_t)
    {
        try
        {
            std::vector<float> env = py_tovect(py_env);

            if (cleared)
                throw std::string("The trait is configured, hence the call of set_trait() is required!");

            evolm::matrix<size_t> shape(2, 1);
            shape = a.shape();

            size_t n_individuals = in_pop.size();
            size_t n_traits = shape[1];

            if (n_individuals < 1)
                throw std::string("Cannot provide observation on empty population!");

            if (env.size() != n_traits)
                throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

            if (base_genome_structure != in_pop.get_genome_structure())
                throw std::string("The genome structure of the base population is not the same as in the population being observed!");

            realloc_traits(in_pop, n_traits);
            calculate_trait(in_pop, env, n_traits);

            evolm::matrix<float> t(n_individuals, n_traits);

            for (size_t i = 0; i < n_traits; i++)
            {
                for (size_t j = 0; j < n_individuals; j++)
                {
                    t(j, i) = ta(j, i) + te(j, i) + t_mean(i, 0);
                }
            }

            // register the observed phenotypes for corresponding individual
            for (size_t i = 0; i < n_individuals; i++)
            {
                std::vector<float> obs;
                for (size_t j = 0; j < n_traits; j++)
                {
                    obs.push_back(t(i, j));
                }
                in_pop.phenotype_at_cpp(i, obs);
            }

            t.printf(out_t, false); // not in append mode

            ta.clear();
            te.clear();
            t.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, pybind11::array_t<float>, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, pybind11::array_t<float>, const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, pybind11::array_t<float>, const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Trait::get_observations(Population &in_pop, pybind11::array_t<float> py_env, const std::string &out_t, const std::string &out_g)
    {
        try
        {
            std::vector<float> env = py_tovect(py_env);

            if (cleared)
                throw std::string("The trait is configured, hence the call of set_trait() is required!");

            evolm::matrix<size_t> shape(2, 1);
            shape = a.shape();

            size_t n_individuals = in_pop.size();
            size_t n_traits = shape[1];

            if (n_individuals < 1)
                throw std::string("Cannot provide observation on empty population!");

            if (env.size() != n_traits)
                throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

            if (base_genome_structure != in_pop.get_genome_structure())
                throw std::string("The genome structure of the base population is not the same as in the population being observed!");

            realloc_traits(in_pop, n_traits);
            calculate_trait(in_pop, env, n_traits);

            evolm::matrix<float> t(n_individuals, n_traits);

            for (size_t i = 0; i < n_traits; i++)
            {
                for (size_t j = 0; j < n_individuals; j++)
                {
                    t(j, i) = ta(j, i) + te(j, i) + t_mean(i, 0);
                }
            }

            // register the observed phenotypes for corresponding individual
            for (size_t i = 0; i < n_individuals; i++)
            {
                std::vector<float> obs;
                for (size_t j = 0; j < n_traits; j++)
                {
                    obs.push_back(t(i, j));
                }
                in_pop.phenotype_at_cpp(i, obs);
            }

            t.printf(out_t, false); // not in append mode

            ta.clear();
            te.clear();
            t.clear();

            Utilites u;
            u.fremove(out_g);

            in_pop.get_all_genotypes(out_g);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, pybind11::array_t<float>, const std::string &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, pybind11::array_t<float>, const std::string &, const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, pybind11::array_t<float>, const std::string &, const std::string &)" << '\n';
            throw;
        }
    }

#endif

    //===============================================================================================================

    void Trait::get_observations(Population &in_pop, std::vector<float> &env)
    {
        try
        {
            if (cleared)
                throw std::string("The trait is configured, hence the call of set_trait() is required!");

            evolm::matrix<size_t> shape(2, 1);
            shape = a.shape();

            size_t n_individuals = in_pop.size();
            size_t n_traits = shape[1];

            if (n_individuals < 1)
                throw std::string("Cannot provide observation on empty population!");

            if (env.size() != n_traits)
                throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

            if (base_genome_structure != in_pop.get_genome_structure())
                throw std::string("The genome structure of the base population is not the same as in the population being observed!");

            realloc_traits(in_pop, n_traits);
            calculate_trait(in_pop, env, n_traits);

            evolm::matrix<float> t(n_individuals, n_traits);

            for (size_t i = 0; i < n_traits; i++)
            {
                for (size_t j = 0; j < n_individuals; j++)
                {
                    t(j, i) = ta(j, i) + te(j, i) + t_mean(i, 0);
                }
            }

            // register the observed phenotypes for corresponding individual
            for (size_t i = 0; i < n_individuals; i++)
            {
                std::vector<float> obs;
                for (size_t j = 0; j < n_traits; j++)
                {
                    obs.push_back(t(i, j));
                }
                in_pop.phenotype_at_cpp(i, obs);
            }

            ta.clear();
            te.clear();
            t.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, std::vector<float> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, std::vector<float> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, std::vector<float> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Trait::get_observations(Population &in_pop, std::vector<float> &env, std::vector<std::vector<float>> &out_t, std::vector<std::vector<short>> &out_g)
    {
        try
        {
            if (cleared)
                throw std::string("The trait is configured, hence the call of set_trait() is required!");

            evolm::matrix<size_t> shape(2, 1);
            shape = a.shape();

            size_t n_individuals = in_pop.size();
            size_t n_traits = shape[1];

            if (n_individuals < 1)
                throw std::string("Cannot provide observation on empty population!");

            if (env.size() != n_traits)
                throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

            if (base_genome_structure != in_pop.get_genome_structure())
                throw std::string("The genome structure of the base population is not the same as in the population being observed!");

            realloc_traits(in_pop, n_traits);
            calculate_trait(in_pop, env, n_traits);

            evolm::matrix<float> t(n_individuals, n_traits);

            for (size_t i = 0; i < n_traits; i++)
            {
                for (size_t j = 0; j < n_individuals; j++)
                {
                    t(j, i) = ta(j, i) + te(j, i) + t_mean(i, 0);
                }
            }

            // register the observed phenotypes for corresponding individual
            for (size_t i = 0; i < n_individuals; i++)
            {
                std::vector<float> obs;
                for (size_t j = 0; j < n_traits; j++)
                {
                    obs.push_back(t(i, j));
                }
                in_pop.phenotype_at_cpp(i, obs);
            }

            t.to_vector2d(out_t);

            ta.clear();
            te.clear();
            t.clear();

            in_pop.get_all_genotypes(out_g);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, std::vector<float> &, std::vector<std::vector<float>> &, std::vector<std::vector<short>> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, std::vector<float> &, std::vector<std::vector<float>> &, std::vector<std::vector<short>> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, std::vector<float> &, std::vector<std::vector<float>> &, std::vector<std::vector<short>> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Trait::get_observations(Population &in_pop, std::vector<float> &env, const std::string &out_t, const std::string &out_g)
    {
        try
        {
            if (cleared)
                throw std::string("The trait is configured, hence the call of set_trait() is required!");

            evolm::matrix<size_t> shape(2, 1);
            shape = a.shape();

            size_t n_individuals = in_pop.size();
            size_t n_traits = shape[1];

            if (n_individuals < 1)
                throw std::string("Cannot provide observation on empty population!");

            if (env.size() != n_traits)
                throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

            if (base_genome_structure != in_pop.get_genome_structure())
                throw std::string("The genome structure of the base population is not the same as in the population being observed!");

            realloc_traits(in_pop, n_traits);
            calculate_trait(in_pop, env, n_traits);

            evolm::matrix<float> t(n_individuals, n_traits);

            for (size_t i = 0; i < n_traits; i++)
            {
                for (size_t j = 0; j < n_individuals; j++)
                {
                    t(j, i) = ta(j, i) + te(j, i) + t_mean(i, 0);
                }
            }

            // register the observed phenotypes for corresponding individual
            for (size_t i = 0; i < n_individuals; i++)
            {
                std::vector<float> obs;
                for (size_t j = 0; j < n_traits; j++)
                {
                    obs.push_back(t(i, j));
                }
                in_pop.phenotype_at_cpp(i, obs);
            }

            t.printf(out_t, false); // not in append mode

            ta.clear();
            te.clear();
            t.clear();

            Utilites u;
            u.fremove(out_g);

            in_pop.get_all_genotypes(out_g);
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, std::vector<float> &, const std::string &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, std::vector<float> &, const std::string &, const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, std::vector<float> &, const std::string &, const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Trait::get_observations(Population &in_pop, std::vector<float> &env, const std::string &out_t)
    {
        try
        {
            if (cleared)
                throw std::string("The trait is configured, hence the call of set_trait() is required!");

            evolm::matrix<size_t> shape(2, 1);
            shape = a.shape();

            size_t n_individuals = in_pop.size();
            size_t n_traits = shape[1];

            if (n_individuals < 1)
                throw std::string("Cannot provide observation on empty population!");

            if (env.size() != n_traits)
                throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

            if (base_genome_structure != in_pop.get_genome_structure())
                throw std::string("The genome structure of the base population is not the same as in the population being observed!");

            realloc_traits(in_pop, n_traits);
            calculate_trait(in_pop, env, n_traits);

            evolm::matrix<float> t(n_individuals, n_traits);

            for (size_t i = 0; i < n_traits; i++)
            {
                for (size_t j = 0; j < n_individuals; j++)
                {
                    t(j, i) = ta(j, i) + te(j, i) + t_mean(i, 0);
                }
            }

            // register the observed phenotypes for corresponding individual
            for (size_t i = 0; i < n_individuals; i++)
            {
                std::vector<float> obs;
                for (size_t j = 0; j < n_traits; j++)
                {
                    obs.push_back(t(i, j));
                }
                in_pop.phenotype_at_cpp(i, obs);
            }

            t.printf(out_t, false); // not in append mode

            ta.clear();
            te.clear();
            t.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, std::vector<float> &, const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, std::vector<float> &, const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, std::vector<float> &, const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Trait::get_observations(Population &in_pop, std::vector<float> &env, std::vector<std::vector<float>> &out_t)
    {
        try
        {
            if (cleared)
                throw std::string("The trait is configured, hence the call of set_trait() is required!");

            evolm::matrix<size_t> shape(2, 1);
            shape = a.shape();

            size_t n_individuals = in_pop.size();
            size_t n_traits = shape[1];

            if (n_individuals < 1)
                throw std::string("Cannot provide observation on empty population!");

            if (env.size() != n_traits)
                throw std::string("The demension of the array ENV array does not correspond to the number of traits!");

            if (base_genome_structure != in_pop.get_genome_structure())
                throw std::string("The genome structure of the base population is not the same as in the population being observed!");

            realloc_traits(in_pop, n_traits);
            calculate_trait(in_pop, env, n_traits);

            evolm::matrix<float> t(n_individuals, n_traits);

            for (size_t i = 0; i < n_traits; i++)
            {
                for (size_t j = 0; j < n_individuals; j++)
                {
                    t(j, i) = ta(j, i) + te(j, i) + t_mean(i, 0);
                }
            }

            // register the observed phenotypes for corresponding individual
            for (size_t i = 0; i < n_individuals; i++)
            {
                std::vector<float> obs;
                for (size_t j = 0; j < n_traits; j++)
                {
                    obs.push_back(t(i, j));
                }
                in_pop.phenotype_at_cpp(i, obs);
            }

            t.to_vector2d(out_t);

            ta.clear();
            te.clear();
            t.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, std::vector<float> &, std::vector<std::vector<float>> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, std::vector<float> &, std::vector<std::vector<float>> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Trait::get_observations(Population &, std::vector<float> &, std::vector<std::vector<float>> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    /*
                Trait::?()
                {
                    try
                    {
                    }
                    catch (const std::exception &e)
                    {
                        std::cerr << "Exception in Trait::?()" << '\n';
                        std::cerr << e.what() << '\n';
                        throw;
                    }
                    catch (const std::string &e)
                    {
                        std::cerr << "Exception in Trait::?()" << '\n';
                        std::cerr <<"Reason: "<< e << '\n';
                        throw;
                    }
                    catch (...)
                    {
                        std::cerr << "Exception in Trait::?()" << '\n';
                        throw;
                    }
                }

                //===============================================================================================================
            */
} // end of the namespace evogen