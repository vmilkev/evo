#include <iostream>
#include <fstream>
#include <chrono>
#include<thread>

#include <memory>

#include "Amat.hpp"
#include "Gmat.hpp"
#include "Hmat.hpp"

int main()
{
    bool test_A = false;
    bool test_G = false;
    bool test_H = false;
    bool test_H2 = true;

    evoped::Amat<float> ap;
    
    //std::string pedigree_file = "/Users/au383883/Documents/MY/codebase/evo/evolm/tests/data/model_2/pedigree.dat";
    //ap.make_matrix(pedigree_file, true);
    //ap.save_matrix("iA", "A1");

    evolm::smatrix<float> iA_sf;
    std::vector<std::int64_t> a_id;
    evolm::smatrix<float> irA_sf;
    std::vector<std::int64_t> ra_id;
    evolm::smatrix<float> iA22_sf;
    evolm::matrix<float> iA22;
    std::vector<std::int64_t> ai22_id;
    evolm::matrix<float> A22f;
    evolm::smatrix<float> A22_sf;
    std::vector<std::int64_t> a22_id;

    std::string snp_file = "../evolm/tests/data/large_data/SNPS/YY_plink";
    //std::string snp_file = "YY_plink";

    std::string ped_file = "../evolm/tests/data/large_data/DMU/data/dmu_pedigree_yy_20240125.txt";
    //std::string ped_file = "dmu_pedigree_yy_20240125.txt";

    std::string gid_file = "../evolm/tests/data/large_data/DMU/YY/yy.grm.id";

std::cout<<"Start testing:"<<"\n";

    if (test_A)
    {
        auto start = std::chrono::high_resolution_clock::now();

        ap.make_matrix_forgenotyped(ped_file, gid_file, true); // A22(-1)

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout <<"A22(-1) (seconds): "<< (double)duration.count() / 1000.0 << std::endl;

        ap.get_matrix("iA22", iA22, ai22_id);
        std::cout<<"num ids iA22: "<<ai22_id.size()<<"\n";

        start = std::chrono::high_resolution_clock::now();

        ap.make_matrix(ped_file, true); // full A(-1)

        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout <<"full A(-1) (seconds): "<< (double)duration.count() / 1000.0 << std::endl;

        start = std::chrono::high_resolution_clock::now();

        ap.make_matrix(ped_file, gid_file, true); // red A(-1)

        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout <<"red A(-1) (seconds): "<< (double)duration.count() / 1000.0 << std::endl;

        ap.get_matrix("iA", iA_sf, a_id);
        ap.get_matrix("irA", irA_sf, ra_id);

        std::cout<<"num ids: "<<a_id.size()<<" "<<ra_id.size()<<"\n";

        // ------------------------------------------------------------
        
        std::cout<<"making all:"<<"\n";

        start = std::chrono::high_resolution_clock::now();

        ap.make_all(ped_file, gid_file);

        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout <<"completed (milliseconds): "<< duration.count() << std::endl;

        //ap.make_all("tests/data/sstep_050/data/id4trace.PED", "tests/data/sstep_050/data/typed_050.dat"); // relatively big data
        //ap.make_all("tests/data/ped_bkg.dat", "tests/data/typed2"); // reduced A(-1)

        ap.get_matrix("iA", iA_sf, a_id);
        ap.get_matrix("irA", irA_sf, ra_id);
        ap.get_matrix("iA22", iA22_sf, ai22_id);
        ap.get_matrix("A22", A22f, a22_id);

        iA_sf.clean();
        irA_sf.clean();
        iA22_sf.clean();
        A22_sf.clean();
    }

    evoped::Gmat<float> gmat;
    evolm::matrix<float> iG;
    std::vector<std::int64_t> g_id;

    if ( test_G )
    {
        std::cout<<"Making big G, sizeof(int): "<<sizeof(int)<<"\n";
        
        gmat.make_matrix(snp_file);
        
        std::cout<<"Completed making big G"<<"\n";
        
        gmat.get_ids(g_id);        
        std::cout<<"num of ids in G: "<<g_id.size()<<"\n";
        std::cout<<"expected size of G: "<<(double)g_id.size()*g_id.size()*sizeof(float)<<"\n";

        std::cout<<"Constructing A22"<<"\n";
        evoped::Amat<float> ap2;
        evolm::matrix<float> A22;
        std::vector<int64_t> a22_ids;
        
        ap2.make_matrix_forgenotyped(ped_file,g_id,false);
        
        ap2.get_matrix("A22", A22, a22_id);
        A22.fread();
        
        std::cout<<"expected size of A22: "<<(double)a22_id.size()*a22_id.size()*sizeof(float)<<"\n";
        std::cout<<"Done"<<"\n";

        std::cout<<"Scalling G by A22:"<<"\n";
        gmat.scale_matrix(A22,0.25);
        std::cout<<"Done"<<"\n";

        std::cout<<"Inverting big G"<<"\n";
        gmat.scale_diag(0.05);        
        gmat.invert_matrix(true); // inverting as a full-store
        std::cout<<"Completed inverting big G"<<"\n";

        std::cout<<"Clear A22 and ap2 class"<<"\n";

        A22.clear();
        ap2.clear();

        std::cout<<"Saving iG:"<<"\n";

        gmat.save_matrix("iG.dmat","g_id.vect");

        std::cout<<"Clear gmat class and g_id"<<"\n";

        gmat.clear();
        g_id.clear();

        std::cout<<"num of IDs in iG: "<<g_id.size()<<"\n";

#ifdef UTEST
        float a, b, a_diag, a_ofd, g_diag, g_ofd;

        gmat.get_alpha_beta(a,b,a_diag,a_ofd,g_diag,g_ofd);

        std::cout<<"alpha: "<<a<<"\n";
        std::cout<<"beta: "<<b<<"\n";
        std::cout<<"a_diag: "<<a_diag<<"\n";
        std::cout<<"a_ofd: "<<a_ofd<<"\n";
        std::cout<<"g_diag: "<<g_diag<<"\n";
        std::cout<<"g_ofd: "<<g_ofd<<"\n";
#endif

        std::cout<<"n IDs in G: "<<g_id.size()<<"\n";
        std::cout<<"size of G: "<<iG.size()*sizeof(float)<<"\n";
        std::cout<<"\n";
    }

    evoped::Hmat<float> h;
    //evolm::smatrix<float> H;
    std::vector<std::int64_t> h_id;

    if ( test_H )
    {
        auto start = std::chrono::high_resolution_clock::now();

        h.make_matrix(snp_file, ped_file, "h_matrix_big");

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout <<"time spent for H matrix (seconds): "<< (double)duration.count() / 1000.0 << std::endl;
    }

    if (test_H2)
    {
        std::string g_matr_file = "tests/data/sstep_050/data/gmat_050.gmatrix";
        std::string ped_file2 = "tests/data/sstep_050/data/id4trace.PED";
        //std::string g_matr_file = "tests/data/big/gmat_050.gmatrix";
        //std::string ped_file2 = "tests/data/big/id4trace.PED";

        auto start = std::chrono::high_resolution_clock::now();

        h.make_matrix(g_matr_file, ped_file2, "h_matrix");

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout <<"time spent for H matrix (seconds): "<< (double)duration.count() / 1000.0 << std::endl;
    }

    return 0;
}
