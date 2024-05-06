#include <iostream>
#include <chrono>
#include<thread>

#include "Amat.hpp"
#include "Gmat.hpp"
#include "Hmat.hpp"

int main()
{
    evoped::Amat<float> ap;

    evolm::matrix<double> iA;
    evolm::smatrix<double> iA_s;
    evolm::smatrix<float> iA_sf;
    std::vector<std::int64_t> a_id;
    evolm::matrix<double> irA;
    evolm::smatrix<double> irA_s;
    evolm::smatrix<float> irA_sf;
    std::vector<std::int64_t> ra_id;
    evolm::matrix<double> iA22;
    evolm::matrix<float> iA22f;
    std::vector<std::int64_t> ai22_id;
    evolm::matrix<double> A22;
    evolm::matrix<float> A22f;
    std::vector<std::int64_t> a22_id;

    evolm::matrix<double> A;
    std::vector<std::int64_t> aa_id;
    evolm::matrix<double> rA;
    std::vector<std::int64_t> ara_id;

std::cout<<"Start testing:"<<"\n";

    auto start = std::chrono::high_resolution_clock::now();

    //ap.make_matrix("/Users/au383883/Documents/MY/codebase/evo/evolm/tests/data/large_data/DMU/data/dmu_pedigree_yy_20240125.txt", true, true); // full A(-1)

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout <<"full A(-1) (milliseconds): "<< duration.count() << std::endl;

        //ap.make_matrix("tests/data/ped_632.dat", true); // full A(-1)
        //ap.make_matrix("tests/data/sstep_050/data/id4trace.PED", true); // full A(-1)
    //ap.get_matrix("iA", iA_s, a_id, false);
        //iA.symtorec();
    //iA_s.print("Full A(-1)");

    std::cout<<"n IDs in A(-1): "<<a_id.size()<<"\n";
    //for (auto &e: a_id)
    //    std::cout<<e<<" ";
    std::cout<<"\n";

    /*    std::cout<<"making reduced A(-1)"<<"\n";ll
        
        ap.make_matrix("tests/data/ped_bkg.dat", "tests/data/typed2", true); // reduced A(-1)
        //ap.make_matrix("tests/data/sstep_050/data/id4trace.PED", "tests/data/sstep_050/data/typed_050.dat", true); // reduced A(-1)
        ap.get_matrix("irA", irA, ra_id);
        irA.print("Reduced A(-1)");
    */
std::cout<<"making all:"<<"\n";

    start = std::chrono::high_resolution_clock::now();

    /*ap.make_all("/Users/au383883/Documents/MY/codebase/evo/evolm/tests/data/large_data/DMU/data/dmu_pedigree_yy_20240125.txt",
                "/Users/au383883/Documents/MY/codebase/evo/evolm/tests/data/large_data/DMU/YY/yy.grm.id",
                true);*/
    ap.make_all("../evolm/tests/data/large_data/DMU/data/dmu_pedigree_yy_20240125.txt",
                "../evolm/tests/data/large_data/DMU/YY/yy.grm.id",
                true);

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout <<"All (milliseconds): "<< duration.count() << std::endl;

    //ap.make_all("tests/data/sstep_050/data/id4trace.PED", "tests/data/sstep_050/data/typed_050.dat", false); // relatively big data
    //ap.make_all("tests/data/ped_bkg.dat", "tests/data/typed2", false); // reduced A(-1)

std::cout<<"getting iA"<<"\n";
    ap.get_matrix("iA", iA_sf, a_id, false);
std::cout<<"getting irA"<<"\n";
    ap.get_matrix("irA", irA_sf, ra_id, false);
std::cout<<"getting iA22"<<"\n";
    ap.get_matrix("iA22", iA22f, ai22_id, false);
std::cout<<"getting A22"<<"\n";
    ap.get_matrix("A22", A22f, a22_id, false);

    std::cout<<"n IDs in A(-1): "<<a_id.size()<<"\n";
    //for (auto &e: a_id)
    //    std::cout<<e<<" ";
    std::cout<<"\n";

    std::cout<<"n IDs in red A(-1): "<<ra_id.size()<<"\n";
    //for (auto &e: ra_id)
    //    std::cout<<e<<" ";
    std::cout<<"\n";

    std::cout<<"n IDs in A22(-1): "<<ai22_id.size()<<"\n";
    //for (auto &e: ai22_id)
    //    std::cout<<e<<" ";
    std::cout<<"\n";

    std::cout<<"n IDs in A22: "<<a22_id.size()<<"\n";
    //for (auto &e: a22_id)
    //    std::cout<<e<<" ";
    std::cout<<"\n";

exit(1);
    evoped::Gmat gmat;
    evolm::matrix<double> iG;
    std::vector<std::int64_t> g_id;

    std::cout<<"reading G"<<"\n";

    //gmat.make_matrix("tests/data/allele.dat");
    //gmat.read_matrix("tests/data/g_mat");
    gmat.read_matrix("tests/data/sstep_050/data/gmat_050.dat");

    std::cout<<"Scaling G matrix"<<"\n";

    gmat.scale_matrix(A22, 0.25);

#ifdef UTEST
    double a, b, a_diag, a_ofd, g_diag, g_ofd;

    gmat.get_alpha_beta(a,b,a_diag,a_ofd,g_diag,g_ofd);

    std::cout<<"alpha: "<<a<<"\n";
    std::cout<<"beta: "<<b<<"\n";
    std::cout<<"a_diag: "<<a_diag<<"\n";
    std::cout<<"a_ofd: "<<a_ofd<<"\n";
    std::cout<<"g_diag: "<<g_diag<<"\n";
    std::cout<<"g_ofd: "<<g_ofd<<"\n";
#endif

    std::cout<<"Inverting G matrix"<<"\n";

    gmat.invert_matrix(true);
    // sparse inverse:
    //std::vector<std::int64_t> core_id{ 18, 20, 22, 25 };
    //gmat.invert_matrix(core_id);

    gmat.get_matrix(iG,g_id);
    //iG.print("sparse G(-1)");

    std::cout<<"n IDs in G: "<<g_id.size()<<"\n";
    std::cout<<"size of G: "<<iG.size()<<"\n";
    //for (auto &e: g_id)
    //    std::cout<<e<<" ";
    std::cout<<"\n";


    evoped::Hmat h;
    evolm::matrix<double> H;
    std::vector<std::int64_t> h_id;
    h.make_matrix(iA, a_id, iA22, ai22_id, iG, g_id);
    
    h.get_matrix(H,h_id,false);
    //H.print("H(-1) from run");

    std::cout<<"n IDs in H(-1): "<<h_id.size()<<"\n";
    //for (auto &e: h_id)
    //    std::cout<<e<<" ";
    std::cout<<"\n";

    std::cout<<"\n";
    std::cout<<"Done with Hmat."<<"\n";
    std::cout<<"\n";

    return 0;
}