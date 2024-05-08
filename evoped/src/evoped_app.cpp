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
    evolm::smatrix<float> iA22_sf;
    std::vector<std::int64_t> ai22_id;
    evolm::matrix<double> A22;
    evolm::matrix<float> A22f;
    evolm::smatrix<float> A22_sf;
    std::vector<std::int64_t> a22_id;

    evolm::matrix<double> A;
    std::vector<std::int64_t> aa_id;
    evolm::matrix<double> rA;
    std::vector<std::int64_t> ara_id;

std::cout<<"Start testing:"<<"\n";

    auto start = std::chrono::high_resolution_clock::now();

    //ap.make_matrix("/Users/au383883/Documents/MY/codebase/evo/evolm/tests/data/large_data/DMU/data/dmu_pedigree_yy_20240125.txt", true); // full A(-1)

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout <<"full A(-1) (milliseconds): "<< duration.count() << std::endl;

        //ap.make_matrix("tests/data/ped_632.dat", true); // full A(-1)
        //ap.make_matrix("tests/data/sstep_050/data/id4trace.PED", true); // full A(-1)
    //ap.get_matrix("iA", iA_s, a_id);
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
                "/Users/au383883/Documents/MY/codebase/evo/evolm/tests/data/large_data/DMU/YY/yy.grm.id");*/
    ap.make_all("../evolm/tests/data/large_data/DMU/data/dmu_pedigree_yy_20240125.txt",
                "../evolm/tests/data/large_data/DMU/YY/yy.grm.id");

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout <<"completed (milliseconds): "<< duration.count() << std::endl;

    //ap.make_all("tests/data/sstep_050/data/id4trace.PED", "tests/data/sstep_050/data/typed_050.dat"); // relatively big data
    //ap.make_all("tests/data/ped_bkg.dat", "tests/data/typed2"); // reduced A(-1)

std::cout<<"getting iA"<<"\n";
    ap.get_matrix("iA", iA_sf, a_id);
std::cout<<"getting irA"<<"\n";
    ap.get_matrix("irA", irA_sf, ra_id);
std::cout<<"getting iA22"<<"\n";
    ap.get_matrix("iA22", iA22_sf, ai22_id);
std::cout<<"getting A22"<<"\n";
    ap.get_matrix("A22", A22f, a22_id);

    std::cout<<"reading iA_sf"<<"\n";
    iA_sf.fread();
    std::cout<<"n IDs in A(-1): "<<a_id.size()<<", values in matrix: "<<iA_sf.size()<<"used memory of sparse A22: "<<(double)(iA_sf.size()*sizeof(float))/(1024*1024)<<" MB + "<<(double)(iA_sf.size()*sizeof(size_t))/(1024*1024)<<" MB"<<"\n";
    std::cout<<"other calcul. of memory usage: "<<(double)iA_sf.get_memory_usage()/(1024*1024)<<" MB"<<"\n";
    //for (auto &e: a_id)
    //    std::cout<<e<<" ";
    std::cout<<"\n";

    std::cout<<"reading irA_sf"<<"\n";
    irA_sf.fread();
    std::cout<<"n IDs in red A(-1): "<<ra_id.size()<<", values in matrix: "<<irA_sf.size()<<"used memory of sparse A22: "<<(double)(irA_sf.size()*sizeof(float))/(1024*1024)<<" MB + "<<(double)(irA_sf.size()*sizeof(size_t))/(1024*1024)<<" MB"<<"\n";
    std::cout<<"other calcul. of memory usage: "<<(double)irA_sf.get_memory_usage()/(1024*1024)<<" MB"<<"\n";
    //for (auto &e: ra_id)
    //    std::cout<<e<<" ";
    std::cout<<"\n";

    std::cout<<"reading iA22_sf"<<"\n";
    iA22_sf.fread();
    std::cout<<"n IDs in A22(-1): "<<ai22_id.size()<<", values in matrix: "<<iA22_sf.size()<<"used memory of sparse A22: "<<(double)(iA22_sf.size()*sizeof(float))/(1024*1024)<<" MB + "<<(double)(iA22_sf.size()*sizeof(size_t))/(1024*1024)<<" MB"<<"\n";
    std::cout<<"other calcul. of memory usage: "<<(double)iA22_sf.get_memory_usage()/(1024*1024)<<" MB"<<"\n";
    //for (auto &e: ai22_id)
    //    std::cout<<e<<" ";
    std::cout<<"\n";

    std::cout<<"reading (dense) A22f"<<"\n";
    A22f.fread();
    std::cout<<"n IDs in A22: "<<a22_id.size()<<"\n";
    std::cout<<"other calcul. of memory usage: "<<(double)a22_id.size()*a22_id.size()/(1024*1024)<<" MB"<<"\n";
    //for (auto &e: a22_id)
    //    std::cout<<e<<" ";
    std::cout<<"\n";

    iA_sf.clean();
    irA_sf.clean();
    iA22_sf.clean();
    A22_sf.clean();

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