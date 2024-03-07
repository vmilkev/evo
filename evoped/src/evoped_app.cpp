#include <iostream>

#include "Amat.hpp"
#include "Gmat.hpp"
#include "Hmat.hpp"

int main()
{

    evoped::Amat ap;
    evolm::matrix<double> iA;
    std::vector<std::int64_t> a_id;
    evolm::matrix<double> irA;
    std::vector<std::int64_t> ra_id;
    evolm::matrix<double> iA22;
    std::vector<std::int64_t> ai22_id;
    evolm::matrix<double> A22;
    std::vector<std::int64_t> a22_id;

    evolm::matrix<double> A;
    std::vector<std::int64_t> aa_id;
    evolm::matrix<double> rA;
    std::vector<std::int64_t> ara_id;


    ap.make_matrix("tests/data/ped_632.dat", true); // full A(-1)
    //ap.make_matrix("tests/data/sstep_050/data/id4trace.PED", true); // full A(-1)
    ap.get_matrix("iA", iA, a_id, false);
    iA.symtorec();
    iA.print("Full A(-1)");

    std::cout<<"n IDs in A(-1): "<<a_id.size()<<"\n";
    for (auto &e: a_id)
        std::cout<<e<<" ";
    std::cout<<"\n";
exit(1);
/*    std::cout<<"making reduced A(-1)"<<"\n";ll
    
    ap.make_matrix("tests/data/ped_bkg.dat", "tests/data/typed2", true); // reduced A(-1)
    //ap.make_matrix("tests/data/sstep_050/data/id4trace.PED", "tests/data/sstep_050/data/typed_050.dat", true); // reduced A(-1)
    ap.get_matrix("irA", irA, ra_id);
    irA.print("Reduced A(-1)");
*/
std::cout<<"making all:"<<"\n";

    ap.make_all("tests/data/sstep_050/data/id4trace.PED", "tests/data/sstep_050/data/typed_050.dat"); // relatively big data
    //ap.make_all("tests/data/ped_bkg.dat", "tests/data/typed2"); // reduced A(-1)

std::cout<<"getting iA"<<"\n";
    ap.get_matrix("iA", iA, a_id, false);
std::cout<<"getting irA"<<"\n";
    ap.get_matrix("irA", irA, ra_id, false);
std::cout<<"getting iA22"<<"\n";
    ap.get_matrix("iA22", iA22, ai22_id, false);
std::cout<<"getting A22"<<"\n";
    ap.get_matrix("A22", A22, a22_id, false);

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


    evoped::Gmat gmat;
    evolm::matrix<double> iG;
    std::vector<std::int64_t> g_id;

    std::cout<<"reading G"<<"\n";

    //gmat.make_matrix("tests/data/allele.dat");
    //gmat.read_matrix("tests/data/g_mat");
    gmat.read_matrix("tests/data/sstep_050/data/gmat_050.dat");

    std::cout<<"Scaling G matrix"<<"\n";

    gmat.scale_matrix(A22, 0.25);

    double a, b, a_diag, a_ofd, g_diag, g_ofd;

    gmat.get_alpha_beta(a,b,a_diag,a_ofd,g_diag,g_ofd);

    std::cout<<"alpha: "<<a<<"\n";
    std::cout<<"beta: "<<b<<"\n";
    std::cout<<"a_diag: "<<a_diag<<"\n";
    std::cout<<"a_ofd: "<<a_ofd<<"\n";
    std::cout<<"g_diag: "<<g_diag<<"\n";
    std::cout<<"g_ofd: "<<g_ofd<<"\n";

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