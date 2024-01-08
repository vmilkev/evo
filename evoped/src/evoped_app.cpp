#include <iostream>

#include "Amat.hpp"
#include "Gmat.hpp"
#include "Hmat.hpp"

int main()
{
    std::cout<<"Here main."<<"\n";
    
    evoped::Amat p;
    evoped::Amat p2;
    std::cout<<"reading pedigree & genotyped files:"<<"\n";

    evolm::matrix <double> A;
    
    p2.make_matrix("tests/data/ped.dat", true);
    p2.get_matrix(A);
    A.symtorec();
    A.print("A1(-1)");

    p2.make_matrix("tests/data/ped.dat", false);
    p2.get_matrix(A);
    A.symtorec();
    A.print("A1");

    p2.make_matrix("tests/data/ped_bkg.dat", "tests/data/typed0", true);
    p2.get_matrix(A);
    A.symtorec();
    A.print("A2(-1)");

    p2.make_matrix("tests/data/ped_bkg.dat", "tests/data/typed0", false);
    p2.get_matrix(A);
    A.symtorec();
    A.print("A2");

    p2.make_matrix("tests/data/ped2.dat", true);  
    p2.get_matrix(A);
    A.symtorec();
    A.print("A3(-1)");

    p2.make_matrix("tests/data/ped2.dat", false);  
    p2.get_matrix(A);
    A.symtorec();
    A.print("A3");

    std::cout<<"\n";
    std::cout<<"Done with Amat."<<"\n";
    std::cout<<"\n";

    evoped::Gmat g;

    //g.read_matrix("tests/data/g_mat2");
    //g.make_matrix("tests/data/obs_1000_snp_1000.txt"); 
    g.make_matrix("tests/data/snps.dat","tests/data/ids.dat");
    
    evolm::matrix <double> G;
    evolm::matrix <double> G2;
    evolm::matrix <double> I;
    std::vector<std::int64_t> gids;

    g.get_matrix(G);
    G.symtorec();
    g.get_ids(gids);

    std::cout<<"address of outside G: "<<std::addressof(G)<<", capacity of G: "<<G.capacity()<<", size of G: "<<G.size()<<", is empty(): "<<G.empty()<<"\n";
    //std::cout<<"address of outside G2: "<<std::addressof(G2)<<", capacity of G2: "<<G2.capacity()<<", size of G2: "<<G2.size()<<", is empty(): "<<G2.empty()<<"\n";

    G.print("G-matrix");

    /*std::cout<<"\n";
    std::cout<<"IDs of the matrix:"<<"\n";
    for (size_t i = 0; i < gids.size(); i++)
    {
        std::cout<<gids[i]<<"\n";
    }
    std::cout<<"\n";*/

    g.scale_matrix(1.05);
    g.invert_matrix(true);
    G.clear();
    g.get_matrix(G);

    G.print("Inverted full-stored G-matrix with scalling");

    std::cout<<"after inv. address of outside G: "<<std::addressof(G)<<", after inv. capacity of G: "<<G.capacity()<<", size of G: "<<G.size()<<"\n";
    // --------------------------------------------
    //G.clear();
    g.make_matrix("tests/data/snps.dat","tests/data/ids.dat");
    g.get_matrix(G);

    G.print("New G-matrix");

    g.scale_matrix(1.05);

    g.get_matrix(G);

    G.print("Scaled G-matrix");
    
    G2 = G;
    G2.print("G2 == G");

    g.invert_matrix(false);
    //G.clear();
    g.get_matrix(G);

    G.print("Inverted scaled G-matrix");

    G = G * (-1.0);

    G.print("G = G * (-1.0)");

    G2.symtorec();
    G.symtorec();

    I = G2 * G;

    I.print("expected the identity matrix");

    std::cout<<"\n";
    std::cout<<"Done with Gmat."<<"\n";
    std::cout<<"\n";

    evoped::Amat ap;
    evolm::matrix<double> iA;
    std::vector<std::int64_t> a_id;
    evolm::matrix<double> irA;
    std::vector<std::int64_t> ra_id;

    ap.make_matrix("tests/data/ped_bkg.dat", true); // full A(-1)
    ap.get_matrix(iA);
    ap.get_ids(a_id);
    iA.print("A(-1)");

    ap.make_matrix("tests/data/ped_bkg.dat", "tests/data/typed2", true); // reduced A(-1)
    ap.get_matrix(irA);
    ap.get_ids(ra_id);
    irA.print("Ared(-1)");

    evoped::Gmat gmat;
    evolm::matrix<double> iG;
    std::vector<std::int64_t> g_id;

    gmat.make_matrix("tests/data/allele.dat");

    gmat.get_ids(g_id);

    // Extract A22 for scaling G matrix
    evoped::Hmat h;
    evolm::matrix<double> A22;
    h.get_submatrix(iA,a_id,g_id,false);
    h.get_matrix(A22);
    A22.print("A22");

    gmat.scale_matrix(A22, 0.25);   
    gmat.get_matrix(iG);
    iG.print("G(-1)");

    evolm::matrix<double> H;
    h.make_matrix(iA, a_id, irA, ra_id, iG, g_id);
    h.get_matrix(H);
    H.print("H(-1)");

    iA.print("A(-1)");
    irA.print("Ared(-1)");
    iG.print("G(-1)");
    H.print("H(-1)");

    std::cout<<"\n";
    std::cout<<"Done with Hmat."<<"\n";
    std::cout<<"\n";

    return 0;
}