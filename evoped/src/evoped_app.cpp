#include <iostream>

#include "Amat.hpp"
#include "Gmat.hpp"

int main()
{
    std::cout<<"Here main."<<"\n";
    
    evoped::Amat p("tests/data/ped_bkg.dat", "tests/data/typed2");
    evoped::Amat p2;
    std::cout<<"reading pedigree & genotyped files:"<<"\n";
    
    p.get_ainv();
    p2.get_ainv("tests/data/ped.dat");
    p2.get_ainv("tests/data/ped_bkg.dat", "tests/data/typed2");
    p2.get_ainv();

    //evoped::Gmat g;
    evoped::Gmat *g = new evoped::Gmat;

    g->read_matrix("tests/data/g_mat");    
    //g.make_matrix("tests/data/allele.dat");
    
    evolm::matrix <double> G;
    evolm::matrix <double> G2;

    g->get_matrix(G);

    std::cout<<"address of outside G: "<<std::addressof(G)<<", capacity of G: "<<G.capacity()<<", size of G: "<<G.size()<<", is empty(): "<<G.empty()<<"\n";
    //std::cout<<"address of outside G2: "<<std::addressof(G2)<<", capacity of G2: "<<G2.capacity()<<", size of G2: "<<G2.size()<<", is empty(): "<<G2.empty()<<"\n";

    std::cout<<"\n";
    std::cout<<"G matrix:"<<"\n";
    for (size_t i = 0; i < 9; i++)
    {
        for (size_t j = 0; j < 9; j++)
        {
            std::cout<<G(i,j)<<" ";
            //G2(i,j) = G(i,j);
        }
        std::cout<<"\n";
    }
    std::cout<<"\n";

    g->invert_matrix();
    G.clear();
    g->get_matrix(G);

    std::cout<<"after inv. address of outside G: "<<std::addressof(G)<<", after inv. capacity of G: "<<G.capacity()<<", size of G: "<<G.size()<<"\n";
    std::cout<<"\n";
    std::cout<<"inverse of G matrix:"<<"\n";
    for (size_t i = 0; i < 9; i++)
    {
        for (size_t j = 0; j < 9; j++)
            std::cout<<G(i,j)<<" ";
        std::cout<<"\n";
    }
    std::cout<<"\n";

exit(1);

    double m[] = {
     1.0447,    -0.53818,     0.63471,    -0.32742,     0.26935,    -0.39388,    -0.83894,    -0.48137,     0.55758,
    -0.53818,       1.547,    -0.97653,    -0.28472,     -0.9009,      0.2452,     0.73787,    0.062788,     0.19598,
     0.63471,    -0.97653,      1.1909,     0.11579,     0.54964,    -0.45469,     -0.7264,    -0.21647,    -0.18195,
    -0.32742,    -0.28472,     0.11579,     0.63878,    0.066364,    0.019099,     0.41955,     0.16539,    -0.81458,
     0.26935,     -0.9009,     0.54964,    0.066364,      1.2477,    -0.28404,    -0.44515,   -0.057478,    -0.51011,
    -0.39388,      0.2452,    -0.45469,    0.019099,    -0.28404,     0.66253,      0.4145,     0.13643,    -0.28844,
    -0.83894,     0.73787,     -0.7264,     0.41955,    -0.44515,      0.4145,      1.0772,     0.16539,    -0.74151,
    -0.48137,    0.062788,    -0.21647,     0.16539,   -0.057478,     0.13643,     0.16539,     0.65316,    -0.37649,
     0.55758,     0.19598,    -0.18195,    -0.81458,    -0.51011,    -0.28844,    -0.74151,    -0.37649,      2.1165
    };

    double m2[81];
    G2.resize(9,9);
    for (size_t i = 0; i < 9*9; i++)
    {
        G2[i] = m[i];
        //m2[i] = G[i];
        //G2[i] = m2[i];
    }

    std::cout<<"G2 matrix:"<<"\n";
    for (size_t i = 0; i < 9; i++)
    {
        for (size_t j = 0; j < 9; j++)
        {
            std::cout<<G2(i,j)<<" ";
        }
        std::cout<<"\n";
    }
    std::cout<<"\n";
    
    G2.invert();

    std::cout<<"after the inverse of G2 matrix:"<<"\n";
    for (size_t i = 0; i < 9; i++)
    {
        for (size_t j = 0; j < 9; j++)
        {
            std::cout<<G2(i,j)<<" ";
        }
        std::cout<<"\n";
    }
    std::cout<<"\n";

    /*for (size_t i = 0; i < g_row.size(); i++)
    {
        std::cout<<g_row[i]<<" "<<g_col[i]<<" "<<val[i]<<"\n";
    }*/
    

    /*
    std::cout<<"Full pedigree."<<"\n";
    for (auto const &e: p.pedigree)
        std::cout<<e.first.val_1<<" "<<e.first.val_2<<" "<<e.second.val_1<<" "<<e.second.val_2<<"\n";

    std::cout<<"Reduced pedigree."<<"\n";
    for (auto const &e: p.r_pedigree)
        std::cout<<e.first.val_1<<" "<<e.first.val_2<<" "<<e.second.val_1<<" "<<e.second.val_2<<"\n";

    std::cout<<"A(-1)."<<"\n";

    std::cout<<"Full A(-1)."<<"\n";
    for (auto const &e: p.ainv)
        std::cout<<e.first.val_1<<" "<<e.first.val_2<<" "<<e.second<<"\n";

    std::cout<<"Reduced A(-1)."<<"\n";
    for (auto const &e: p.r_ainv)
        std::cout<<e.first.val_1<<" "<<e.first.val_2<<" "<<e.second<<"\n";
*/
    std::cout<<"Done."<<"\n";

    return 0;
}