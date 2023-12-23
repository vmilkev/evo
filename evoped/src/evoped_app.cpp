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

    evoped::Gmat g;

    std::vector<std::int64_t> g_row, g_col;
    std::vector<double> val;
    g.read_gmatrix("tests/data/g_mat", g_row, g_col, val);

    for (size_t i = 0; i < g_row.size(); i++)
    {
        std::cout<<g_row[i]<<" "<<g_col[i]<<" "<<val[i]<<"\n";
    }
    
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