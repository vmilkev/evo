#include "Amat.hpp"

namespace evoped
{
    //===============================================================================================================

    Amat::Amat()
    {
        try
        {
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::Amat()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::Amat()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::Amat()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    Amat::Amat(const std::string &ped_file)
    {
        try
        {
            pedigree_file = ped_file;

            fread_pedigree(ped_file);
            trace_pedigree();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::Amat(std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::Amat(std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::Amat(std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    Amat::Amat(const std::string &ped_file, const std::string &g_file)
    {
        try
        {
            pedigree_file = ped_file;
            genotyped_file = g_file;

            fread_pedigree(ped_file);
            fread_genotyped_id(g_file);
            trace_pedigree();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::Amat(std::string &, std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::Amat(std::string &, std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::Amat(std::string &, std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::get_ainv()
    {
        try
        {
            if ( !pedigree.empty() )
            {
                get_ainv(pedigree, ainv, true);
                pedigree.clear();
            }
            else if ( !pedigree_file.empty() )
            {
                if ( !genotyped_file.empty() )
                {
                    fread_pedigree(pedigree_file);
                    fread_genotyped_id(genotyped_file);
                    trace_pedigree();
                    get_ainv(pedigree, ainv, true);
                    get_ainv(r_pedigree, r_ainv, true);
                    pedigree.clear();
                    r_pedigree.clear();
                    return;
                }
                else
                {
                    fread_pedigree(pedigree_file);
                    trace_pedigree();
                    get_ainv(pedigree, ainv, true);
                    pedigree.clear();
                }
            }
            else
                throw std::string("The pedigree is empty!");

            if ( !r_pedigree.empty() )
            {
                get_ainv(r_pedigree, r_ainv, true);
                r_pedigree.clear();
            }
            else if ( !genotyped_file.empty() )
            {
                fread_genotyped_id(genotyped_file);
                if ( pedigree.empty() && !pedigree_file.empty() )
                    fread_pedigree(pedigree_file);
                else
                    return;
                
                trace_pedigree();
                get_ainv(r_pedigree, r_ainv, true);
                r_pedigree.clear();
            }

        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_ainv()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_ainv()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_ainv()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::get_ainv(const std::string &ped_file)
    {
        try
        {
            pedigree_file = ped_file;

            fread_pedigree(ped_file);

            trace_pedigree();

            get_ainv(pedigree, ainv, true);

            pedigree.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_ainv(std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_ainv(std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_ainv(std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::get_ainv(const std::string &ped_file, const std::string &g_file)
    {
        try
        {
            pedigree_file = ped_file;
            genotyped_file = g_file;

            fread_pedigree(ped_file);
            fread_genotyped_id(g_file);

            trace_pedigree();

            get_ainv(pedigree, ainv, true);
            get_ainv(r_pedigree, r_ainv, true);

            pedigree.clear();
            r_pedigree.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_ainv(std::string &, std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_ainv(std::string &, std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_ainv(std::string &, std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::fread_pedigree(const std::string &ped_file)
    {
        std::ifstream ped;

        try
        {
            std::string line;
            PedPair key;
            PedPair id_pair;

            char *end;
            const char *p;
            std::vector<double> tmp_list;

            ped.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            ped.open(ped_file, std::fstream::in);

            if (ped)
            {
                while (getline(ped, line))
                {
                    p = line.c_str();
                    for (double f = std::strtod(p, &end); p != end; f = std::strtod(p, &end))
                    {
                        p = end;
                        if (errno == ERANGE)
                        {
                            errno = 0;
                            throw std::string("Range error during reading pedigree file");
                        }
                        tmp_list.push_back(f);
                    }

                    if (tmp_list.size() > 0)
                    {
                        key.val_1 = static_cast<std::int64_t>(tmp_list[3]);
                        key.val_2 = static_cast<std::int64_t>(tmp_list[0]);
                        id_pair.val_1 = static_cast<std::int64_t>(tmp_list[1]);
                        id_pair.val_2 = static_cast<std::int64_t>(tmp_list[2]);
                        pedigree_from_file[key] = id_pair;
                        birth_id_map[static_cast<std::int64_t>(tmp_list[0])] = static_cast<std::int64_t>(tmp_list[3]);
                        pedID.push_back(static_cast<std::int64_t>(tmp_list[0]));

                        tmp_list.erase(tmp_list.begin(), tmp_list.end());
                    }

                    if (ped.eof())
                    {
                        ped.close();
                        if (!is_unique(pedID))
                            pedID.erase(unique(pedID.begin(), pedID.end()), pedID.end()); // here the vector should be sorted and unique
                        return;
                    }
                }

                ped.close();
                if (!is_unique(pedID))
                    pedID.erase(unique(pedID.begin(), pedID.end()), pedID.end()); // here the vector should be sorted and unique
            }
            else
                throw std::string("Cannot open pedigree file!");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::fread_pedigree(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            if (ped.eof())
                std::cerr << "Check for empty line(s) at the endd of the file!"
                          << "\n";
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::fread_pedigree(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::fread_pedigree(const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    bool Amat::is_unique(std::vector<std::int64_t> &x)
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
            std::cerr << "Exception in Amat::is_unique(std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::is_unique(std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::is_unique(std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::fread_genotyped_id(const std::string &g_file)
    {
        std::ifstream ped;

        try
        {
            /* reads genotyped and core IDs from typed file into the vectors: genotyped, core */

            std::string line;
            std::int64_t t_gtyp, t_gcor;
            t_gtyp = t_gcor = 0;

            char *end;
            const char *p;
            std::vector<double> tmp_list;

            ped.exceptions(std::ifstream::failbit | std::ifstream::badbit);
            ped.open(g_file, std::fstream::in);

            if (ped)
            {
                while (std::getline(ped, line))
                {
                    p = line.c_str();

                    for (double f = std::strtod(p, &end); p != end; f = std::strtod(p, &end))
                    {
                        p = end;
                        if (errno == ERANGE)
                        {
                            errno = 0;
                            throw std::string("Range error during reading genotyped IDs file");
                        }
                        tmp_list.push_back(f);
                    }

                    if (tmp_list.size() > 0)
                    {
                        if (static_cast<std::int64_t>(tmp_list[0]) != 0 && static_cast<std::int64_t>(tmp_list[0]) != t_gtyp)
                            genotypedID.push_back(static_cast<std::int64_t>(tmp_list[0]));

                        if (static_cast<std::int64_t>(tmp_list[1]) != 0 && static_cast<std::int64_t>(tmp_list[0]) != t_gtyp)
                            coreID.push_back(static_cast<std::int64_t>(tmp_list[0]));

                        t_gtyp = static_cast<std::int64_t>(tmp_list[0]);
                        t_gcor = static_cast<std::int64_t>(tmp_list[1]);

                        tmp_list.erase(tmp_list.begin(), tmp_list.end());
                    }

                    if (ped.eof())
                    {
                        ped.close();

                        if (!is_unique(genotypedID))
                            genotypedID.erase(unique(genotypedID.begin(), genotypedID.end()), genotypedID.end());

                        if (!is_unique(coreID))
                            coreID.erase(unique(coreID.begin(), coreID.end()), coreID.end());

                        return;
                    }
                }

                ped.close();

                if (!is_unique(genotypedID))
                    genotypedID.erase(unique(genotypedID.begin(), genotypedID.end()), genotypedID.end());

                if (!is_unique(coreID))
                    coreID.erase(unique(coreID.begin(), coreID.end()), coreID.end());

                return;
            }
            else
                throw std::string("Cannot open genotyped ids file!");
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::fread_genotyped_id(const std::string &)" << '\n';
            std::cerr << e.what() << '\n';
            if (ped.eof())
                std::cerr << "Check for empty line(s) at the endd of the file!"
                          << "\n";
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::fread_genotyped_id(const std::string &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::fread_genotyped_id(const std::string &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    bool Amat::is_invect(std::vector<std::int64_t> &where, std::int64_t what)
    {
        try
        {
            std::vector<std::int64_t>::iterator it;
            it = find(where.begin(), where.end(), what);
            if (it != where.end())
                return true;
            else
                return false;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::is_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::is_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::is_invect(std::vector<std::int64_t> &, std::int64_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::trace_pedigree(std::map<PedPair, PedPair> &in_ped, std::map<PedPair, PedPair> &out_ped, std::vector<std::int64_t> &traced_id)
    {
        try
        {
            PedPair id_pair;
            PedPair key;
            std::vector<std::int64_t> gen_pedID(traced_id);
            std::int64_t elem_v;
            size_t iter_v = 0;
            bool exists;

            std::vector<std::int64_t> days;
            std::vector<std::int64_t> ids;
            std::vector<std::int64_t> sire;
            std::vector<std::int64_t> dame;

            for (auto const &elem_m : in_ped)
            {
                days.push_back(elem_m.first.val_1);
                ids.push_back(elem_m.first.val_2);
                sire.push_back(elem_m.second.val_1);
                dame.push_back(elem_m.second.val_2);
            }

            do
            {
                elem_v = gen_pedID[iter_v];

                exists = false;

                for (size_t i = 0; i < days.size(); i++)
                {

                    if (ids[i] == elem_v)
                    {

                        exists = true;
                        id_pair.val_1 = sire[i];
                        id_pair.val_2 = dame[i];

                        if ((birth_id_map[sire[i]] >= days[i]) || (birth_id_map[dame[i]] >= days[i]))
                            throw std::string("Amat is not correct: parents born before offspring!");

                        key.val_1 = days[i];
                        key.val_2 = elem_v;
                        out_ped[key] = id_pair;

                        if (id_pair.val_1 != 0)
                        {
                            if (!is_invect(gen_pedID, id_pair.val_1))
                                gen_pedID.push_back(id_pair.val_1);
                        }
                        if (id_pair.val_2 != 0)
                        {
                            if (!is_invect(gen_pedID, id_pair.val_2))
                                gen_pedID.push_back(id_pair.val_2);
                        }

                        break;
                    }
                }

                if (!exists)
                {
                    key.val_1 = -1;
                    key.val_2 = elem_v;
                    id_pair.val_1 = 0;
                    id_pair.val_2 = 0;
                    out_ped[key] = id_pair;
                    exists = false;
                }

                iter_v++;

            } while (iter_v < gen_pedID.size());

            gen_pedID.clear();
            days.clear();
            ids.clear();
            sire.clear();
            dame.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::trace_pedigree(std::map<PedPair, PedPair> &, std::map<PedPair, PedPair> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::trace_pedigree(std::map<PedPair, PedPair> &, std::map<PedPair, PedPair> &, std::vector<std::int64_t> &)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::trace_pedigree(std::map<PedPair, PedPair> &, std::map<PedPair, PedPair> &, std::vector<std::int64_t> &)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::trace_pedigree()
    {
        // requires fread_pedigree(file) call first, and, possibly, fread_genotyped_id(file) call as well
        try
        {
            if (!pedigree_from_file.empty())
            {
                if (pedID.empty())
                    throw std::string("Empty pedigree IDs.");

                if (!pedigree.empty())
                    pedigree.clear();

                trace_pedigree(pedigree_from_file, pedigree, pedID); // tracing full pedigree

                if (!genotypedID.empty())
                {
                    if (!r_pedigree.empty())
                        r_pedigree.clear();

                    trace_pedigree(pedigree_from_file, r_pedigree, genotypedID); // tracing reduced pedigree for genotyped individuals
                }
            }
            else
                throw std::string("There is no base pedigree which needs to be traced.");

            pedigree_from_file.clear();
            birth_id_map.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::trace_pedigree()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::trace_pedigree()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::trace_pedigree()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::int64_t Amat::pos_inped(std::map<std::int64_t, std::int64_t> &codemap, std::int64_t id)
    {
        try
        {
            auto pos = codemap.find(id);

            if (pos == codemap.end())
                return 0;
            else
                return pos->second;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::pos_inped(std::map<std::int64_t, std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::pos_inped(std::map<std::int64_t, std::int64_t> &, std::int64_t)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::pos_inped(std::map<std::int64_t, std::int64_t> &, std::int64_t)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::get_dinv(std::map<PedPair, PedPair> &ped, std::vector<double> &dinv, bool inbreed)
    {
        // diagonal elements of A(-1)
        // Algorithm: M.Sargolzaei et al. 'A fast algorithm for computing inbreeding ...' (works correct!)
        try
        {
            if (inbreed)
            {
                size_t n = ped.size();
                size_t m = n + 1;

                std::vector<std::vector<std::int64_t>> Ped((m), std::vector<std::int64_t>(2, 0)); // main pedigree, defaults to zero initial value

                std::map<std::int64_t, std::int64_t> code_map;
                size_t code_id = 1;
                inbrF.push_back(0.0); // we shall be able to get the coeff. starting from index '1'

                for (auto const &elem : ped)
                {

                    code_map[elem.first.val_2] = code_id;

                    Ped[code_id][0] = pos_inped(code_map, elem.second.val_1);
                    Ped[code_id][1] = pos_inped(code_map, elem.second.val_2);

                    inbrF.push_back(0.0); // initialize vector of inbreeding coeff.

                    code_id++;
                }

                size_t i, j, k, rN, rS, S, D, MIP;

                std::vector<std::vector<std::int64_t>> rPed((m + 1), std::vector<std::int64_t>(2, 0)); // reduced pedigree
                std::vector<std::int64_t> SId(m, 0);                                                   // will contain the sorted animals ID based on the ID of their sires
                std::vector<std::int64_t> Link(m, 0);                                                  // will contain new ID of ancestors at position of their original ID
                std::vector<std::int64_t> MaxIdP(m + 1, 0);                                            // will contain maximum new ID of parents for each paternal group at position of the new ID of each sire

                std::vector<double> F((n + 1), 0.0); // inbreeding coefficients
                std::vector<double> B((m + 1), 0.0); // within family segregation variances
                std::vector<double> x((m + 1), 0.0); // x arrays

                // set values for the unknown parent
                F[0] = -1.0;
                x[0] = 0.0;
                Link[0] = 0;
                rN = 1;
                i = 1;

                for (; i <= n; ++i)
                { // extract and recode ancestors

                    SId[i] = i;
                    Link[i] = 0;
                    if (i <= m)
                        x[i] = 0.0;

                    S = Ped[i][0];
                    D = Ped[i][1];

                    if (S && !Link[S])
                    {
                        MaxIdP[rN] = Link[S] = rN;
                        rPed[rN][0] = Link[Ped[S][0]];
                        rPed[rN++][1] = Link[Ped[S][1]];
                    }
                    if (D && !Link[D])
                    {
                        /*MaxIdP[rN] =*/Link[D] = rN;
                        rPed[rN][0] = Link[Ped[D][0]];
                        rPed[rN++][1] = Link[Ped[D][1]];
                    }
                    if (MaxIdP[Link[S]] < Link[D])
                        MaxIdP[Link[S]] = Link[D]; // determine maximum ID of parents for each paternal group
                }

                // sort animals according to ID of their sires into SId
                // The pedigree is already sorted so the oldest animal IDs appears first
                for (size_t i = 1; i <= n; i++)
                {
                    SId[i] = i;
                }

                k = 1;
                i = 1;
                for (; i <= n;)
                { // do from the oldest sire to the youngest sire

                    int t = Ped[SId[i]][0];
                    if (!t)
                    {
                        F[SId[i++]] = 0.0; // sire is not known
                    }
                    else
                    {

                        S = Ped[SId[i]][0];
                        rS = Link[S];
                        MIP = MaxIdP[rS];
                        x[rS] = 1.0;

                        for (; k <= S; ++k)
                        { // compute within family segregation variances
                            if (Link[k])
                                B[Link[k]] = 0.5 - 0.25 * (F[Ped[k][0]] + F[Ped[k][1]]);
                        }

                        for (j = rS; j; --j)
                        { // trace back the reduced pedigree
                            if (x[j])
                            { // consider only ancestors of the sire

                                if (rPed[j][0])
                                    x[rPed[j][0]] += x[j] * 0.5;
                                if (rPed[j][1])
                                    x[rPed[j][1]] += x[j] * 0.5;
                                x[j] *= B[j];
                            }
                        }

                        for (j = 1; j <= MIP; ++j)
                        { // trace forth the reduced pedigree
                            x[j] += (x[rPed[j][0]] + x[rPed[j][1]]) * 0.5;
                        }

                        for (; i <= n; ++i)
                        { // obtain F for progeny of the current sire
                            if ((int)S != Ped[SId[i]][0])
                                break;
                            else
                                F[SId[i]] = x[Link[Ped[SId[i]][1]]] * 0.5;
                        }

                        for (j = 1; j <= MIP; ++j)
                        {
                            x[j] = 0.0; // set to 0 for next evaluation of sire
                        }
                    }
                }

                for (size_t i = 1; i <= n; i++)
                {

                    auto s = Ped[i][0];
                    auto d = Ped[i][1];

                    inbrF[i] = F[i]; // fill global vector of inbr. coeff.

                    if (s && d)
                    {
                        double g = 1.0 / (0.5 - 0.25 * (F[s] + F[d]));
                        dinv.push_back(g);
                    }
                    else if (!s && !d)
                    {
                        double g = 1.0;
                        dinv.push_back(g);
                    }
                    else
                    {
                        double g;
                        if (s)
                            g = 1.0 / (0.75 - 0.25 * (F[s]));
                        else
                            g = 1.0 / (0.75 - 0.25 * (F[d]));
                        dinv.push_back(g);
                    }
                }
            }
            else
            {
                for (size_t i = 0; i <= ped.size(); i++)
                    inbrF.push_back(0.0);

                std::int64_t s, d;
                for (auto const &elem : ped)
                {
                    s = elem.second.val_1;
                    d = elem.second.val_2;
                    if (s && d)
                        dinv.push_back(2.0);
                    else if (!s && !d)
                        dinv.push_back(1.0);
                    else
                        dinv.push_back(4.0 / 3.0);
                }
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_dinv(std::map<PedPair, PedPair> &, std::vector<double> &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_dinv(std::map<PedPair, PedPair> &, std::vector<double> &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_dinv(std::map<PedPair, PedPair> &, std::vector<double> &, bool)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    void Amat::get_ainv(std::map<PedPair, PedPair> &ped, std::map<PedPair, double> &ai, bool inbreed)
    {
        // Calculates A(-1) matrix (as map representation)
        try
        {
            PedPair akey;
            std::vector<double> di; // store D(-1)

            get_dinv(ped, di, inbreed);

            std::int64_t s, d, id;
            double dinv;

            for (auto const &elem : ped)
            {
                dinv = di.front();
                s = elem.second.val_1;
                d = elem.second.val_2;
                id = elem.first.val_2;

                if (s && d)
                {
                    akey.val_1 = id;
                    akey.val_2 = id;
                    ai[akey] = ai[akey] + dinv;

                    akey.val_1 = id;
                    akey.val_2 = s;
                    if (s > id)
                    {
                        akey.val_1 = s;
                        akey.val_2 = id;
                    }
                    ai[akey] = ai[akey] - dinv / 2.0;

                    akey.val_1 = id;
                    akey.val_2 = d;
                    if (d > id)
                    {
                        akey.val_1 = d;
                        akey.val_2 = id;
                    }
                    ai[akey] = ai[akey] - dinv / 2.0;

                    akey.val_1 = s;
                    akey.val_2 = s;
                    ai[akey] = ai[akey] + dinv / 4.0;

                    akey.val_1 = d;
                    akey.val_2 = d;
                    ai[akey] = ai[akey] + dinv / 4.0;

                    if (s >= d)
                    {
                        akey.val_1 = s;
                        akey.val_2 = d;
                    }
                    else
                    {
                        akey.val_1 = d;
                        akey.val_2 = s;
                    }
                    ai[akey] = ai[akey] + dinv / 4.0;
                }
                else if (!s && !d)
                {
                    akey.val_1 = id;
                    akey.val_2 = id;
                    ai[akey] = ai[akey] + dinv;
                }
                else
                {
                    akey.val_1 = id;
                    akey.val_2 = id;
                    ai[akey] = ai[akey] + dinv;

                    akey.val_1 = id;
                    if (s)
                    {
                        akey.val_2 = s;
                        if (s > id)
                        {
                            akey.val_1 = s;
                            akey.val_2 = id;
                        }
                    }
                    else
                    {
                        akey.val_2 = d;
                        if (d > id)
                        {
                            akey.val_1 = d;
                            akey.val_2 = id;
                        }
                    }
                    ai[akey] = ai[akey] - dinv / 2.0;

                    if (s)
                    {
                        akey.val_1 = s;
                        akey.val_2 = s;
                    }
                    else
                    {
                        akey.val_1 = d;
                        akey.val_2 = d;
                    }
                    ai[akey] = ai[akey] + dinv / 4.0;
                }
                di.erase(di.begin());
            }
            di.clear();
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_ainv(std::map<PedPair, PedPair> &, std::map<PedPair, double> &, bool)" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_ainv(std::map<PedPair, PedPair> &, std::map<PedPair, double> &, bool)" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_ainv(std::map<PedPair, PedPair> &, std::map<PedPair, double> &, bool)" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<std::int64_t> Amat::get_genotyped_ids()
    {
        try
        {
            return genotypedID;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_genotyped_ids()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_genotyped_ids()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_genotyped_ids()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<std::int64_t> Amat::get_core_ids()
    {
        try
        {
            return coreID;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_core_ids()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_core_ids()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_core_ids()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<std::int64_t> Amat::get_pedigree_ids()
    {
        try
        {
            return pedID;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_pedigree_ids()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_pedigree_ids()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_pedigree_ids()" << '\n';
            throw;
        }
    }

    //===============================================================================================================

    std::vector<double> Amat::get_inbreeding()
    {
        try
        {
            return inbrF;
        }
        catch (const std::exception &e)
        {
            std::cerr << "Exception in Amat::get_inbreeding()" << '\n';
            std::cerr << e.what() << '\n';
            throw;
        }
        catch (const std::string &e)
        {
            std::cerr << "Exception in Amat::get_inbreeding()" << '\n';
            std::cerr << "Reason: " << e << '\n';
            throw;
        }
        catch (...)
        {
            std::cerr << "Exception in Amat::get_inbreeding()" << '\n';
            throw;
        }
    }

    //===============================================================================================================
    
        void Amat::clear()
        {
            try
            {
                pedigree_from_file.clear();
                birth_id_map.clear();
                pedigree.clear();
                r_pedigree.clear();
            }
            catch (const std::exception &e)
            {
                std::cerr << "Exception in Amat::clear()" << '\n';
                std::cerr << e.what() << '\n';
                throw;
            }
            catch (const std::string &e)
            {
                std::cerr << "Exception in Amat::clear()" << '\n';
                std::cerr << "Reason: " << e << '\n';
                throw;
            }
            catch (...)
            {
                std::cerr << "Exception in Amat::clear()" << '\n';
                throw;
            }
        }
    
    //===============================================================================================================

} // end of namespace evoped