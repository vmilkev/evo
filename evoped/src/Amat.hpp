#ifndef Amat_hpp__
#define Amat_hpp__

#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "cs_matrix.hpp"

namespace evoped
{
    class Amat
    {
    public:

        Amat();
        ~Amat();

        void make_matrix(const std::string &ped_file, bool use_ainv);
        void make_matrix(const std::string &ped_file, const std::string &g_file, bool use_ainv);
        void get_ids(std::vector<std::int64_t> &out);
        void get_inbreeding(std::vector<double> &out);
        void get_matrix(evolm::matrix<double>& arr);
        void get_matrix(std::vector<double>& arr);
        void bin_write();
        void bin_read();
        void clear();
        
    private:

        struct PedPair
        {
            std::int64_t val_1; // KEY: day, PAIR: sire
            std::int64_t val_2; // KEY: id, PAIR: dame

            bool operator<(const PedPair &rhs) const
            {
                if (this->val_1 != rhs.val_1)
                    return this->val_1 < rhs.val_1;
                else
                    return this->val_2 < rhs.val_2;
            }
        };

        std::map<std::int64_t, std::int64_t> birth_id_map; // the map which holds <animal_id, birth_day>, used when check correctness of pedigree
        evolm::matrix<double> A;                           // A or A(-1) matrix container
        std::vector<std::int64_t> traced_pedID;            // traced pedigree IDs
        std::vector<double> inbrF;                         // inbreeding coeffecients

        void fread_pedigree(const std::string &ped_file,
                            std::map<PedPair, PedPair> &out_ped,
                            std::vector<std::int64_t> &out_ids);
        void fread_genotyped_id(const std::string &g_file,
                                std::vector<std::int64_t> &out_ids);        
        void get_ainv(std::map<PedPair, PedPair> &ped,
                      std::map<PedPair, double> &ai,
                      bool inbreed);

        void get_a(std::map<PedPair, PedPair> &in_ped,
                   std::map<PedPair, double> &out_a);

        void trace_pedigree(std::map<PedPair, PedPair> &in_ped,
                            std::map<PedPair, PedPair> &out_ped,
                            std::vector<std::int64_t> &traced_id);
        void get_dinv(std::map<PedPair, PedPair> &ped,
                      std::vector<double> &dinv,
                      bool inbreed);
        bool is_invect(std::vector<std::int64_t> &where,
                       std::int64_t what);               // was find_invect
        std::int64_t pos_inped(std::map<std::int64_t,std::int64_t> &codemap,
                               std::int64_t id);
        bool is_unique(std::vector<std::int64_t> &x);
        void get_RecodedIdMap(std::map<std::int64_t,std::int64_t> &id_map,
                              std::vector<std::int64_t> &idVect);
        void map_to_matr(std::map<PedPair, double> &amap,
                         std::vector<std::int64_t> &ids,
                         bool use_ainv);
    };

} // end of namespace evoped

#endif // Amat_hpp__