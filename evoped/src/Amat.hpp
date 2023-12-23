#ifndef Amat_hpp__
#define Amat_hpp__

#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

namespace evoped
{
    class Amat
    {
    public:

        Amat();
        Amat(const std::string &ped_file);
        Amat(const std::string &ped_file, const std::string &g_file);

        void get_ainv(); // format of the output ainv ???
        void get_ainv(const std::string &ped_file); // format of the output ainv ???
        void get_ainv(const std::string &ped_file, const std::string &g_file); // format of the output ainv ???

        std::vector<std::int64_t> get_genotyped_ids();
        std::vector<std::int64_t> get_core_ids();
        std::vector<std::int64_t> get_pedigree_ids();
        std::vector<double> get_inbreeding();

        void clear();

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
        
        //int get_A22(std::map<PedPair, PedPair> &ped, float *A22);        
        //int getA22vector(std::vector<float> &w, std::vector<std::int64_t> &v, std::vector<std::vector<std::int64_t>> &Ped);        
        //std::map<PedPair, float> a; // A matrix container

    private:

        std::string pedigree_file;
        std::string genotyped_file;

        // ------- temporal pedigree storages -------------
        std::map<PedPair, PedPair> pedigree_from_file;        /* original pedigree from red file. will be deleted after use */
        std::map<std::int64_t, std::int64_t> birth_id_map; /* the map which holds <animal_id, birth_day>, used when check correctness of pedigree */

        std::map<PedPair, PedPair> pedigree;   /* traced back original pedigree */
        std::map<PedPair, PedPair> r_pedigree; /* reduced pedigree container, sorted by day order during initialization */
        // ------------------------------------------------
        
        // ------ A(-1) matrix container ------------------
        std::map<PedPair, double> ainv;
        std::map<PedPair, double> r_ainv;
        // ------------------------------------------------

        std::vector<std::int64_t> genotypedID; /* container for genotyped IDs (from 'typed' file) */
        std::vector<std::int64_t> coreID;      /* container for core IDs (from 'typed' file) */
        std::vector<std::int64_t> pedID;       /* container for the list of pedigree IDs */
        std::vector<double> inbrF;             /* inbreeding coeffecients */

        void fread_pedigree(const std::string &ped_file);
        void fread_genotyped_id(const std::string &g_file);
        void trace_pedigree();
        
        void get_ainv(std::map<PedPair, PedPair> &ped,
                      std::map<PedPair, double> &ai,
                      bool inbreed);

        void trace_pedigree(std::map<PedPair, PedPair> &in_ped,
                            std::map<PedPair, PedPair> &out_ped,
                            std::vector<std::int64_t> &traced_id);

        void get_dinv(std::map<PedPair, PedPair> &ped,
                      std::vector<double> &dinv,
                      bool inbreed);

        bool is_invect(std::vector<std::int64_t> &where, std::int64_t what); // was find_invect
        std::int64_t pos_inped(std::map<std::int64_t, std::int64_t> &codemap, std::int64_t id);
        bool is_unique(std::vector<std::int64_t> &x);
    };

} // end of namespace evoped

#endif // Amat_hpp__