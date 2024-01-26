#ifndef Amat_hpp__
#define Amat_hpp__

#include <fstream>
#include <algorithm>

#include "cs_matrix.hpp"
#include "Utilities2.hpp"

namespace evoped
{
    class Amat
    {
    public:

        Amat();
        ~Amat();

        void make_matrix(const std::string &ped_file,
                         bool use_ainv);
        void make_matrix(const std::string &ped_file,
                         const std::string &g_file,
                         bool use_ainv);
        void make_all(const std::string &ped_file,
                      const std::string &g_file);
        void get_inbreeding(std::vector<double> &out);        
        void get_matrix(const std::string &name,
                        evolm::matrix<double>& arr,
                        std::vector<std::int64_t> &out,
                        bool keep_ondisk);
        void get_matrix(const std::string &name,
                        std::vector<double>& arr,
                        std::vector<std::int64_t> &out,
                        bool keep_ondisk);
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

        evolm::matrix<double> iA;
        evolm::matrix<double> irA;
        evolm::matrix<double> iA22;
        evolm::matrix<double> A22;
        std::vector<std::int64_t> id_iA;
        std::vector<std::int64_t> id_irA;
        std::vector<std::int64_t> id_A22;

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
        std::int64_t pos_inped(std::map<std::int64_t,std::int64_t> &codemap,
                               std::int64_t id);
        void map_to_matr(std::map<PedPair, double> &amap,
                         std::vector<std::int64_t> &ids,
                         bool use_ainv);
        void get_A22(std::map <PedPair, PedPair> &ped,
                           std::vector<std::int64_t> &genotypedID);
        void getA22vector(std::vector <double> &w,
                          std::vector <std::int64_t> &v,
                          std::vector<std::vector<std::int64_t> > &Ped);
        void get_iA22(evolm::matrix<double>& full_matr,
                      std::vector<std::int64_t>& matr_ids,
                      std::vector<std::int64_t>& selected_ids);
    };

} // end of namespace evoped

#endif // Amat_hpp__