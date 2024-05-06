#ifndef Amat_hpp__
#define Amat_hpp__

// #include <fstream>
// #include <algorithm>
#include "sparse_matrix.hpp"
#include "cs_matrix.hpp"
#include "Utilities2.hpp"

namespace evoped
{
    template <typename T>
    class Amat
    {
    public:
        Amat();
        ~Amat();

        void make_matrix(const std::string &ped_file,
                         bool use_ainv,
                         bool use_large);
        void make_matrix(const std::string &ped_file,
                         const std::string &g_file,
                         bool use_ainv,
                         bool use_large);
        void make_all(const std::string &ped_file,
                      const std::string &g_file,
                      bool use_large);
        void get_inbreeding(std::vector<T> &out);
        void get_matrix(const std::string &name,
                        evolm::matrix<T> &arr,
                        std::vector<std::int64_t> &out,
                        bool keep_ondisk);
        void get_matrix(const std::string &name,
                        evolm::smatrix<T> &arr,
                        std::vector<std::int64_t> &out,
                        bool keep_ondisk);
        void get_matrix(const std::string &name,
                        std::vector<T> &arr,
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
        std::vector<std::int64_t> traced_pedID;            // traced pedigree IDs
        std::vector<T> inbrF;                              // inbreeding coeffecients

        // dense matrices
        evolm::matrix<T> A;    // A or A(-1) matrix container (only for internal purposes)
        evolm::matrix<T> iA;   // full A(-1) container
        evolm::matrix<T> irA;  // reduced A(-1) container
        evolm::matrix<T> iA22; // reduced A(-1) container for genotypeed individuals (this is always dense ?)
        evolm::matrix<T> A22;  // reduced A container for genotypeed individuals (this is always dense ?)

        // sparse matrices
        evolm::smatrix<T> A_s;
        evolm::smatrix<T> iA_s;
        evolm::smatrix<T> irA_s;
        evolm::smatrix<T> iA22_s;
        // evolm::matrix<T> A22_s;

        std::vector<std::int64_t> id_iA;
        std::vector<std::int64_t> id_irA;
        std::vector<std::int64_t> id_A22;

        void fread_pedigree(const std::string &ped_file,
                            std::map<PedPair, PedPair> &out_ped,
                            std::vector<std::int64_t> &out_ids);
        void fread_genotyped_id(const std::string &g_file,
                                std::vector<std::int64_t> &out_ids);
        void get_ainv(std::map<PedPair, PedPair> &ped,
                      std::map<PedPair, T> &ai,
                      bool inbreed);

        void get_a(std::map<PedPair, PedPair> &in_ped,
                   std::map<PedPair, T> &out_a);

        void trace_pedigree(std::map<PedPair, PedPair> &in_ped,
                            std::map<PedPair, PedPair> &out_ped,
                            std::vector<std::int64_t> &traced_id);
        void get_dinv(std::map<PedPair, PedPair> &ped,
                      std::vector<T> &dinv,
                      bool inbreed);
        std::int64_t pos_inped(std::map<std::int64_t, std::int64_t> &codemap,
                               std::int64_t id);
        void map_to_matr(std::map<PedPair, T> &amap,
                         std::vector<std::int64_t> &ids,
                         bool use_ainv,
                         bool use_large);
        void get_A22(std::map<PedPair, PedPair> &ped,
                     std::vector<std::int64_t> &genotypedID);
        void getA22vector(std::vector<T> &w,
                          std::vector<std::int64_t> &v,
                          std::vector<std::vector<std::int64_t>> &Ped);
        void get_iA22(evolm::matrix<T> &full_matr,
                      std::vector<std::int64_t> &matr_ids,
                      std::vector<std::int64_t> &selected_ids);
        void get_iA22(evolm::smatrix<T> &full_matr,
                      std::vector<std::int64_t> &matr_ids,
                      std::vector<std::int64_t> &selected_ids);

        // multithreading
        void thread_loads(std::vector<std::int64_t> &in,
                          std::vector<size_t> &out);
        void trace_operation(std::vector<std::int64_t> &days,
                             std::vector<std::int64_t> &ids,
                             std::vector<std::int64_t> &sire,
                             std::vector<std::int64_t> &dame,
                             std::vector<std::int64_t> &in_traced_id,
                             std::vector<std::int64_t> &out_traced_id,
                             std::map<PedPair, PedPair> &out_ped,
                             std::vector<size_t> &loads_vect,
                             size_t thr_id);
    };

} // end of namespace evoped

#endif // Amat_hpp__