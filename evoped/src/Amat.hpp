#ifndef Amat_hpp__
#define Amat_hpp__

#include "Utilities2.hpp"

namespace evoped
{
    template <typename T>
    class Amat
    {
    public:
        Amat();
        Amat( double threshold );
        ~Amat();

        void make_matrix(const std::string &ped_file,
                         bool use_ainv);
        void make_matrix(const std::string &ped_file,
                         const std::string &g_file,
                         bool use_ainv);
        void make_matrix_forgenotyped(
                         const std::string &ped_file,
                         const std::string &g_file,
                         bool use_ainv);        
        void make_matrix_forgenotyped(
                         const std::string &ped_file,
                         std::vector<std::int64_t> &genotyped_ids,
                         bool use_ainv);        
        void make_all(const std::string &ped_file,
                      const std::string &g_file);
        void make_all(const std::string &ped_file,
                      std::vector<std::int64_t> &g_ids);
        void get_inbreeding(std::vector<T> &out);
        
        void get_matrix(const std::string &name,
                        evolm::matrix<T> &arr,
                        std::vector<std::int64_t> &out);

        void get_matrix(const std::string &name,
                        evolm::smatrix<T> &arr,
                        std::vector<std::int64_t> &out);
        void get_matrix(const std::string &name,
                        std::vector<T> &arr,
                        std::vector<std::int64_t> &out);
        
        void clear();
        void set_sparsiity_threshold( double threshold );

    private:

        double data_sparsity;
        bool use_sparse; // the default is TRUE
        double sparsity_threshold; // the default is 90.0

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

        struct EmptyDataStatus
        {
            bool A;
            bool iA;
            bool irA;
            bool iA22;
            bool A22;

            bool A_s;
            bool iA_s;
            bool irA_s;
            bool iA22_s;
            bool A22_s;
        };
        EmptyDataStatus IsEmpty;

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
        evolm::smatrix<T> A22_s;

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
        void map_to_matr(std::map<PedPair, T> &in_amap,
                         std::vector<std::int64_t> &in_ids,
                         evolm::smatrix<T> &out_matr);
        void map_to_matr(std::map<PedPair, T> &in_amap,
                         std::vector<std::int64_t> &in_ids,
                         evolm::matrix<T> &out_matr);
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