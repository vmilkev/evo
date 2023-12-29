#ifndef Gmat_hpp__
#define Gmat_hpp__

#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#include "cs_matrix.hpp"

#define workload 100000

namespace evoped
{
    class Gmat
    {
    public:
        Gmat();
        ~Gmat();

        void read_matrix(const std::string &gmat_file);
        void make_matrix(const std::string &fname);
        void invert_matrix();
        void invert_matrix(std::vector<std::int64_t>& core_id); // sparse inverse (APY)
        void bin_write();
        void bin_read();
        void get_matrix(evolm::matrix<double>& arr);
        void clear();

    private:
        evolm::matrix<double> G;          // G or inv. of G matrix container
        //evolm::matrix<double> invG;       // inversed G matrix container
        
        evolm::matrix<double> Z;          // Z (snp) matrix, temporal
        std::vector<std::int64_t> gmatID; // container for the list of G matrix IDs, initiated while reading pre-built G-matrix from file
        double freq;
        std::map<std::int64_t, std::string> snp_map; // initial markers data, temporal
        std::map<size_t, std::int64_t> anim_id_map;  // key: consecutive index, value: animal ID

        void find_RecodedIdMap(std::map<std::int64_t, std::int64_t> &id_map,
                               std::vector<std::int64_t> &whereIidList,
                               std::vector<std::int64_t> &whatIdList);      // not sure if this is needed
        void get_RecodedIdMap(std::map<std::int64_t, std::int64_t> &id_map,
                              std::vector<std::int64_t> &idVect);           // not sure if this is needed
        size_t find_invect(std::vector<std::int64_t> &where,
                           std::int64_t what);                              // was find_invect2
        bool is_unique(std::vector<std::int64_t> &x);
        template <typename T>
        void get_gvalues(std::vector<std::int64_t> &row,
                         std::vector<std::int64_t> &col,
                         std::vector<T> &val,
                         double diag_val);             // expected to be usefull in the Hmat class
        void get_gids(std::vector<std::int64_t> &ids); // expected to be usefull in the Hmat class
        void read_snp(const std::string &snp_file);
        void parse_string(std::string &snp_str, std::vector<int> &markers);
        void make_zmatrix();
        void make_matrix();
        template <typename T>
        void read_matrix(const std::string &gmat_file, std::vector<std::int64_t> &g_row, std::vector<std::int64_t> &g_col, std::vector<T> &g_val);
    };

} // end of namespace evoped

#endif // Gmat_hpp__