#ifndef effects_storage_hpp__
#define effects_storage_hpp__

#include "compact_storage.hpp"

namespace evolm
{
    class effects_storage
    {
    private:
        compact_storage<int> i_effect;    // type 1
        compact_storage<float> f_effect;  // type 2
        compact_storage<double> d_effect; // type 3

    public:
        
        bool is_diagonal;
        int type = 0;

        effects_storage();

        void set(compact_storage<int> &in_effect);
        void set(compact_storage<float> &in_effect);
        void set(compact_storage<double> &in_effect);

        void get(compact_storage<int> &out_effect);
        void get(compact_storage<float> &out_effect);
        void get(compact_storage<double> &out_effect);

        void fread();
        void fwrite();

        void clear();
        std::vector<size_t> shape();
        size_t size();

        void get_fcast(matrix<float> &out, size_t *irow, size_t *icol);
        void get_fcast(smatrix<float> &out, size_t *irow, size_t *icol);
        void get_fcast(float **out, size_t *irow, size_t *icol);
        void get_fcast(std::vector<float> &vals_out, std::vector<size_t> &keys_out, size_t *irow, size_t *icol);

        void row_dot_float_vect(float **in_vect, size_t irow, size_t icol, float &out_res);

        void extend_by(effects_storage &rhs_matr);
        void element_wise_dot(effects_storage &rhs_matr);
        void optimize();

        void print(std::string &message);
    };
}

#endif // effects_hpp__
