/*
    dense_matrix.cpp

    Explicit instantiation declaration of matrix class.

*/

#include "dense_matrix.hpp"

namespace evolm
{
    extern template class matrix<float>;
    extern template class matrix<double>;
}
