/*
    Explicit instantiation declaration of the compact_storage class.
*/

#include "compact_storage.hpp"

namespace evolm
{
    extern template class compact_storage<float>;
    extern template class compact_storage<double>;
    extern template class compact_storage<int>;
}