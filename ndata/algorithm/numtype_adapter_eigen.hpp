#ifndef NUMTYPE_ADAPTERS_EIGEN_HPP_G9DJMXVS
#define NUMTYPE_ADAPTERS_EIGEN_HPP_G9DJMXVS

#include <algorithm/numtype_adapters_fundamental.hpp>

namespace ndata {
namespace helpers {

template <typename EigenT>
using DerivedOf = std::result_of<EigenT.derived>::type
template <typename EigenT>
using DenseBaseParent = Eigen::DenseBase<DerivedOf<EigenT>>

template <typename T>
struct numtype_adapter<
    T,
    typename std::enable_if<
        std::is_base_of<
            DenseBaseParent<T>,
            T
        >::value
    >::type
>
{

    static constexpr T ZERO = T::ZERO;

};

} //end namespace ndata
} //end namespace helpers




#endif /* end of include guard: NUMTYPE_ADAPTERS_EIGEN_HPP_G9DJMXVS */ 
