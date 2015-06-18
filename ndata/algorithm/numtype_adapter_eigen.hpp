#ifndef NUMTYPE_ADAPTERS_EIGEN_HPP_G9DJMXVS
#define NUMTYPE_ADAPTERS_EIGEN_HPP_G9DJMXVS

#include <ndata/algorithm/numtype_adapter_fundamental.hpp>

//TODO test this with interp

namespace ndata {
namespace helpers {

template <typename EigenT>
using DerivedOf = typename std::result_of<typename EigenT::derived>::type;

template <typename EigenT>
using DenseBaseParent = Eigen::DenseBase<DerivedOf<EigenT>>;

template <typename T>
struct numtype_adapter<
    T,
    typename std::enable_if<
        std::is_base_of<
            Eigen::EigenBase<T>,
            T
        >::value
    >::type
>
{

    static const
    T ZERO;

};

template <typename T>
const T
numtype_adapter<
    T,
    typename std::enable_if<
        std::is_base_of<
            Eigen::EigenBase<T>,
            T
        >::value
    >::type
>::ZERO = T::Zero();

} //end namespace ndata
} //end namespace helpers




#endif /* end of include guard: NUMTYPE_ADAPTERS_EIGEN_HPP_G9DJMXVS */ 
