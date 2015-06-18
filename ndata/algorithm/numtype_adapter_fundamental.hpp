#ifndef NUMTYPE_ADAPTERS_FUNDAMENTAL_HPP_H0TGUOAB
#define NUMTYPE_ADAPTERS_FUNDAMENTAL_HPP_H0TGUOAB

#include <type_traits>

//TODO document

namespace ndata {
namespace helpers {

template<typename T, typename Enable = void>
struct numtype_adapter {

    static_assert(std::is_arithmetic<T>::value == false, "");
    static_assert(
            std::is_base_of<
                Eigen::EigenBase<T>,
                T
            >::value == false,
            "is base"
        );
    //static_assert(
    //        false,
    //        "Numerical type unsupported. Did you include a relevant numtype adapter definition?"
    //        );
};

template<typename T>
struct numtype_adapter<
        T,
        typename std::enable_if<std::is_arithmetic<T>::value>::type
        >
{

    static constexpr T ZERO = T(0);

};

} //end namespace ndata
} //end namespace helpers

#endif /* end of include guard: NUMTYPE_ADAPTERS_FUNDAMENTAL_HPP_H0TGUOAB */
