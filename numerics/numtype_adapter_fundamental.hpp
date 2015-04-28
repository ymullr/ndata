#ifndef NUMTYPE_ADAPTERS_FUNDAMENTAL_HPP_H0TGUOAB
#define NUMTYPE_ADAPTERS_FUNDAMENTAL_HPP_H0TGUOAB

#include <type_traits>

//TODO document

namespace ndata {

template<typename T, typename Enable = void>
struct numtype_adapter {

    static_assert(std::is_arithmetic<T>::value == false, "");

};

template<typename T>
struct numtype_adapter<
        T,
        typename std::enable_if<std::is_arithmetic<T>::value>::type
        >
{

    static constexpr T ZERO = T(0);

};

}

#endif /* end of include guard: NUMTYPE_ADAPTERS_FUNDAMENTAL_HPP_H0TGUOAB */ 
