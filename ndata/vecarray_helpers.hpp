#ifndef VECARRAY_HELPERS_HPP_POIQTBTA
#define VECARRAY_HELPERS_HPP_POIQTBTA

namespace ndata {
namespace helpers {

    /**
     * Compile time "function" with recursive templated struct
     */
    template <long ... vals>
    struct static_max_or_dynamic {
    };

    template <long val1, long val2, long ... vals>
    struct static_max_or_dynamic<val1, val2, vals...> {

        static constexpr long value =
                (val1==DYNAMICALLY_SIZED or val2==DYNAMICALLY_SIZED)?
                        DYNAMICALLY_SIZED:
                    (val1 >= val2)?
                        static_max_or_dynamic<val1, vals...>::value:
                            static_max_or_dynamic<val2, vals...>::value;
    };

    template <long val1>
    struct static_max_or_dynamic<val1> {
        static constexpr long value = val1;
    };


    /**
     * Compile time "function" with recursive templated struct
     */
    template <long ... vals>
    struct is_any_dynamically_sized {
    };

    template <long val1, long ... vals>
    struct is_any_dynamically_sized<val1, vals...> {
        static constexpr bool value = is_any_dynamically_sized<val1>::value or is_any_dynamically_sized<vals...>::value;
    };

    template <long val1>
    struct is_any_dynamically_sized<val1> {
        static constexpr bool value = val1 == DYNAMICALLY_SIZED;
    };


    //make vecarray like biggest with overloads to return a static if both operands are static,
    //or a dynamic one if any or both operands are dynamic
    template<
        long StatSize1,
        long StatSize2
        >
    auto
    make_vecarray_like_biggest(
                vecarray<
                    long,
                    std::enable_if<StatSize1!=DYNAMICALLY_SIZED, std::integral_constant<long, StatSize1>>::type::value
                    > v1,
                vecarray<
                    long,
                    std::enable_if<StatSize2!=DYNAMICALLY_SIZED, std::integral_constant<long, StatSize2>>::type::value
                    > v2
            )
    {
        return vecarray<long,
                (v1.STATIC_SIZE_OR_DYNAMIC >= v2.STATIC_SIZE_OR_DYNAMIC)?
                                v1.STATIC_SIZE_OR_DYNAMIC :
                                    v2.STATIC_SIZE_OR_DYNAMIC
                       >();
    };


    //make vecarray like biggest with overloads to return a static vecarray if both operands are static,
    //or a dynamic one if any or both operands are dynamic
    //specialization 1
    template<
        long StatSize1,
        long StatSize2
        >
    auto
    make_vecarray_like_biggest_dyn(
                vecarray<
                    long,
                    StatSize1//std::enable_if<StatSize1!=DYNAMICALLY_SIZED, std::integral_value<StatSize1>>::type::value
                    > v1,
                vecarray<
                    long,
                    StatSize2 //std::enable_if<StatSize2!=DYNAMICALLY_SIZED, std::integral_value<StatSize2>>::type::value
                    > v2
            )
    {
        return vecarray<long, DYNAMICALLY_SIZED> (std::max(v1.size(), v2.size()));
    };


    //specialization 2
    template<
        long StatSize1,
        long StatSize2
        >
    auto
    make_vecarray_like_biggest(
                vecarray<
                    long,
                    std::enable_if<StatSize1==DYNAMICALLY_SIZED, std::integral_constant<long, StatSize1>>::type::value
                    > v1,
                vecarray<
                    long,
                    StatSize2 //std::enable_if<StatSize2!=DYNAMICALLY_SIZED, std::integral_value<StatSize2>>::type::value
                    > v2
            )
    {
        return make_vecarray_like_biggest_dyn(v1, v2);
    };

    //specialization 3
    template<
        long StatSize1,
        long StatSize2
        >
    auto
    make_vecarray_like_biggest(
                vecarray<
                    long,
                    StatSize1
                    > v1,
                vecarray<
                    long,
                    std::enable_if<StatSize2==DYNAMICALLY_SIZED, std::integral_constant<long, StatSize2>>::type::value
                    > v2
            )
    {
        return make_vecarray_like_biggest_dyn(v1, v2);
    };


} //end namespace helpers
} //end namespace ndata

#endif /* end of include guard: VECARRAY_HELPERS_HPP_POIQTBTA */
