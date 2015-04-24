#ifndef HELPERS_HPP_YIOQCJSG
#define HELPERS_HPP_YIOQCJSG

#include <cstdlib>
#include <utility>
#include "tuple_utility.hpp"
#include <exception>
#include <type_traits>

#include "vecarray.hpp"

namespace ndata {

namespace helpers {


    using ShapeStridePair = std::pair<size_t, long>;

    // ( shape[i], strides[i] )
    template<size_t N> using SliceAcc = vecarray<ShapeStridePair, N>;


    template <long ... vals>
    struct static_max_or_dynamic {
    };

    template <long val1, long val2, long ... vals>
    struct static_max_or_dynamic<val1, val2, vals...> {

        static constexpr long value =
                (val1==DYNAMIC_SIZE or val2==DYNAMIC_SIZE)?
                        DYNAMIC_SIZE:
                    (val1 >= val2)?
                        static_max_or_dynamic<val1, vals...>::value:
                            static_max_or_dynamic<val2, vals...>::value;
    };

    template <long val1>
    struct static_max_or_dynamic<val1> {
        static constexpr long value = val1;
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
                    size_t,
                    std::enable_if<StatSize1!=DYNAMIC_SIZE, std::integral_constant<long, StatSize1>>::type::value
                    > v1,
                vecarray<
                    size_t,
                    std::enable_if<StatSize2!=DYNAMIC_SIZE, std::integral_constant<long, StatSize2>>::type::value
                    > v2
            )
    {
        return vecarray<size_t,
                (v1.STATIC_SIZE_OR_DYNAMIC >= v2.STATIC_SIZE_OR_DYNAMIC)?
                                v1.STATIC_SIZE_OR_DYNAMIC :
                                    v2.STATIC_SIZE_OR_DYNAMIC
                       >();
    };


    //make vecarray like biggest with overloads to return a static if both operands are static,
    //or a dynamic one if any or both operands are dynamic
    template<
        long StatSize1,
        long StatSize2
        >
    auto
    make_vecarray_like_biggest_dyn(
                vecarray<
                    size_t,
                    StatSize1//std::enable_if<StatSize1!=DYNAMIC_SIZE, std::integral_value<StatSize1>>::type::value
                    > v1,
                vecarray<
                    size_t,
                    StatSize2 //std::enable_if<StatSize2!=DYNAMIC_SIZE, std::integral_value<StatSize2>>::type::value
                    > v2
            )
    {
        return vecarray<size_t, ndata::DYNAMIC_SIZE> (std::max(v1.size(), v2.size()));
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
                    size_t,
                    std::enable_if<StatSize1==DYNAMIC_SIZE, std::integral_constant<long, StatSize1>>::type::value
                    > v1,
                vecarray<
                    size_t,
                    StatSize2 //std::enable_if<StatSize2!=DYNAMIC_SIZE, std::integral_value<StatSize2>>::type::value
                    > v2
            )
    {
        return make_vecarray_like_biggest_dyn(v1, v2);
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
                    size_t,
                    StatSize1
                    > v1,
                vecarray<
                    size_t,
                    std::enable_if<StatSize2==DYNAMIC_SIZE, std::integral_constant<long, StatSize2>>::type::value
                    > v2
            )
    {
        return make_vecarray_like_biggest_dyn(v1, v2);
    };

    template<typename ContentT, size_t idim>
    vecarray<ContentT, idim>
    array_from_argpack(vecarray<ContentT, idim> acc) {
        return acc;
    }


    /**
     * Accumulates the ndimensional indices passed as argument in an array
     * */
    template<typename ContentT, size_t idim, typename... ContentTPackT>
    auto
    array_from_argpack(
            vecarray<ContentT, idim> acc,
            ContentT i,
            ContentTPackT... rest
            ) {
        vecarray<ContentT, idim+1> new_acc = acc.append(i);
        return array_from_argpack<ContentT, idim+1>(new_acc, rest...);
    }

    template <typename ... Ts>
    struct static_check_valid_indice_types {
    };

    template <typename T, typename ... Ts>
    struct static_check_valid_indice_types<T, Ts...>{
        static constexpr bool value =
                static_check_valid_indice_types<T>::value
                and
                static_check_valid_indice_types<Ts...>::value ;
    };

    //termination
    template <typename T>
    struct static_check_valid_indice_types<T> {
        static constexpr bool value = std::is_integral<T>::value;
    };

    //template <typename Indexer, typename... Indexers>
    //auto // vecarray<Indexer, N>
    //broadcast(Indexers... others) {
    //    constexpr size_t arrsize = sizeof...(others);

    //    vecarray<Indexer, arrsize> arr_args = array_from_argpack(vecarray<Indexer, 0>(), others...);

    //    vecarray<Indexer, arrsize> ret;

    //    for (size_t i = 0; i < ret.size(); ++i) {
    //        vecarray<Indexer, arrsize> reord;
    //        reord[0] = arr_args[i];
    //        size_t io = 1;
    //        for (size_t i2 = 0; i2 < arrsize; ++i2) {
    //            if (io != i) {
    //                reord[io] = arr_args[i2];
    //                io++;
    //            }
    //        }
    //        ret[i] = broadcast_left_list(reord);
    //    }

    //    return ret;
    //}

    //termination
    template <typename Indexer>
    auto
    broadcast_left(Indexer ind, std::tuple<> t) {
        t=t;//silence warning
        return ind;
    }

    //broadcasts first element of the tuple on the other elements in the tuple
    template <typename Indexer, typename ... Indexers>
    auto
    broadcast_left(Indexer v1, std::tuple<Indexers...> t) {
        auto h_t = tuple_utility::split_ht(std::move(t));
        auto v2 = h_t.first;
        auto tuple_rest = h_t.second;

        auto shape_v1 = v1.get_shape();
        auto shape_v2 = v2.get_shape();

        size_t new_size = max(
                    shape_v1.size(),
                    shape_v2.size()
                    );

        auto strides_v1 = v1.get_strides();
        //auto strides_v2 = v2.get_strides();

        //new shape can be a static vecaray or a dynamic vecarray if any
        //of v1 or v2 is dynamic
        auto new_shape = make_vecarray_like_biggest<
                decltype(shape_v1)::STATIC_SIZE_OR_DYNAMIC,
                decltype(shape_v2)::STATIC_SIZE_OR_DYNAMIC
                >(v1.get_shape(), v2.get_shape());

        static_assert(
                    new_shape.STATIC_SIZE_OR_DYNAMIC ==  decltype(shape_v1)::STATIC_SIZE_OR_DYNAMIC
                    or
                    new_shape.STATIC_SIZE_OR_DYNAMIC ==  decltype(shape_v2)::STATIC_SIZE_OR_DYNAMIC,
                    ""
                    );

        //static_assert(
        //            new_shape.STATIC_SIZE_OR_DYNAMIC == std::max(
        //                decltype(shape_v1)::STATIC_SIZE_OR_DYNAMIC,
        //                decltype(shape_v2)::STATIC_SIZE_OR_DYNAMIC
        //                )
        //            or
        //            (
        //                (
        //                    decltype(shape_v1)::STATIC_SIZE_OR_DYNAMIC == DYNAMIC_SIZE
        //                    or
        //                    decltype(shape_v2)::STATIC_SIZE_OR_DYNAMIC == DYNAMIC_SIZE
        //                )
        //                and
        //                new_shape.STATIC_SIZE_OR_DYNAMIC == DYNAMIC_SIZE
        //            ),
        //            ""
        //            );

        vecarray<long, new_shape.STATIC_SIZE_OR_DYNAMIC>
            new_strides (new_shape.dynsize());
            //new_strides_v2 = shape;

        for (size_t is = 0; is < new_size; ++is) {
            size_t i_v2 = v2.get_shape().size()-is-1,
                   i_v1 = v1.get_shape().size()-is-1,
                   i_new = new_size - is -1;

            if (i_v1 < shape_v1.size() && i_v2 < shape_v2.size()) {

                if (shape_v1[i_v1] == shape_v2[i_v2]) {
                    new_shape[i_new] = shape_v1[i_v1];
                    new_strides[i_new] = strides_v1[i_v1];
                    //new_strides_v2[i_new] = strides_v2[i_v2];
                } else if (shape_v1[i_v1] == 1) {
                    new_shape[i_new] = shape_v2[i_v2];
                    new_strides[i_new] = 0;
                    //new_strides_v2[i_new] = strides_v2[i_v2];
                } else if (shape_v2[i_v2] == 1) {
                    new_shape[i_new] = shape_v1[i_v1];
                    new_strides[i_new] = strides_v1[i_v1];
                    //new_strides_v2[i_new] = 0;
                } else {
                    throw(std::domain_error("Dimension size must match or be equal to 1"));
                }

            } else if (i_v2 >= v2.get_shape().size()) {
                new_shape[i_new] = shape_v1[i_v1];
                new_strides[i_new] = strides_v1[i_v1];
                //new_strides_v2[i_new] = 0;
            } else if (i_v1 >= v1.get_shape().size()) {
                new_shape[i_new] = shape_v2[i_v2];
                new_strides[i_new] = 0;
                //new_strides_v2[i_new] = strides_v2[i_v2];
            } else {
                assert("false");
            }
        }
        auto broadcasted_tuple = v1.template reshape<new_shape.STATIC_SIZE_OR_DYNAMIC>(new_shape, new_strides);
        return broadcast_left(
                    broadcasted_tuple,
                    tuple_rest
                    );
    }


    ////recursion termination, when toproc is empty
    //template <
    //          typename... IndexersAcc,
    //          typename... IndexersFull
    //        >
    //auto //std::tuple<IndexersFull...>
    //broadcast_rec(
    //        //the accumulated result
    //        std::tuple<IndexersAcc...> acc,
    //        //this is the pack from which elements to be processed are taken, function returns when empty
    //        std::tuple<> toproc
    //    )
    //{
    //    toproc = toproc; //silence warnings
    //    return acc;
    //}

    //template <
    //          typename... IndexersAcc,
    //          typename... IndexersToProcess
    //        >
    //auto //std::tuple<IndexersFull...>
    //broadcast_rec(
    //        //the accumulated result
    //        std::tuple<IndexersAcc...> acc,
    //        //this is the pack from which elements to be processed are taken, function returns when empty
    //        std::tuple<IndexersToProcess...> toproc
    //        ) {
    //    auto proc_head = tuple_utility::head(toproc);
    //    auto proc_tail = tuple_utility::tail(toproc);

    //    auto processed_indexer =
    //            broadcast_left(
    //                proc_head,
    //                full_head,
    //                full_tail
    //                );
    //    return broadcast_rec(
    //        //new acc
    //        std::tuple_cat(
    //            acc,
    //            make_tuple(processed_indexer)
    //        ),
    //        proc_tail
    //    );
    //}

    template <
            long ndims,
            typename ... Indexers
            >
    auto
    broadcast_on_target(
            indexer<ndims> target,
            tuple<Indexers...> all
            )
    {
        return tuple_utility::tuple_transform(
                [target] (auto&& tup) {
                    return broadcast_left(
                                tup,
                                std::make_tuple(target)
                                );
                },
                all
        );
    }


    template <typename ... Indexers>
    auto //std::tuple<Indexers...> tuple of indexers
    broadcast(std::tuple<Indexers...> tup_ind) {
        auto h_t = tuple_utility::split_ht(std::move(tup_ind));

        //indexer<helpers::static_max_or_dynamic<
        //        decltype(Indexers::shape_)::STATIC_SIZE_OR_DYNAMIC...>::value>
        auto
            indexer_target = broadcast_left(
                //tuple_utility::head(tup_ind),
                //tuple_utility::tail(tup_ind)
                h_t.first,
                h_t.second
                );

        //std::tuple<
        //        ndatacontainer<
        //            typename Indexers::type_ContainerT,
        //            typename Indexers::type_T,
        //            helpers::static_max_or_dynamic<Indexers::STATIC_SIZE_OR_DYNAMIC...>::value
        //        >...
        //    >
        auto
        ret = broadcast_on_target(
            make_indexer(indexer_target),
            //toproc
            std::tuple_cat(
                std::make_tuple(h_t.first),
                h_t.second
                )
            );
        return ret;
    }

    template <typename TupIndexers>
    auto
    broadcast_views(
            TupIndexers & tup_ind
            )
    {
        auto tupind_views = tuple_utility::tuple_transform_ptr(
                    [] (auto& tupelt) {
                        return tupelt.to_view();
                    },
                    tup_ind
                );
        return broadcast(tupind_views);
    }

} //end .namespace helpers

} //end namespace ndata

#endif /* end of include guard: HELPERS_HPP_YIOQCJSG */
