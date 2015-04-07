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
                (v1.static_size_or_dynamic >= v2.static_size_or_dynamic)?
                                v1.static_size_or_dynamic :
                                    v2.static_size_or_dynamic
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

    //broadcasts v1 to match v2
    template <typename Indexer1, typename Indexer2>
    auto
    broadcast_left(Indexer1 v1 , Indexer2 v2, std::tuple<> empty_tup) {
        empty_tup = empty_tup; //silence a warning

        auto shape_v1 = v1.get_shape(),
             shape_v2 = v2.get_shape();

        size_t new_size = max(
                    shape_v1.size(),
                    shape_v2.size()
                    );

        auto strides_v1 = v1.get_strides();
        //auto strides_v2 = v2.get_strides();

        //new shape can be a static vecaray or a dynamic vecarray if any
        //of v1 or v2 is dynamic
        auto new_shape = make_vecarray_like_biggest<
                decltype(v1.get_shape())::static_size_or_dynamic,
                decltype(v2.get_shape())::static_size_or_dynamic
                >(v1.get_shape(), v2.get_shape());

        vecarray<long, new_shape.static_size_or_dynamic>
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
        return v1.template reshape<new_shape.static_size_or_dynamic>(new_shape, new_strides);
    }

    template <typename Indexer1, typename Indexer2, typename... Indexers>
    auto
    broadcast_left(Indexer1 arg0, Indexer2 arg1, std::tuple<Indexers...> argN) {
        auto ret = broadcast_left(
                    broadcast_left(arg0, arg1, std::tuple<>()),
                    tuple_utility::head(argN),
                    tuple_utility::tail(argN)
                    );
        return ret;
    }

    //recursion termination, when toproc is empty
    template <
              typename... IndexersAcc,
              typename... IndexersFull
            >
    auto //std::tuple<IndexersFull...>
    broadcast_rec(
            //the accumulated result
            std::tuple<IndexersAcc...> acc,
            //this is the pack from which elements to be processed are taken, function returns when empty
            std::tuple<> toproc,
            //full list of indexers that is kept around
            std::tuple<IndexersFull...> full
        )
    {
        toproc = toproc; //silence warnings
        full=full;
        return acc;
    }

    template <
              typename... IndexersAcc,
              typename... IndexersToProcess,
              typename... IndexersFull
            >
    auto //std::tuple<IndexersFull...>
    broadcast_rec(
            //the accumulated result
            std::tuple<IndexersAcc...> acc,
            //this is the pack from which elements to be processed are taken, function returns when empty
            std::tuple<IndexersToProcess...> toproc,
            //full list of indexers that is kept around
            std::tuple<IndexersFull...> full
            ) {
        auto proc_head = tuple_utility::head(toproc);
        auto proc_tail = tuple_utility::tail(toproc);
        auto full_head = tuple_utility::head(full);
        auto full_tail = tuple_utility::tail(full);

        auto processed_indexer =
                broadcast_left(
                    proc_head,
                    full_head,
                    full_tail
                    );
        return broadcast_rec(
            //new acc
            std::tuple_cat(
                acc,
                make_tuple(processed_indexer)
            ),
            proc_tail,
            full
        );
    }


    template <typename TupIndexers>
    auto //std::tuple<Indexers...> tuple of broadcast indexers
    broadcast(TupIndexers iovs) {
        auto ret = broadcast_rec(
                    //empty accumulator
                    std::tuple<>(),
                    //toproc
                    iovs,
                    //full
                    iovs
                    );
        return ret;
    }

} //end .namespace helpers

} //end namespace ndata

#endif /* end of include guard: HELPERS_HPP_YIOQCJSG */
