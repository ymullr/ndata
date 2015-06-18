#ifndef HELPERS_HPP_YIOQCJSG
#define HELPERS_HPP_YIOQCJSG

#include <cstdlib>
#include <utility>
#include "tuple_utilities.hpp"
#include <exception>
#include <type_traits>

#include "vecarray.hpp"

namespace ndata {
namespace helpers {


    using shape_stride_pair = std::pair<long, long>;

    // ( shape[i], strides[i] )
    template<long N> using SliceAcc = vecarray<shape_stride_pair, N>;

    template <typename ... Ts>
    struct static_check_valid_indice_types {
    };

    //struct used to implement "type level function"
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

    //termination
    template <typename Indexer>
    auto
    broadcast_left(Indexer ind, std::tuple<> t) {
        t=t;//silence warning
        return ind;
    }

    /**
     * @brief broadcasts v1 on all elements of the tuple t (all must have compatible dimensions)
     */
    template <typename Indexer, typename ... Indexers>
    auto
    broadcast_left(Indexer v1, std::tuple<Indexers...> t) {
        auto h_t = tuple_utilities::split_ht(std::move(t));
        auto v2 = h_t.first;
        auto tuple_rest = h_t.second;

        auto shape_v1 = v1.get_shape();
        auto shape_v2 = v2.get_shape();

        size_t new_size = std::max(
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
                    throw(std::domain_error("Dimension sizes must match or be equal to 1"));
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

    /**
     * @brief Broadcasts all elements of "all" on the target "target" and returns them.
     *  Fails if dimensions are not compatible.
     */
    template <
            long ndims,
            typename ... Indexers
            >
    auto
    broadcast_on_target(
            indexer<ndims> target,
            std::tuple<Indexers...> all
            )
    {
        return tuple_utilities::tuple_transform(
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
        auto h_t = tuple_utilities::split_ht(std::move(tup_ind));

        //indexer<helpers::static_max_or_dynamic<
        //        decltype(Indexers::shape_)::STATIC_SIZE_OR_DYNAMIC...>::value>
        auto
            indexer_target = broadcast_left(
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
            helpers::make_indexer(indexer_target),
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
        auto tupind_views = tuple_utilities::tuple_transform(
                    [] (auto& tupelt) {
                        return tupelt.as_view();
                    },
                    tup_ind
                );
        return broadcast(tupind_views);
    }

} //end .namespace helpers

} //end namespace ndata

#endif /* end of include guard: HELPERS_HPP_YIOQCJSG */
