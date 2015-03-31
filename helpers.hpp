#ifndef HELPERS_HPP_YIOQCJSG
#define HELPERS_HPP_YIOQCJSG

#include "ndata.hpp"
#include <cstdlib>
#include <exception>

namespace ndata {

namespace helpers {


    using ShapeStridePair = std::pair<size_t, long>;

    // ( shape[i], strides[i] )
    template<size_t N> using SliceAcc = vecarray<ShapeStridePair, N>;

    /**
     *
     *
     * tuple : size_t start_index, vecarray<array<long, 2>, ndims_slice> slices):
     *
     */
    template <long ndimslices>
    indexer<ndimslices> make_indexer_helper(std::pair<size_t, SliceAcc<ndimslices>> pr) {

        vecarray<size_t, ndimslices> shape (pr.second.dynsize());
        vecarray<long, ndimslices> strides (pr.second.dynsize());

        for (size_t i = 0; i < pr.second.size(); ++i) {
            ShapeStridePair sh_st_p = pr.second[i];
            shape[i] = sh_st_p.first;
            strides[i] = sh_st_p.second;
        }

        return indexer<ndimslices>(pr.first, shape, strides);

    }

    template<typename T, size_t static_size_1, size_t static_size_2>
    auto
    make_vecarray_like_biggest(
            vecarray<T, static_size_1> v1,
            vecarray<T, static_size_2> v2
            )
    {
        if (
                static_size_1 != ndata::DYNAMIC_SIZE
            and static_size_2 != ndata::DYNAMIC_SIZE
                )
        {
            return vecarray<T, std::max(static_size_1, static_size_2)>();
        } else {
            return vecarray<T, ndata::DYNAMIC_SIZE> (std::max(v1.size(), v2.size()));
        }
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

    template <typename IndexerOrView, typename... IndexerOrViews>
    auto // vecarray<IndexerOrView, N>
    broadcast(IndexerOrViews... others) {
        constexpr size_t arrsize = sizeof...(others);

        vecarray<IndexerOrView, arrsize> arr_args = array_from_argpack(vecarray<IndexerOrView, 0>(), others...);

        vecarray<IndexerOrView, arrsize> ret;

        for (size_t i = 0; i < ret.size(); ++i) {
            vecarray<IndexerOrView, arrsize> reord;
            reord[0] = arr_args[i];
            size_t io = 1;
            for (size_t i2 = 0; i2 < arrsize; ++i2) {
                if (io != i) {
                    reord[io] = arr_args[i2];
                    io++;
                }
            }
            ret[i] = broadcast_left_list(reord);
        }

        return ret;
    }

    template <typename IndexerOrView, size_t nargs>
    IndexerOrView
    broadcast_left_list(vecarray<IndexerOrView,  nargs> args) {
        vecarray<IndexerOrView, nargs-2> rest;
        for (size_t i = 0; i < rest.size(); ++i) {
            rest[i] = args[i+2];
        }
        return broadcast_left_list(broadcast_left(args[0], args[1]), rest);
    }

    template <typename IndexerOrView>
    IndexerOrView
    broadcast_left_list(vecarray<IndexerOrView,  2> args) {
        return broadcast_left(args[0], args[1]);
    }

    //broadcasts v1 to match v2
    template <typename IndexerOrView>
    IndexerOrView
    broadcast_left(IndexerOrView v1 , IndexerOrView v2) {

        auto shape_v1 = v1.get_shape(),
             shape_v2 = v2.get_shape();

        size_t new_size = max(
                    shape_v1.size(),
                    shape_v1.size()
                    );

        auto strides_v1 = v1.get_strides();
        //auto strides_v2 = v2.get_strides();

        auto new_shape = make_vecarray_like_biggest(v1, v2);

        auto
            new_strides_v1 = new_shape;
            //new_strides_v2 = shape;

        for (size_t is = 0; is < new_size; ++is) {
            size_t i_v2 = v2.get_shape().size()-is-1,
                   i_v1 = v1.get_shape().size()-is-1,
                   i_new = new_size - is -1;

            if (i_v1 < shape_v1.size() && i_v2 < shape_v2.size()) {

                if (shape_v1[i_v1] == shape_v2[i_v2]) {
                    new_shape[i_new] = shape_v1[i_v1];
                    new_strides_v1[i_new] = strides_v1[i_v1];
                    //new_strides_v2[i_new] = strides_v2[i_v2];
                } else if (shape_v1[i_v1] == 1) {
                    new_shape[i_new] = shape_v2[i_v2];
                    new_strides_v1[i_new] = 0;
                    //new_strides_v2[i_new] = strides_v2[i_v2];
                } else if (shape_v2[i_v2] == 1) {
                    new_shape[i_new] = shape_v1[i_v1];
                    new_strides_v1[i_new] = strides_v1[i_v1];
                    //new_strides_v2[i_new] = 0;
                } else {
                    throw(std::domain_error("Dimension size must match or be equal to 1"));
                }

            } else if (i_v2 >= v2.get_shape().size()) {
                new_shape[i_new] = shape_v1[i_v1];
                new_strides_v1[i_new] = strides_v1[i_v1];
                //new_strides_v2[i_new] = 0;
            } else if (i_v1 >= v1.get_shape().size()) {
                new_shape[i_new] = shape_v2[i_v2];
                new_strides_v1[i_new] = 0;
                //new_strides_v2[i_new] = strides_v2[i_v2];
            } else {
                assert("false");
            }

            //if (is < shape_.size()) {
            //    size_t i_this = shape_.size() - is - 1;

            //    assert(
            //                targ_shape[i_v2] == shape_[i_this]
            //             or shape_[i_this] == 1
            //        );
            //    new_shape[i_v2] = shape_[i_this];
            //    if (shape_[i_this] == 1) {
            //        //that's how this dimension is "extended" without creating new data
            //        new_strides[i_v2] = 0;
            //    } else {
            //        new_strides[i_v2] = strides_[i_this];
            //    }
            //} else {
            //    new_shape[i_v2] = targ_shape[i_v2];
            //    new_strides[i_v2] = 0;
            //}

        }

        return v1.reindex(new_shape, new_strides_v1);
    }
}
}

#endif /* end of include guard: HELPERS_HPP_YIOQCJSG */
