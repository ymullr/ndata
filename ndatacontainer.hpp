#ifndef NDVIEW_HPP_1LWJBSCE
#define NDVIEW_HPP_1LWJBSCE

#include <functional>

#include "vecarray.hpp"
#include "ndata.hpp"
#include "indexer.hpp"
#include "helpers.hpp"

namespace ndata {

/**
 * An special ndindexer holding its data (or a pointer to it).
 * inherits from indexer, so it as all of its methods,
 *
 * in addition to [], val,  and transform.
 *
 * The method slice and broadcast return nviews instead of indexer
 **/
template<typename ContainerT, typename T, size_t ndims>
struct ndatacontainer: indexer<ndims> {

    ContainerT data_;

    //nview(ContainerT data, vecarray<size_t, ndims> shape):
    //    indexer<ndims>(shape),
    //    data_(data)
    //{ }

    ndatacontainer(ContainerT data, indexer<ndims> idxr):
        indexer<ndims>(idxr),
        data_(data)
    { }


    template<typename... SizeT>
    ndatacontainer(ContainerT data, size_t shape0, SizeT... shape):
        indexer<ndims>(shape0, shape...),
        data_(data)
    {  }

    T&
    operator[](size_t i) {
        return data_[i];
    }

    template <typename... SizeT>
    T&
    val(SizeT... indices) {
        return data_[indexer<ndims>::index(indices...)];
    }

    T&
    val(vecarray<size_t, ndims> indices) {
        return data_[indexer<ndims>::index(indices)];
    }

    template <size_t ndims_broad>
    ndatacontainer<ContainerT, T, ndims_broad>
    broadcast(ndatacontainer<ContainerT, T, ndims_broad> target) {
        return broadcast_helper(this, target);
    }

    template <typename... NdViews, typename ReturnT, typename... Ts>
    ndatacontainer<ContainerT, T, ndims>
    transform(NdViews... ndv, std::function<ReturnT(Ts... args)> func) {
        ndatacontainer<ContainerT, T, ndims> retdat = this;
        for (size_t ielt = 0; ielt < indexer<ndims>::size(); ++ielt) {
            retdat[ielt] = func(ndv[ielt]...);
        }
        return retdat;
    }

    /**
     * Returns a new ndindexer with eventually a smaller number of dimensions,
     * The new ndindexer computes indices matching the requested slice of the array.
     */
    template <typename... DimSliceT>
    auto
    slice(DimSliceT ... slice_or_index) {
        //use slice method of the parent class and use it to own_data this.data_
        //and return a new slice
        return  indexer<ndims>::slice(
                    slice_or_index...
                ).template own_data<ContainerT, T>(data_);
    }

    /**
     * Returns a new ndindexer with eventually a smaller number of dimensions,
     * The new ndindexer computes indices matching the requested slice of the array.
     *
     * Unlike the slice(...) method, slice_view only returns an ndataview which only
     * holds a pointer to its data.
     */
    template <typename... DimSliceT>
    auto
    slice_view(DimSliceT ... slice_or_index) {
        //use slice method of the parent class and use it to own_data this.data_
        //and return a new slice
        ndataview<T, ndims> ret (&data_[0], indexer<ndims>::slice(
                    slice_or_index...));
        return ret;
    }

    ndataview<T, ndims>
    to_view() {
        //use slice method of the parent class and use it to own_data this.data_
        //and return a new slice
        return ndataview<T, ndims>(&data_[0], *this);
    }

    //
    // * Returns a new ndindexer with eventually a smaller number of dimensions,
    // * The new ndindexer computes indices matching the requested slice of the array.
    // */
    //auto
    //to_view() {
    //    //use slice method of the parent class and use it to own_data this.data_
    //    //and return a new slice
    //    return this->template view_data<T>(data_);
    //}


    /**
     * used by broadcast
     */
    template<size_t new_ndims>
    ndatacontainer<ContainerT, T, new_ndims>
    reshape(vecarray<size_t, new_ndims> new_shape, vecarray<long, new_ndims> new_strides) {
        return ndatacontainer<ContainerT, T, new_ndims>(data_, indexer<new_ndims>(this->start_index_, new_shape, new_strides));
    }
};



}

#endif /* end of include guard: NDVIEW_HPP_1LWJBSCE */

