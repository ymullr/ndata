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
struct nview: indexer<ndims> {

    ContainerT data_;

    //nview(ContainerT data, vecarray<size_t, ndims> shape):
    //    indexer<ndims>(shape),
    //    data_(data)
    //{ }

    nview(ContainerT data, indexer<ndims> idxr):
        indexer<ndims>(idxr),
        data_(data)
    { }


    template<typename... SizeT>
    nview(ContainerT data, size_t shape0, SizeT... shape):
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
    nview<ContainerT, T, ndims_broad>
    broadcast(nview<ContainerT, T, ndims_broad> target) {
        return broadcast_helper(this, target);
    }

    template <typename... NdViews, typename ReturnT, typename... Ts>
    nview<ContainerT, T, ndims>
    transform(NdViews... ndv, std::function<ReturnT(Ts... args)> func) {
        nview<ContainerT, T, ndims> retdat = this;
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
        //use slice method of the parent class and use it to view this.data_
        //and return a new slice
        return  indexer<ndims>::slice(
                    slice_or_index...
                ).view(data_);
    }

    /**
     * used by broadcast
     */
    template<size_t new_ndims>
    nview<ContainerT, T, new_ndims>
    reshape(vecarray<size_t, ndims> new_shape, vecarray<long, ndims> new_strides) {
        return nview(data_, indexer<new_ndims>(indexer<ndims>::start_index_, new_shape, new_strides));
    }
};



}

#endif /* end of include guard: NDVIEW_HPP_1LWJBSCE */

