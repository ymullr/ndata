#ifndef NDVIEW_HPP_1LWJBSCE
#define NDVIEW_HPP_1LWJBSCE

#include <functional>

#include "vecarray.hpp"
#include "ndata.hpp"
#include "nindexer.hpp"
#include "nhelpers.hpp"

namespace ndata {

namespace helpers {

//returns sliced indexer, general case
template<typename ContainerT, typename T, long ndimslices, long ndims_base>
struct return_sliced_view_or_value {

    static
    ndataview<T, ndimslices>
    do_it(
            std::pair<size_t, helpers::SliceAcc<ndimslices>> pr,
            ContainerT&& data
            )
    {
        return ndataview<T, ndimslices>(
                    make_indexer_from_slices_helper(pr),
                    std::move(data)
                    );
    }
};

//specialisation when pr contains no slice (only straight indexes were given
//in that case returns a reference to a numerical value
template<typename ContainerT, typename T, long ndims_base>
struct return_sliced_view_or_value<ContainerT, T, 0, ndims_base> {

    static
    T&
    do_it(
            std::pair<size_t, helpers::SliceAcc<0>> pr,
            ContainerT& data
            )
    {
        return data[pr.first];
    }
};



}

/**
 * An special ndindexer holding its data (or a pointer to it).
 * inherits from indexer, so it as all of its methods,
 *
 * in addition to [], val,  and transform.
 *
 * The method slice and broadcast return nviews instead of indexer
 **/
template<typename ContainerT, typename T, long ndims>
struct ndatacontainer: indexer<ndims> {

    ContainerT data_;

    //nview(ContainerT data, vecarray<size_t, ndims> shape):
    //    indexer<ndims>(shape),
    //    data_(data)
    //{ }

    /**
     * @brief Construct from indexer and data
     * @param idxr
     * @param data
     */
    ndatacontainer(indexer<ndims> idxr, ContainerT data):
        indexer<ndims>(idxr),
        data_(data)
    { }

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


    //
    ///**
    // * Returns a new ndindexer with eventually a smaller number of dimensions,
    // * The new ndindexer computes indices matching the requested slice of the array.
    // */
    //template <typename... DimSliceT>
    //auto
    //slice(DimSliceT ... slice_or_index) {
    //    //use slice method of the parent class and use it to own_data this.data_
    //    //and return a new slice
    //    return  indexer<ndims>::slice(
    //                slice_or_index...
    //            ).template own_data<ContainerT, T>(data_);
    //}

    template <typename... IndexOrRangeT>
    auto
    operator()(IndexOrRangeT ... index_or_range) {
        helpers::SliceAcc<0> sa;

        //type = std::pair<size_t start_index, helpers::SliceAcc<STATIC_SIZE_OR_DYNAMIC>>
        auto pr = indexer<ndims>::template slice_rec<0, 0>(
            indexer<ndims>::start_index_,
            sa,
            index_or_range...
        );

        return helpers::return_sliced_view_or_value<ContainerT, T, pr.second.STATIC_SIZE_OR_DYNAMIC, ndims>::do_it(
                    pr,
                    std::move(data_)
        );

    }

    ///**
    // *
    // * Returns a new ndindexer with eventually a smaller number of dimensions,
    // * The new ndindexer computes indices matching the requested slice of the array.
    // *
    // */
    //template <typename... DimSliceT>
    //auto
    //slice(Rng index_range,  DimSliceT ... slice_or_index)  {
    //    return helpers::make_indexer_from_slices_helper(
    //            slice_rec<0, 0>(
    //                start_index_,
    //                helpers::SliceAcc<0>(),
    //                index_range,
    //                slice_or_index...
    //                )
    //            );
    //}



    //
    //**
    // * Returns a new ndindexer with eventually a smaller number of dimensions,
    // * The new ndindexer computes indices matching the requested slice of the array.
    // *
    // * Unlike the slice(...) method, slice_view only returns an ndataview which only
    // * holds a pointer to its data.
    // */
    //template <typename... DimSliceT>
    //auto
    //slice_view(DimSliceT ... slice_or_index) {
    //    //use slice method of the parent class and use it to own_data this.data_
    //    //and return a new slice
    //    ndataview<T, ndims> ret (
    //                indexer<ndims>::slice(slice_or_index...),
    //                &data_[0]
    //    );
    //    return ret;
    //}

    //default conversion operator to a view
    operator ndataview<T, ndims>() {
        //use slice method of the parent class and use it to own_data this.data_
        //and return a new slice
        return ndataview<T, ndims>(*this, &data_[0]);
    }

    //convenience function to call default conversion without specifying template parameters
    auto to_view() {
        return ndataview<T, ndims>(*this);
    }

    void fill(T val) {
        for (size_t i = 0; i < this->size(); ++i) {
            operator[](i) = val;
        }
    }

    /**
     * used by broadcast
     */
    template<size_t new_ndims>
    ndatacontainer<ContainerT, T, new_ndims>
    reshape(vecarray<size_t, new_ndims> new_shape, vecarray<long, new_ndims> new_strides) {
        return ndatacontainer<ContainerT, T, new_ndims>(indexer<new_ndims>(this->start_index_, new_shape, new_strides), data_);
    }
};



}

#endif /* end of include guard: NDVIEW_HPP_1LWJBSCE */

