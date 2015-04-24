#ifndef NDVIEW_HPP_1LWJBSCE
#define NDVIEW_HPP_1LWJBSCE

#include <functional>

#include "vecarray.hpp"
#include "ndata.hpp"
#include "nindexer.hpp"
#include "nhelpers.hpp"

namespace ndata {

template <typename ContainerT, typename T, long ndims>
ndatacontainer<ContainerT, T, ndims>
make_ndatacontainer(
        indexer<ndims> idx,
        ContainerT data
        )
{
    return ndatacontainer<ContainerT, T, ndims>(
                idx,
                data
                );
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

    //contained data
    ContainerT data_;

    //static members
    typedef ContainerT type_ContainerT;

    typedef T type_T;

    //static constexpr long STATIC_SIZE_OR_DYNAMIC = ndims;


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

    //template <typename... SizeT>
    //T&
    //val(SizeT... indices) {
    //    return data_[indexer<ndims>::index(indices...)];
    //}

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


public :

    /**
     * Returns a view
     */
    template <typename... IndexOrRangeT>
    auto
    slice(IndexOrRangeT ... index_or_range) {
        return  make_ndatacontainer<T*, T>(
                    this->index_slice(index_or_range...),
                    &data_[0]
                );
    }

    template <typename ... Long>
    T&
    operator()(Long... indices) {
        size_t index = indexer<ndims>::index(indices...);
        return data_[index];
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

    template <typename ContainerT_rhs, typename T_rhs, long ndims_rhs>
    void
    assign(ndatacontainer<ContainerT_rhs, T_rhs, ndims_rhs> rhs) {

        static_assert(ndims_rhs == ndims or ndims_rhs == DYNAMIC_SIZE, "");

        auto ndi = this->ndindex(0);
        for (size_t i = 0; i < this->size(); ++i) {
            this->operator()(ndi) = rhs(ndi);
            this->increment_ndindex(ndi);
        }
    }


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

