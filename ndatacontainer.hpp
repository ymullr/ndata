#ifndef NDVIEW_HPP_1LWJBSCE
#define NDVIEW_HPP_1LWJBSCE

#include <functional>

#include "vecarray.hpp"
#include "ndata.hpp"
#include "nindexer.hpp"
#include "ndata/helpers.hpp"
#include "ndata/loops.hpp"

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
        data_(std::move(data))
    { }

    /**
     * @brief empty constructor for later assignment
     */
    ndatacontainer() { }


    T&
    operator[](size_t i) {
        return data_[i];
    }

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


    /**
     * @brief elementwise copy of the values of rhs to the internal data. Doesn't perform broadcasting,
     * use assign_transform if you want to benefit from broadcasting
     */
    template <int loop_type = SERIAL, typename ContainerT_rhs, typename T_rhs, long ndims_rhs>
    void
    assign(ndatacontainer<ContainerT_rhs, T_rhs, ndims_rhs> rhs) {
        static_assert(ndims_rhs == ndims or ndims_rhs == DYNAMICALLY_SIZED, "");

        for (size_t i = 0; i < this->get_shape().size(); ++i) {
            assert(rhs.get_shape()[i] == this->get_shape()[i]);
        }

        nforeach<loop_type>(
                    std::tie(*this, rhs),
                    [] (T& this_val, T_rhs rhs_val) {
                        this_val = rhs_val;
                    }
            );
    }

    /**
     * @brief equivalent to calling .assign(ntransform(...)) but skips the extra temporary
     */
    template <int loop_type = SERIAL, typename... Ndatacontainer, typename FuncT>
    void
    assign_transform(std::tuple<Ndatacontainer...> ndata_tup, FuncT func)  {

        //static_assert(ndims_rhs == ndims or ndims_rhs == DYNAMICALLY_SIZED, "");

        auto ndata_tuple_bcviews = helpers::broadcast_views(ndata_tup);

        nforeach_base<loop_type>(
                    std::tuple_cat(
                        std::make_tuple(this->to_view()),
                        ndata_tuple_bcviews
                        ),
                    [func] (T& this_val, auto ... rhs_vals) {
                        this_val = func(rhs_vals...);
                    }
            );
    }


    template <typename... Ndatacontainer, typename FuncT>
    void
    assign_transform_parallel(std::tuple<Ndatacontainer...> ndata_tup, FuncT func)  {
        assign_transform<PARALLEL>(ndata_tup, func);
    }

    //operator ndataview<T, ndims>() {
    //    //use slice method of the parent class and use it to own_data this.data_
    //    //and return a new slice
    //    return ndataview<T, ndims>(*this, &data_[0]);
    //}

    //convenience function to call default conversion without specifying template parameters
    auto to_view() {
        return ndataview<T, ndims>(*this, &data_[0]);
    }

    //void fill(T val) {
    //    for (size_t i = 0; i < this->size(); ++i) {
    //        operator[](i) = val;
    //    }
    //}

    /**
     * Used by broadcast
     *
     * TODO
     *
     */
    template<long new_ndims>
    ndataview<T, new_ndims>
    reshape(vecarray<long, new_ndims> new_shape, vecarray<long, new_ndims> new_strides) {
        return ndataview<T, new_ndims>(indexer<new_ndims>(this->start_index_, new_shape, new_strides), &data_[0]);
    }

};



}

#endif /* end of include guard: NDVIEW_HPP_1LWJBSCE */

