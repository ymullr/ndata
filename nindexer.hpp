#ifndef NDINDEXER_HPP_9FICI4GD
#define NDINDEXER_HPP_9FICI4GD

#include <vector>

#include <initializer_list>
#include <cstdlib>
#include <type_traits>

#include "vecarray.hpp"
#include "ndata.hpp"
#include "nhelpers.hpp"


namespace ndata {

/**
 * it is possible to represent the last element of a dimension in a slice
 * by using or offset from the end of the array by using a negative value.
 * 
 * This is a convenience alias
 */
static const long END = -1;


/**
 * A range object to pass to ndindexer.slice
 *
 * CamelCased to make it stand out!
 */
struct Rng {
    long start;
    long stop;
    long step;

    Rng ():
    start(0),
    stop(END),
    step(1)
    {};

    explicit
    Rng (long index):
    start(index),
    stop(index+1),
    step(1)
    { };

    Rng (long start, long stop):
    start(start),
    stop(stop),
    step(1)
    { };

    Rng (long start, long stop, long step):
    start(start),
    stop(stop),
    step(step)
    { }
};

/**
 * Helper function allowing the user to create an indexer, dimensionality is infered
 * from the number of arguments.
 */
template <typename... SizeT>
auto
make_indexer(size_t shape0, SizeT ... shapen) {
    return indexer<sizeof...(shapen)+1>(shape0, shapen...);
}

///**
// * Helper function allowing the user to create an indexer, create indexer from its members
// */
//template <size_t ndims>
//auto
//make_indexer(vecarray<size_t, ndims> shape, vecarray<long, ndims> strides, size_t start_index) {
//    return indexer<ndims>(start_index, shape, strides);
//}

template<long ndims>
auto
make_indexer(indexer<ndims> idx) {
    return idx;
}


namespace helpers {

    //forward decl
    /**
     * tuple : size_t start_index, vecarray<array<long, 2>, ndims_slice> slices):
     */
    template <long ndimslices>
    indexer<ndimslices> make_indexer_from_slices_helper(std::pair<size_t, SliceAcc<ndimslices>> pr);

}


template<long ndims>
struct indexer {

    template<typename... SizeT>
    indexer(size_t shape0, SizeT... shape):
        start_index_(0)
    {
        static_assert(
            ndims != 0,
            "This overload is only available when the number of dimensions"
            " is a compile time constant"
            );
        static_assert(
                    ndims == sizeof...(shape)+1,
            "Number of parameters doesn't match the number of dimensions."
        );

        assert(shape0>=0);

        //vecarray size 0 for dynamic version is going to be a problem
        shape_=helpers::array_from_argpack<size_t, 0>(vecarray<size_t, 0> (), shape0, shape...);
        strides_ = calc_strides_from_shape(shape_);
    }

    //initialize from shape with default strides
    //you can also pass a vector or array instead, it will autoconvert
    indexer(vecarray<size_t, ndims> shape):
        start_index_(0),
        shape_(shape),
        strides_(calc_strides_from_shape(shape)) {
    }

    //construct from members
    indexer(size_t start_index, vecarray<size_t, ndims> shape, vecarray<long, ndims> strides):
        start_index_(start_index),
        shape_(shape),
        strides_(strides) { };


    //empty constructor for later assignment
    indexer() { };


    /**
     * ndindex may here be passed as a vecarray (custom type), but also work with an 
     * a std::array or a std::vector (when the number of dimensions is
     * known at compile time), or an std::vector when the number of dimensions
     * is known at runtime. 
     */
    size_t index(vecarray<size_t, ndims> ndindex){
        assert(ndindex.size() == shape_.size());

        size_t indexacc=start_index_;

        //looping on array dimensions
        for (size_t i = 0; i < ndindex.size(); ++i) {
            //bound checking when debugging
            assert(ndindex[i] < shape_[i]);

            indexacc+=ndindex[i]*strides_[i];
        }

        assert(indexacc < this->unsliced_size());

        return indexacc;
    }

    size_t index(vecarray<long, ndims> ndindex) {
        vecarray<size_t, ndims> rev_index (ndindex.dynsize());
        for (size_t i = 0; i < ndindex.size(); ++i) {
            rev_index[i] = reverse_negative_index(ndindex[i]);
        }
        return index(rev_index);
    }

    template <typename ... IndexT>
    size_t index(IndexT ... indices) {
        static_assert(sizeof...(indices) == ndims, "Number of indices doesn't match dimensionality");
        static_assert(ndims != DYNAMIC_SIZE, "Please use an overload taking a vecarray with dynamically dimensioned arrays");

        static_assert(
                    helpers::static_check_valid_indice_types<IndexT...>::value,
                    "Indices passed as argument must be of integral type"
                    );

        auto pr = slice_rec<0, 0>(
            indexer<ndims>::start_index_,
            helpers::SliceAcc<0>(),
            indices...
        );

        static_assert(
            decltype(pr.second)::STATIC_SIZE_OR_DYNAMIC == 0
            //decltype(pr.second)::STATIC_SIZE_OR_DYNAMIC == DYNAMIC_SIZE
                    , "");

        assert(pr.first < this->unsliced_size());

        return pr.first;
    }


    /**
     * index_or_range may be a long or Rng
     */
    template <typename... IndexOrRangeT>
    auto
    index_slice(IndexOrRangeT ... index_or_range) {
        //type = std::pair<size_t start_index, helpers::SliceAcc<STATIC_SIZE_OR_DYNAMIC>>
        auto pr = slice_rec<0, 0>(
            start_index_,
            helpers::SliceAcc<0>(),
            index_or_range...
        );

        return  helpers::make_indexer_from_slices_helper(pr);
    }



    //template<typename... SizeT >
    //size_t index(size_t i0, SizeT... in){
    //    static_assert(
    //        ndims != DYNAMIC_SIZE,
    //        "This overload is only available when the number of dimensions"
    //        " is known at compile time"
    //        );
    //    static_assert(
    //        ndims == sizeof...(in)+1,
    //        "Number of parameters doesn't match the number of dimensions."
    //    );

    //    return index_rec<0>(start_index_, i0, in...);
    //}


    size_t get_start_index() {
        return start_index_;
    }

    vecarray<size_t, ndims> get_shape() {
        return shape_;
    }

    /**
     * Return the stride for the current dimension
     */
    vecarray<long, ndims> get_strides() {
        return strides_;
    }

    //size_t stride(size_t i) {
    //    assert(i<shape_.size());
    //    return strides_[i];
    //}

    //size_t shape(size_t i) {
    //    assert(i<shape_.size());
    //    return shape_[i];
    //}


    size_t size(){
        size_t acc=1;
        for(size_t i=0;i<shape_.size();i++){
            acc*=shape_[i];
        }
        return acc;
    }

    /**
     * Underlying container size, bigger than size() if sliced
     */
    size_t unsliced_size(){
        size_t acc=0;
        for(size_t i=0;i<shape_.size();i++){
            acc+=shape_[i]*std::abs(strides_[i]);
        }
        return acc;
    }

    /**
     * For a given ndarray, this will increments the dimensional indices for
     * each dimension in the correct order.
     *
     * This will work with stl containers such as std::array, a std::vector or a vecarray.
     *
     * Useful to avoid deeply nested loop with a large number of dimension by
     * folding everything into one. Also makes things easier when the number of
     * dimensions is only known at runtime.
     *
     * Use at the end of a for loop when looping on all the elements of an ndarray, like :
     *
     * ndindexer<3> ndindexer (N1, N2, N3);
     *
     * vector<float> ndarray (ndindexer.size());
     *
     * //your multidimensional index
     * //all dimensional indices are 0 initialized
     * vector<size_t> ndi (3, 0); //will also work with a std::array or vecarray
     *
     * for (size_t i = 0; i < ndarray.size(); ++i) * {
     *
     *     //here ndindex.index(ndi) == i
     *
     *     ndindex.increment_ndindex(ndi)
     * }
     */
    template<class ContainerT>
    void
    increment_ndindex(ContainerT & ndindex){
        assert(ndindex.size() == shape_.size());
        increment_ndindex(&ndindex[0]);
    }

    /**
     * Raw pointer version
     */
    template<typename SizeT>
    void
    increment_ndindex(SizeT * ndindex){
        //would also be doable with a bunch of modulo (%) ops but that might be slower/less explicit
        for (long idim = shape_.size()-1; idim >= 0; --idim) {
            //update.back-most dimensional index unless its maxxed
            if (ndindex[idim] < shape_[idim]-1) {
                ndindex[idim]++;
                return;  //were done
            } else if (ndindex[idim] == shape_[idim]-1) { //if maxxed
                ndindex[idim]=0; //reset this dimensional index to 0 
                //continue loop, skip to next outer dimension
                //if this is already the outermost dimension then the index is simply reset
                //no error since this is probably going to happen during the
                //last iteration of an eventual for loop
            }
        }
    }

    /**
     * Probably quite slow, use the increment_ndindex function for more speed
     *
     * Should be implemented with some modulo ops
     */
    vecarray<size_t, ndims> 
    ndindex(size_t flatindex) {

        vecarray<size_t, ndims> ndindex (shape_);

        ndindex.fill(0);

        for (size_t i = 0; i <= flatindex; ++i) {
            increment_ndindex(ndindex);
        }

        return ndindex;
    }


    ///**
    // * Returns a new ndindexer with eventually a smaller number of dimensions,
    // * The new ndindexer computes indices matching the requested slice of the array.
    // */
    //template <typename... DimSliceT>
    //auto
    //index_slice(long index, DimSliceT ... slice_or_index) {
    //    return helpers::make_indexer_from_slices(
    //            slice_rec<0, 0>(
    //                start_index_,
    //                helpers::SliceAcc<0>(),
    //                index,
    //                slice_or_index...
    //                )
    //            );
    //}

    ///**
    // * Returns a new ndindexer with eventually a smaller number of dimensions,
    // * The new ndindexer computes indices matching the requested slice of the array.
    // *
    // */
    //template <typename... DimSliceT>
    //auto
    //index_slice(Rng index_range,  DimSliceT ... slice_or_index)  {
    //    return helpers::make_indexer_from_slices(
    //            slice_rec<0, 0>(
    //                start_index_,
    //                helpers::SliceAcc<0>(),
    //                index_range,
    //                slice_or_index...
    //                )
    //            );
    //}

    //TODO fortran_order

    ///**
    // * produces a ndatacontainer. see ndatacontainer.
    // * To use in for loops
    // *
    // * auto ind_it = my_ndindexes.slice(3, {2, 6}, {}); //{} is equivalent to {0:-1}
    // *
    // * for(size_t ind_it.begin(); ind_it < ind_it.end();ind_it++) {
    // *    do_stuff(my_vector[ind_it]);
    // * }
    // *
    // * //Or even with a for range loop :
    // *
    // * for(size_t ind : ind_it) {
    // *     do_stuff(my_vector[ind_it]);
    // * }
    // *
    // */
    //template<typename ContainerT, typename T>
    //ndatacontainer<ContainerT, T, ndims>
    //own_data(ContainerT data) {
    //    return ndatacontainer<ContainerT, T, ndims>(*this, data);
    //}

    //**
    // * produces a ndataview, which is like an ndatacontainer but only
    // * holds a pointer to its data. See ndatacontainer.
    // */
    //template<typename ContainerT, typename T>
    //ndataview<T, ndims>
    //view_data(ContainerT data) {
    //    return ndataview<T, ndims>(&data[0], *this);
    //}

    /**
     * used by broadcast
     */
    template<size_t new_ndims>
    indexer<new_ndims>
    reshape(vecarray<size_t, ndims> new_shape, vecarray<long, ndims> new_strides) {
        return indexer<new_ndims>(start_index_, new_shape, new_strides);
    }

    protected:

    size_t start_index_;

    //kept for bound checking
    vecarray<size_t, ndims> shape_; 

    vecarray<long, ndims> strides_;

    static
    vecarray<long, ndims> 
    calc_strides_from_shape(vecarray<size_t, ndims> shape) {

        vecarray<long, ndims> ret (shape.dynsize());

        for (size_t i = 0; i < shape.size(); ++i) {

            long this_stride=1;
            for(size_t j=i+1;j<shape.size();j++){
                this_stride*=shape[j];
            }


            ret[i] = this_stride;
        }

        return ret;
    }


protected:


    //template<size_t idim, typename... SizeT>
    //size_t index_rec(size_t acc, long i_idim, SizeT... in){
    //    assert(acc<size());
    //    long i_idim_rev = reverse_negative_index(idim, i_idim);
    //    assert(i_idim_rev < long(shape_[idim]));

    //    size_t new_acc = acc+strides_[idim]*i_idim_rev;

    //    return index_rec<idim+1>(new_acc, in...);
    //}

    //template<size_t idim>
    //size_t index_rec(size_t acc){
    //    static_assert(idim == ndims, "");
    //    return acc;
    //}

    //now handling various arities and slice length/shorthand possible
    //with template recursion

    /**
     * An actual index is given instead of a slice, only one element is thus
     * included and the ndindexer will have one less dimension.
     *
     * the current shape/stride element is discarded and only the start_index is
     * repositioned to the relevant position
     */
    template <
            size_t idim,
            long ndimslices,
            typename IntegerType,
            typename std::enable_if<std::is_integral<IntegerType>::value,
                IntegerType>::type* = nullptr,
            typename... SliceIndex
            >
    auto
    slice_rec(
            size_t start_ind,
            helpers::SliceAcc<ndimslices> slices,
            IntegerType index,
            SliceIndex... slice_or_index
            )
    {
        start_ind += strides_[idim]*reverse_negative_index(idim, index);

        return slice_rec<idim+1, ndimslices>(start_ind, slices, slice_or_index...);
    }

    //full slice {start, stop, step}
    template <size_t idim, long ndimslices, typename... SliceIndex>
    auto
    slice_rec(
            size_t start_ind,
            helpers::SliceAcc<ndimslices> slices,
            Rng range,
            SliceIndex... slice_or_index
            ) //-> decltype(slice_rec<idim+1, ndimsslices>(size_t, SlicesT<ndimslices>, SliceIndex...))
    {

        //would be trouble
        assert(range.step!=0); 

        start_ind += strides_[idim]*reverse_negative_index(idim, range.start);

        long new_stride = strides_[idim]*range.step;

        helpers::ShapeStridePair full_slice = make_pair(
            labs(reverse_negative_index(idim, range.stop) - reverse_negative_index(idim, range.start))
                /labs(range.step),
            new_stride
        );

        //DOES increment dimension
        return slice_rec<idim+1, ndimslices+1>(start_ind, slices.append(full_slice), slice_or_index...); 
    }

    //termination
    template <size_t idim, long ndimslices>
    std::pair<size_t, helpers::SliceAcc<ndimslices>>
    slice_rec(size_t start_ind, helpers::SliceAcc<ndimslices> slices) {

        static_assert(
            idim == ndims,
            ""
        );

        static_assert(ndims == -1 or ndimslices<=ndims, "");

        return make_pair(start_ind, slices);
    }


    //TODO template specialization based on signedness
    long reverse_negative_index(size_t idim, long ind) { 
        return (ind>=0)? ind : long(shape_[idim])+(ind)+1;
    }

    //done handling slice signatures


};


namespace helpers {

    /**
     * tuple : size_t start_index, vecarray<array<long, 2>, ndims_slice> slices):
     */
    template <long ndimslices>
    indexer<ndimslices> make_indexer_from_slices_helper(std::pair<size_t, SliceAcc<ndimslices>> pr) {

        vecarray<size_t, ndimslices> shape (pr.second.dynsize());
        vecarray<long, ndimslices> strides (pr.second.dynsize());

        for (size_t i = 0; i < pr.second.size(); ++i) {
            ShapeStridePair sh_st_p = pr.second[i];
            shape[i] = sh_st_p.first;
            strides[i] = sh_st_p.second;
        }

        return indexer<ndimslices>(pr.first, shape, strides);

    }

}

}
#endif /* end of include guard: NDINDEXER_HPP_9FICI4GD */
