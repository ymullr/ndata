#ifndef NDINDEXER_HPP_9FICI4GD
#define NDINDEXER_HPP_9FICI4GD

#include <vector>

#include <initializer_list>
#include <cstdlib> //math stuff
#include "vecarray.hpp"


namespace ndata {

/**
 * it is possible to represent the last element of a dimension in a slice
 * by using or offset from the end of the array by using a negative value.
 * 
 * This is a convenience alias
 */
static const long END = -1;


// ( shape[i], strides[i] )
template<size_t N> using SlicesT = vecarray<pair<size_t, long>, N>;


//some forward declarations

template<size_t ndims>
struct indexer;


//please ignore
namespace helpers {

template <long ndimslices>
indexer<ndimslices>
make_indexer_helper(pair<size_t, SlicesT<ndimslices>> pr);

}


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
 * Helper function allowing the user to create an indexer, dimensionalty is infered
 * from the number of arguments.
 */
template <typename... Shape>
auto
make_indexer(Shape... shape) {
    return indexer<sizeof...(shape)>(shape...);
}

template<size_t ndims>
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
        static_assert(ndims == sizeof...(shape)+1,
            "Number of parameters doesn't match the number of dimensions."
        );

        assert(shape0>=0);

        //vecarray size 0 for dynamic version is going to be a problem
        shape_=array_from_argpack<size_t, 0>(vecarray<size_t, 0> (), shape0, shape...);
        strides_ = calc_strides_from_shape(shape_);
    }

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

    //ndindexer(vector<size_t> shape):
    //    shape_(shape)
    //{ 
    //    strides_=calc_strides_from_shape(shape);
    //    //static_assert(ndims == 0, 
    //    //   "This constructor always produce a dynamic vecarray,"
    //    //   " please set ndims template argument to 0 to use this constructor."
    //    //);
    //    assert(shape.size() == ndims or ndims == 0);
    //}

    //ndindexer(array<size_t, ndims> shape):
    //    shape_(shape),
    //    strides_(calc_strides_from_shape(shape)) { }

    //empty constructor for later assignment
    indexer() { };

    /**
     * ndindex may here be passed as a vecarray (custom type), but also work with an 
     * a std::array or an initializer_list (when the number of dimensions is
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

        return indexacc;
    }

    template<typename... SizeT >
    size_t index(size_t i0, SizeT... in){
        static_assert(
            ndims != 0,
            "This overload is only available when the number of dimensions"
            " is known at compile time"
            );
        static_assert(
            ndims == sizeof...(in)+1,
            "Number of parameters doesn't match the number of dimensions."
        );

        return index_rec<0>(start_index_, i0, in...);
    }


    /**
     * Return the stride for the current dimension
     */
    vecarray<size_t, ndims> get_strides() {
        return strides_;
    }

    ///**
    // * Just returns the shape used to construct the class (if any).
    // */
    vecarray<size_t, ndims> get_shape() {
        return shape_;
    }

    size_t stride(size_t i) {
        assert(i<shape_.size());
        return strides_[i];
    }

    size_t shape(size_t i) {
        assert(i<shape_.size());
        return shape_[i];
    }


    size_t size(){
        size_t acc=1;
        for(size_t i=0;i<shape_.size();i++){
            acc*=shape_[i];
        }
        return acc;
    }


    /**
     * For a given ndarray, this will increments the dimensional indices for
     * each dimension in the correct order.
     *
     * This will work with any decent container like std::array, a std::vector or a vecarray.
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

    /**
     * Returns a new ndindexer with eventually a smaller number of dimensions,
     * The new ndindexer computes indices matching the requested slice of the array.
     */
    //template <typename... RngOrIndexT>
    //auto
    //slice(RngOrIndexT ... slice_or_index) 
    //{
    //    return make_ndindexer(
    //            slice_rec<0, 0>(
    //                start_index_,
    //                SlicesT<0>(),
    //                slice_or_index...
    //                )
    //            );
    //}

    /**
     * Returns a new ndindexer with eventually a smaller number of dimensions,
     * The new ndindexer computes indices matching the requested slice of the array.
     */
    template <typename... DimSliceT>
    auto
    slice(long index, DimSliceT ... slice_or_index) 
    {
        return indexer(
                slice_rec<0, 0>(
                    start_index_,
                    SlicesT<0>(),
                    index,
                    slice_or_index...
                    )
                );
    }

    /**
     * Returns a new ndindexer with eventually a smaller number of dimensions,
     * The new ndindexer computes indices matching the requested slice of the array.
     *
     */
    template <typename... DimSliceT>
    auto
    slice(Rng index_range,  DimSliceT ... slice_or_index) 
    {
        return helpers::make_indexer_helper(
                slice_rec<0, 0>(
                    start_index_,
                    SlicesT<0>(),
                    index_range,
                    slice_or_index...
                    )
                );
    }

    ///**
    // * Returns a new ndindexer with eventually a smaller number of dimensions,
    // * The new ndindexer computes indices matching the requested slice of the array.
    // *
    // * array size max 3
    // *
    // */
    //template <typename... DimSliceT>
    //auto
    //slice(array<long, 2> arrslice,  DimSliceT ... slice_or_index) 
    //{
    //    return make_ndindexer(
    //            slice_rec<0, 0>(
    //                start_index_,
    //                SlicesT<0>(),
    //                arrslice,
    //                std::forward<DimSliceT...>(slice_or_index)...
    //                )
    //            );
    //}


    //fortran_order

    //broadcast

    /**
     * produces a ndview. see ndview.
     * To use in for loops
     *
     * auto ind_it = my_ndindexes.slice(3, {2, 6}, {}); //{} is equivalent to {0:-1}
     *
     * for(size_t ind_it.begin(); ind_it < ind_it.end();ind_it++) {
     *    do_stuff(my_vector[ind_it]);
     * }
     *
     * //Or even with a for range loop :
     *
     * for(size_t ind : ind_it) {
     *     do_stuff(my_vector[ind_it]);
     * }
     *
     * Also good for reshaping
     */
    //view

    //own

    private:

    size_t start_index_;

    //for bound checking
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

        //static_assert(
        //    sizeof...(rest) == ndims-idim-1,
        //    ""
        //);
        vecarray<ContentT, idim+1> new_acc = acc.append(i);
        return array_from_argpack<ContentT, idim+1>(new_acc, rest...);
    }

    template<typename ContentT, size_t idim>
    vecarray<ContentT, idim> array_from_argpack(vecarray<ContentT, idim> acc) {
        static_assert(
                idim == ndims,
                ""
                );
        return acc;
    }


    template<size_t idim, typename... SizeT>
    size_t index_rec(size_t acc, long i_idim, SizeT... in){
        assert(acc<size());
        long i_idim_rev = reverse_negative_index(idim, i_idim);
        assert(i_idim_rev < long(shape_[idim]));

        size_t new_acc = acc+strides_[idim]*i_idim_rev;

        assert(new_acc < size());
        return index_rec<idim+1>(new_acc, in...);
    }

    template<size_t idim>
    size_t index_rec(size_t acc){
        static_assert(idim == ndims, "");
        assert(acc<size());
        return acc;
    }

    //now handling various arities and slice length/shorthand possible
    //with template recursion



    /**
     * An actual index is given instead of a slice, only one element is thus
     * included and the ndindexer will have one less dimension.
     *
     * the current shape/stride element is discarded and only the start_index is
     * repositioned to the relevant position
     */
    template <size_t idim, size_t ndimslices, typename... SliceIndex>
    auto
    slice_rec(
            size_t start_ind,
            SlicesT<ndimslices> slices,
            long index,
            SliceIndex... slice_or_index
            ) //-> decltype(slice_rec<idim+1, ndimslices>(size_t, SlicesT<ndimslices>, SliceIndex...))
    {
        start_ind += strides_[idim]*reverse_negative_index(idim, index);

        return slice_rec<idim+1, ndimslices>(start_ind, slices, slice_or_index...);
    }

    //full slice {start, stop, step}
    template <size_t idim, size_t ndimslices, typename... SliceIndex>
    auto
    slice_rec(
            size_t start_ind,
            SlicesT<ndimslices> slices,
            Rng range,
            SliceIndex... slice_or_index
            ) //-> decltype(slice_rec<idim+1, ndimsslices>(size_t, SlicesT<ndimslices>, SliceIndex...))
    {

        //would be trouble
        assert(range.step!=0); 

        start_ind += strides_[idim]*reverse_negative_index(idim, range.start);

        long new_stride = strides_[idim]*range.step;

        array<long, 2> full_slice = {
            labs(reverse_negative_index(idim, range.stop) - reverse_negative_index(idim, range.start))
                /labs(range.step),
            new_stride
        };

        //DOES increment dimension
        return slice_rec<idim+1, ndimslices+1>(start_ind, slices.append(full_slice), slice_or_index...); 
    }

    //termination
    template <size_t idim, size_t ndimslices>
    pair<size_t, SlicesT<ndimslices>>
    slice_rec(size_t start_ind, SlicesT<ndimslices> slices) {

        static_assert(
            idim == ndims,
            ""
        );

        static_assert(ndimslices<=ndims, "");

        return make_pair(start_ind, slices);
    }

    long reverse_negative_index(size_t idim, long ind) { 
        return (ind>=0)? ind : long(shape_[idim])+(ind)+1;
    }

    //done handling slice signatures


};

namespace helpers {

/**
 *
 *
 * tuple : size_t start_index, vecarray<array<long, 2>, ndims_slice> slices):
 *
 */
    template <long ndimslices>
    indexer<ndimslices> make_indexer_helper(pair<size_t, SlicesT<ndimslices>> pr) {

        vecarray<size_t, ndimslices> shape (pr.second.dynsize());
        vecarray<long, ndimslices> strides (pr.second.dynsize());

        for (size_t i = 0; i < pr.second.size(); ++i) {
            shape[i] = pr.second[i].first;
            strides[i] = pr.second[i].second;
        }

        return indexer<ndimslices>(pr.first, shape, strides);

    }

}


}
#endif /* end of include guard: NDINDEXER_HPP_9FICI4GD */
