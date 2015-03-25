#include <vector>

#include <initializer_list>
#include "vecarray.hpp"

#ifndef NDINDEXER_HPP_9FICI4GD
#define NDINDEXER_HPP_9FICI4GD

template<size_t ndims>
struct ndindexer {

    template<typename... SizeT>
    ndindexer(size_t shape0, SizeT... shape):
        start_index_(0)
    {

        static_assert(
            ndims != 0,
            "This overload is only available when the number of dimensions is a compile time constant"
            );
        static_assert(ndims == sizeof...(shape)+1,
            "Number of parameters doesn't match the number of dimensions."
        );

        shape_=array_from_argpack<size_t, ndims>(shape0, shape...);
        strides_ = calc_strides_from_shape(shape_);
    }

    ndindexer(vecarray<size_t, ndims> shape):
        start_index_(0),
        shape_(shape),
        strides_(calc_strides_from_shape(shape)) { 
        }

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
    ndindexer() { };

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
            "This overload is only available when the number of dimensions is known at compile time"
            );
        static_assert(
            ndims == sizeof...(in)+1,
            "Number of parameters doesn't match the number of dimensions."
        );

        return index_rec<0>(0ul, i0, in...);
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
    vecarray<size_t, ndims> ndindex(size_t flatindex) {

        vecarray<size_t, ndims> ndindex (shape_);

        ndindex.fill(0);

        for (size_t i = 0; i <= flatindex; ++i) {
            increment_ndindex(ndindex);
        }

        return ndindex;
    }



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

    //private:
    public:

    // ( shape[i], strides[i] )
    template<size_t N> using SlicesT = vecarray<pair<size_t, long>, N>;

    size_t start_index_;

    //for bound checking
    vecarray<size_t, ndims> shape_; 

    vecarray<long, ndims> strides_;

    /**
     *
     * Private constructor for slicing
     *
     * tuple : size_t start_index, vecarray<array<long, 2>, ndims_slice> slices): 
     *
     */
    template <size_t ndimslices>
    ndindexer(pair<size_t, SlicesT<ndimslices>> pr) :
        start_index_(pr.first),
        shape_(pr.second.dynsize()),
        strides_(pr.second.dynsize()) 
    {
        for (size_t i = 0; i < pr.second.size(); ++i) {
            shape_[i] = pr.second[i][0];
            strides_[i] = pr.second[i][1];
        }
    }

    static
    vecarray<long, ndims> 
    calc_strides_from_shape(vecarray<size_t, ndims> shape) {

        vecarray<long, ndims> ret (shape);

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
    template<size_t ndims_variadic, typename ContentT, typename... ContentTPackT>
    vecarray<ContentT, ndims_variadic> array_from_argpack(ContentT i, ContentTPackT... rest) {

        static_assert(ndims_variadic != 0, "Cannot be used when ndims = 0 (Dynamic)");

        static_assert(
            sizeof...(rest) == ndims_variadic-1,
            ""
        );

        vecarray<ContentT, ndims_variadic> ret;
        ret[0] = i;

        vecarray<ContentT, ndims_variadic-1> rest_arr = array_from_argpack<ContentT, ndims_variadic-1>(rest...);

        for (size_t irest = 0; irest < rest_arr.size(); ++irest) {
            ret[irest+1] = rest_arr[irest];
        }

        return ret;
    }

    template<size_t ndims_variadic, typename ContentT>
    array<ContentT, 0> array_from_argpack() {
        //static_assert(ndims_variadic == 0, 
        //        );
        array<ContentT, 0> ret;
        return ret;
    }


    template<size_t idim, typename... SizeT>
    size_t index_rec(size_t acc, size_t i_idim, SizeT... in){
        return index_rec<idim+1>(acc+strides_[idim]*i_idim, in...);
    }

    template<size_t idim>
    size_t index_rec(size_t acc){
        static_assert(idim == ndims, "");
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
        start_ind += strides_[idim]*reverse_negative_index(index);

        return slice_rec<idim+1, ndimslices>(start_ind, slices, slice_or_index...);
    }


    // {} type of slice, where all elements on a dimension are included
    template <size_t idim, size_t ndimslices, typename... SliceIndex>
    auto
    slice_rec(
            size_t start_ind,
            SlicesT<ndimslices> slices,
            array<long, 0> slice,
            SliceIndex... slice_or_index
            ) //-> decltype(slice_rec<idim+1, ndimslices+1>(size_t, SlicesT<ndimslices>, SliceIndex...))
    {

        //taken as is
        array<long, 2> full_slice {
            shape_[idim], 
            strides_[idim]
        };

        //doesn't increment dimension, only deals with conversion to full slice
        return slice_rec<idim+1, ndimslices+1>(start_ind, slices.cons(full_slice), slice_or_index...); 
    }

    
    /**
     * {start, stop} type of slice where all elements btw bounds are included
     * this function doesnt advance recursion but only processes its arguments
     */
    template <size_t idim, size_t ndimslices, typename... SliceIndex>
    auto
    slice_rec(
            size_t start_ind,
            SlicesT<ndimslices> slices,
            array<long, 2> slice,
            SliceIndex... slice_or_index
            ) //-> decltype(slice_rec<idim, ndimsslices>(size_t, SlicesT<ndimslices>, SliceIndex...))
    {

        array<long, 3> full_slice = {
            slice[0],
            slice[1],
            1
        };

        //doesn't increment dimension, only deals with conversion to full slice
        return slice_rec<idim, ndimslices>(start_ind, slices, slice_or_index...); 
    }

    //full slice {start, stop, step}
    template <size_t idim, size_t ndimslices, typename... SliceIndex>
    auto
    slice_rec(
            size_t start_ind,
            SlicesT<ndimslices> slices,
            array<long, 3> slice,
            SliceIndex... slice_or_index
            ) //-> decltype(slice_rec<idim+1, ndimsslices>(size_t, SlicesT<ndimslices>, SliceIndex...))
    {

        start_ind += strides_[idim]*reverse_negative_index(slice[0]);

        //would be trouble
        assert(slice[2]!=0); 

        long new_stride = strides_[idim]*slice[2];

        array<long, 2> full_slice = {
            abs(reverse_negative_index(slice[1]) - reverse_negative_index(slice[0]))
                /abs(slice[2]),
            new_stride
        };

        //DOES increment dimension
        return slice_rec<idim+1, ndimslices>(start_ind, slices.cons(full_slice), slice_or_index...); 
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
        return (ind>=0)? ind : shape_[idim]+(ind+1)*stride(idim);
    }

    //done handling slice signatures
    //
    public :

    /**
     * Returns a new ndindexer with eventually a smaller number of dimensions,
     * The new ndindexer computes indices matching the requested slice of the array.
     */
    template <typename... DimSliceT>
    auto
    slice(DimSliceT... slice_or_index) 
        //-> decltype(ndindexer(decltype(slice_rec<0, 0>(size_t, SlicesT<0>, DimSliceT...))))
    {
        return ndindexer(slice_rec<0, 0>(start_index_, SlicesT<0>(), slice_or_index...));
    }


};


/**
 *
 * An special ndindexer linked to its data.
 *
 * Can produce iterators in a way similar to STL containers
 *
 * See ndview::iterator
 *
 * Becomes invalidated when an iterator would be
 *
 **/
template<typename ContainerT, typename T, size_t ndims>
struct ndview: ndindexer<ndims> {

    ContainerT v;

    ndview(ndindexer<ndims> ndi, ContainerT ndarray)
        : ndindexer<ndims>(ndi), v(ndarray) { };  

    struct iterator;
};

//nested iterator class
//an iterator is produced by calling ndindexer.slice()
template<typename ContainerT, typename T, size_t ndims>
struct ndview<ContainerT, T, ndims>::iterator {

    //vecarray<array<long, 3>, ndims> slices_

    vecarray<size_t, ndims> ndindexer;

    iterator(ndview<ContainerT, T, ndims> ndv):
        ndview<ContainerT, T, ndims>(ndv)
    { }

    T* begin();

    size_t end();

};



#endif /* end of include guard: NDINDEXER_HPP_9FICI4GD */
