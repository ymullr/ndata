#include <vector>

#include <initializer_list>
#include "vecarray.hpp"

#ifndef NDINDEXER_HPP_9FICI4GD
#define NDINDEXER_HPP_9FICI4GD

template<size_t ndims>
struct ndindexer {

    template<typename... SizeT>
    ndindexer(SizeT... shape)
    {
        static_assert(
            ndims != 0,
            "This overload is only available when the number of dimensions is a compile time constant"
            );
        static_assert(ndims == sizeof...(shape),
            "Number of parameters doesn't match the number of dimensions."
        );
        shape_=array_from_variadic_argument<ndims>(shape...);
        strides_ = calc_strides_(shape_);
    }

    ndindexer(vecarray<size_t, ndims> shape):
        shape_(shape),
        strides_(calc_strides_(shape)) { }


    ndindexer(vector<size_t> shape):
        shape_(shape)
    { 
        strides_=calc_strides_(shape);
        //static_assert(ndims == 0, 
        //   "This constructor always produce a dynamic vecarray,"
        //   " please set ndims template argument to 0 to use this constructor."
        //);
        assert(shape.size() == ndims or ndims == 0);
    }

    ndindexer(array<size_t, ndims> shape):
        shape_(shape),
        strides_(calc_strides_(shape)) { }

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


        size_t indexacc=0;

        //looping on array dimensions
        for (size_t i = 0; i < ndindex.size(); ++i) {
            //bound checking when debugging
            assert(ndindex[i] < shape_[i]);

            indexacc+=ndindex[i]*strides_[i];
        }

        return indexacc;
    }

    template<typename... SizeT>
    size_t index(SizeT... ndindex){
        static_assert(
            ndims != 0,
            "This overload is only available when the number of dimensions is a compile time constant"
            );
        static_assert(
            ndims == sizeof...(ndindex),
            "Number of parameters doesn't match the number of dimensions."
        );
        array<size_t, ndims> ind = array_from_variadic_argument<ndims>(ndindex...);

        return index(vecarray<size_t, ndims>(ind));
    }

    /**
     * Return the stride for the current dimension
     */
    vecarray<size_t, ndims> get_strides() {
        return strides_;
    }

    /**
     * Just returns the shape used to construct the class (if any).
     */
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
     *
     *
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
     *
     */
    vecarray<size_t, ndims> ndindex(size_t flatindex) {

        vecarray<size_t, ndims> ndindex (shape_);

        ndindex.fill(0);

        for (size_t i = 0; i <= flatindex; ++i) {
            increment_ndindex(ndindex);
        }

        return ndindex;
    }

    private:

    vecarray<size_t, ndims> shape_;

    vecarray<size_t, ndims> strides_;

    static
    vecarray<size_t, ndims> 
    calc_strides_(vecarray<size_t, ndims> shape) {

        vecarray<size_t, ndims> ret (shape);

        for (size_t i = 0; i < shape.size(); ++i) {

            size_t this_stride=1;
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
    template<size_t ndims_variadic, typename... SizeT>
    array<size_t, ndims_variadic> array_from_variadic_argument(size_t i, SizeT... rest) {
        static_assert(
            sizeof...(rest) == ndims_variadic-1,
            ""
        );
        array<size_t, ndims_variadic> ret;
        ret[0] = i;

        array<size_t, ndims_variadic-1> rest_arr = array_from_variadic_argument<ndims_variadic-1>(rest...);

        for (size_t irest = 0; irest < rest_arr.size(); ++irest) {
            ret[irest+1] = rest_arr[irest];
        }

        return ret;
    }

    template<size_t ndims_variadic>
    array<size_t, 0> array_from_variadic_argument() {
        //static_assert(ndims_variadic == 0, 
        //        );
        array<size_t, 0> ret;
        return ret;
    }

};

#endif /* end of include guard: NDINDEXER_HPP_9FICI4GD */
