#include <vector>

#include <initializer_list>
#include "vecarray.hpp"

#ifndef NDINDEXER_HPP_9FICI4GD
#define NDINDEXER_HPP_9FICI4GD


template<size_t ndims>
struct ndindexer {

    vecarray<size_t, ndims> shape;

    vecarray<size_t, ndims> strides;

    static const size_t ndims_ = ndims;

    ndindexer(initializer_list<size_t> shape_):
        shape(shape_),
        strides(calc_strides())
    {
        //array<size_t> sh = 
        assert(shape.size() == shape_.size());
    }


    ndindexer(vector<size_t> shape):
        shape(shape),
        strides(calc_strides)
    { 
        static_assert(ndims == 0, "This constructor always produce a dynamic vecarray, please set ndims template argument to 0.");
    }

    ndindexer(array<size_t, ndims> shape):
        shape(shape),
        strides(calc_strides()) { }

    ndindexer(vecarray<size_t, ndims> shape):
        shape(shape),
        strides(calc_strides()) { }

    size_t i(initializer_list<size_t> ndindex){
        assert(ndindex.size() == shape.size());

        size_t i=0;
        size_t indexacc=0;

        //looping on array dimensions
        for(size_t dim_index: ndindex){

                //bound checking when NOT compiled with -DNDEBUG 
                #ifndef NDEBUG
                assert(dim_index < shape[i]);
                #endif

                indexacc+=dim_index*strides[i];

                i++;
        }

        return indexacc;
    }

    size_t stride(size_t dim){
        size_t stride=1;
        for(size_t j=dim+1;j<shape.size();j++){
            stride*=shape[j];
        }
        return stride;
    }

    size_t size(){
        size_t acc=1;
        for(size_t i=0;i<shape.size();i++){
            acc*=shape[i];
        }
        return acc;
    }

    private:
    vecarray<size_t, ndims> calc_strides() {

        vecarray<size_t, ndims> strides_ (shape);

        for (size_t i = 0; i < shape.size(); ++i) {
            strides_[i] = stride(i);
        }

        return strides_;
    }

};


#endif /* end of include guard: NDINDEXER_HPP_9FICI4GD */
