#ifndef NVECTOR_HPP_B0EOQJGU
#define NVECTOR_HPP_B0EOQJGU

#include <vector>

namespace ndata {


//type alias for uniformity with nvector
template <typename T, size_t ndims>
using ndataview = ndatacontainer<(&std::raw vector<T>), T, ndims>;


//forward declaration
template<typename T, size_t StatSize>
auto //nvector<T, somedim>
make_nvector(ndatacontainer<std::vector<T>, T, StatSize> idxr);

template<typename T, size_t ndims>
struct nvector: public ndatacontainer<std::vector<T>, T, ndims>  {

    //construct from ndatacontainer
    nvector(ndatacontainer<std::vector<T>, T, ndims> ndv):
        ndatacontainer<std::vector<T>, T, ndims>(
            ndv
            )
    { }

    //nvector(vecarray<size_t, ndims> shape):
    //    ndatacontainer<ndims, std::vector<T>, T>(
    //            std::vector<T>(make_indexer(shape).size()),
    //            shape
    //        )
    //{  }

    //construct from indexer or shape
    nvector(indexer<ndims> idxr):
        ndatacontainer<std::vector<T>, T, ndims>(
            std::vector<T>(idxr.size()),
            idxr
            )
    {  }

    //construct from shape, variadic
    template<typename... SizeT>
    nvector(size_t shape0, SizeT... shapeN):
             ndatacontainer<std::vector<T>, T, ndims>(
             std::vector<T>(make_indexer(shape0, shapeN...).size()),
             shape0,
             shapeN...
             )
    {  }


    //construct from data and shape (default strides = 1)
    //nvector(std::vector<T> data, vecarray<size_t, ndims> shape):
    //    ndatacontainer<ndims, std::vector<T>, T>(
    //            data,
    //            shape
    //        )
    //{  assert(data.size() == indexer<ndims>::size()); }

    //construct from data and its indexer shape
    nvector(std::vector<T> data, indexer<ndims> idxr):
        ndatacontainer<std::vector<T>, T, ndims>(
                data,
                idxr
            )
    { assert(data.size() >= idxr.size() ); }

    /**
     * used by broadcast
     */
    //template<size_t new_ndims>
    //ndatacontainer<ContainerT, T, new_ndims>
    //reshape(vecarray<size_t, ndims> new_shape, vecarray<long, ndims> new_strides) {
    //    return ndatacontainer(data_, indexer<new_ndims>(indexer<ndims>::start_index_, new_shape, new_strides));
    //}

    /**
     * Returns a new ndindexer with eventually a smaller number of dimensions,
     * The new ndindexer computes indices matching the requested slice of the array.
     */
    template <typename... DimSliceT>
    auto
    slice(DimSliceT ... slice_or_index) {
        //use slice method of the parent class and use it to view this.data_
        //and return a new slice
        return  make_vector(
                    ndatacontainer<std::vector<T>, T, ndims>::slice(
                    slice_or_index...
                )
                    );
    }

};

//make new nvector from shape
template<typename T, typename... SizeTs>
auto //nvector<T, somedim>
make_nvector(size_t shape0, SizeTs... shapeN) {
    nvector<T, sizeof...(shapeN)+1> ret (shape0, shapeN...);
    return ret;
}

//make new nvector from shape
template<typename T, typename... SizeTs>
auto //nvector<T, somedim>
make_nvector(std::vector<T> data, size_t shape0, SizeTs... shapeN) {
    nvector<T, sizeof...(shapeN)+1> ret (data, shape0, shapeN...);
    return ret;
}

//make new nvector from indexer or shape with 0 initialized data
template<typename T, size_t StatSize>
auto //nvector<T, somedim>
make_nvector(indexer<StatSize> idxr) {
    nvector<T, StatSize> ret (idxr);
    return ret;
}

//make new nvector from shape and data
template<typename T, size_t StatSize>
auto //nvector<T, somedim>
make_nvector(std::vector<T> data, indexer<StatSize> idxr) {
    nvector<T, StatSize> ret (data, idxr);
    return ret;
}


//make new nvector from ndatacontainer
template<typename T, size_t StatSize>
auto //nvector<T, somedim>
make_nvector(ndatacontainer<std::vector<T>, T, StatSize> idxr) {
    nvector<T, StatSize> ret (idxr);
    return ret;
}

}

#endif /* end of include guard: NVECTOR_HPP_B0EOQJGU */ 
