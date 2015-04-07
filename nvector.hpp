#ifndef NVECTOR_HPP_B0EOQJGU
#define NVECTOR_HPP_B0EOQJGU

#include <vector>
#include "ndatacontainer.hpp"

namespace ndata {



//forward declaration
template<typename T, size_t StatSize>
auto //nvector<T, somedim>
make_nvector(ndatacontainer<std::vector<T>, T, StatSize> idxr);

template<typename T, size_t ndims>
struct nvector: ndatacontainer<std::vector<T>, T, ndims>  {

    //construct from ndatacontainer
    nvector(ndatacontainer<std::vector<T>, T, ndims> ndv):
        ndatacontainer<std::vector<T>, T, ndims>(
            ndv
            )
    { }

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

    //construct from data and its indexer shape
    nvector(std::vector<T> data, indexer<ndims> idxr):
        ndatacontainer<std::vector<T>, T, ndims>(
                data,
                idxr
            )
    { assert(data.size() >= idxr.size() ); }

    //TODO resize() function which would reallocate data and remove unreachable elements
    //to be used after some slicing

    /**
     * Returns a new ndindexer with eventually a smaller number of dimensions,
     * The new ndindexer computes indices matching the requested slice of the array.
     */
    template <typename... DimSliceT>
    auto
    slice(DimSliceT ... slice_or_index) {
        //use slice method of the parent class and use it to view this.data_
        //and return a new slice
        return  make_nvector(
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
