#ifndef NVECTOR_HPP_B0EOQJGU
#define NVECTOR_HPP_B0EOQJGU

#include <vector>
#include "ndatacontainer.hpp"

namespace ndata {



//forward declaration
template<typename T, long ndims>
auto //nvector<T, somedim>
make_nvector(ndatacontainer<std::vector<T>, T, ndims> idxr);

template<typename T, long ndims>
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
            idxr,
            std::vector<T>(idxr.size())
            )
    {  }

    //construct from data and its indexer shape
    nvector(indexer<ndims> idxr, std::vector<T> data):
        ndatacontainer<std::vector<T>, T, ndims>(
                idxr,
                data
            )
    { assert(data.size() >= idxr.size() ); }

    //TODO? resize() function which would reallocate data and remove unreachable elements
    //to be used after some slicing

    ///**
    // * Returns a new ndindexer with eventually a smaller number of dimensions,
    // * The new ndindexer computes indices matching the requested slice of the array.
    // */
    //template <typename... DimSliceT>
    //auto
    //slice(DimSliceT ... slice_or_index) {
    //    //use slice method of the parent class and use it to view this.data_
    //    //and return a new slice
    //    return make_nvector(
    //        ndatacontainer<std::vector<T>, T, ndims>::slice(
    //        slice_or_index...
    //        )
    //    );
    //}

};


//make new nvector from indexer or shape with 0 initialized data
template<typename T, long ndims>
auto //nvector<T, somedim>
make_nvector(indexer<ndims> idxr) {
    nvector<T, ndims> ret (idxr);
    return ret;
}

//make new nvector from shape and data
template<typename T, long ndims>
auto //nvector<T, somedim>
make_nvector(indexer<ndims> idxr, std::vector<T> data) {
    nvector<T, ndims> ret (data, idxr);
    return ret;
}


//make new nvector from ndatacontainer
template<typename T, long ndims>
auto //nvector<T, somedim>
make_nvector(ndatacontainer<std::vector<T>, T, ndims> idxr) {
    nvector<T, ndims> ret (idxr);
    return ret;
}

}

#endif /* end of include guard: NVECTOR_HPP_B0EOQJGU */ 
