#ifndef NVECTOR_HPP_B0EOQJGU
#define NVECTOR_HPP_B0EOQJGU

#include <vector>

namespace ndata {

template<size_t ndims, typename T>
struct nvector: ndview<ndims, std::vector<T>, T>  {

    //construct from shape with default initialized vector
    nvector(vecarray<size_t, ndims> shape):
        ndview<ndims, std::vector<T>, T>(
                std::vector<T>(make_indexer(shape).size()),
                shape
            )
    {  }

    //nvector(indexer<ndims> shape):
    //    ndview<ndims, std::vector, T>(
    //            std::vector<T>(make_indexer(shape).size()),
    //            shape
    //        )
    //{  }

    //construct from data and shape (default strides = 1)
    nvector(std::vector<T> data, vecarray<size_t, ndims> shape):
        ndview<ndims, std::vector<T>, T>(
                data,
                shape
            )
    {  assert(data.size() == indexer<ndims>::size()); }

    //construct from data and its indexer
    nvector(std::vector<T> data, indexer<ndims> idxr):
        ndview<ndims, std::vector<T>, T>(
                data,
                idxr
            )
    {  assert(data.size() >= idxr.size() );}

    //construct from ndview
    nvector(ndview<ndims, std::vector<T>, T> ndv):
        ndview<ndims, std::vector<T>, T>(
            ndv
            )
    { }
};

}

#endif /* end of include guard: NVECTOR_HPP_B0EOQJGU */ 
