#ifndef NVECTOR_HPP_B0EOQJGU
#define NVECTOR_HPP_B0EOQJGU

#include <vector>

namespace ndata {

template<size_t ndims, typename T>
struct nvector: ndview<ndims, std::vector<T>, T>  {

    //construct from ndview
    nvector(ndview<ndims, std::vector<T>, T> ndv):
        ndview<ndims, std::vector<T>, T>(
            ndv
            )
    { }

    //nvector(vecarray<size_t, ndims> shape):
    //    ndview<ndims, std::vector<T>, T>(
    //            std::vector<T>(make_indexer(shape).size()),
    //            shape
    //        )
    //{  }

    nvector(indexer<ndims> idxr):
        ndview<ndims, std::vector<T>, T>(
            std::vector<T>(idxr.size()),
            idxr
            )
    {  }

    //template<typename... SizeT>
    //nvector(size_t shape0, SizeT... shapeN):
    //         ndview<ndims, std::vector<T>, T>(
    //         std::vector<T>(make_indexer(shape0, shapeN...).size()),
    //         shape0,
    //         shapeN...
    //         )
    //{  }


    //construct from data and shape (default strides = 1)
    //nvector(std::vector<T> data, vecarray<size_t, ndims> shape):
    //    ndview<ndims, std::vector<T>, T>(
    //            data,
    //            shape
    //        )
    //{  assert(data.size() == indexer<ndims>::size()); }

    //construct from data and its indexer
    nvector(std::vector<T> data, indexer<ndims> idxr):
        ndview<ndims, std::vector<T>, T>(
                data,
                idxr
            )
    {  assert(data.size() >= idxr.size() );}

};

}

#endif /* end of include guard: NVECTOR_HPP_B0EOQJGU */ 
