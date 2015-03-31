#ifndef NVECTOR_HPP_B0EOQJGU
#define NVECTOR_HPP_B0EOQJGU

#include <vector>

namespace ndata {

template<typename T, size_t ndims>
struct nvector: public ndview<std::vector<T>, T, ndims>  {

    //construct from ndview
    nvector(ndview<std::vector<T>, T, ndims> ndv):
        ndview<std::vector<T>, T, ndims>(
            ndv
            )
    { }

    //nvector(vecarray<size_t, ndims> shape):
    //    ndview<ndims, std::vector<T>, T>(
    //            std::vector<T>(make_indexer(shape).size()),
    //            shape
    //        )
    //{  }

    //construct from indexer or shape
    nvector(indexer<ndims> idxr):
        ndview<std::vector<T>, T, ndims>(
            std::vector<T>(idxr.size()),
            idxr
            )
    {  }

    //construct from shape, variadic
    template<typename... SizeT>
    nvector(size_t shape0, SizeT... shapeN):
             ndview<std::vector<T>, T, ndims>(
             std::vector<T>(make_indexer(shape0, shapeN...).size()),
             shape0,
             shapeN...
             )
    {  }


    //construct from data and shape (default strides = 1)
    //nvector(std::vector<T> data, vecarray<size_t, ndims> shape):
    //    ndview<ndims, std::vector<T>, T>(
    //            data,
    //            shape
    //        )
    //{  assert(data.size() == indexer<ndims>::size()); }

    //construct from data and its indexer shape
    nvector(std::vector<T> data, indexer<ndims> idxr):
        ndview<std::vector<T>, T, ndims>(
                data,
                idxr
            )
    {  assert(data.size() >= idxr.size() );}

    /**
     * used by broadcast
     */
    //template<size_t new_ndims>
    //ndview<ContainerT, T, new_ndims>
    //reshape(vecarray<size_t, ndims> new_shape, vecarray<long, ndims> new_strides) {
    //    return ndview(data_, indexer<new_ndims>(indexer<ndims>::start_index_, new_shape, new_strides));
    //}

};

//make new nvector from shape
template<typename T, typename... SizeTs>
auto //nvector<T, somedim>
make_nvector(SizeTs... args) {
    nvector<T, sizeof...(args)> ret (args...);
    return ret;
}


}

#endif /* end of include guard: NVECTOR_HPP_B0EOQJGU */ 
