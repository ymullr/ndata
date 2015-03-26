#ifndef NDVIEW_HPP_1LWJBSCE
#define NDVIEW_HPP_1LWJBSCE

#include "vecarray.hpp"
#include "ndindexer.hpp"

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




#endif /* end of include guard: NDVIEW_HPP_1LWJBSCE */

