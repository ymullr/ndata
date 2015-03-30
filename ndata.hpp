#ifndef NDATA_HPP_VFUXJBDN
#define NDATA_HPP_VFUXJBDN

#include <vector>


namespace ndata {

    struct Rng;

    template <size_t ndims>
    struct indexer;

    template<size_t ndims, typename ContainerT, typename T>
    struct ndview;

    template<size_t ndims, typename T>
    struct nvector;
}


#include "indexer.hpp"
#include "ndview.hpp"
#include "nvector.hpp"
//#include "ndview_impl.hpp"


#endif /* end of include guard: NDATA_HPP_VFUXJBDN */ 
