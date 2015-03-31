#ifndef NDATA_HPP_VFUXJBDN
#define NDATA_HPP_VFUXJBDN

#include <vector>
#include <functional>
#include "helpers.hpp"


namespace ndata {

    struct Rng;

    template <size_t ndims>
    struct indexer;

    template<typename ContainerT, typename T, size_t ndims>
    struct ndatacontainer;

    template<typename T, size_t ndims>
    struct nvector;

}

#include "indexer.hpp"
#include "ndatacontainer.hpp"
#include "nvector.hpp"
#include "ndata_functions.hpp"


#endif /* end of include guard: NDATA_HPP_VFUXJBDN */ 
