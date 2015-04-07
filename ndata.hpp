#ifndef NDATA_HPP_VFUXJBDN
#define NDATA_HPP_VFUXJBDN

#include <vector>
#include <functional>
#include "nhelpers.hpp"

namespace ndata {

    struct Rng;

    template <size_t ndims>
    struct indexer;

    /**
     * inherits from indexer
     */
    template<typename ContainerT, typename T, size_t ndims>
    struct ndatacontainer;

    //an ndatacontainer that doesn't own its data
    template <typename T, size_t ndims>
    using ndataview = ndatacontainer<T*, T, ndims>;

    /**
     * Inherits from ndatacontainer, essentially an ndatacontainer based on a
     * std::vector with some convenience constructors.
     */
    template<typename T, size_t ndims>
    struct nvector;

}

#include "nindexer.hpp"
#include "ndatacontainer.hpp"
#include "nvector.hpp"
#include "ndata_functions.hpp"


#endif /* end of include guard: NDATA_HPP_VFUXJBDN */ 
