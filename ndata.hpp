#ifndef NDATA_HPP_VFUXJBDN
#define NDATA_HPP_VFUXJBDN

#include <vector>
#include <functional>
#include <memory>
#include "nhelpers.hpp"

namespace ndata {

    struct Rng;

    template <long ndims>
    struct indexer;

    /**
     * inherits from indexer
     */
    template<typename ContainerT, typename T, long ndims>
    struct ndatacontainer;

    //an ndatacontainer that doesn't own its data
    //template <typename T, long ndims>
    //using ndataview = ndatacontainer<T*, T, ndims>;
    template <typename T, long ndims>
    using ndataview = ndatacontainer<std::unique_ptr<T>, T, ndims>;

    /**
     * Inherits from ndatacontainer, essentially an ndatacontainer based on a
     * std::vector with some convenience constructors.
     */
    template<typename T, long ndims>
    struct nvector;

}

#include "nindexer.hpp"
#include "ndatacontainer.hpp"
#include "nvector.hpp"
#include "ndata_functions.hpp"


#endif /* end of include guard: NDATA_HPP_VFUXJBDN */ 
