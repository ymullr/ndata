#ifndef NDATA_HPP_VFUXJBDN
#define NDATA_HPP_VFUXJBDN

#include <vector>
#include <functional>
#include <memory>

namespace ndata {

    struct range;

    template <long ndims>
    struct indexer;

    /**
     * inherits from indexer
     */
    template<typename ContainerT, typename T, long ndims>
    struct ndatacontainer;

    //an ndatacontainer that doesn't own its data
    template <typename T, long ndims>
    using ndataview = ndatacontainer<T*, T, ndims>;

    /**
     * Inherits from ndatacontainer, essentially an ndatacontainer based on a
     * std::vector with some convenience constructors.
     */
    template<typename T, long ndims>
    struct nvector;

}

#include "ndata/forward_declarations.hpp"
#include "ndata/helpers.hpp"
#include "ndata/indexer.hpp"
#include "ndata/ndatacontainer.hpp"
#include "ndata/nvector.hpp"
#include "ndata/loops.hpp"


#endif /* end of include guard: NDATA_HPP_VFUXJBDN */ 
