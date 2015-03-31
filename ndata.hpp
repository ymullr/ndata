#ifndef NDATA_HPP_VFUXJBDN
#define NDATA_HPP_VFUXJBDN

#include <vector>
#include <functional>


namespace ndata {

    struct Rng;

    template <size_t ndims>
    struct indexer;

    template<size_t ndims, typename ContainerT, typename T>
    struct ndview;

    template<size_t ndims, typename T>
    struct nvector;

    template <typename... NVectorTs, typename... Ts, typename Tret>
    ntransform(NVectorTs... vN, std::function<Tret(Ts...)> func)  {
        tuple<nvector<Ts...>> vec_bc = helpers::broadcast(v1, vN);
        nvector
        for (int i = 0; i < vec_bc[0].size(); ++i) {

        }

    }
}


#include "indexer.hpp"
#include "ndview.hpp"
#include "nvector.hpp"
//#include "ndview_impl.hpp"


#endif /* end of include guard: NDATA_HPP_VFUXJBDN */ 
