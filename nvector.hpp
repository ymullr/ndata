#ifndef NVECTOR_HPP_B0EOQJGU
#define NVECTOR_HPP_B0EOQJGU

#include <vector>
#include "ndatacontainer.hpp"

namespace ndata {



//forward declaration
template<typename T, long ndims>
auto //nvector<T, somedim>
make_nvector(ndatacontainer<std::vector<T>, T, ndims> idxr);

template<typename T, long ndims>
struct nvector: ndatacontainer<std::vector<T>, T, ndims>  {


    /**
     * @brief construct from an indexer and associated (pointer to) data
     * @param idxr
     * @param data
     */
    template <typename T_rhs>
    nvector(
            indexer<ndims> idxr,
            T_rhs* data
            ):
        nvector(idxr)
    {
        ndataview<T_rhs, ndims> arg_data_view (idxr, data);

        this->assign(
                    arg_data_view
                    );
    }

    /**
     * @brief construct from an indexer and initial value to which the elements of the vector are initialized
     * @param idxr
     * @param data
     */
    nvector(
            indexer<ndims> idxr,
            T initial_value = T()
            ):
        ndatacontainer<std::vector<T>, T, ndims>(idxr.get_shape(), std::vector<T>(idxr.size(), initial_value))
    { }


    /**
     * @brief construct from indexer and data perform elementwise copy of data passed as argument
     *  to the inner vector container only the elements reachable vie the indexer slice are
     *  copied the inner container may be smaller than the copied from data
     * @param idxr
     * @param data
     */
    template <typename ContainerT_rhs, typename T_rhs>
    nvector(
            ndatacontainer<ContainerT_rhs, T_rhs, ndims> ndv
            ):
        nvector(
            ndv,
            &ndv.data_[0]
            )
    { }

    //TODO? resize() function which would reallocate data and remove unreachable elements
    //to be used after some slicing

    ///**
    // * Returns a new ndindexer with eventually a smaller number of dimensions,
    // * The new ndindexer computes indices matching the requested slice of the array.
    // */
    //template <typename... DimSliceT>
    //auto
    //slice(DimSliceT ... slice_or_index) {
    //    //use slice method of the parent class and use it to view this.data_
    //    //and return a new slice
    //    return make_nvector(
    //        ndatacontainer<std::vector<T>, T, ndims>::slice(
    //        slice_or_index...
    //        )
    //    );
    //}

};


//make new nvector from indexer or shape with 0 initialized data
template<typename T, long ndims>
auto //nvector<T, somedim>
make_nvector(indexer<ndims> idxr, T initial_value = T()) {
    nvector<T, ndims> ret (idxr, initial_value);
    return ret;
}

//make new nvector from shape and data
template<typename T, long ndims>
auto //nvector<T, somedim>
make_nvector(indexer<ndims> idxr, std::vector<T> data) {
    nvector<T, ndims> ret (data, idxr);
    return ret;
}


//make new nvector from ndatacontainer
template<typename T, long ndims>
auto //nvector<T, somedim>
make_nvector(ndatacontainer<std::vector<T>, T, ndims> idxr) {
    nvector<T, ndims> ret (idxr);
    return ret;
}

}

#endif /* end of include guard: NVECTOR_HPP_B0EOQJGU */ 
