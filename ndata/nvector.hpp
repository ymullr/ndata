/*! \file */
#ifndef NVECTOR_HPP_B0EOQJGU
#define NVECTOR_HPP_B0EOQJGU

#include <vector>
#include "ndata/ndatacontainer.hpp"

namespace ndata {

enum UNINITIALIZED { };

//constexpr enum ValueInitialization UNINITIALIZED = ValueInitialization::UNINITIALIZED;

//forward declaration
template<typename T, long ndims>
auto //nvector<T, somedim>
make_nvector(ndatacontainer<std::vector<T>, T, ndims> idxr);

/**
 * @brief A n-dimensional container using std::vector from the STL.
 *
 * In addition to the methods inherited from ndatacontainer<std::vector<T>, T, ndims>, it has the property
 * of re-adjusting the size of the underlying std::vector when copy constructed.
 */
template<typename T, long ndims>
struct nvector: ndatacontainer<std::vector<T>, T, ndims>  {

    /**
     * @brief construct from an indexer and initial value to which the elements of the vector are initialized
     * @param idxr
     * @param initial_value must be specified to avoid ambiguity with nvector(ndatacontainer<...> ).
     *  See nvector(indexer<...>, ValueInitialization) to leave internal data uninitialized
     *  (actually makes no difference with std::vector as ContainerT)
     */
    nvector(
            indexer<ndims> idxr,
            T initial_value = T()
            ):
        ndatacontainer<std::vector<T>, T, ndims>(idxr.get_shape(), std::vector<T>(idxr.size(), initial_value))
    { }

    /**
     * @brief construct from an indexer and initial value to which the elements of the vector are initialized
     * @param idxr
     * @param data
     */
    nvector(
            indexer<ndims> idxr,
            UNINITIALIZED
            ):
        ndatacontainer<std::vector<T>, T, ndims>(idxr.get_shape(), std::vector<T>(idxr.size()))
    { }


    /**
     * @brief Construct from indexer and data as a vector
     * @param idxr
     * @param data
     */
    nvector(indexer<ndims> idxr, std::vector<T> data):
        nvector(idxr, &data[0])
    { }

    /**
     * @brief construct from an indexer and associated (pointer to) data.
     * @param idxr
     * @param data
     */
    template <typename T_rhs>
    nvector(
            indexer<ndims> idxr,
            T_rhs* data
            ):
        nvector(idxr)
        //nvector(idxr, 0)
    {
        ndataview<T_rhs, ndims> arg_data_view (idxr, data);

        this->assign(
                    arg_data_view
                    );
    }

    /**
     * @brief construct from an ndatacontainer, copying its data.
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

    //empty constructor for later assignment
    nvector(){};

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
    nvector<T, ndims> ret (idxr, data);
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
