#ifndef NDATA_FUNCTIONS_HPP_J3F2DFOA
#define NDATA_FUNCTIONS_HPP_J3F2DFOA


#include <tuple>
//#include "nvector.hpp"
#include "nhelpers.hpp"
#include <utility>
#include "tuple_utility.hpp"

namespace ndata {

    template <typename Tret, typename FuncT, typename... Ndatacontainer>
    void
    nforeach(std::tuple<Ndatacontainer&...> ndata_tup_refs, FuncT func)  {

        auto tuple_params_bcv = helpers::broadcast_views(ndata_tup_refs);

        //broadcast arguments against each others
        //all broadcasted indexers should have the same shape
        //let's get the first
        auto& idxr = std::get<0>(tuple_params_bcv);

        //multidimensional index common to all indexers
        auto ndindex = idxr.ndindex(0);

        for (size_t i = 0; i < idxr.size(); ++i) {
            //transform tuple of Ndatacontainer to tuple of refs to scalar values for current ndindex
            auto tuple_params_scalar = tuple_utility::tuple_transform([ndindex] (auto & A) {
                return &A(ndindex);
            }, tuple_params_bcv);

            tuple_utility::apply(func, tuple_params_scalar);

            idxr.increment_ndindex(ndindex);
        }
    }


    template <typename Tret, typename FuncT, typename... Ndatacontainer>
    auto //nvector<Tret, some_ndims>
    ntransform(std::tuple<Ndatacontainer...> ndata_tup, FuncT func)  {

        auto ndata_tuple_bcviews = helpers::broadcast_views(ndata_tup);


        //broadcast arguments against each others
        //auto ndata_tuple_bc = helpers::broadcast(ndata_tuple_views);

        //init return nvector with correct shape
        auto retshape = std::get<0>(ndata_tuple_bcviews).get_shape();
        nvector<Tret, retshape.STATIC_SIZE_OR_DYNAMIC> ret (std::get<0>(ndata_tuple_bcviews), 0);

        //our multidimensional index
        auto ndindex = ret.ndindex(0);

        for (size_t i = 0; i < ret.size(); ++i) {
            //spread the tuple to pass the i-indexed as scalar argument to func
            //ret[i] = tuple_utility::invoke_helper(func, ndata_tuple_bcviews, i_tup(), i); //func(get<i_tup>(ndata_tuple_bcviews)[i]...);
            //std::tuple<typename Ndatacontainer::type_T...>
            auto
                tuple_params_scalar = tuple_utility::tuple_transform(
                        [ndindex] (auto A) {
                            return A(ndindex);
                            },
                        ndata_tuple_bcviews);
            ret(ndindex) = tuple_utility::apply(func, tuple_params_scalar); //func(get<i_tup>(ndata_tuple_bcviews)[i]...);
            ret.increment_ndindex(ndindex);
        }

        return ret;
    }

    //alternative overload which avoids the creation of a tuple but seems less readable
    //template <typename Tret, typename FuncT, typename... Ndatacontainer>
    //auto //nvector<Tret, some_ndims>
    //ntransform(FuncT func, Ndatacontainer... ndata_tup)  {
    //    return ntransform<Tret, FuncT>(make_tuple(ndata_tup...), func);
    //}


}

#endif /* end of include guard: NDATA_FUNCTIONS_HPP_J3F2DFOA */
