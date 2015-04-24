#ifndef NDATA_FUNCTIONS_HPP_J3F2DFOA
#define NDATA_FUNCTIONS_HPP_J3F2DFOA


#include <tuple>
//#include "nvector.hpp"
#include "nhelpers.hpp"
#include <utility>
#include "tuple_utility.hpp"

namespace ndata {

    template <typename Tret, typename FuncT, typename... Ndatacontainers>
    void
    nforeach(std::tuple<Ndatacontainers*...> ndata_tup_ptrs, FuncT func)  {

        //ensure all elements inside ndata_tup_ptrs are converted to views
        //necessing for returning pointers to elements in ndata later
        //transform all pointers to ndatacontainers to ndviews
        auto ndata_tup = tuple_utility::tuple_transform_ptr([] (auto * elt) {
            return elt->to_view();
        }, ndata_tup_ptrs);

        auto tuple_params_bc = helpers::broadcast(ndata_tup);

        //broadcast arguments against each others
        //all broadcasted indexers should have the same shape
        //let's get the first
        auto idxr = std::get<0>(tuple_params_bc);

        //multidimensional index common to all indexers
        auto ndindex = idxr.get_shape();
        ndindex.fill(0);

        for (size_t i = 0; i < idxr.size(); ++i) {
            //transform tuple of ndatacontainers to tuple of refs to scalar values for current ndindex
            auto tuple_params_scalar = tuple_utility::tuple_transform_ptr([ndindex] (auto A) {
                return &A(ndindex);
            }, tuple_params_bc);

            tuple_utility::apply(func, tuple_params_scalar);

            idxr.increment_ndindex(ndindex);
        }
    }

    //template <typename Tret, typename FuncT, typename... Ndatacontainers>
    //auto //nvector<Tret, some_ndims>
    //nforeach(FuncT func, Ndatacontainers... ndata_tup)  {
    //    nforeach(make_tuple(ndata_tup...), func);
    //}

    //disabled bc FuncT can't be properly inferred with this signature which leads to awkward
    //call syntax
    //template <typename Tret, typename FuncT, typename... Ndatacontainers>
    //auto //nvector<Tret, some_ndims>
    //ntransform(Ndatacontainers... ndata_tup, std::function<FuncT> func)  {
    //    return ntransform<Tret, FuncT>(make_tuple(ndata_tup...), func);
    //}

    template <typename Tret, typename FuncT, typename... Ndatacontainers>
    auto //nvector<Tret, some_ndims>
    ntransform(std::tuple<Ndatacontainers...> ndata_tup, FuncT func)  {

        //convert to views before broadcasting (broadcasting needs to be able to copy,
        //ndatacontainer<unique_pointer, ...> is not copyable but ndataviews are)
        //std::tuple<
        //            ndataview<
        //                typename Ndatacontainers::type_T,
        //                //Ndatacontainers::STATIC_SIZE_OR_DYNAMIC
        //                helpers::static_max_or_dynamic<Ndatacontainers::STATIC_SIZE_OR_DYNAMIC...>::value
        //            > ...
        //        >
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
            //std::tuple<typename Ndatacontainers::type_T...>
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
    //template <typename Tret, typename FuncT, typename... Ndatacontainers>
    //auto //nvector<Tret, some_ndims>
    //ntransform(FuncT func, Ndatacontainers... ndata_tup)  {
    //    return ntransform<Tret, FuncT>(make_tuple(ndata_tup...), func);
    //}


}

#endif /* end of include guard: NDATA_FUNCTIONS_HPP_J3F2DFOA */
