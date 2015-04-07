#ifndef NDATA_FUNCTIONS_HPP_J3F2DFOA
#define NDATA_FUNCTIONS_HPP_J3F2DFOA


#include <tuple>
//#include "nvector.hpp"
#include "helpers.hpp"
#include <utility>
#include "tuple_utility.hpp"

namespace ndata {

    template <typename Tret, typename FuncT, typename... Ndatacontainers>
    void
    nforeach(std::tuple<Ndatacontainers...> vN_containers, FuncT func)  {

        //ensure all elements inside vN_containers are converted to views
        //necessing for returning pointers to elements in ndata later
        auto vN = tuple_utility::tuple_transform_ptr([] (auto & elt) {
            return elt.to_ndataview();
        }, vN_containers);

        auto tuple_params_bc = helpers::broadcast(vN);

        //broadcast arguments against each others
        //all broadcasted indexers should have the same shape
        //let's get the first
        auto idxr = std::get<0>(tuple_params_bc);

        //multidimensional index common to all indexers
        auto ndindex = idxr.get_shape();
        ndindex.fill(0);

        for (int i = 0; i < idxr.size(); ++i) {
            //transform tuple of ndatacontainers to tuple of refs to scalar values for current ndindex
            auto tuple_params_scalar = tuple_utility::tuple_transform([ndindex] (auto A) {
                return &A.val(ndindex);
            }, tuple_params_bc);

            tuple_utility::apply(func, tuple_params_scalar);

            idxr.increment_ndindex(ndindex);
        }
    }

    //template <typename Tret, typename FuncT, typename... Ndatacontainers>
    //auto //nvector<Tret, some_ndims>
    //nforeach(FuncT func, Ndatacontainers... vN)  {
    //    nforeach(make_tuple(vN...), func);
    //}

    //disabled bc FuncT can't be properly inferred with this signature which leads to awkward
    //call syntax
    //template <typename Tret, typename FuncT, typename... Ndatacontainers>
    //auto //nvector<Tret, some_ndims>
    //ntransform(Ndatacontainers... vN, std::function<FuncT> func)  {
    //    return ntransform<Tret, FuncT>(make_tuple(vN...), func);
    //}

    template <typename Tret, typename FuncT, typename... Ndatacontainers>
    auto //nvector<Tret, some_ndims>
    ntransform(std::tuple<Ndatacontainers...> vN, FuncT func)  {

        //broadcast arguments against each others
        auto tuple_params_bc = helpers::broadcast(vN);

        //init return nvector with correct shape
        auto retshape = std::get<0>(tuple_params_bc).get_shape();
        nvector<Tret, retshape.static_size_or_dynamic> ret (retshape);

        //our multidimensional index
        auto ndindex = retshape;
        ndindex.fill(0);
        for (size_t i = 0; i < ret.size(); ++i) {
            //spread the tuple to pass the i-indexed as scalar argument to func
            //ret[i] = tuple_utility::invoke_helper(func, tuple_params_bc, i_tup(), i); //func(get<i_tup>(tuple_params_bc)[i]...);
            auto tuple_params_scalar = tuple_utility::tuple_transform([ndindex] (auto A) {
                return A.val(ndindex);
            }, tuple_params_bc);
            ret.val(ndindex) = tuple_utility::apply(func, tuple_params_scalar); //func(get<i_tup>(tuple_params_bc)[i]...);
            ret.increment_ndindex(ndindex);
        }

        return ret;
    }

    //alternative overload which avoids the creation of a tuple but seems less readable
    //template <typename Tret, typename FuncT, typename... Ndatacontainers>
    //auto //nvector<Tret, some_ndims>
    //ntransform(FuncT func, Ndatacontainers... vN)  {
    //    return ntransform<Tret, FuncT>(make_tuple(vN...), func);
    //}


}

#endif /* end of include guard: NDATA_FUNCTIONS_HPP_J3F2DFOA */
