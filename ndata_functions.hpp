#ifndef NDATA_FUNCTIONS_HPP_J3F2DFOA
#define NDATA_FUNCTIONS_HPP_J3F2DFOA


#include <tuple>
//#include "nvector.hpp"
#include "helpers.hpp"
#include <utility>
#include "tuple_utility.hpp"

namespace ndata {

    template <typename Tret, typename FuncT, typename... Ts, size_t... Ndims>
    auto //nvector<some_ndims, Tret>
    ntransform_helper(std::tuple<nvector<Ts, Ndims>...> vN, FuncT func)  {

        static_assert(std::tuple_size<decltype(vN)>::value == sizeof...(Ts), "function arity doesn't match tuple length");

        //broadcast arguments against each others
        auto tuple_params_bc = helpers::broadcast(vN);
        //tuple_utility::apply<decltype(helpers::broadcast<decltype(vN)>), decltype(vN)>(
        //    helpers::broadcast<decltype(vN)>,
        //    vN
        //);

        //init return nvector with correct shape

        auto retshape = std::get<0>(tuple_params_bc).get_shape();

        nvector<Tret, retshape.static_size_or_dynamic> ret (retshape);

        constexpr size_t tup_size = std::tuple_size<decltype(tuple_params_bc)>::value;

        auto ndindex = retshape;
        ndindex.fill(0);
        for (int i = 0; i < ret.size(); ++i) {
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
}

#endif /* end of include guard: NDATA_FUNCTIONS_HPP_J3F2DFOA */
