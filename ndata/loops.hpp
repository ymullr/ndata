#ifndef NDATA_FUNCTIONS_HPP_J3F2DFOA
#define NDATA_FUNCTIONS_HPP_J3F2DFOA


#include <tuple>
#include "ndata/helpers.hpp"
#include <utility>
#include "tuple_utility.hpp"

#ifdef _OPENMP
   #include <omp.h>
   #define NDATA_OMP_GET_NUM_THREADS() omp_get_num_threads()
#else
   #define NDATA_OMP_GET_NUM_THREADS() 1
#endif

namespace ndata {

    constexpr int
        SERIAL=0,
        PARALLEL=1;


    template <
            int loop_type = SERIAL,
            long ndims,
            typename ... Ts,
            typename FuncT
            >
    void
    nforeach_views(std::tuple<ndataview<Ts, ndims>...> ndata_views, FuncT func)  {

        //broadcast arguments against each others
        //all broadcasted indexers should have the same shape
        //let's get the first
        auto& idxr = std::get<0>(ndata_views);

        if (loop_type == SERIAL or NDATA_OMP_GET_NUM_THREADS() == 1) {

            auto ndindex = idxr.ndindex(0);

            for (size_t i = 0; i < idxr.size(); ++i) {
                //TODO auto infer tuple<T&...> and get rid of pointers in apply signature (aliasing?)
                //transform tuple of Ndatacontainer to tuple of refs to scalar values for current ndindex
                auto tuple_params_scalar = tuple_utility::tuple_transform([ndindex] (auto & A) {
                    return &A(ndindex);
                }, ndata_views);

                tuple_utility::apply(func, tuple_params_scalar);

                idxr.increment_ndindex(ndindex);
            }

        } else {

#pragma omp parallel for schedule(static)
            for (size_t i_start=0; i_start<idxr.size(); i_start+= idxr.size()/NDATA_OMP_GET_NUM_THREADS()) {

                //our multidimensional index
                auto ndindex = idxr.ndindex(i_start);

                for (size_t i = 0; i < idxr.size()/NDATA_OMP_GET_NUM_THREADS(); ++i) {
                        //TODO
                        //transform tuple of Ndatacontainer to tuple of refs to scalar values for current ndindex
                        //infer ref types tuple<T&...> and get rid of raw pointers in apply signature (also : what of issues with pointer aliasing?)
                        //here auto is std::tuple<T*...>
                        //problem is T& decay to T with auto
                        auto tuple_params_scalar = tuple_utility::tuple_transform([ndindex] (auto & A) {
                            return &A(ndindex);
                        }, ndata_views);

                        tuple_utility::apply(func, tuple_params_scalar);

                        idxr.increment_ndindex(ndindex);
                 }
            }
        }

    }

    template <int loop_type = SERIAL , typename FuncT, typename... Ndatacontainer>
    void
    nforeach(std::tuple<Ndatacontainer&...> ndata_tup_refs, FuncT func)  {
        auto ndata_views = helpers::broadcast_views(ndata_tup_refs);
        nforeach_views<loop_type>(ndata_views, func);
    }

    template <typename FuncT, typename... Ndatacontainer>
    void
    nforeach_parallel(std::tuple<Ndatacontainer&...> ndata_tup_refs, FuncT func)  {
        auto ndata_views = helpers::broadcast_views(ndata_tup_refs);
        nforeach_views<PARALLEL>(ndata_views, func);
    }


    //TODO make this able to infer Tret
    template <typename Tret, int loop_type = SERIAL, typename FuncT, typename... Ndatacontainer>
    auto //nvector<Tret, ndims_broadcasted>
    ntransform(std::tuple<Ndatacontainer...> ndata_tup, FuncT func)  {

        auto ndata_tuple_bcviews = helpers::broadcast_views(ndata_tup);

        //broadcast arguments against each others
        //auto ndata_tuple_bc = helpers::broadcast(ndata_tuple_views);

        //init return nvector with correct shape
        auto retshape = std::get<0>(ndata_tuple_bcviews).get_shape();
        nvector<Tret, retshape.STATIC_SIZE_OR_DYNAMIC> ret (std::get<0>(ndata_tuple_bcviews), 0);

        nforeach_views<loop_type>(
                    std::tuple_cat(
                        std::make_tuple(ret.to_view()),
                        ndata_tuple_bcviews
                        )
                    ,
                    [func] (Tret & ret_val, auto ... param_vals) {
                        ret_val=func(param_vals...);
                    }
            );

        return ret;
    }


    template <typename Tret, typename FuncT, typename... Ndatacontainer>
    auto //nvector<Tret, ndims_broadcasted>
    ntransform_parallel(std::tuple<Ndatacontainer...> ndata_tup, FuncT func)  {
        ntransform<Tret, PARALLEL>(ndata_tup, func);
    }

}

#endif /* end of include guard: NDATA_FUNCTIONS_HPP_J3F2DFOA */
