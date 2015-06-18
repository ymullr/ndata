/*! \file Contains looping constructs for ndatacontainers */
#ifndef NDATA_FUNCTIONS_HPP_J3F2DFOA
#define NDATA_FUNCTIONS_HPP_J3F2DFOA


#include <tuple>
#include "ndata/helpers.hpp"
#include <utility>
#include "tuple_utilities.hpp"

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

    namespace helpers {

        template <
            int loop_type,
            long idim,
            long ndims
        >
        struct dim_loop_recur {

            template <
                typename FuncT,
                typename ... Ts,
                typename ... VecarrayLong
            >
            static
            void
            do_it(
                    std::tuple<Ts*...> tup_ndata_ptrs, //pointers to data
                    FuncT  func,
                    vecarray<long, ndims-idim>  shape, //only "remaining" dimensions
                    std::tuple<
                        VecarrayLong...
                        > arr_strides
                    )
            {
                static_assert(std::tuple_size<decltype(tup_ndata_ptrs)>() == std::tuple_size<decltype(arr_strides)>(), "");
                static_assert(idim<ndims, "");
                static_assert(ndims!=DYNAMICALLY_SIZED, "not implemented");

                std::function<void(size_t)> loop_inner_block (
                            [&] (size_t i) {
                                auto tup_current_stride = tuple_utilities::tuple_transform(
                                            [] (vecarray<long, ndims-idim> strides) {return strides[0];},
                                            arr_strides
                                        );

                                auto tup_strides_tail = tuple_utilities::tuple_transform(
                                            [] (vecarray<long, ndims-idim> strides) {return strides.drop_front();},
                                            arr_strides
                                        );

                                auto tup_tup_ndata_ptr_stride = tuple_utilities::zip(
                                            tup_ndata_ptrs,
                                            tup_current_stride
                                        );

                                //add current index and stride to the data pointers
                                std::tuple<Ts*...> new_ndata_ptrs = tuple_utilities::tuple_transform(
                                            [i] (auto tup_ndata_ptr_stride) {
                                    auto ndata_ptr = std::get<0>(tup_ndata_ptr_stride);
                                    long strd = std::get<1>(tup_ndata_ptr_stride);
                                    return ndata_ptr+i*strd; },
                                tup_tup_ndata_ptr_stride
                                );

                                //run loop on next dimensions
                                dim_loop_recur<SERIAL, idim+1, ndims>::do_it(
                                            new_ndata_ptrs,
                                            func,
                                            shape.drop_front(),
                                            tup_strides_tail
                                            );
                            }
                );

                if (loop_type == SERIAL) {
                    for (size_t i = 0; i < size_t(shape[0]); ++i) {
                        loop_inner_block(i);
                    }
                } else {
#pragma omp parallel for schedule(static)
                    for (size_t i = 0; i < size_t(shape[0]); ++i) {
                        loop_inner_block(i);
                    }
                }
            }
        };

        //recursion termination idim == ndims
        template <
            int loop_type,
            long ndims
        >
        struct dim_loop_recur<loop_type, ndims, ndims> {


            template <
                    typename FuncT,
                    typename ... Ts,
                    typename ... VecarrayLong
                    >
            static
            void
            do_it(
                    std::tuple<Ts*...> tup_ndata_ptrs, //pointers to data
                    FuncT func,
                    vecarray<long, 0> , //shape, no remaining dimensions (last reached)
                    std::tuple<
                    VecarrayLong... //no remaining dimensions (last reached)
                    >
                    //arr_strides
                    )
            {
                //static_assert(idim==ndims, "");
                //apply the function to the scalars pointed by tup_ndata_ptrs
                tuple_utilities::apply(func, tup_ndata_ptrs);
            };
        };



    }


    template <
            int loop_type = SERIAL,
            long ... ndims,
            typename ... Ts,
            typename FuncT
            >
    void
    nforeach_base(std::tuple<ndataview<Ts, ndims>...> ndata_views, FuncT func)  {

        //get pointers to first element of data
        std::tuple<Ts*...> tup_ndata_ptrs = tuple_utilities::tuple_transform(
                    [] (auto ndv) {return ndv.data_+ndv.get_start_index();},
                    ndata_views
                );

        //all containers have been broadcasted and have the same shape at this point
        auto shape = std::get<0>(ndata_views).get_shape();

        //std::tuple<
        //        vecarray<long, ndims>
        //    >
        auto
        strides = tuple_utilities::tuple_transform(
                    [] (auto ndv) {return ndv.get_strides();},
                    ndata_views
                );

        helpers::dim_loop_recur<loop_type, 0, shape.STATIC_SIZE_OR_DYNAMIC>::do_it(
                    tup_ndata_ptrs,
                    func,
                    shape,
                    strides
                    );

    }

    template <int loop_type = SERIAL , typename FuncT, typename... Ndatacontainer>
    void
    nforeach(std::tuple<Ndatacontainer&...> ndata_tup_refs, FuncT func)  {
        auto ndata_views = helpers::broadcast_views(ndata_tup_refs);
        nforeach_base<loop_type>(ndata_views, func);
    }

    template <int loop_type = SERIAL , typename FuncT, typename... Ts, long ndims>
    void
    nforeach(std::tuple<ndataview<Ts, ndims>...> ndata_tup_refs, FuncT func)  {
        auto ndata_views = helpers::broadcast(ndata_tup_refs);
        nforeach_base<loop_type>(ndata_views, func);
    }


    template <typename FuncT, typename... Ndatacontainer>
    void
    nforeach_parallel(std::tuple<Ndatacontainer&...> ndata_tup_refs, FuncT func)  {
        auto ndata_views = helpers::broadcast_views(ndata_tup_refs);
        nforeach_base<PARALLEL>(ndata_views, func);
    }

    template <typename FuncT, typename... Ts, long ndims>
    void
    nforeach_parallel(std::tuple<ndataview<Ts, ndims>...> ndata_tup_refs, FuncT func)  {
        auto ndata_views = helpers::broadcast(ndata_tup_refs);
        nforeach_base<PARALLEL>(ndata_views, func);
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

        nforeach_base<loop_type>(
                    std::tuple_cat(
                        std::make_tuple(ret.as_view()),
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



    //TODO concatenate

}

#endif /* end of include guard: NDATA_FUNCTIONS_HPP_J3F2DFOA */
