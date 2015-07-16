/*! \file Contains looping constructs for ndatacontainers */
#ifndef NDINTERP_HPP_CUZ3R9SV
#define NDINTERP_HPP_CUZ3R9SV

#ifndef NDEBUG
#define NDINTERP_DEBUG
#endif

#include <vector>
#include <array>
#include <math.h>
#include <cassert>
#include "ndata.hpp"
#include "ndata/algorithm/numtype_adapter_fundamental.hpp"
#include "ndata/algorithm/sequences.hpp"

namespace ndata {
namespace interp {

//-----------------------------------------------------------------------------
//	HELPERS
//-----------------------------------------------------------------------------

template<class NumT>
NumT clamp(NumT x, NumT min, NumT max)
{
    if(x>max) {
        x = max;
    }

    if (x<min) {
        x = min;
    }

    return x;
}


//-----------------------------------------------------------------------------
//	CONSTANTS
//-----------------------------------------------------------------------------

namespace overflow_behaviour {

    struct zero {
        static
        void handle_istart_istop(long & i_start, long & i_stop, size_t size, float) {

            i_start = clamp(i_start, 0l, long(size));
            i_stop = std::min(i_stop, long(size));

            //also make sure we got floating points error/rounding right
            assert(i_start >= 0 and i_stop <= long(size)+1); //this is not true for other overflow behaviours
            return;
        }

        static
        void handle_overflow(long &, size_t) {
            //do nothing
            return;
        }

    };

    struct stretch {
        static
        void handle_istart_istop(long &, long &, size_t, float) {
            //nothing
            return;
        }

        static
        void handle_overflow(long & i_uold, size_t size) {
            i_uold = clamp(long(i_uold), 0l, long(size)-1l);
        }
    };

    struct cyclic {
        static
        void handle_istart_istop(long &, long &, size_t, float) {
            //nothing
            return;
        }

        static
        void handle_overflow(long & i_uold, size_t size) {
                i_uold = i_uold%long(size);
                if(i_uold < 0) {
                    i_uold += size;
                }
        }
    };

    struct throw_ {
        static
        void handle_istart_istop(long & i_start, long & i_stop, size_t size, float index_frac) {
            if (index_frac < 0 or index_frac > long(size-1)) {
                throw(std::out_of_range(""));
            }
            assert(not (i_start < 0 or i_stop > long(size)));
            return;
        }

        static
        void handle_overflow(long &, size_t) {
            //do nothing
            return;
        }

    };

    struct assert_ {
        static
        void handle_istart_istop(long & i_start, long & i_stop, size_t size, float index_frac) {
            assert(not (index_frac < 0 or index_frac > long(size-1)));

            //also make sure we also got floating points error/rounding right
            assert(not (i_start < 0 or i_stop > long(size)));
            return;
        }

        static
        void handle_overflow(long &, size_t) {
            //do nothing
            return;
        }
    };
}

namespace helpers {

    template <long ndims>
    struct run_overflow_handlers {

        template <typename ... OverflowBehaviour>
        static
        void
        do_it(
            std::tuple<OverflowBehaviour ...> overflow_behaviours,
            long * i_starts, //total hack
            long * i_stops,
            vecarray<long, ndims> shape,
            vecarray<float, ndims> index_frac
            )
        {
            static_assert(
                std::tuple_size<decltype(overflow_behaviours)>() == ndims,
                "Number of overflow behaviours in the tupledoesn't match the number of dimensions that need to be interpolated"
                );
            auto tup_ht = tuple_utilities::split_ht(std::move(overflow_behaviours));

            tup_ht.first.handle_istart_istop(*i_starts, *i_stops, shape[0], index_frac[0]);
            run_overflow_handlers<ndims-1>::do_it(
                tup_ht.second,
                i_starts++,
                i_stops++,
                shape.drop_front(),
                index_frac.drop_front()
                );
            return;
        }
    };


    template <>
    struct run_overflow_handlers<0> {

        template <typename ... OverflowBehaviour>
        static
        void do_it(
            std::tuple<OverflowBehaviour ...> overflow_behaviours,
            long *,// i_starts,
            long *,// i_stops,
            vecarray<long, 0>,// shape,
            vecarray<float, 0>// index_frac
            )
        {
            static_assert(
                std::tuple_size<decltype(overflow_behaviours)>() == 0,
                ""
                );
            return;
        }
    };


    /**
     * Process the copying the correct slice of u_old in u_new. It is necessary to use recursion
     * to deal with cyclic overflow behaviour appropriately.
     */
    template <long ndims_fold, long ndims_to_keep, typename T>
    struct copy_values_from_uold_to_unew {
        template <typename ... OverflowBehaviours>
        static
        void
        do_it(
            std::tuple<OverflowBehaviours...> overflow_behaviours,
            vecarray<size_t, ndims_fold> axis_to_fold,
            vecarray<size_t, ndims_to_keep> axis_to_keep,
            vecarray<long, ndims_fold> i_starts,
            vecarray<long, ndims_fold> i_stops, //useless?
            ndataview<T, ndims_fold+ndims_to_keep> u_slice,
            ndataview<T, ndims_fold+ndims_to_keep> u_new_slice
        )
        {
            static_assert(std::tuple_size<decltype(overflow_behaviours)>() == ndims_fold, "");
            auto tup_overfl_behav_ht = tuple_utilities::split_ht(std::move(overflow_behaviours));

            vecarray<size_t, ndims_to_keep+ndims_fold-1> axis_ranges (STATICALLY_SIZED);
            for (size_t i = 0; i < ndims_fold-1; ++i) {
                axis_ranges[i] = axis_to_fold[i+1];
            }
            for (size_t i = 0; i < ndims_to_keep; ++i) {
                axis_ranges[ndims_fold-1+i] = axis_to_keep[i];
            }

            size_t current_axis = axis_to_fold[0];

            //adjust axis numbers for the subslices
            for (size_t iax = 1; iax < axis_to_fold.size(); ++iax) {
                if(axis_to_fold[iax] > current_axis) {
                    axis_to_fold[iax]--;
                }
            }

            for (size_t iax = 0; iax < axis_to_keep.size(); ++iax) {
                if(axis_to_keep[iax] > current_axis) {
                    axis_to_keep[iax]--;
                }
            }

            long current_axis_size = u_new_slice.get_shape()[current_axis];
            //recurse with one less dimension on the subslices associated
            //with each element of the current dimension
            for (long i = 0; i < current_axis_size ; ++i) {

                long u_old_index = i_starts[0]+i;

                tup_overfl_behav_ht.first.handle_overflow(u_old_index, u_slice.get_shape()[current_axis]);

                auto u_subslice = u_slice.slice_alt(
                            vecarray<range, ndims_to_keep+ndims_fold-1>(STATICALLY_SIZED, range()),
                            axis_ranges,
                            make_vecarray(u_old_index),
                            make_vecarray(current_axis)
                            );
                auto u_new_subslice = u_new_slice.slice_alt(
                            vecarray<range, ndims_to_keep+ndims_fold-1>(STATICALLY_SIZED, range()),
                            axis_ranges,
                            make_vecarray(long(i)),
                            make_vecarray(current_axis)
                            );

                //recursive call
                copy_values_from_uold_to_unew<ndims_fold-1, ndims_to_keep, T>::do_it(
                            tup_overfl_behav_ht.second,
                            axis_to_fold.drop_front(),
                            axis_to_keep,
                            i_starts.drop_front(),
                            i_stops.drop_front(),
                            u_subslice,
                            u_new_subslice
                            );

            }
        }
    };

    template <long ndims_to_keep, typename T>
    struct copy_values_from_uold_to_unew<0, ndims_to_keep, T> {
        template <typename ... OverflowBehaviours>
        static
        void
        do_it(
            std::tuple<>,// overflow_behaviours,
            vecarray<size_t, 0>,// axis_to_fold,
            vecarray<size_t, ndims_to_keep>,// axis_to_keep,
            vecarray<long, 0>,// i_starts,
            vecarray<long, 0>,// i_stops,
            ndataview<T, ndims_to_keep> u_slice,
            ndataview<T, ndims_to_keep> u_new_slice
        )
        {
            u_new_slice.assign(u_slice);
            return;
        }
    };


}

/**
 * Linear interpolation kernel
 */
struct kern_linear {
    float kern(float x) {
        assert(fabs(x)<=ONE_SIDED_WIDTH);
        return float(1)-fabs(x);
    };

    static constexpr long ONE_SIDED_WIDTH = 1; //from zero (included) to upper bound
};

//triangular kernel reproducing linear interpolation
/**
 * nearest-neighbor kernel
 */
struct kern_nearest_neighbor {

    float kern(float x) {
        assert(fabs(x)<=ONE_SIDED_WIDTH);
        return (x >= -0.5 and x < 0.5)? 1 : 0;
    };

    static constexpr long ONE_SIDED_WIDTH = 1; //from zero (included) to upper bound
};


/**
 * Approximated tricubic interpolation with Keys' convolution kernel
 * (cf wikipedia bicubic interpolation, convolution version)
 * R. Keys, (1981). "Cubic convolution interpolation for digital image processing".
 * IEEE Transactions on Acoustics, Speech, and Signal Processing 29 (6): 1153–1160.
 */
struct kern_cubic {
    float kern(float x) {

        assert(fabs(x)<=ONE_SIDED_WIDTH);

        float dx = fabs(x);
        float dx2 = dx*dx;
        float dx3 = dx2*dx;

        const float a = -0.5;

        float convCoeff = 0;

        if (dx <= 1) {
            convCoeff = (a+2)*dx3-(a+3)*dx2+1;
        } else if (dx < 2){
            convCoeff = a*dx3-5*a*dx2+8*a*dx-4*a;
        } else {
            convCoeff = 0;
        }

        return convCoeff;
    }

    static constexpr long ONE_SIDED_WIDTH = 2; //kernel width from zero to upper bound
};

//TODO Lanczos kernel

template<class KernT, long ndims, long ndims_fold, class ContainerT, class T>
struct interpolate_inner {

    static
    nvector<T, ndims-ndims_fold>
    do_it(
            ndatacontainer<ContainerT, T, ndims> u,
            vecarray<float, ndims_fold> index_frac,
            vecarray<size_t, ndims_fold> axis
            )
    {
        static_assert(ndims != DYNAMICALLY_SIZED, "Dynamic case not implemented");
        static_assert(ndims >= ndims_fold or ndims_fold == DYNAMICALLY_SIZED or ndims == DYNAMICALLY_SIZED,
                      "The number of fractional indices must not exceed the number of dimensions");
        static_assert(ndims_fold > 0 or ndims_fold == DYNAMICALLY_SIZED, "");


        //prepare vecarrays of ranges and indices for the reduced u array
        //we fold the array by performing the interpolation on the last dimension contained in axis

        auto conv_coeffs_shape = u.get_shape();

        //conv_coeffs_shape is in effect 1D array which is going to be broadcasted along all axis (where shape[i]=1)
        //but the one on which the interpolation is performed
        //exemple shape (1, 1, u.get_shape()[axis.back()], 1, 1 ...)
        for (size_t i = 0; i < conv_coeffs_shape.size(); ++i) {
            if(i!=axis.back()) {
                conv_coeffs_shape[i]=1;
            }
        }

        nvector<float, conv_coeffs_shape.STATIC_SIZE_OR_DYNAMIC>
            conv_coeffs (conv_coeffs_shape);

        vecarray<long, ndims-1>
                indices_ (STATICALLY_SIZED, 0ul);

        //all axes but the folded-on dimension
        std::vector<size_t> all_axes_but_the_folded_one_vec = numrange(0ul, u.get_shape().size()).data_;
        all_axes_but_the_folded_one_vec.erase(
                    all_axes_but_the_folded_one_vec.begin()+axis.back()
                    );

        vecarray<size_t, ndims-1> all_axes_but_the_folded_one (all_axes_but_the_folded_one_vec);

        //all indices_[i] are 0, as there's only one element along those dimensions anyway
        //so we are in fact assigning to the whole conv_coeffs array
        conv_coeffs.slice_alt(
                make_vecarray(range()),
                make_vecarray(axis.back()),
                indices_,
                all_axes_but_the_folded_one
            ).assign_transform(
                std::make_tuple(numrange(u.get_shape()[axis.back()])),
            [&index_frac] (long ix) {
                float dx = float(ix) - index_frac.back();
                KernT kern;
                return kern.kern(dx);
            }
        );

        auto new_shape = u.get_shape();
        new_shape[axis.back()] = 1; //array is going to be reduced (by summation) along this axis

        //zero initialized
        nvector<float, new_shape.STATIC_SIZE_OR_DYNAMIC>
                u_new (new_shape, 0.f);

        //perform reduction through broadcasting with nforeach
        nforeach(
            std::tie(u_new, u, conv_coeffs),
            [=] (auto & vu_new, auto vu, auto vcc) {
                vu_new += vu * vcc;
            });

        auto new_axis_predrop = axis.drop_back();
        auto new_axis = new_axis_predrop;
        //decrement axis values to account for the folded axis if necessary
        for (size_t i = 0; i < new_axis.size(); ++i) {
            if(new_axis[i]>axis.back()) {
                new_axis[i] -= 1;
            }
        }


        return interpolate_inner<KernT, ndims-1, ndims_fold-1, T*, T>::do_it(
            u_new.slice_alt(
                vecarray<range, ndims-1>(STATICALLY_SIZED, range()),
                all_axes_but_the_folded_one,
                make_vecarray(0l),
                make_vecarray(axis.back())
                ),
                index_frac.drop_back(),
                axis.drop_back()
            );
    }
};

template <class KernT, long ndims, typename ContainerT, typename T>
struct interpolate_inner<KernT, ndims, 0, ContainerT, T> {

    static 
    nvector<T, ndims>
    do_it(
            ndatacontainer<ContainerT, T, ndims> u,
            vecarray<float, 0>,//index_frac,
            vecarray<size_t, 0>// axis
            )
    {
        return make_nvector(u);
    }
};


/**
 * Interpolate one value among a regularly sampled grid of data. The position must be passed as a
 * fraction of an index on each dimension.
 */
template<class KernT, long ndims, long ndims_fold, typename ContainerT, typename T, typename ... OverflowBehaviours>
nvector<T, ndims-ndims_fold>
interpolate (
        ndatacontainer<ContainerT, T, ndims> u,
        //sampling is "normalized" by delta so it is expressed in terms of a fraction of
        //the sourceField indices instead of a real position
        //[ndims]
        vecarray<float, ndims_fold> index_frac,
        vecarray<size_t, ndims_fold> axis,
        //[from 1 to ndims]
        std::tuple<OverflowBehaviours...> overflow_behaviours
        //std::tuple<OverflowBehaviour...> overflow_behaviours
        )
{
    vecarray<long, ndims_fold> shape (axis.dynsize());

    for (size_t i = 0; i < axis.size(); ++i) {
        shape[i] = u.get_shape()[axis[i]];
    }

    static_assert(
        sizeof...(OverflowBehaviours) == ndims_fold,
        "The number of elements in the overflow_behaviours tuple doesn't match the number of interpolation axes"
        );

    vecarray<long, ndims_fold> i_starts, i_stops;
    i_starts = i_stops = vecarray<long, ndims_fold>(axis.dynsize());

    for (size_t i = 0; i < shape.size(); ++i) {
        i_starts[i] = floor(index_frac[i])-long(KernT::ONE_SIDED_WIDTH)+1;
        i_stops[i] = ceil(index_frac[i])+long(KernT::ONE_SIDED_WIDTH);
    }

    //depending on the overflow behaviour, i_starts and i_stops may be processed by appropriate handlers at this point
    //(mostly to do some clamping and bound checking)
    helpers::run_overflow_handlers<ndims_fold>::do_it(
            overflow_behaviours,
            &i_starts[0],
            &i_stops[0],
            shape,
            index_frac
            );

    for (size_t i = 0; i < shape.size(); ++i) {
        //zero or one value for each dimension
        assert(i_starts[i] <= i_stops[i]);
    }

    //now extracting the hypercube from u where the interpolation kernel is non-zero
    //however there's a twist: for the cyclic overflow behaviour we need to be able to "wrap around" during slice,
    //so in fact we still need to iterate all the elements of the reachable dimension and apply the wrap around

    //find out the shape of the sliced hypercube
    vecarray<long, ndims> unew_shape (u.get_shape().dynsize());
    for (size_t idim = 0; idim < ndims; ++idim) {
        bool found = false;
        long i_axis;
        for (size_t iax = 0; iax < axis.size(); ++iax) {
            if(idim == axis[iax]) {
                found = true;
                i_axis = iax;
            }
        }
        if(found == true) {
            unew_shape[idim] = i_stops[i_axis] - i_starts[i_axis];
        } else {
            unew_shape[idim] = u.get_shape()[idim];
        }
    }

    nvector<T, ndims> unew (unew_shape, ndata::helpers::numtype_adapter<T>::ZERO);

    //vecarray<size_t, ndims> ndi_unew = unew.ndindex(0);
    /*
    //1D broadcastable indice_ranges in each foldable dimension
    vecarray<
            nvector<T, ndims>,
            ndims_fold
            >
            indice_ranges (u.get_shape().dynsize);

    for (size_t i = 0; i < indice_ranges.size(); ++i) {
        vecarray<long, ndims> indice_range_shape (u.get_shape().dynsize());
        for (size_t i_shape = 0; i_shape < indice_range_shape.size(); ++i_shape) {
            if(i_shape == axis[i]) {
                indice_range_shape[i_shape] = numrange(u.get_shape()[axis[i]]);
            } else {
                indice_range_shape[i_shape] = 1;
            }
        }

        indice_ranges[i] = nvector<T, ndims>();
    }
    */


    vecarray<size_t, ndims-ndims_fold> axis_to_keep (STATICALLY_SIZED);
    size_t i_rest = 0;
    for (size_t i = 0; i < u.get_shape().size(); ++i) {
        bool found_in_axis = false;
        for (size_t i_fold = 0; i_fold < axis.size(); ++i_fold) {
            if(i == axis[i_fold]) {
                found_in_axis = true;
                break;
            }
        }
        if(not found_in_axis) {
            axis_to_keep[i_rest] = i;
            i_rest++;
        }
    }

    //filling all the values of the new (reduced) array u_new from the values of the old u
    helpers::copy_values_from_uold_to_unew<ndims_fold, ndims-ndims_fold, T>::do_it(
            overflow_behaviours,
            axis,
            axis_to_keep,
            i_starts,
            i_stops,
            u.as_view(),
            unew.as_view()
            );

    /*
    for (size_t i = 0; i < unew.size(); ++i) {

        //compute the multidimensional index for the current value in the old u array
        //ndi_uold is the multidimensional index in u, which corresponds to i in unew
        //the overflowBehaviour logic comes here in play
        vecarray<size_t, ndims> ndi_uold (unew_shape.dynsize());

        for (size_t idim = 0; idim < ndims; ++idim) {
            long iUThisDim=(i_starts[idim]+long(ndi_unew[idim]));

            switch (overflow_behaviours[idim]) {
                case CYCLIC:
                    iUThisDim = iUThisDim%shape[idim];
                    if(iUThisDim < 0) {
                        iUThisDim += shape[idim];
                    }
                    break;
                case STRETCH:
                    iUThisDim = clamp(long(iUThisDim), 0l, long(shape[idim])-1l);
                    break;
                case ZERO:
                    //nothing to do
                    break;
                case THROW:
                    //nothing to do
                    break;
                case ASSERT:
                    //nothing to do
                    break;
                default:
                    abort();//shouldn't happen
                    break;
            }

            assert(iUThisDim < long(shape[idim]));
            assert(iUThisDim >= 0);
            ndi_uold[idim] = size_t(iUThisDim);
        }

        unew(ndi_unew) = u(ndi_uold); //no stride in uOld either since we grab the last dimension

        //update the multidimensional index to match the flat one for the next iteration
        //needed to locate matching values in the old u array
        unew.increment_ndindex(ndi_unew);
    }
    */

    //adjust indexfrac for new reduced array size
    for (size_t i = 0; i < index_frac.size(); ++i) {
        index_frac[i] = index_frac[i] - i_starts[i];
    }

    assert(index_frac.size() == unew_shape.size() );

    return interpolate_inner<KernT, ndims, ndims_fold, ContainerT, T>::do_it(
            unew,
            index_frac,
            axis
    );
}

//-----------------------------------------------------------------------------
//	NOW A BUNCH OF OVERLOADS FOR DEALING WITH SIMPLER CASES
//-----------------------------------------------------------------------------

//Perform on all axes, one overflow_behaviour specified for all dimensions
template<class KernT, typename OverflowBehaviour = overflow_behaviour::throw_, long ndims, typename ContainerT, class T>
T interpolate (
        ndatacontainer<ContainerT, T, ndims> u,
        vecarray<float, ndims> index_frac
        ) {
    return interpolate<KernT>(
        u,
        index_frac,
        vecarray<size_t, ndims>(numrange(size_t(ndims)).data_),
        tuple_utilities::make_uniform_tuple<ndims>(OverflowBehaviour())
        ).data_[0];//unwrap the scalar from the dim 0 nvector
}


//Perform on all axes
template<class KernT, long ndims, typename ContainerT, class T, typename ... OverflowBehaviours>
T interpolate (
        ndatacontainer<ContainerT, T, ndims> u,
        vecarray<float, ndims> index_frac,
        std::tuple<OverflowBehaviours...> overflow_behaviours
        ) {
    return interpolate<KernT>(
        u,
        index_frac,
        vecarray<size_t, ndims>(numrange(size_t(ndims)).data_),
        overflow_behaviours
        ).data_[0];//unwrap the scalar from the dim 0 nvector;
}


//one dimensional case
template<class KernT, typename OverflowBehaviour = overflow_behaviour::throw_, typename ContainerT, class T>
T interpolate (
        ndatacontainer<ContainerT, T, 1> u,
        float index_frac
        ) {
    return interpolate<KernT>(
        u,
        vecarray<float, 1>({index_frac}),
        std::make_tuple(OverflowBehaviour())
        );
}


/**
 * nD case
 */
template<long ndims>
vecarray<float, ndims> position_to_ifrac(
        vecarray<float, ndims> position,
        vecarray<float, ndims> origin,
        vecarray<float, ndims> steps
        ) {

    assert(
            position.size()==steps.size()
        and origin.size()==steps.size()
        );

    vecarray<float, ndims> index_frac = position;

    for (size_t i = 0; i < index_frac.size(); ++i) {
        index_frac[i] = (position[i] - origin[i]) / steps[i];
    }

    return index_frac;
}


/**
 * 1D case
 */
float position_to_ifrac(
        float position,
        float origin,
        float step
        )
{
    return position_to_ifrac(
                make_vecarray(position),
                make_vecarray(origin),
                make_vecarray(step)
                )[0];
}

}//end namespace interp
}//end namespace ndata

#endif /* end of include guard: NDINTERP_HPP_CUZ3R9SV */
