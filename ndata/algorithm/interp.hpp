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


enum overflow_behaviour {
    ZERO,
    STRETCH,
    CYCLIC,
    THROW
};

/**
 * Linear interpolation kernel
 */
struct kern_linear {
    float kern(float x) {
        assert(fabs(x)<=ONESIDEDWIDTH);
        return float(1)-fabs(x);
    };

    static constexpr long ONESIDEDWIDTH = 1; //from zero (included) to upper bound
};

//triangular kernel reproducing linear interpolation
/**
 * nearest-neighbor kernel
 */
struct kern_nearest_neighbor {

    float kern(float x) {
        assert(fabs(x)<=ONESIDEDWIDTH);
        return (x >= -0.5 and x < 0.5)? 1 : 0;
    };

    static constexpr long ONESIDEDWIDTH = 1; //from zero (included) to upper bound
};


/**
 * Approximated tricubic interpolation with Keys' convolution kernel
 * (cf wikipedia bicubic interpolation, convolution version)
 * R. Keys, (1981). "Cubic convolution interpolation for digital image processing".
 * IEEE Transactions on Acoustics, Speech, and Signal Processing 29 (6): 1153–1160.
 */
struct kern_cubic {
    float kern(float x) {

        assert(fabs(x)<=ONESIDEDWIDTH);

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
    };

    static constexpr long ONESIDEDWIDTH = 2; //kernel width from zero to upper bound
};

//TODO Lanczos kernel

template<class KernT, long ndims, class ContainerT, class T>
struct interpolate_inner {

    static
    T do_it(
            ndatacontainer<ContainerT, T, ndims> u,
            vecarray<float, ndims> index_frac
            )
    {
        static_assert(ndims > 1, "Dynamic case not implemented");
        
        //u_new loses the.back dimension compared to u
        //the.back dimension is the one on which the interpolation is performed
        nvector<T, ndims-1> u_new (u.get_shape().drop_back(), 0);

        auto ndi_unew = u_new.ndindex(0);

        //iterating on all the values of the new array u_new
        //also means iterating on all dimensions (but the last) of the old array
        //that is we iterate on all the 1D columns of the old array
        //
        //iunew_flat is the flattened multidimensional index (with VectorT stride)
        for (size_t iunew_flat = 0; iunew_flat < u_new.size(); ++iunew_flat) {

            //setting the matching flatened multidimensional index from the
            //start of the current column in the old u array
            //size_t iu0 = 0;
            auto ndi_u = ndi_unew.append(0);

            //get a slice from u matching the missing dimension of u_new
            nvector<T, 1> u_col (make_indexer(u.get_shape().back()), 0);
            for (size_t iCol = 0; iCol < u_col.size(); ++iCol) {
                u_col(iCol) = u(ndi_u); //no stride in uOld either since we grab the.back dimension
                ndi_u.back()++;
            }

            //update u_new
            //recursive call, this one doesn't have subrecursion since ndims == 1
            u_new(ndi_unew) = interpolate_inner<KernT, 1, ContainerT, T>::do_it(
                u_col,
                make_vecarray(index_frac.back())
            );

            u_new.increment_ndindex(ndi_unew);
        }

        //recurse with one less dimension
        //will terminate when ndims = 1 (see partially specialized template under)
        return interpolate_inner<KernT, ndims-1, std::vector<T>, T>::do_it(u_new, index_frac.drop_back());
    }
};

template <class KernT, typename ContainerT, typename T>
struct interpolate_inner<KernT, 1, ContainerT, T> {

    static 
    T
    do_it(
            ndatacontainer<ContainerT, T, 1> u,
            vecarray<float, 1> index_frac
            )
    {
        assert(
                u.get_shape().size() == 1
            and index_frac.size() == 1
              );

        //termination condition
        //interpolate in 1D
        //... easy
        KernT kern;

        T ret = helpers::numtype_adapter<T>::ZERO;

        for (long i = 0; i < u.get_shape()[0]; ++i) {
            float x = float(i)-index_frac[0];
            float convCoeff = kern.kern(x);
            ret +=  convCoeff * u[i];
        }

        return ret;
    }
};

/**
 * Interpolate one value among a regularly sampled grid of data. The position must be passed as a
 * fraction of an index on each dimension.
 */
template<class KernT, long ndims, typename ContainerT, typename T>
T interpolate (
        //flattened array of VectorT
        //TODO Note on VectorT packing/padding/casting
        ndatacontainer<ContainerT, T, ndims> u,
        //sampling is "normalized" by delta so it is expressed in terms of a fraction of
        //the sourceField indices instead of a real position
        //[ndims]
        vecarray<float, ndims> index_frac,
        //[from 1 to ndims]
        vecarray<overflow_behaviour, ndims> overflow_behaviours
        )
{
    auto shape = u.get_shape();

    assert(
            index_frac.size() == shape.size()
        and overflow_behaviours.size() == shape.size()
        );

    vecarray<long, ndims> i_starts, i_stops;
    i_starts = i_stops = vecarray<long, ndims>(shape.dynsize());

    for (size_t i = 0; i < shape.size(); ++i) {

        overflow_behaviour ob = overflow_behaviours[i];

        if (ob == overflow_behaviour::STRETCH) {
            index_frac[i] = clamp(index_frac[i], float(0), float(shape[i])-1);
        }

        long kernWidth = KernT::ONESIDEDWIDTH*2;

        i_starts[i] = floor(index_frac[i])-long(KernT::ONESIDEDWIDTH)+1;
        i_stops[i] = i_starts[i] + kernWidth;

        switch (ob) {
            case ZERO:
                i_starts[i] = std::max(0l, i_starts[i]);
                i_stops[i] = std::min(i_stops[i], long(shape[i]));

                //early return in case the intersect btw the kernel and the field is empty
                if(i_starts[i] >= i_stops[i]) return helpers::numtype_adapter<T>::ZERO;

                assert(i_starts[i] >= 0 and i_stops[i] <= long(shape[i])); //this is not true for other overflow behaviours
                break;
            case THROW:
                if (i_starts[i] < 0 or i_stops[i] >= long(shape[i])) {
                    throw(std::out_of_range(""));
                }
            default:
                //nothing
                break;
        }

        //at least one value in each dimension
        assert(i_starts[i] < i_stops[i]);
    }

    //now extracting the hypercube from u where the kernel is non-zero

    vecarray<long, ndims> unew_shape = shape;

    for (size_t idim = 0; idim < ndims; ++idim) {
        unew_shape[idim] = i_stops[idim] - i_starts[idim];
    }

    nvector<T, ndims> unew (unew_shape, helpers::numtype_adapter<T>::ZERO);

    vecarray<size_t, ndims> ndi_unew = unew.ndindex(0);

    //iterating on all the values of the new array u_new
    //and their matching values from the old array
    //
    //i is the flattened index in unew
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

    //adjust indexfrac for new reduced array size
    for (size_t i = 0; i < index_frac.size(); ++i) {
        index_frac[i] = index_frac[i] - i_starts[i];
    }

    assert(index_frac.size() == unew_shape.size() );

    return interpolate_inner<KernT, ndims, std::vector<T>, T>::do_it(
            unew,
            index_frac
    );
}

//-----------------------------------------------------------------------------
//	NOW A BUNCH OF OVERLOADS FOR DEALING WITH SIMPLER CASES
//-----------------------------------------------------------------------------

//Only one overflow_behaviour specified for all dimensions
template<class KernT, long ndims, typename ContainerT, class T>
T interpolate (
        ndatacontainer<ContainerT, T, ndims> u,
        vecarray<float, ndims> index_frac,
        overflow_behaviour overflowBehaviour = overflow_behaviour::STRETCH
        ) {
    vecarray<overflow_behaviour, ndims> overflow_behaviours (index_frac.dynsize());
    overflow_behaviours.fill(overflowBehaviour);
    return interpolate<KernT>(
        u,
        index_frac,
        overflow_behaviours
        );
}

//one dimensional case
template<class KernT, typename ContainerT, class T>
T interpolate (
        ndatacontainer<ContainerT, T, 1> u,
        float index_frac,
        overflow_behaviour overflowBehaviour = overflow_behaviour::STRETCH
        ) {
    return interpolate<KernT>(
        u,
        vecarray<float, 1>({index_frac}),
        vecarray<overflow_behaviour, 1>({overflowBehaviour})
        );
}

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

}//end namespace interp
}//end namespace ndata

#endif /* end of include guard: NDINTERP_HPP_CUZ3R9SV */
