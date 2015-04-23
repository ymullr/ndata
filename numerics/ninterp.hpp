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
#include "numerics/numtype_adapter_fundamental.hpp"

namespace ndata {

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


enum OverflowBehaviour {
    ZERO,
    STRETCH,
    CYCLIC
};

//triangular kernel reproducing linear interpolation
struct KernLinear {
    float kern(float x) {
        assert(fabs(x)<=ONESIDEDWIDTH);
        return float(1)-fabs(x);
    };

    static constexpr long ONESIDEDWIDTH = 1; //from zero (included) to upper bound
};

//triangular kernel reproducing linear interpolation
struct KernNearestNeighbor {

    float kern(float x) {
        assert(fabs(x)<=ONESIDEDWIDTH);
        return (x >= -0.5 and x < 0.5)? 1 : 0;
    };

    static constexpr long ONESIDEDWIDTH = 1; //from zero (included) to upper bound
};


//Approximated tricubic interpolation with Keys' convolution kernel
//(cf wikipedia bicubic interpolation, convolution version)
//R. Keys, (1981). "Cubic convolution interpolation for digital image processing".
//IEEE Transactions on Acoustics, Speech, and Signal Processing 29 (6): 1153–1160.
struct KernCubic {
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


template<class KernT, class T, long ndims>
struct InterpolateInner {

    static
    T interpolateInner(
            nvector<T, ndims> u,
            vecarray<float, ndims> indexFrac
            )
    {
        static_assert(ndims > 1, "Dynamic case not implemented");
        
        //uNew loses the.back dimension compared to u
        //the.back dimension is the one on which the interpolation is performed
        //
        //the other dimensions are reduced to a number of values matching the kernel's width
        auto unewShape = u.get_shape().drop_back();

        indexer<ndims-1> unewIndexer (unewShape);

        size_t unew_nvecs = unewIndexer.size();

        auto uNew = nvector<T, ndims-1>(unew_nvecs, {0});

        //multidimensional index (VectorT stride ignored) in uNew
        auto ndindex_unew = uNew.make_zero_ndindex();

        //iterating on all the vectors values of the new array uNew
        //also means iterating on all dimensions (but the last) of the old array
        //that is we iterate on all the 1D columns of the old array
        //
        //iunew_flat is the flattened multidimensional index (with VectorT stride)
        for (size_t iunew_flat = 0; iunew_flat < unew_nvecs; ++iunew_flat) {

            //compute the matching flatened multidimensional index from the
            //start of the current column in the old u array
            size_t iu0 = 0;
            for (size_t idim = 0; idim < ndims-1; ++idim) {
                iu0 += ndindex_unew[idim]*u.shape.back();
            }

            nvector<T, 1> uCol (
                        {u.shape.back()}
                    );

            for (size_t iCol = 0; iCol < uCol.size(); ++iCol) {
                uCol[iCol] = u[iu0+iCol]; //no stride in uOld either since we grab the.back dimension
            }

            //update uNew
            //recursive call, this one doesn't have subrecursion since ndims == 1
            uNew[iunew_flat] = InterpolateInner<KernT, T, 1>::interpolateInner(
                uCol,
                {indexFrac.back()}
            );

            unewIndexer.increment_ndindex(ndindex_unew);
        }

        //recurse with one less dimension
        //will terminate when ndims = 1 (see partially specialized template under)
        indexFrac.pop_back();
        return InterpolateInner<KernT, T, ndims-1>::interpolateInner(unewShape, uNew, indexFrac);
    }
};

template <class KernT, typename T>
struct InterpolateInner<KernT, T, 1> {

    static 
    T
    interpolateInner(
            ndataview<T, 1> u,
            vecarray<float, 1> indexFrac
            )
    {
        assert(
                u.get_shape().size() == 1
            and indexFrac.size() == 1
              );

        //termination condition
        //interpolate in 1D
        //... easy
        KernT kern;

        T ret = numtype_adapter<T>::ZERO;

        for (size_t i = 0; i < u.shape(0); ++i) {
            float x = float(i)-indexFrac[0];
            float convCoeff = kern.kern(x);
            ret +=  convCoeff * u[i];
        }

        return ret;
    }
};

template<class KernT, typename ContainerT, typename T, long ndims>
T interpolate (
        //flattened array of VectorT
        //TODO Note on VectorT packing/padding/casting
        ndatacontainer<ContainerT, T, ndims> u,
        //sampling is "normalized" by delta so it is expressed in terms of a fraction of
        //the sourceField indices instead of a real position
        //[ndims]
        vecarray<float, ndims> indexFrac,
        //[from 1 to ndims]
        vecarray<OverflowBehaviour, ndims> overflowBehaviours
        )
{

    auto shape = u.get_shape();

    assert(
            indexFrac.size() == shape.size()
        and overflowBehaviours.size() == shape.size()
        );

    vecarray<long, ndims> iStarts, iStops;
    iStarts = iStops = vecarray<long, ndims>(shape.dynsize());

    for (size_t i = 0; i < shape.size(); ++i) {

        OverflowBehaviour ob = overflowBehaviours[i];

        if (ob == OverflowBehaviour::STRETCH) {
            indexFrac[i] = clamp(indexFrac[i], float(0), float(shape[i])-1);
        }

        long kernWidth = KernT::ONESIDEDWIDTH*2;

        iStarts[i] = floor(indexFrac[i])-long(KernT::ONESIDEDWIDTH)+1;
        iStops[i] = iStarts[i] + kernWidth;

        switch (ob) {
            case ZERO:
                iStarts[i] = std::max(0l, iStarts[i]);
                iStops[i] = std::min(iStops[i], long(shape[i]));

                //early return in case the intersect btw the kernel and the field is empty
                if(iStarts[i] >= iStops[i]) return numtype_adapter<T>::ZERO;

                assert(iStarts[i] >= 0 and iStops[i] <= long(shape[i])); //this is not true for other overflow behaviours
                break;
            default:
                //nothing
                break;
        }

        //at least one value in each dimension
        assert(iStarts[i] < iStops[i]); 
    }

    //now extracting the hypercube from u where the kernel is non-zero

    vecarray<size_t, ndims> unewShape = shape;

    for (size_t idim = 0; idim < ndims; ++idim) {
        unewShape[idim] = iStops[idim] - iStarts[idim];
    }

    nvector<T, ndims> unew (unewShape);
    unew.fill(numtype_adapter<T>::ZERO);

    vecarray<size_t, ndims> ndindex_unew = unew.ndindex(0);

    //iterating on all the values of the new array uNew
    //and their atching values from the old array
    //
    //i is the flattened index in unew
    for (size_t i = 0; i < unew.size(); ++i) {

        //compute the multidimensional index for the current column in the u array

        //ndi_uold is the multidimensional index in u, which corresponds to i in unew
        //the overflowBehaviour logic comes here in play
        vecarray<size_t, ndims> ndi_uold = shape;

        for (size_t idim = 0; idim < ndims; ++idim) {
            long iUThisDim=(iStarts[idim]+ndindex_unew[idim]);

            switch (overflowBehaviours[idim]) {
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
                default:
                    abort();//shouldn't happen
                    break;
            }

            assert(iUThisDim < long(shape[idim]));
            assert(iUThisDim >= 0);
            ndi_uold[idim] = iUThisDim;
        }

        unew[i] = u.ndindex(ndi_uold); //no stride in uOld either since we grab the last dimension

        //update the multidimensional index to match the flat one for the next iteration
        //needed to locate matching values in the old u array
        unew.increment_ndindex(ndindex_unew);
    }

    //adjust indexfrac for new reduced array size
    for (size_t i = 0; i < indexFrac.size(); ++i) {
        indexFrac[i] = indexFrac[i] - iStarts[i];
    }

    assert(indexFrac.size() == unewShape.size() );

    return InterpolateInner<KernT, T, ndims>::interpolateInner(
            unew,
            indexFrac
    );
}

//-----------------------------------------------------------------------------
//	NOW A BUNCH OF OVERLOADS FOR DEALING WITH SIMPLER CASES
//-----------------------------------------------------------------------------

//Only one OverflowBehaviour specified for all dimensions
template<class KernT, typename ContainerT, class T, long ndims>
T interpolate (
        ndatacontainer<ContainerT, T, ndims> u,
        vecarray<float, ndims> indexFrac,
        OverflowBehaviour overflowBehaviour = OverflowBehaviour::STRETCH
        ) {
    vecarray<OverflowBehaviour, ndims> overflowBehaviours (indexFrac.dynsize());
    overflowBehaviours.fill(overflowBehaviour);
    return interpolate<KernT>(
        u,
        indexFrac,
        overflowBehaviours
        );
}

//one dimensional case
template<class KernT, typename ContainerT, class T>
T interpolate (
        ndatacontainer<ContainerT, T, 1> u,
        float indexFrac,
        OverflowBehaviour overflowBehaviour = OverflowBehaviour::STRETCH
        ) {
    return interpolate<KernT>(
        u,
        {indexFrac},
        {overflowBehaviour}
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

    vecarray<float, ndims> indexFrac = position;

    for (size_t i = 0; i < indexFrac.size(); ++i) {
        indexFrac[i] = (position[i] - origin[i]) / steps[i];
    }

    return indexFrac;
}

}

#endif /* end of include guard: NDINTERP_HPP_CUZ3R9SV */
