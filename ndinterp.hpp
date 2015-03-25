#ifndef NDINTERP_HPP_CUZ3R9SV
#define NDINTERP_HPP_CUZ3R9SV

#ifndef NDEBUG
#define NDINTERP_DEBUG
#endif

#include <vector>
#include <array>
#include <math.h>
#include <cassert>
#include "ndindexer.hpp"

using namespace std;

namespace ndinterp {

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
//R. Keys, (1981). "Cubic convolution interpolation for digital image processing". IEEE Transactions on Acoustics, Speech, and Signal Processing 29 (6): 1153–1160.
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


template<class KernT, size_t ndims, size_t vectorSize>
struct InterpolateInner {

    static
    array<float, vectorSize> interpolateInner(
            vector<size_t> shape,
            vector<array<float, vectorSize>> u, 
            vector<float> indexFrac
            )
    {
        assert(ndims > 1);
        assert(shape.size() == ndims);
        
        //uNew loses the.back dimension compared to u
        //the.back dimension is the one on which the interpolation is performed
        //
        //the other dimensions are reduced to a number of values matching the kernel's width
        auto unewShape = shape;
        unewShape.pop_back();

        ndindexer<ndims-1> unewIndexer (unewShape);

        size_t unew_nvecs = unewIndexer.size();

        auto uNew = vector<array<float, vectorSize>>(unew_nvecs, {0});

        //multidimensional index (VectorT stride ignored) in uNew
        auto ndindex_unew = vector<size_t>(unewShape.size(),0);

        //
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
                iu0 += ndindex_unew[idim]*shape.back();
            }

            vector<array<float, vectorSize>> uCol (
                    shape.back(),
                    {0}
                    );

            for (size_t iCol = 0; iCol < uCol.size(); ++iCol) {
                uCol[iCol] = u[iu0+iCol]; //no stride in uOld either since we grab the.back dimension
            }

            //update uNew
            //recursive call, this one doesn't have subrecursion since ndims == 1
            uNew[iunew_flat] = InterpolateInner<KernT, 1, vectorSize>::interpolateInner(
                vector<size_t>(1, shape.back()),
                uCol,
                vector<float>(1, indexFrac.back())
            );

            unewIndexer.increment_ndindex(ndindex_unew);
        }

        //recurse with one less dimension
        //will terminate when ndims = 1 (see partially specialized template under)
        indexFrac.pop_back();
        return InterpolateInner<KernT, ndims-1, vectorSize>::interpolateInner(unewShape, uNew, indexFrac);
    }
};

template <class KernT, size_t vectorSize>
struct InterpolateInner<KernT, 1, vectorSize> {

    static 
    array<float, vectorSize> 
    interpolateInner(
            vector<size_t> shape,
            vector<array<float, vectorSize>> u, 
            vector<float> indexFrac
            )
    {
        assert(
                shape.size() == 1
            and indexFrac.size() == 1
              );

        //termination condition
        //interpolate in 1D
        //... easy
        KernT kern;

        array<float, vectorSize> ret {0};

        for (size_t i = 0; i < shape[0]; ++i) {
            float x = float(i)-indexFrac[0];
            float convCoeff = kern.kern(x);
            for (size_t iscal = 0; iscal < vectorSize; ++iscal) {
                ret[iscal] +=  convCoeff * u[i][iscal];
            }
        }

        return ret;
    }
};

template<class KernT, size_t ndims, size_t vectorSize>
array<float, vectorSize> interpolate (
        //lists the size of each dimension of the array of VectorT
        //(the size of VectorT must not included)
        //[ndims]
        vector<size_t> shape,
        //flattened array of VectorT
        //TODO Note on VectorT packing/padding/casting
        float * u, 
        //sampling is "normalized" by delta so it is expressed in terms of a fraction of
        //the sourceField indices instead of a real position
        //[ndims]
        vector<float> indexFrac,
        //[from 1 to ndims]
        vector<OverflowBehaviour> overflowBehaviours
        )
{
    assert(shape.size()>0);
    assert(
            indexFrac.size() == ndims
        and shape.size() == ndims
        and indexFrac.size() == ndims
        );
    assert(overflowBehaviours.size() <= shape.size() and overflowBehaviours.size() > 0);

    vector<long> iStarts, iStops;
    iStarts = iStops = vector<long>(ndims);

    array<float, vectorSize> zeroVec;
    zeroVec.fill(0);

    vecarray<size_t, ndims+1> shapeWithVectorSize (0, shape);

    shapeWithVectorSize[ndims]=vectorSize;

    ndindexer<ndims+1> uIndexer (shapeWithVectorSize);

    for (size_t i = 0; i < ndims; ++i) {
        if (i >= overflowBehaviours.size()) {
            //copy last element to fill missing ones
            overflowBehaviours.push_back(
                    overflowBehaviours.back()
                    );
        }

        OverflowBehaviour ob = overflowBehaviours[i];

        if (ob == OverflowBehaviour::STRETCH) {
            indexFrac[i] = clamp(indexFrac[i], float(0), float(shape[i])-1);
        }

        long kernWidth = KernT::ONESIDEDWIDTH*2;

        iStarts[i] = floor(indexFrac[i])-long(KernT::ONESIDEDWIDTH)+1;
        iStops[i] = iStarts[i] + kernWidth;

        switch (ob) {
            case ZERO:
                iStarts[i] = max(0l, iStarts[i]);
                iStops[i] = min(iStops[i], long(shape[i]));

                //early return in case the intersect btw the kernel and the field is empty
                if(iStarts[i] >= iStops[i]) return zeroVec;

                assert(iStarts[i] >= 0 and iStops[i] <= shape[i]); //this is not true for other overflow behaviours
                break;
            default:
                //nothing
                break;
        }

        //at least one value in each dimension
        assert(iStarts[i] < iStops[i]); 
    }

    //extract the hypercube for u where the kernel is non-zero
    //the overflowBehaviour logic comes here in play
    array<size_t, ndims> unewShape;

    for (size_t idim = 0; idim < ndims; ++idim) {
        unewShape[idim] = iStops[idim] - iStarts[idim];
    }

    ndindexer<ndims> unewIndexer (unewShape);

    //unlike the field u passed as argument, we keep each vector in a distinct
    //fixed size array (safer). The packed representation of u is only better for
    //interoperation with external code.
    //
    //This is still stored continuously (fixed size arrays are stack
    //allocated), but may or may not be binary compatible with the packed
    //storage format depending on the presence or not of extra padding of the
    //std::array, which is compiler/architecture dependent.
    vector<array<float, vectorSize> > unew (unewIndexer.size(), array<float, vectorSize>());

    //multidimensional index (VectorT dimension ignored) in uNew
    auto ndindex_unew = vector<size_t>(ndims,0);

    //iterating on all the vectors values of the new array uNew
    //also means iterating on all dimensions (but the.back) of the old array
    //that is we iterate on all the 1D columns of the old array
    //
    //iunew_flat is the flattened multidimensional index (with VectorT stride)
    for (size_t iunew_flat = 0; iunew_flat < unew.size(); ++iunew_flat) {


        //compute the multidimensional index from the
        //start of the current column in the old u array
        array<size_t, ndims+1> ndi_uold;

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

            assert(iUThisDim < shape[idim]);
            assert(iUThisDim >= 0);
            ndi_uold[idim] = iUThisDim;
        }

        ndi_uold[ndims] = 0; //first element of the vector
        size_t iu0 = uIndexer.index(ndi_uold);
        for (size_t iScalar = 0; iScalar < vectorSize; ++iScalar) {
            unew[iunew_flat][iScalar] = u[iu0+iScalar]; //no stride in uOld either since we grab the last dimension
        }

        //update the multidimensional index to match the flat one for the next iteration
        //needed to locate matching values in the old u array
        unewIndexer.increment_ndindex(ndindex_unew);
    }

    //adjust indexfrac for new reduced array size
    for (size_t i = 0; i < indexFrac.size(); ++i) {
        indexFrac[i] = indexFrac[i] - iStarts[i];
    }

    assert(indexFrac.size() == unewShape.size() );

    auto retArray = InterpolateInner<KernT, ndims, vectorSize>::interpolateInner(
            vector<size_t>(unewShape.begin(), unewShape.end()),
            unew,
            indexFrac
            );

    array<float, vectorSize> retVec;
    //convert from array to VectorT
    for (size_t i = 0; i < retVec.size(); ++i) {
        retVec[i] = retArray[i];
    }

    return retVec;
}

//-----------------------------------------------------------------------------
//	NOW A BUNCH OF OVERLOADS FOR DEALING WITH SIMPLER CASES
//-----------------------------------------------------------------------------

/**
 * The main function to perform interpolation
 */
template<class KernT, size_t ndims, size_t vectorSize>
array<float, vectorSize> interpolate (
        //lists the size of each dimension of the array of VectorT
        //(the size of VectorT must not included)
        //[ndims]
        vector<size_t> shape,
        //flattened array of VectorT
        float * u, 
        //sampling is "normalized" by delta so it is expressed in terms of a fraction of
        //the sourceField indices instead of a real position
        //[ndims]
        vector<float> indexFrac,
        //[from 1 to ndims]
        OverflowBehaviour overflowBehaviour = OverflowBehaviour::STRETCH
        ) {
    return interpolate<KernT, ndims, vectorSize>(
        shape,
        u,
        indexFrac,
        vector<OverflowBehaviour>(1, overflowBehaviour)
        );
}

template<class KernT, size_t ndims>
float interpolate (
        //[ndims]
        vector<size_t> shape,
        //vector of float
        float * u, 
        //sampling is "normalized" by delta so it is expressed in terms of a fraction of
        //the sourceField indices instead of a real position
        //[ndims]
        vector<float> indexFrac,
        //[from 1 to ndims]
        OverflowBehaviour overflowBehaviour = OverflowBehaviour::STRETCH
        ) {
    return interpolate<KernT, ndims, 1>(
        shape,
        u,
        indexFrac,
        {overflowBehaviour}
        )[0];
}

template<class KernT, size_t ndims>
float interpolate (
        //[ndims]
        vector<size_t> shape,
        //vector of float
        float * u, 
        //sampling is "normalized" by delta so it is expressed in terms of a fraction of
        //the sourceField indices instead of a real position
        //[ndims]
        vector<float> indexFrac,
        //[from 1 to ndims]
        vector<OverflowBehaviour> overflowBehaviours
        ) {
    return interpolate<KernT, ndims, 1>(
        shape,
        u,
        indexFrac,
        overflowBehaviours
        )[0];
}

template<class KernT>
float interpolate (
        size_t shape,
        //array of VectorT
        float * u, 
        //sampling is "normalized" by delta so it is expressed in terms of a fraction of
        //the sourceField indices instead of a real position
        float indexFrac,
        //[from 1 to ndims]
        OverflowBehaviour overflowBehaviour = OverflowBehaviour::STRETCH
        ) {
    return interpolate<KernT, 1, 1>(
        {shape},
        u,
        {indexFrac},
        {overflowBehaviour}
        )[0];
}

template<class KernT, size_t vectorSize>
array<float, vectorSize> interpolate (
        size_t shape,
        //array of VectorT
        float * u, 
        //sampling is "normalized" by delta so it is expressed in terms of a fraction of
        //the sourceField indices instead of a real position
        float indexFrac,
        //[from 1 to ndims]
        OverflowBehaviour overflowBehaviour = OverflowBehaviour::STRETCH
        ) {
    return interpolate<KernT, 1, vectorSize>(
        {shape},
        u,
        {indexFrac},
        {overflowBehaviour}
        );
}


vector<float> positionToIfrac(
        vector<float> position,
        vector<float> origin,
        vector<float> steps
        ) {

    assert(
            position.size()==steps.size()
        and origin.size()==steps.size()
        );

    vector<float> indexFrac (position.size());

    for (size_t i = 0; i < indexFrac.size(); ++i) {
        indexFrac[i] = (position[i] - origin[i]) / steps[i];
    }

    return indexFrac;
}

}

#endif /* end of include guard: NDINTERP_HPP_CUZ3R9SV */
