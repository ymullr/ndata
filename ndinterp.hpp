
#ifndef NDINTERP_HPP_CUZ3R9SV
#define NDINTERP_HPP_CUZ3R9SV

#include <vector>
#include <math.h>
#include <cassert>
#include "ndarray.h"

using namespace std;

namespace Ndinterp {

//-----------------------------------------------------------------------------
//	HELPERS
//-----------------------------------------------------------------------------

template<typename NumT>
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

enum OverflowBehaviour {
    ZERO,
    INFINITE,
    CYCLIC
};

auto defaultOverflowBehaviour = vector<OverflowBehaviour>(1, OverflowBehaviour::ZERO);

/*
class NdInterp<typename KernT, typename VectorT, size_t ndims>{
    public :
        NdInterp(array<ndims, float> ){};


        Array3f position_to_ifrac(Array3f position, Array3f origin, Array3f steps, Arr3ui sizes) {

            Array3f i_frac;

            for (size_t i = 0; i < 3; ++i) {
                i_frac[i] = (position[i] - origin[i]) / steps[i];
            }

            //wrap around x,y
            for (size_t i = 0; i < 2; ++i) {
                i_frac[i] = fmod(i_frac[i], sizes.v[i]);//i_frac always > 0 at that point
            }

            
            if (i_frac[2]<0 || i_frac[2]>sizes.v[2]-1) {
#ifdef NOAUTOCONV
                //bound check on z, but doesn't play well with autoconv
                cerr<<"Error : Position index is out of bouu_shape, can't interpolate value. Aborting."
                    <<endl
                    <<"Index (fractional for interpolation) = ("
                    <<endl
                    <<i_frac[0]<<", "
                    <<i_frac[1]<<", "
                    <<i_frac[2]<<") "
                    <<" vs number of points = "
                    <<"Nx="<<sizes.v[0]<<", "
                    <<"Ny="<<sizes.v[1]<<", "
                    <<"Nz="<<sizes.v[2]<<", ";
                throw("");
#else 
                //no checking when autoconv is enabled bc points could be carried away from domain
                //here we just clamp them inside the domain so that it is prolonged to infty
                //should use 0 pading instead? //TODO
                i_frac[2] = clamp(i_frac[2], 0, sizes.v[2]);
#endif
            };

            return i_frac;
        }

}
*/

//triangular kernel reproducing linear interpolation
class KernLinear {
    public: 
        float kern(float x) {

            float dx = abs(x);


            float conv_coeff = 0;

            if (dx < 1) {
                conv_coeff = 1-dx;
            } else {
                conv_coeff = 0;
            }

            return conv_coeff;
        };

        static constexpr long ONESIDEDWIDTH = 1; //from zero (included) to upper bound
};

//Approximated tricubic interpolation with Keys' convolution kernel
//(cf wikipedia bicubic interpolation, convolution version)
//R. Keys, (1981). "Cubic convolution interpolation for digital image processing". IEEE Transactions on Acoustics, Speech, and Signal Processing 29 (6): 1153–1160.
class KernCubic {
    public: 
        float kern(float x) {

            float dx = abs(x);
            float dx2 = dx*dx;
            float dx3 = dx2*dx;

            const float a = -0.5;

            float conv_coeff = 0;

            if (dx <= 1) {
                conv_coeff = (a+2)*dx3-(a+3)*dx2+1;
            } else if (dx < 2){
                conv_coeff = a*dx3-5*a*dx2+8*a*dx-4*a;
            } else {
                conv_coeff = 0;
            }

            return conv_coeff;
        };

        static constexpr long ONESIDEDWIDTH = 2; //kernel width from zero to upper bound
};


template<typename KernT, typename VectorT, size_t ndims>
VectorT interpolate (
        //lists the size of each dimension of the array of VectorT
        //(the size of VectorT must not included)
        //[ndims]
        vector<size_t> shape,
        //flattened vector of VectorT
        float * u, 
        //sampling is "normalized" by delta so it is expressed in terms of a fraction of
        //the sourceField indices instead of a real position
        //[ndims]
        vector<float> index_frac,
        //[from 1 to ndims]
        vector<OverflowBehaviour> overflowBehaviours = defaultOverflowBehaviour
        )
{
    assert(shape.size()>0);
    assert(index_frac.size() == shape.size() == ndims);
    assert(overflowBehaviours.size() <= shape.size() and overflowBehaviours.size() > 0);

    vector<long> i_starts, i_stops, u_strides;
    i_starts = i_stops = u_strides = vector<long>(ndims);

    VectorT zeroVec ();
    for (size_t i = 0; i < zeroVec.size(); ++i) {
        zeroVec[i] = 0;
    }

    NdShape u_shape = {ndims+1, {}};

    assert(ndims <= NDAMAXRANK-1);
    for (size_t i = 0; i < shape.size(); ++i) {
        u_shape.shape[i] = shape[i];
    }

    //including the VectorT dimension
    u_shape.shape[shape.size()] = zeroVec.size();

    for (size_t i = 0; i < ndims; ++i) {
        OverflowBehaviour ob;

        if (i < overflowBehaviours.size()) ob = overflowBehaviours[i];
        else ob = overflowBehaviours.back();

        long kernWidth = KernT::ONESIDEDWIDTH*2;

        if (ob == OverflowBehaviour::INFINITE) {
            index_frac[i] = clamp(index_frac[i], 0, shape[i]-1);
        }

        i_starts[i] = floor(index_frac[i])-1;
        i_stops[i] = i_starts[i] + kernWidth;

        switch (ob) {
            case ZERO:
                i_starts[i] = max(0l, i_starts[i]);
                i_stops[i] = min(i_stops[i], long(shape[i]));

                //early return in case the intersect btw the kernel and the field is empty
                if(i_starts[i] >= i_stops[i]) return zeroVec;

                assert(i_starts[i] >= 0 and i_stops[i] <= shape[i]); //this is not true for other overflow behaviours
                break;
            case default:
                //nothing
                break;
        }

        //at least one value in each dimension
        assert(i_starts[i] < i_stops[i]); 

        u_strides[i] = nd_stride(u_shape, i);
    }

    //extract the hypercube for u where the kernel is non-zero
    //the overflowBehaviour logic comes here in play
    vector<size_t> unew_shape (ndims);
    //NdShape unew_shape_scalar = {ndims+1, {}};
    //vector<size_t> unew_strides (ndims);

    size_t unew_size = 1;
    for (size_t idim = 0; idim < ndims; ++idim) {
        unew_shape[idim] = i_stops[idim] - i_starts[idim];
        //unew_shape_scalar.shape[idim] = unew_shape[idim];
        unew_size *= unew_shape[idim];
    }

    //including the VectorT dimension
    //unew_shape_scalar.shape[ndims] = zeroVec.size();
    unew_size *= zeroVec.size();

    //for (size_t idim = 0; idim < ndims; ++idim) {
    //    unew_strides[idim] = nd_stride(unew_shape_scalar, idim);
    //}

    vector<float> unew (unew_size, 0);

    //multidimensional index (VectorT dimension ignored) in u_new
    auto i_unew = vector<size_t>(ndims-1,0);

    //iterating on all the vectors values of the new array u_new
    //also means iterating on all dimensions (but the.back) of the old array
    //that is we iterate on all the 1D columns of the old array
    //
    //i_unew_flat is the flattened multidimensional index (with VectorT stride)
    for (size_t i_unew_flat = 0; i_unew_flat < u_new_size; i_unew_flat +=zeroVec.size()) {

        //update the multidimensional index to match the flat one
        //needed to locate matching values in the old u array
        //would also be doable with a bunch of modulo (%) ops but that might be slower/less explicit
        for (size_t idim = i_unew.size()-1; idim >= 0; --idim) {
            //update.back-most dimensional index unless its maxxed
            if (i_unew[idim] < shape[idim]-1) {
                i_unew[idim]++;
                break;  //were done
            } else { //if maxxed
                i_unew[idim]=0; //reset this dimensional index to 0 
                //continue loop, skip to next outer dimension
            }
        }

        //compute the matching flatened multidimensional index from the
        //start of the current column in the old u array
        size_t i_u0 = 0;
        for (size_t idim = 0; idim < ndims-1; ++idim) {
            size_t i_u_this_dim=(i_starts[idim]+i_unew[idim]);

            switch (ob) {
                case CYCLIC:
                    i_u_this_dim = %shape[idim];
                    break;
                case INFINITE:
                    i_u_this_dim = clamp(i_u_this_dim, 0, shape[idim]-1);
                    break;
                case ZERO:
                    //nothing
                    break;
                case default:
                    throw("");//shouldn't happen
                    break;
            }

            assert(i_u_this_dim < shape[dim] and i_u_this_dim >= 0);

            i_u0 += i_u_this_dim*u_strides[idim];
        }

        for (size_t i_scalar = 0; i_scalar < zeroVec.size(); ++i_scalar) {
            u_new[i_unew_flat+i_scalar] = u[i_u0+i_scalar]; //no stride in u_old either since we grab the.back dimension
        }
    }

    assert(...)

    interpolate_inner(...);
}

void increment_multidimensional_index(){};

template<>
interpolate_inner()
{
    ...

    //termination condition
    if (shape.size() == 1) {

    } 
    
    //size_t u_old_nvecs = 1;
    //for (size_t i = 0; i < shape.length()-1; ++i) u_old_nvecs *= shape[i];

    //u_new loses the.back dimension compared to u
    //the.back dimension is the one on which the interpolation is performed
    //
    //the other dimensions are reduced to a number of values matching the kernel's width
    auto unew_shape = shape;
    unew_shape.pop_back();

    size_t u_new_nvecs = 1;
    for (size_t idim = 0; idim < ndims-1; ++idim) {
        u_new_nvecs *= unew_shape[idim];
    }

    _____________________

    size_t u_new_nscalars = u_new_nvecs*zeroVec.size();

    auto u_new = vector<float>(u_new_nscalars,0);

    //multidimensional index (VectorT stride ignored) in u_new
    auto i_unew = vector<size_t>(ndims-1,0);

    //
    //iterating on all the vectors values of the new array u_new
    //also means iterating on all dimensions (but the.back) of the old array
    //that is we iterate on all the 1D columns of the old array
    //
    //i_unew_flat is the flattened multidimensional index (with VectorT stride)
    for (size_t i_unew_flat = 0; i_unew_flat < u_new_nscalars; i_unew_flat +=zeroVec.size()) {

        //update the multidimensional index to match the flat one
        //needed to locate matching values in the old u array
        //would also be doable with a bunch of modulo (%) ops but that might be slower/less explicit
        for (size_t idim = i_unew.size()-1; idim >= 0; --idim) {
            //update.back-most dimensional index unless its maxxed
            if (i_unew[idim] < shape[idim]-1) {
                i_unew[idim]++;
                break;  //were done
            } else { //if maxxed
                i_unew[idim]=0; //reset this dimensional index to 0 
                //continue loop, skip to next outer dimension
            }
        }

        //compute the matching flatened multidimensional index from the
        //start of the current column in the old u array
        size_t i_u0 = 0;
        for (size_t idim = 0; idim < ndims-1; ++idim) {
            i_u0 += (i_starts[idim]+i_unew[idim])*u_strides[idim];
        }

        auto u_col = vector<float>(
                (i_stops.back() - i_starts.back()) * zeroVec.size(),
                0
                );

        for (size_t i_col = 0; i_col < u_col.size(); ++i_col) {
            u_col[i_col] = u[i_u0+i_col]; //no stride in u_old either since we grab the.back dimension
        }

    _____________________

        auto index_frac_col = vector<float>(1, index_frac.back());
        auto ob_col = vector<OverflowBehaviour>(1, overflowBehaviours.back());

        auto shape_col = vector<size_t>(1, shape.back());

        //recursive call, this one doesn't have subrecursion since ndims == 1
        VectorT u_new_vec = interpolate<KernT, VectorT, 1>(shape_col, u_col, index_frac_col, ob_col);

        //update u_new
        for (size_t i = 0; i < u_new_vec.size(); ++i) {
            u_new[i_unew_flat+i] = u_new_vec[i];
        }

#ifdef DEBUG
        for (size_t idim = 0; idim < ndims-1; ++idim) {
            assert(i_unew[idim] == unew_shape[idim]-1);
        }
#endif

    }

    //recurse with one less dimension
    index_frac.pop_back();

    if (overflowBehaviours.size() == ndims) {
        overflowBehaviours.pop_back();
    }

    return interpolate_inner<KernT, VectorT, ndims-1>(unew_shape, u_new, index_frac, overflowBehaviours);
}

/*
template<typename KernT, typename VectorT, size_t ndims>
VectorT interpolate (
        //sampling is "normalized" by delta so it is expressed in terms of a fraction of
        //the sourceField indices instead of a real position
        //[ndims]
        vector<size_t> sizes,
        //flattened vector of VectorT
        vector<float> u, 
        //[ndims]
        vector<float> index_frac,
        //[ndims]
        vector<OverflowBehaviour> OverflowBehaviours = defaultOverflowBehaviour
        )
{
    assert(index_frac.size() == sizes.size());
    assert(OverflowBehaviours.size() <= sizes.size());

    long ixfloor=floor(index_frac[0]),
        iyfloor=floor(index_frac[1]), 
        izfloor=floor(index_frac[2]);

    //indices used for the convolution
    //ixyzfloor is >= 0
    //i_start is >= -1 (cyclic/wrap around index)
    //except for z where i_start is >= 0
    long
        i_start[3] = {
            (ixfloor-1),
            (iyfloor-1),
            (izfloor <= 0)? 0 : (izfloor-1) //handle special case at boundary
        },
        i_stop[3];
   
    for (long i = 0; i < 3; ++i) {
        //i_shift max = 4 it will be the high bound of the half open interval [i_start: i_stop[


        if (i == 0 or i == 1) { //special treatment on x and y for wrap around
            i_stop[i] = i_start[i]+4;
        } else {

            //TODO
            //autoconv may bring points outside the field, maybe boundary check should be
            //disabled and we should use zero padding and issue a warning?
            //this is to enable infinite zero padding at the boundaries
            //on z only
            for (long i_shift = 1; i_shift < 5; ++i_shift) {
                if (i_start[i] + i_shift == sizes.v[i] or i_shift == 4) {
                    //reached boundary
                    i_stop[i] = i_start[i]+ i_shift;
                    assert(i_stop[2] <= sizes.v[2]);
                    break;
                }
            }
            assert(i_start[2] >= 0 && i_start[2] <= sizes.v[2]);
        }
        assert(i_stop[i] <= i_start[i]+4);
        assert(i_stop[i] > i_start[i]);// && i_stop[i] <= sizes.v[i]);
    }

    VectorT ret = VectorT();
    ret.fill(0);

    size_t 
        stridez = ret.size(),
        stridey = sizes.v[2] * stridez,
        stridex = stridey * sizes.v[1];

    ////DEBUG
    ////nearest neighbor
    //for (int i = 0; i < array_len; ++i)
    //{
    //    ret[i] = u[stridex * ixfloor + stridey * iyfloor  * izfloor + i];
    //}
    ////DEBUG
    //float min_dx = 100,
    //      max_dx = -100;
    //return ret;

    #define POW2(x) (x)*(x)


    float xx2[4],
          yy2[4],
          zz2[4];

    //precomputing the cartesian distances to neighboring nodes (for speed)
    for (long i = 0; i < i_stop[0]-i_start[0]; ++i) {
        xx2[i] = POW2(index_frac[0] - (i_start[0] + i));
    }
    for (long i = 0; i < i_stop[1]-i_start[1]; ++i) {
        yy2[i] = POW2(index_frac[1] - (i_start[1] + i));
    }
    for (long i = 0; i < i_stop[2]-i_start[2]; ++i) {
        zz2[i] = POW2(index_frac[2] - (i_start[2] + i));
    }

    //DEBUG
    //for (int i = 0; i < 4; ++i)
    //{
    //        cout<<" | xx2: "<<xx2[i]<<endl;
    //        cout<<" | yy2: "<<yy2[i]<<endl;
    //        cout<<" | zz2: "<<zz2[i]<<endl;
    //}

#ifdef DEBUG
    float coeffsum = 0;
#endif
    

    for (long i0 = 0; i0 < i_stop[0]-i_start[0]; ++i0) {
        //modulo number of points in direction to handle wrap around on x and y
        size_t ix = (i_start[0] + i0)%sizes.v[0];
        assert(ix < sizes.v[0]);

        for (long i1 = 0; i1 < i_stop[1]-i_start[1]; ++i1) {
            //modulo number of points in direction to handle wrap around on x and y
            size_t iy = (i_start[1] + i1)%sizes.v[1];
            assert(iy < sizes.v[1]);

            for (long i2 = 0 ; i2 < i_stop[2]-i_start[2]; ++i2) {
                size_t iz = i_start[2] + i2;

                float dx2 = xx2[i0] + yy2[i1] + zz2[i2];
                float dx = sqrt(dx2);
                float dx3 = dx2*dx;

                

#ifdef DEBUG
                assert(conv_coeff <= 1 && conv_coeff >= -1);
                coeffsum += conv_coeff;
#endif

                for (int i = 0; i < ret.size(); ++i)
                {
                    ret[i] += 
                        conv_coeff * 
#ifdef DEBUG            
                        u.at(stridex * ix + stridey * iy + iz * stridez + i) //bound checked
#else
                        u[stridex * ix + stridey * iy + iz * stridez + i]
#endif
                        );
                    
                }

                assert(stridex * ix + stridey * iy + iz * stridez 
                        <= sizes.v[0]*sizes.v[1]*sizes.v[2]*3);
            }
        }
    }


#ifdef DEBUG
    cout<<"coeffsum: "<<coeffsum<<endl;
#endif
    //DEBUG
    //cout<<" | max_dx: "<<max_dx<<endl;
    //cout<<" | min_dx: "<<min_dx<<endl;
    //cout<<" | index_frac[0]: "<<index_frac[0]<<endl;
    //cout<<" | i_start[0]: "<<i_start[0]<<endl;
    //cout<<" | i_stop[0]: "<<i_stop[0]<<endl;
    
    //cout<<" | index_frac: "<<
    //    endl<<index_frac<<endl;
    
    return ret;
}
*/


//turn a numerical function into a look-up table based variant
//function<float(float, float, float)>
//inline function_approx_interp(float (*func) (float), size_t n, float max_val) {
//    float width = max_val;
//    float min_val = 0;
//    vector<float> precalc (n*n*n);
//    NdShape LUT_shp = {3, {n, n, n}};
//    for (unsigned int i_x = 0; i_x < n; ++i_x)
//        for (unsigned int i_y = 0; i_y < n; ++i_y)
//            for (unsigned int i_z = 0; i_z < n; ++i_z)
//                precalc[i_x*n*n + i_y*n + i_z] = func(
//                        sqrt(
//                             pow(i_x*width/(n-1),2)
//                            +pow(i_y*width/(n-1)*squeeze_coeff,2)
//                            +pow(i_z*width/(n-1)*squeeze_coeff,2)
//                            )
//                        ); 
//
//    return [=] (float x1, float x2, float x3) {
//        assert(x1<=max_val);
//        assert(x2<=max_val);
//        assert(x3<=max_val);
//        //abs bc we assume symmetry around 0
//        int i1 = round(abs(x1)/width*(n-1));
//        int i2 = round(abs(x2)/width*(n-1));
//        int i3 = round(abs(x3)/width*(n-1));
//        assert(i1 >= 0);
//        assert(i2 >= 0);
//        assert(i3 >= 0);
//        assert(i1<n);
//        assert(i2<n);
//        assert(i3<n);
//
//        return precalc[n*(i1*n+i2)+i3];
//    };
//}
//
//

}

#endif /* end of include guard: NDINTERP_HPP_CUZ3R9SV */
