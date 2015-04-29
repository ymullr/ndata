#ifndef NARRAY_H_3BHPEXVR
#define NARRAY_H_3BHPEXVR

#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>
#include <stdarg.h>
#include <stddef.h>


#ifndef NDAMAXRANK
	#define NDAMAXRANK 10
#endif

#ifndef NdShape_H
#define NdShape_H

/**
 * structure holding relevant information for determining the indexation of a
 * pseudo multidimensional array. Maximum rank can be adjusted through the
 * macro NDAMAXRANK or corresponding compiler flag -DNDAMAXRANK=... (default
 * 10).
 */
//class NdShape;

struct NdShape {

	//rank of the array (ie number of dimensions)
	size_t rank;

	//shape of each dimension
	size_t shape[NDAMAXRANK];
};
typedef struct NdShape NdShape;

/**
 * This is for accessing a memory contiguous dynamic variable length array
 * as you would if it were declared as a static multidimensional array on the
 * stack (which CANNOT be heap allocated)
 *
 * USAGE:
 *
 * preparing the array:
 *
 *   struct NdShape rshp={size_t rank=n, {size_t i1_max, ..., in_max} };
 *   float *ndarr=malloc(sizeof(float)*nd_length(rshp));
 *
 * accessing the array:
 *
 *   float ndarr_elem=ndarr[nd_i(rshp, i1, ..., in)];
 *
 * don't forget to free memory:
 *
 *   free(ndarr);
 *
 *
 * Memory contiguous array are allocated like a 1D array, not pointer to
 * pointer so you cannot access it with the form array[i][j][k].... Same
 * applies for nested std::vector<> in c++.  Instead you can access it with
 * array[ndarray(size_t rank, size_t rshp, ...size_t i,size_t j,size_t k)] by using
 * this function which is easier than handling the pointer arithmetic yourself.
 *
 * Boost multiarray would be an alternative for multidimensional dynamic
 * contiguous arrays with c++, although it is said to be slow.
 *
 * see also http://stackoverflow.com/questions/4810664/how-do-i-use-arrays-in-c
 */
size_t nd_i(NdShape nd_shape, ...){
        va_list indexes;
        size_t stride,i,j;
        size_t indexacc=0;

        va_start(indexes, nd_shape);

        //looping on array dimensions
        for(i=0;i<nd_shape.rank;i++){

                //bound checking when NOT compiled with -DNDEBUG
                #ifndef NDEBUG
                size_t dim_index = va_arg(indexes,size_t);
                assert(dim_index < nd_shape.shape[i]);
                #define NDARRAY_CURRENT_INDEX dim_index
                #else
                #define NDARRAY_CURRENT_INDEX va_arg(indexes,size_t)
                #endif

                stride=1;
                //looping on "remaining" array dimensions which matter for the
                //stride
                for(j=i+1;j<nd_shape.rank;j++){
                    stride*=nd_shape.shape[j];
                }

                indexacc+=NDARRAY_CURRENT_INDEX*stride;
        }
        va_end(indexes);

        return indexacc;
}

/**
 * returns the inner stride for the given dimension. the stride is the array
 * index offset between consecutive elements.
 */
size_t nd_stride(NdShape nd_shape, size_t dim){
    size_t stride=1;
    size_t j;
    for(j=dim+1;j<nd_shape.rank;j++){
            stride*=nd_shape.shape[j];
    }
    return stride;
}

/**
 * returns the length (indices) of a ndarray of the provided shape
 */
size_t nd_len(NdShape nd_shape){
    size_t acc=1,i;
    for(i=0;i<nd_shape.rank;i++){
        acc*=nd_shape.shape[i];
    }
    return acc;
}

#endif //ifndef NdShape_H

#ifdef __cplusplus
}
#endif

#endif /* end of include guard: NARRAY_H_3BHPEXVR */
