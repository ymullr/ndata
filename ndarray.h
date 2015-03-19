#include <stdarg.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

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
 *   struct NdShape rshp={size_t rank=n, {size_t i1_max, ..., in_max} };
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
size_t nd_i(NdShape nd_shape, ...);

/** 
 * returns the inner stride for the given dimension. the stride is the array
 * index offset between consecutive elements.
 */
size_t nd_stride(NdShape nd_shape, size_t dimension);

/**
 * returns the length (indices) of a ndarray of the provided shape
 */
size_t nd_len(NdShape nd_shape);


#endif //ifndef NdShape_H

#ifdef __cplusplus
}
#endif
