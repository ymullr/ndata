#include <stdarg.h>
#include "ndarray.h"
#include <assert.h>

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

size_t nd_stride(NdShape nd_shape, size_t dim){
    size_t stride=1;
    size_t j;
	for(j=dim+1;j<nd_shape.rank;j++){
			stride*=nd_shape.shape[j];
	}
	return stride;
}

size_t nd_len(NdShape nd_shape){
	size_t acc=1,i;
	for(i=0;i<nd_shape.rank;i++){
		acc*=nd_shape.shape[i];
	}
	return acc;
}
