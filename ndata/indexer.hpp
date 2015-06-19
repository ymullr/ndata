/*! \file */
#ifndef NDINDEXER_HPP_9FICI4GD
#define NDINDEXER_HPP_9FICI4GD

#include <vector>

#include <initializer_list>
#include <cstdlib>
#include <type_traits>

#include "vecarray.hpp"
#include "ndata.hpp"
#include "ndata/helpers.hpp"

namespace ndata {

/**
 * it is possible to represent the last element of a dimension in a slice
 * by using or offset from the end of the array by using a negative value.
 * 
 * This is a convenience alias
 */
constexpr
long END = -1;


/**
 * @brief An index-range object to used for slicing arrays
 */
struct range {
    long start;
    long stop;
    long step;

    range ():
    start(0),
    stop(END),
    step(1)
    {};

    explicit
    range (long index):
    start(index),
    stop(index+1),
    step(1)
    { };

    range (long start, long stop):
    start(start),
    stop(stop),
    step(1)
    { };

    range (long start, long stop, long step):
    start(start),
    stop(stop),
    step(step)
    { }
};

struct newdimT {};
constexpr newdimT NEWDIM = newdimT();

/**
 * @brief Helper function allowing the user to create an indexer, dimensionality is infered
 * from the number of arguments.
 */
template <typename... LongT>
auto
make_indexer(long size_0, LongT ... size_n) {
    return indexer<sizeof...(size_n)+1>(size_0, size_n...);
}

namespace helpers {

    //forward decl
    /**
     * tuple : size_t start_index, vecarray<array<long, 2>, ndims_slice> slices):
     */
    template <long ndimslices>
    indexer<ndimslices> make_indexer_from_slices(std::pair<size_t, SliceAcc<ndimslices>> pr);

    //downcasting helper
    template<long ndims>
    indexer<ndims>
    make_indexer(indexer<ndims> idx) {
        return idx;
    }

    /*
    enum rin_type {
        is_range,
        is_index,
        is_newdim
    };

    struct range_or_index_or_newdim {
        rin_type type_;
        union content {
            std::pair<size_t, range>,
            std::pair
        };
    };
    */

    template <long ndims>
    void check_shape(vecarray<long, ndims> shape) {
        for (size_t i = 0; i < shape.size(); ++i) {
            assert(shape[i] >= 0);
        }
    }
}

/**
 * @brief An indexer object constructed with a given shape can compute the "flat" index
 *  for an element in a contiguous memory array from a set of indexes along each dimension.
 */
template<long ndims>
struct indexer {

    /**
     * @brief Initialize from size along each dimension.
     */
    template<typename... LongT>
    explicit
    indexer(long shape0, LongT... shape):
        start_index_(0),
        shape_(helpers::array_from_argpack<long>(vecarray<long, 0> (), shape0, shape...))
    {
        static_assert(
            ndims != 0,
            "This overload is only available when the number of dimensions"
            " is a compile time constant"
            );
        static_assert(
                    ndims == sizeof...(shape)+1,
            "Number of parameters doesn't match the number of dimensions."
        );

        assert(shape0>=0);

        //vecarray size 0 for dynamic version is going to be a problem
        helpers::check_shape(shape_);

        strides_ = calc_strides_from_shape(shape_);
    }

    /**
     * @brief Initialize from shape with default strides
     * @param shape The shape of the indexed data (the size along along each dimension).
     *  An std::vector, an std::array or an std::initialization_list may be passed instead of a vecarray,
     *  as vecarray is implicitly constructible from either of these.
     */
    indexer(vecarray<long, ndims> shape):
        start_index_(0),
        shape_(shape),
        strides_(calc_strides_from_shape(shape))
    {
        helpers::check_shape(shape_);
    }

    /**
     * @brief construct from members
     */
    indexer(size_t start_index, vecarray<long, ndims> shape, vecarray<long, ndims> strides):
        start_index_(start_index),
        shape_(shape),
        strides_(strides)
    {
        helpers::check_shape(shape_);
    };


    //empty constructor for later assignment
    indexer() { };


    /**
     * ndindex may here be passed as a vecarray (custom type), but also work with an
     * a std::array or a std::vector (when the number of dimensions is
     * known at compile time), or an std::vector when the number of dimensions
     * is known at runtime. 
     */
    size_t index(vecarray<size_t, ndims> ndindex){
        assert(ndindex.size() == shape_.size());

        size_t indexacc=start_index_;

        //looping on array dimensions
        for (size_t i = 0; i < ndindex.size(); ++i) {
            //bound checking when debugging
            assert(long(ndindex[i]) < shape_[i]);

            indexacc+=ndindex[i]*strides_[i];
        }

        assert(indexacc < this->unsliced_size());

        return indexacc;
    }

    size_t index(vecarray<long, ndims> ndindex) {
        vecarray<size_t, ndims> rev_index (ndindex.dynsize());
        for (size_t i = 0; i < ndindex.size(); ++i) {
            rev_index[i] = reverse_negative_index(ndindex[i]);
        }
        return index(rev_index);
    }

    template <typename ... IndexT>
    size_t index(IndexT ... indices) {
        static_assert(sizeof...(indices) == ndims, "Number of indices doesn't match dimensionality");
        static_assert(ndims != DYNAMICALLY_SIZED, "Please use an overload taking a vecarray with dynamically dimensioned arrays");

        static_assert(
                    helpers::static_check_valid_indice_types<IndexT...>::value,
                    "Indices passed as argument must be of integral type"
                    );

        auto pr = slice_rec<0, 0>(
            indexer<ndims>::start_index_,
            helpers::SliceAcc<0>(),
            indices...
        );

        static_assert(
            decltype(pr.second)::STATIC_SIZE_OR_DYNAMIC == 0
            //decltype(pr.second)::STATIC_SIZE_OR_DYNAMIC == DYNAMICALLY_SIZED
                    , "");

        assert(pr.first < this->unsliced_size());

        return pr.first;
    }

    //TODO note on size_t/long indices as doxygen variable
    /**
     * @brief Returns an indexer to a slice of the indexed data.
     *
     * Note that straight indices may be passed as size_t or a long. the size_t overload may be slightly
     * faster, but the long overload provides extra functionality by dealing the negative indices
     * (which indicate an offset to the end of the array).
     *
     * @param index_or_range May be a size_t index, a long "reversible index"
     *  (like the constant ndata::END or -1), the ndata::NEWDIM constant, or an ndata::range.
     */
    template <typename... IndexOrRangeT>
    auto
    slice_indexer(IndexOrRangeT ... index_or_range) {
        //type = std::pair<size_t start_index, helpers::SliceAcc<STATIC_SIZE_OR_DYNAMIC>>
        auto pr = slice_rec<0, 0>(
            start_index_,
            helpers::SliceAcc<0>(),
            index_or_range...
        );

        return  helpers::make_indexer_from_slices(pr);
    }

    /**
     * @brief Returns an indexer to a slice of the indexed data. This overload is more verbose to call
     *  but more amenable to writing dimensionality agnostic functions by virtue of taking ranges and
     *  indices and the dmension to which they apply as pairs in vecarrays.
     * @param ranges A vecarray of std::pair containing the dimension indice and the actual indice range.
     *  Technically idim should be unsigned but using long here might be less bug prone.
     * @param indices A vecarray of std::pair containing the dimension indice and an index ay be a size_t index, a long "reversible index"
     *  (like the constant ndata::END or -1), the ndata::NEWDIM constant, or an ndata::range.
     */
    template <long nranges, long nindices = 0>
    auto
    slice_indexer_alt(
            vecarray<range, nranges>
                ranges,
            vecarray<size_t, nranges>
                axis_ranges,
            vecarray<long, nindices>
                indices = vecarray<long, 0> (ndata::STATICALLY_SIZED),
            vecarray<size_t, nindices>
                axis_indices = vecarray<size_t, 0> (ndata::STATICALLY_SIZED)
        ) -> indexer<
                (nranges == DYNAMICALLY_SIZED)?
                    DYNAMICALLY_SIZED :
                    nranges
             >
    {
        static_assert(
            helpers::is_any_dynamically_sized<nranges, nindices, ndims>::value
         or nranges == ndims-nindices,
            "Statically known number of dimensions does not match"
        );

        constexpr long RET_STATIC_SIZE =
            helpers::is_any_dynamically_sized<nranges, nindices>::value?
                DYNAMICALLY_SIZED:
                nranges;

#ifndef DNDEBUG
        assert(ranges.size() == shape_.size()-indices.size());
        assert(ranges.size() == axis_ranges.size());
        assert(indices.size() == axis_indices.size());
        assert(indices.size() + ranges.size() == shape_.size());

        for (size_t ir = 0; ir < ranges.size(); ++ir) {
            for (size_t i_ind = 0; i_ind < indices.size(); ++i_ind) {
                //make sure that no dimensions is found in ranges and in indices
                assert(axis_ranges[ir] != axis_indices[i_ind]);
                //make sure dimension index doesn't exceed the number of dimensions
                assert(
                           axis_ranges[ir] < shape_.size()
                       and axis_indices[i_ind] < shape_.size()
                      );
            }
        }


        for (size_t i1 = 0; i1 < ranges.size(); ++i1) {

            long idim = axis_ranges[i1];
            assert(idim >= 0 && idim < long(shape_.size()));

            for (size_t i2 = 0; i2 < ranges.size(); ++i2) {
                if (i1 != i2) {
                    //check uniqueness of dimension number in ranges
                    assert(axis_ranges[i1]!= axis_ranges[i2]);
                }
            }
        }

        for (size_t i1 = 0; i1 < indices.size(); ++i1) {

            long idim = axis_indices[i1];
            long ind = indices[i1];
            assert(idim >= 0 && idim < long(shape_.size()));
            assert((ind >= 0)? ind < long(shape_[idim]) : ind >= -long(shape_[idim]));

            for (size_t i2 = 0; i2 < indices.size(); ++i2) {
                if (i1 != i2) {
                    //check uniqueness of dimension number in indices
                    assert(axis_indices[i1] != axis_indices[i2]);
                }
            }
        }
#endif

        size_t ret_start_index = start_index_;

        //this could have been a vecarray of size_t, however unsigned are bug prone (underflow)
        //so a long type is a safer choice
        auto ret_shape_predrop = shape_;
        auto ret_strides_predrop = strides_;

        //updating temporary shape with new ranges
        for (size_t i = 0; i < ranges.size(); ++i) {
            size_t idim = axis_ranges[i];
            range rng = ranges[i];

            assert(rng.step!=0);
           //TODO check index
            ret_start_index += strides_[idim] * reverse_negative_index(idim, rng.start);

            assert(ret_start_index < unsliced_size() or unsliced_size() == 0); //special case unsliced_size() == 0 for 0D vectors
            assert(unsliced_size() != 0 or size() ==0);
            ret_shape_predrop[idim] = (reverse_negative_index(idim, rng.stop)
                                      -reverse_negative_index(idim, rng.start))
                                      /rng.step;
            ret_strides_predrop[idim] = strides_[idim]*rng.step;
        }

        //dropping indexed elements from shape and stride
        vecarray<size_t, indices.STATIC_SIZE_OR_DYNAMIC> droplist (indices.dynsize());

        for (size_t i_ind = 0; i_ind < indices.size(); ++i_ind) {
            size_t idim = axis_indices[i_ind];
            long ind = indices[i_ind];

            ret_start_index += ind * strides_[idim];
            droplist[i_ind] = idim;
        }

        vecarray<long, RET_STATIC_SIZE> ret_shape = ret_shape_predrop.drop(droplist);
        vecarray<long, RET_STATIC_SIZE> ret_strides = ret_strides_predrop.drop(droplist);

        return indexer<RET_STATIC_SIZE>(ret_start_index, ret_shape, ret_strides);
    }

    size_t get_start_index() {
        return start_index_;
    }

    vecarray<long, ndims> get_shape() {
        return shape_;
    }

    /**
     * Return the stride for the current dimension
     */
    vecarray<long, ndims> get_strides() {
        return strides_;
    }

    size_t size(){
        size_t acc=1;
        for(size_t i=0;i<shape_.size();i++){
            acc*=shape_[i];
        }
        return acc;
    }

    /**
     * For a given ndarray, this will increments the dimensional indices for
     * each dimension in the correct order.
     *
     * This will work with stl containers such as std::array, a std::vector or a vecarray.
     *
     * Useful to avoid deeply nested loop with a large number of dimension by
     * folding everything into one. Also makes things easier when the number of
     * dimensions is only known at runtime.
     *
     * Use at the end of a for loop when looping on all the elements of an ndarray, like :
     *
     * ndindexer<3> ndindexer (N1, N2, N3);
     *
     * vector<float> ndarray (ndindexer.size());
     *
     * //your multidimensional index
     * //all dimensional indices are 0 initialized
     * vector<size_t> ndi (3, 0); //will also work with a std::array or vecarray
     *
     * for (size_t i = 0; i < ndarray.size(); ++i) {
     *
     *     //here ndindex.index(ndi) == i
     *
     *     ndindex.increment_ndindex(ndi)
     * }
     */
    template<class ContainerT>
    void
    increment_ndindex(ContainerT & ndindex){
        assert(ndindex.size() == shape_.size());
        increment_ndindex(&ndindex[0]);
    }

    /**
     * Raw pointer version
     */
    template<typename LongT>
    void
    increment_ndindex(LongT * ndindex){
        //would also be doable with a bunch of modulo (%) ops but that might be slower/less explicit
        for (long idim = shape_.size()-1; idim >= 0; --idim) {
            //update.back-most dimensional index unless its maxxed
            if (long(ndindex[idim]) < shape_[idim]-1) {
                ndindex[idim]++;
                return;  //were done
            } else if (long(ndindex[idim]) == shape_[idim]-1) { //if maxxed
                ndindex[idim]=0; //reset this dimensional index to 0 
                //continue loop, skip to next outer dimension
                //if this is already the outermost dimension then the index is simply reset
                //no error since this is probably going to happen during the
                //last iteration of an eventual for loop
            }
        }
    }

    /**
     * Probably quite slow
     *
     * Should be implemented with some modulo ops
     */
    vecarray<size_t, ndims>
    ndindex(size_t flatindex) {

        vecarray<size_t, ndims> ndindex (shape_.dynsize(), 0);

        for (size_t i = 0; i <= flatindex; ++i) {
            increment_ndindex(ndindex);
        }

        return ndindex;
    }

    //TODO fortran_order, swap_dimensions

    /**
     * used by broadcast
     */
    template<long new_ndims>
    indexer<new_ndims>
    reshape(vecarray<long, ndims> new_shape, vecarray<long, ndims> new_strides) {
        return indexer<new_ndims>(start_index_, new_shape, new_strides);
    }

    protected:

    size_t start_index_;

    //kept for bound checking
    //this could have been a vecarray of size_t, however unsigned are bug prone (underflow)
    //so a long type is a safer choice
    vecarray<long, ndims> shape_;

    vecarray<long, ndims> strides_;

    static
    vecarray<long, ndims> 
    calc_strides_from_shape(vecarray<long, ndims> shape) {

        vecarray<long, ndims> ret (shape.dynsize());

        for (size_t i = 0; i < shape.size(); ++i) {

            long this_stride=1;
            for(size_t j=i+1;j<shape.size();j++){
                this_stride*=shape[j];
            }


            ret[i] = this_stride;
        }

        return ret;
    }

    // FIXME: A note: the slice_rec family of functions is redundant with the the vecarray based
    // slice_indexer implementation. Maybe it should go away at some point, and the variadic slice_indexer implementation
    // should rely on the vecarray based one instead. However it seems more info (namely the ordering
    // of ranges and indices and newdims) is known at compile time with the recursive slice_rec implementation.
    // It is unclear if this brings any real benefit speed or safety-wise, but this should be cleared before removal.

    /**
     * An actual index is given instead of a slice, only one element is thus
     * included and the ndindexer will have one less dimension.
     *
     * the current shape/stride element is discarded and only the start_index is
     * repositioned to the relevant position
     */
    template <
            size_t idim,
            long ndimslices,
            typename IntegerType,
            typename std::enable_if<std::is_integral<IntegerType>::value,
                IntegerType>::type* = nullptr,
            typename... SliceIndex
            >
    auto
    slice_rec(
            size_t start_ind,
            helpers::SliceAcc<ndimslices> slices,
            IntegerType index,
            SliceIndex... slice_or_index
            )
    {
        start_ind += strides_[idim]*reverse_negative_index(idim, index);

        return slice_rec<idim+1, ndimslices>(start_ind, slices, slice_or_index...);
    }

    /**
     * A full slice is given {start, stop, step}
     */
    template <size_t idim, long ndimslices, typename... SliceIndex>
    auto
    slice_rec(
            size_t start_ind,
            helpers::SliceAcc<ndimslices> slices,
            range range,
            SliceIndex... slice_or_index
            ) //-> decltype(slice_rec<idim+1, ndimsslices>(size_t, SlicesT<ndimslices>, SliceIndex...))
    {

        //would be trouble
        assert(range.step!=0); 

        start_ind += strides_[idim]*reverse_negative_index(idim, range.start);

        assert(start_ind < unsliced_size());

        long new_stride = strides_[idim]*range.step;

        helpers::shape_stride_pair full_slice = std::make_pair(
            labs(reverse_negative_index(idim, range.stop) - reverse_negative_index(idim, range.start))
                /labs(range.step),
            new_stride
        );

        //DOES increment dimension
        return slice_rec<idim+1, ndimslices+1>(start_ind, slices.append(full_slice), slice_or_index...); 
    }

    /**
     * A newdimT object is given. An new dimension with one element is added to the slice.
     */
    template <size_t idim, long ndimslices, typename... SliceIndex>
    auto
    slice_rec(
            size_t start_ind,
            helpers::SliceAcc<ndimslices> slices,
            newdimT,
            SliceIndex... slice_or_index
            ) //-> decltype(slice_rec<idim+1, ndimsslices>(size_t, SlicesT<ndimslices>, SliceIndex...))
    {
        long new_stride = 0; //shouldn't matter since only one item

        helpers::shape_stride_pair full_slice = std::make_pair(
            1,
            new_stride
        );

        //DOES increment ndimsslice
        //DOESNT increment idim
        //start_ind is untouched
        return slice_rec<idim, ndimslices+1>(start_ind, slices.append(full_slice), slice_or_index...);
    }


    //termination
    template <size_t idim, long ndimslices>
    std::pair<size_t, helpers::SliceAcc<ndimslices>>
    slice_rec(size_t start_ind, helpers::SliceAcc<ndimslices> slices) {

        static_assert(
            idim == ndims,
            ""
        );

        //not true if newdim added
        //static_assert(ndims == -1 or ndimslices<=ndims, "");
        static_assert(ndims == -1 or ndimslices>=0, "");

        return make_pair(start_ind, slices);
    }


    //TODO template specialization based on signedness
    long reverse_negative_index(size_t idim, long ind) { 
        return (ind>=0)? ind : long(shape_[idim])+(ind)+1;
    }


 private:

    /**
     * Underlying container size, bigger than size() if sliced
     *
     * not 100% exhaustive
     */
    size_t unsliced_size(){
        size_t acc=0;
        for(size_t i=0;i<shape_.size();i++){
            acc+=shape_[i]*std::abs(strides_[i]);
        }

        //start_index necessary if the slice has an "outer stride"
        //in the container
        return acc+start_index_;
    }

};


namespace helpers {

    /**
     * tuple : size_t start_index, vecarray<array<long, 2>, ndims_slice> slices):
     */
    template <long ndimslices>
    indexer<ndimslices> make_indexer_from_slices(std::pair<size_t, SliceAcc<ndimslices>> pr) {

        vecarray<long, ndimslices> shape (pr.second.dynsize());
        vecarray<long, ndimslices> strides (pr.second.dynsize());

        for (size_t i = 0; i < pr.second.size(); ++i) {
            shape_stride_pair sh_st_p = pr.second[i];
            shape[i] = sh_st_p.first;
            strides[i] = sh_st_p.second;
        }

        return indexer<ndimslices>(pr.first, shape, strides);

    }

}

}
#endif /* end of include guard: NDINDEXER_HPP_9FICI4GD */
