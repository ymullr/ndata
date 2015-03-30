#ifndef HELPERS_HPP_YIOQCJSG
#define HELPERS_HPP_YIOQCJSG

#include "ndata.hpp"

namespace ndata {

namespace helpers {


    using ShapeStridePair = pair<size_t, long>;

    // ( shape[i], strides[i] )
    template<size_t N> using SliceAcc = vecarray<ShapeStridePair, N>;

    /**
     *
     *
     * tuple : size_t start_index, vecarray<array<long, 2>, ndims_slice> slices):
     *
     */
    template <long ndimslices>
    indexer<ndimslices> make_indexer_helper(pair<size_t, SliceAcc<ndimslices>> pr) {

        vecarray<size_t, ndimslices> shape (pr.second.dynsize());
        vecarray<long, ndimslices> strides (pr.second.dynsize());

        for (size_t i = 0; i < pr.second.size(); ++i) {
            ShapeStridePair sh_st_p = pr.second[i];
            shape[i] = sh_st_p.first;
            strides[i] = sh_st_p.second;
        }

        return indexer<ndimslices>(pr.first, shape, strides);

    }

    template <typename IndexerOrView>
    IndexerOrView
    broadcast_helper(IndexerOrView this_ , IndexerOrView target) {
        auto targ_shape = target.get_shape();

        auto shape_ = this_.get_shape();
        auto strides_ = this_.get_strides();

        decltype(shape_) new_shape;
        decltype(strides_) new_strides;

        assert(shape_.size() <= targ_shape.size());
        size_t new_size = targ_shape.size();


        for (size_t is = 0; is < new_size; ++is) {
            size_t i_target = targ_shape.size()-is-1;

            if (is < shape_.size()) {
                size_t i_this = shape_.size() - is - 1;
                assert(
                            targ_shape[i_target] == shape_[i_this]
                         or shape_[i_this] == 1
                    );
                new_shape[i_target] = shape_[i_this];
                if (shape_[i_this] == 1) {
                    //that's how this dimension is "extended" without creating new data
                    new_strides[i_target] = 0;
                } else {
                    new_strides[i_target] = strides_[i_this];
                }
            } else {
                new_shape[i_target] = targ_shape[i_target];
                new_strides[i_target] = 0;
            }

        }

        return this_.reindex(new_shape, new_strides);
    }

}
}

#endif /* end of include guard: HELPERS_HPP_YIOQCJSG */
