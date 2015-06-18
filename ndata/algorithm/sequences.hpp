#ifndef SEQUENCES_HPP_YKEBA2N3
#define SEQUENCES_HPP_YKEBA2N3

namespace ndata {


template<typename NumT>
nvector<NumT, 1> num_range(NumT start, NumT stop, NumT step=1) {

    auto comp = [step] (NumT v, NumT boundary) {
        if (step < 0) {
            return v > boundary;
        }
        return v < boundary;
    };

    std::vector<NumT> ret_vec (0);

    for (int i = start; comp(i,stop); i+=step) {
        ret_vec.push_back(i);
    }
    return make_nvector(ret_vec);
}

template<typename NumT>
nvector<NumT, 1> linspace(NumT start, NumT stop, size_t size) {

    assert(size > 1);
    NumT step = (stop-start)/(size-1);

    std::vector<NumT> ret_vec (0);

    for (size_t i = 0; i < size ; ++i) {
        ret_vec.push_back(start+long(i)*step);
    }
    return make_nvector(ret_vec);
}


}//end namespace


#endif /* end of include guard: SEQUENCES_HPP_YKEBA2N3 */
 
