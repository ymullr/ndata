#ifndef VECARRAY_HPP_EPBF1RNR
#define VECARRAY_HPP_EPBF1RNR

#include <array>
#include <vector>
#include <cassert>
#include <initializer_list>
#include <exception>

#include <type_traits>

namespace ndata {


//base template
template<class T, long static_size, class Enable=void>
struct vecarray {
    static_assert(static_size>=-1, "Wrong static size!");
};


template <long T> using IsMinusOne = std::is_same<std::integral_constant<long, T>, std::integral_constant<long, -1>>;
template <typename T> using Negate = std::integral_constant<bool, !T::value>;
template <long T> using IsNotMinusOne = Negate<IsMinusOne<T>>;
template <long T> using MinusOne = std::enable_if<IsMinusOne<T>::value, std::integral_constant<long, T>>;
template <long T> using NotMinusOne = std::enable_if<IsNotMinusOne<T>::value, std::integral_constant<long, T>>;

//special value to pass as static_size template argument when a dynamicv size is wanted
constexpr long DYNAMICALLY_SIZED = -1;

//special value to pass as dynamic_size parameter to vecarray constructor when a dynamic vecarray is wanted
constexpr long STATICALLY_SIZED = -2;

//specialization static vecarray
template<class T, long static_size>
struct vecarray<T, static_size, typename std::enable_if<(static_size >= 0)>::type> {

    static constexpr long STATIC_SIZE_OR_DYNAMIC = static_size;

    typedef T type_T;

    std::array<T, static_size> stackStorage;

    explicit
    vecarray(long dynamic_size, T init_val = T()) {
        assert(dynamic_size == STATICALLY_SIZED);

        //data members with content of the vector
        //imitates the behaviour of dataializer lists
        for (size_t i = 0; i < static_size; ++i) {
            stackStorage[i] = init_val;
        }
    }

    vecarray(std::array<T, static_size> iniArr):
        stackStorage(iniArr) {
    }

    vecarray(std::vector<T> ini) {
        //static_assert(static_size == 0, "This constructor always produce a dynamic vecarray");
        assert(ini.size() == static_size);
        for (size_t i = 0; i < ini.size(); ++i) {
            stackStorage[i]=ini[i];
        }
    }

    //overloaded to avoid ambiguity btw vector and array
    vecarray(std::initializer_list<T> iniArr)
    {
        assert(iniArr.size() == static_size);// "Size of dataializer list doesn't match");
        //static_assert(iniArr.size() == static_size, "");// "Size of dataializer list doesn't match");

        size_t i = 0;
        for (T val: iniArr) {
            stackStorage[i] = val;
            i++;
        }
    }

    //empty dataializer for later assignment
    vecarray() {};

    size_t size() {
        return static_size;
    }

    /**
     * For compatibility with dynarray.
     */
    constexpr long dynsize () { return STATICALLY_SIZED;}

    T& operator[](size_t index) {
        assert(index<size());
        return stackStorage[index];
    }

    void fill(T val) {
        for (size_t i = 0; i < size(); ++i) {
            this->operator[](i) = val;
        }
    }

    vecarray<T, static_size>
    reverse () {
        vecarray<T, static_size> new_vecarray;

        for (size_t i = 0; i < static_size; ++i) {
            long rev_ind = static_size-1-i;
            assert(rev_ind > 0);
            new_vecarray[rev_ind] = stackStorage[i];
        }

        return new_vecarray;
    }


    /**
     * Like push_back but returns new vecarray instead of mutating 
     */
    vecarray<T, static_size+1>
    append (T val) {
        vecarray<T, static_size+1> new_vecarray;

        for (size_t i = 0; i < size(); ++i) {
            new_vecarray[i] = operator[](i);
        }

        new_vecarray[static_size] = val;

        return new_vecarray;
    }

    /**
     * returns all elements of the vecarray but the last
     */
    vecarray<T, static_size-1>
    drop_back () {
        static_assert(static_size!=0, "woops already empty");

        vecarray<T, static_size-1> new_vecarray;

        for (size_t i = 0; i < new_vecarray.size(); ++i) {
            new_vecarray[i] = operator[](i);
        }

        return new_vecarray;
    }

    /**
     * returns all elements of the vecarray but the first
     */
    vecarray<T, static_size-1>
    drop_front () {
        static_assert(static_size!=0, "woops already empty");

        vecarray<T, static_size-1> new_vecarray;

        for (size_t i = 0; i < new_vecarray.size(); ++i) {
            new_vecarray[i] = operator[](i+1);
        }

        return new_vecarray;
    }


    T& back() {
        return operator[](size()-1);
    }

};

/**
 * Specialization for a dynamic vecarray
 */
template<class T, long static_size>
struct vecarray<T, static_size, typename std::enable_if<(static_size == DYNAMICALLY_SIZED)>::type > {

    static constexpr long STATIC_SIZE_OR_DYNAMIC = DYNAMICALLY_SIZED;

    std::vector<T> heapStorage;

    vecarray(long dynamic_size) {
        heapStorage = std::vector<T>(dynamic_size);
        assert(dynamic_size>=0);
    }

    vecarray(long dynamic_size, T init_val = T()) {

        heapStorage = std::vector<T>(dynamic_size);

        assert(dynamic_size>=0);

        //data members with content of the vector
        //imitates the behaviour of dataializer lists
        for (size_t i = 0; i < dynamic_size; ++i) {
            heapStorage[i] = init_val;
        }
    }

    vecarray(std::vector<T> ini) {
        heapStorage=ini;
    }

    vecarray(std::array<T, static_size> iniArr): heapStorage(iniArr.begin(), iniArr.end()) { }

    //overloaded to avoid ambiguity btw vector and array
    vecarray(std::initializer_list<T> iniArr) {
        size_t i = 0;
        for (T val: iniArr) {
            heapStorage[i] = val;
            i++;
        }
    }

    //empty dataializer for later assignment
    vecarray() {};

    size_t size() {
        return heapStorage.size();
    }

    long dynsize () { return heapStorage.size();}

    T& operator[](size_t index) {
        assert(index<size());
        return heapStorage[index];
    }

    void fill(T val) {
        for (size_t i = 0; i < size(); ++i) {
            this[i] = val;
        }
    }

    vecarray<T, DYNAMICALLY_SIZED>
    reverse () {
        vecarray<T, 0> new_vecarray (size());
        size_t vecsize = new_vecarray.size();

        for (size_t i = 0; i < vecsize; ++i) {
            long rev_ind = vecsize-1-i;
            assert(rev_ind > 0);
            new_vecarray[rev_ind] = this[i];
        }

        return new_vecarray;
    }


    /**
     * Like push_back but returns new vecarray instead of mutating 
     *
     * Dynamic version
     */
    vecarray<T, DYNAMICALLY_SIZED>
    append (T val) {
        vecarray<T, DYNAMICALLY_SIZED> new_vecarray (size()+1);

        for (size_t i = 0; i < size(); ++i) {
            new_vecarray[i] = this[i];
        }

        new_vecarray[size()] = val;

        return new_vecarray;
    }

    /**
     * returns all elements of the vecarray but the last
     */
    vecarray<T, DYNAMICALLY_SIZED>
    drop_back () {
        if (size()==0){
            throw(std::underflow_error("this.size() is already 0"));
        }

        vecarray<T, DYNAMICALLY_SIZED> new_vecarray (dynsize()-1);

        for (size_t i = 0; i < new_vecarray.size(); ++i) {
            new_vecarray[i] = operator[](i);
        }

        return new_vecarray;
    }

    /**
     * returns all elements of the vecarray but the first
     */
    vecarray<T, DYNAMICALLY_SIZED>
    drop_front () {
        if(size()==0){
            throw(std::underflow_error("this.size() is already 0"));
        }

        vecarray<T, DYNAMICALLY_SIZED> new_vecarray (dynsize()-1);

        for (size_t i = 0; i < new_vecarray.size(); ++i) {
            new_vecarray[i] = operator[](i+1);
        }

        return new_vecarray;
    }

    T back() {
        return operator[](size()-1);
    }

};

namespace helpers {

    template<typename ContentT, long idim>
    vecarray<ContentT, idim>
    array_from_argpack(vecarray<ContentT, idim> acc) {
        return acc;
    }


    template<typename ContentT, long idim, typename... ContentTPackT>
    auto
    array_from_argpack(
            vecarray<ContentT, idim> acc,
            ContentT i,
            ContentTPackT... rest
            ) {
        vecarray<ContentT, idim+1> new_acc = acc.append(i);
        return array_from_argpack<ContentT, idim+1>(new_acc, rest...);
    }
}


template <typename T, typename ... Ts>
auto
make_vecarray(T val1, Ts ... valn) {
    constexpr long static_size = sizeof...(valn)+1;
    vecarray<T, static_size> ret = helpers::array_from_argpack(
                vecarray<T, 0>(),
                val1, valn...
                );
    return ret;
}

} //end namespace


#include "ndata/vecarray_helpers.hpp"

#endif /* end of include guard: VECARRAY_HPP_EPBF1RNR */
