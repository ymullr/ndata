#include <array>
#include <vector>
#include <initializer_list>

#include <type_traits>

using namespace std;

#ifndef VECARRAY_HPP_EPBF1RNR
#define VECARRAY_HPP_EPBF1RNR


template <long T> using IsMinusOne = is_same<integral_constant<long, T>, integral_constant<long, -1>>;
template <typename T> using Negate = integral_constant<bool, !T::value>;
template <long T> using IsNotMinusOne = Negate<IsMinusOne<T>>;
template <long T> using MinusOne = std::enable_if<IsMinusOne<T>::value, long>;
template <long T> using NotMinusOne = enable_if<IsNotMinusOne<T>::value, long>;

template<class T, long static_size, class Enable = void>
struct vecarray {};

//specialization static vecarray
template<class T, long static_size>
struct vecarray<T, static_size, typename NotMinusOne<static_size>::type> {

    array<T, static_size> stackStorage ;

    //provided for compatibility with dynamic version
    vecarray(long dynamic_size) {
        assert(dynamic_size == 0);
    };

    vecarray(long dynamic_size, vector<T> init) {
        assert(dynamic_size == 0);
        assert(init.size() <= static_size);

        //init members with content of the vector
        //imitates the behaviour of initializer lists
        for (size_t i = 0; i < init.size(); ++i) {
            stackStorage[i] = init[i];
        }
    };

    vecarray(array<T, static_size> iniArr):
        stackStorage(iniArr) {
    };

    vecarray(vector<T> ini) {
        //static_assert(static_size == 0, "This constructor always produce a dynamic vecarray");
        assert(ini.size() <= static_size);
        for (size_t i = 0; i < ini.size(); ++i) {
            stackStorage[i]=ini[i];
        }
    };

    //overloaded to avoid ambiguity btw vector and array
    vecarray(initializer_list<T> iniArr)
    {
        static_assert(iniArr.size() == static_size, "Size of initializer list doesn't match");

        size_t i = 0;
        for (T val: iniArr) {
            stackStorage[i] = val;
            i++;
        }
    };

    //empty initializer for later assignment
    vecarray() {};

    size_t size() {
        return static_size;
    };

    /**
     * For compatibility with dynarray.
     */
    constexpr long dynsize () { return -1;}

    T& operator[](size_t index) {
        return stackStorage[index];
    }

    void fill(T val) {
        for (size_t i = 0; i < size(); ++i) {
            this[i] = val;
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
     * Like push_front but returns new vecarray instead of mutating 
     *
     * Static version
     *
     */
    vecarray<T, static_size+1>
    cons (T val) {
        vecarray<T, static_size+1> new_vecarray;
        new_vecarray[0] = val;

        for (size_t i = 0; i < static_size; ++i) {
            new_vecarray[i+1] = this[i];
        }

        return new_vecarray;
    };

};

/**
 * Specialization for a dynamic vecarray
 */
template<class T, long static_size>
struct vecarray<T, static_size, typename MinusOne<static_size>::type> {

    vector<T> heapStorage;

    vecarray(long dynamic_size) {
        heapStorage = vector<T>(dynamic_size); 
        assert(dynamic_size>=0);
    };

    vecarray(long dynamic_size, vector<T> init) {

        heapStorage = vector<T>(dynamic_size);
        size=dynamic_size;

        assert(init.size() <= size);

        assert(dynamic_size>=0);

        //init members with content of the vector
        //imitates the behaviour of initializer lists
        for (size_t i = 0; i < init.size(); ++i) {
            heapStorage[i] = init[i];
        }
    };

    vecarray(vector<T> ini) {
        heapStorage=ini;
    };

    vecarray(array<T, static_size> iniArr): heapStorage(iniArr.begin(), iniArr.end()) { }

    //overloaded to avoid ambiguity btw vector and array
    vecarray(initializer_list<T> iniArr) {
        size_t i = 0;
        for (T val: iniArr) {
            heapStorage[i] = val;
            i++;
        }
    };

    //empty initializer for later assignment
    vecarray() {};

    size_t size() {
        return heapStorage.size();
    };

    long dynsize () { return heapStorage.size();}

    T& operator[](size_t index) {
        return heapStorage[index];
    }

    void fill(T val) {
        for (size_t i = 0; i < size(); ++i) {
            this[i] = val;
        }
    }

    vecarray<T, static_size>
    reverse () {
        vecarray<T, 0> new_vecarray (this.size());
        size_t vecsize = new_vecarray.size();

        for (size_t i = 0; i < vecsize; ++i) {
            long rev_ind = vecsize-1-i;
            assert(rev_ind > 0);
            new_vecarray[rev_ind] = this[i];
        }

        return new_vecarray;
    }


    /**
     * Like push_front but returns new vecarray instead of mutating 
     *
     * Dynamic version
     */
    vecarray<T, -1>
    cons (T val) {
        vecarray<T, -1> new_vecarray (this.size()+1);
        new_vecarray[0] = val;

        for (size_t i = 0; i < new_vecarray.size(); ++i) {
            new_vecarray[i+1] = this[i];
        }

        return new_vecarray;
    } 


};


#endif /* end of include guard: VECARRAY_HPP_EPBF1RNR */
