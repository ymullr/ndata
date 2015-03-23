#include <array>
#include <vector>
#include <initializer_list>

using namespace std;

#ifndef VECARRAY_HPP_EPBF1RNR
#define VECARRAY_HPP_EPBF1RNR




template<class T, size_t static_size>
struct vecarray {

    array<T, static_size> stackStorage ;

    vector<T> heapStorage;

    vecarray(size_t dynamic_size) {
        if (static_size != 0) {
            assert(dynamic_size == 0);
        } else {
           heapStorage = vector<T>(dynamic_size); 
        }
    };

    vecarray(size_t dynamic_size, vector<T> init) {

        T* dest = &stackStorage[0];
        size_t size = static_size;

        if (static_size != 0) {
            assert(dynamic_size == 0);
        } else {
            heapStorage = vector<T>(dynamic_size);
            dest = &heapStorage[0];
            size=dynamic_size;
        }

        assert(init.size() <= size);

        //init members with content of the vector
        //imitates the behaviour of initializer lists
        for (size_t i = 0; i < init.size(); ++i) {
            dest[i] = init[i];
            assert(i<size);
        }
    };

    vecarray(vector<T> ini) {
        //static_assert(static_size == 0, "This constructor always produce a dynamic vecarray");
        assert(ini.size() == static_size or static_size == 0);
        if (static_size == 0) {
            heapStorage=ini;
        } else {
            for (size_t i = 0; i < ini.size(); ++i) {
                stackStorage[i]=ini[i];
            }
        }
    };

    vecarray(array<T, static_size> iniArr):
        stackStorage(iniArr) {
    };

    //overloaded to avoid ambiguity btw vector and array
    vecarray(initializer_list<T> iniArr)
    {
        size_t i = 0;
        for (T val: iniArr) {
            stackStorage[i] = val;
            i++;
        }

        static_assert(iniArr.size() == static_size, "Size of initializer list doesn't match");
    };

    //empty initializer for later assignment
    vecarray() {};

    size_t size() {
        if (static_size == 0) {
            return heapStorage.size();
        } else {
            return static_size;
        }
        abort();
        return 0;
    };

    T& operator[](size_t index) {
        if (static_size == 0) {
            return heapStorage[index];
        } else {
            return stackStorage[index];
        }
        abort();
        return stackStorage[index];
    }

    void fill(T val) {
        for (size_t i = 0; i < size(); ++i) {
            this[i] = val;
        }
    }

};

#endif /* end of include guard: VECARRAY_HPP_EPBF1RNR */
