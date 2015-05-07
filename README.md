# Ndata
##### Ndata is a C++ template library for working with n-dimensional arrays of arbitrary data.

It is aimed at being lightweight, fast, and easy to integrate into your own projects. It offers an interface that will be familiar to users of Python/numpy or Matlab for indexing and slicing multidimensional arrays with any number of dimensions.

Here is a short feature list:
* Header only.
* Low buy-in cost. You can use the indexing logic with your own containers and data type, and do not have to design around the provided data structures.
* Efficient memory contiguous storage, compatible with most popular numerical libraries.
* Bound checking when compiled in debug mode. Makes debugging a lot easier.
* Compile time checking of arrays dimensionality. Makes debugging a lot less necessary.
* N-dimensional loop constructs for iterating over multidimensional arrays. Also available in multithreaded flavor with OpenMP.
* Array shape broadcasting, in the same way that numpy does it. More details [here](http://wiki.scipy.org/EricsBroadcastingDoc).

Not a feature : linear algebra. If you are looking for a good linear algebra library, have a look at [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). 

If, however, you would need to make computations with n-dimensional arrays of matrices and vectors, then you can use Eigen AND Ndata, as they play very well together.

This library should be considered under development, the API may change and the documentation may be lacking. It is however quite useable already.

## How to use

1. Clone this reposiory with `git clone https://github.com/ymullr/ndata`
2. Add the new ndata directory the the include paths of your project. All the code is header only, except the tests.
3. Make sure you have the option -std=c++14 in your compiler flags. If you are using GCC, version >= 4.9 is required.
4. Start working.

If you want to create a new project and use cmake as your build tool, the CMakeLists.cmake file provided to compile the 
tests may be a good starting point.

## Exemples

Import the library and namespace

~~~
#include "ndata.hpp"
using namespace std;
using namespace ndata;
~~~

The basic multidimensional container is called an nvector. it can be built from an indexer object and some data.

Let's create a 2D nvector of long integers initialized at 0 called u1. This one has a shape equal to (5, 2) :

~~~
//The size of our arrays on x and y axis
size_t Nx = 5, Ny = 2;
nvector<long, 2> u1 (indexer<2>(Nx, Ny), 0l);
~~~

This syntax requires a bit of repetition of the template parameters, so some helper functions for creating nvectors are provided, in the same spirit as std::make_tuple. Here is an alternative way to declare and initialize u1 by using those functions :

~~~
auto u1 = make_nvector<long>(make_indexer(Nx, Ny), 0l);
~~~

Now let's create another 2D nvector of std::pair<long, long> and set the value of the elements of both arrays by using nested for loops :

~~~
auto u2 = make_nvector<std::pair<long, long>>(make_indexer(Nx, Ny));

//two dimensions, two loops
for (size_t ix = 0; ix < Nx ; ++ix) {
    for (size_t iy = 0; iy < Ny ; ++iy) {
        //to reference an element of an ndimensional array,
        //we use the paren operator "()", or function call syntax,
        //and pass ix the x-axis index and iy the y-axis index
        u1(ix ,iy) = ix*Ny + iy;
        u2(ix, iy) = make_pair(1,1);
    }
}
~~~

Using nested for loops is OK, but it gets quite cumbersome when the number of dimensions increases. To alleviate this problem, we can use the nforeach function.

The nforeach function takes a tuple of nvector references and a function operating on the elements of the nvectors. The function is then called once for each set of matching elements in the nvectors. 

~~~
nforeach(
    std::tie(u1, u2),
    [] (long & s1, pair<long, long> s2) { //C++11 lambda function syntax
        s1=s1+s2.first;
    }
);
~~~

After this loop, all the values in u1 have increased by one. Note that the lambda function takes its first argument by reference, otherwise it wouldn't have been able to update the original values in u1.

You can look at the tests for more exemples.

Other interesting functionalities are : array slicing, reshaping, elementwise assignment, ntransform loops and broadcasting.

## Numerical algorithms

Some implementations of numerical algorithms can be found in the ndata/numerics/ directory. At the moment there are only functions to perform n-dimensional interpolations of various kinds : nearest neighbor, n-linear, and cubic. 

## API Reference

Doxygen generated API documentation [may be found here](http://ymullr.github.io/ndata/doc/index.html)
