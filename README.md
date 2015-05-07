# Ndata
##### Ndata is a C++ template library for working with n-dimensional arrays of arbitrary data.

It is aimed at being lightweight, fast, and easy to integrate into your own projects.

Here is a short feature list:
* Header only.
* Low buy-in cost. You can use the indexing logic with your own containers and data type, and do not have to design around the provided data structures.
* Efficient memory contiguous storage, compatible with most popular numerical libraries.
* An interface that will be familiar to users of Python/numpy and Matlab for indexing and slicing multidimensional arrays with an arbitrary number of dimensions.
* Bound checking when compiled in debug mode. Makes debugging a lot easier.
* Compile time checking of arrays dimensionality. Makes debugging a lot less necessary.
* N-dimensional loop constructs for iterating over multidimensional arrays. Also available in multithreaded flavor with OpenMP.
* Array shape broadcasting, in the same way numpy does. More details [here](http://wiki.scipy.org/EricsBroadcastingDoc).

Not a feature : linear algebra. If you are looking for a good linear algebra library, have a look at [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). 

If, however, you need to make computations with n-dimensional arrays of matrices and vectors, then you can use Eigen AND Ndata, as they play very well together.

## Exemples



## How to use

1. Clone this reposiory with "git clone https://github.com/ymullr/ndata".
2. Add the new ndata directory the the include paths of your project. All the code is header only, except the tests.
3. Make sure you have the option -std=C++14 in your compiler flags. If you are using GCC, version >= 4.9 is required.
4. Start working.

## API Reference

Doxygen generated API documentation [may be found here](http://ymullr.github.io/ndata/doc/index.html)
