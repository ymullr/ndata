cmake_minimum_required(VERSION 2.6)

#Declaring projet
project(turbgen)

set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_VERBOSE_MAKEFILE ON)


#including script for handling find_package on FFTW3 OpenMP and TINYCC
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules/")

find_package(OpenMP) #not available with clang

#comment one
set(GCC_FLAGS "-std=c++14")#libc++ must be installed (clang standard lib)
set(CMAKE_CXX_COMPILER "/usr/bin/g++-4.9")
#or
#set(CMAKE_CXX_COMPILER "/usr/bin/clang++-3.5")
#set(CLANG_FLAGS "-std=c++14 -stdlib=libc++")#libc++ must be installed (clang standard lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_FLAGS} ${CLANG_FLAGS} ${OpenMP_CXX_FLAGS}")

#sometimes useful: -Wsign-conversion
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -ggdb3 -O0 -fno-inline -DDEBUG  -Wall -Wextra  -Wno-missing-braces ")#-ggdb3 plutot que -g genere et integre les info de debug pour les macros
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")# -ffast-math

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")


include_directories(
	./
    ./c89/includes
	)

find_library(M_LIB m)
find_library(DL_LIB dl)

#add_executable(
#    narray_test
#    ./tests/narray_test.cpp
#	)
#
add_executable(
        indexer_test
        ./tests/indexer_test.cpp
        )

add_executable(
        interp_test
        ./tests/interp_test.cpp
        )

add_executable(
        vecarray_test
        ./tests/vecarray_test.cpp
        )

#add_executable(
#    debug_test
#    ./tests/debug_test.cpp
#)


#fake target
FILE(GLOB_RECURSE INC_ALL "./*.hpp" "./tests/*.hpp")
add_custom_target(
    make_qtcreator_happy_and_able_to_find_my_headers
    SOURCES ${INC_ALL}
)
