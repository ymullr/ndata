#include "ndata.hpp"
#include "ndata/algorithm/interp.hpp"
#include "ndata/algorithm/sequences.hpp"
#include <ndata/debug_helpers.hpp>

//#include <iostream>
#include <fstream>

#include <string>


using namespace ndata;
using namespace ndata::interp;
using namespace std;


/// \brief 1D case
template <typename ContainerT, typename T>
string nvector2csv(ndatacontainer<ContainerT, T, 1> u) {
    string ret ("");

    for (long i = 0; i < u.get_shape()[0]; ++i) {
        ret.append( MakeString() << u(i) << ", ");
    }
    return ret;
}


/// \brief 2D case
template <typename ContainerT, typename T>
string nvector2csv(ndatacontainer<ContainerT, T, 2> u) {
    string ret ("");
    for (long i = 0; i < u.get_shape()[0]; ++i) {
        ret.append(
            MakeString() 
                << nvector2csv(u.slice(i, range())) //recursive call
                << "\n");
    }
    return ret;
}



int main(int /*argc*/, char** /*argv */)
{
    float x_start = -5.f;
    nvector<float, 1> x = numrange(x_start, 6.f, 1.f);
    nvector<float, 1> y = ntransform<float>(
        make_tuple(x),
        [](float xs) {
            if(xs >= 0) {
                return 1.f;
            } else {
                return 0.f;
            }
        });


    auto index_frac = ndata::numrange(0.f, 10.f, 0.1f);
    nvector<float, 1> y_interp = 
        ntransform<float>(
                make_tuple(index_frac),
                [y] (float index_frac_s) {
                    return interpolate<kern_lanczos<2>, overflow_behaviour::zero>(y, make_vecarray(index_frac_s));
                });

    auto x_interp = ntransform<float>(
            make_tuple(index_frac),
            [x_start](float index_frac_s) {
                return index_frac_s-x_start;
            });

    //concatenation of x and y
    nvector<float, 2> xy (make_indexer(2, x_interp.get_shape()[0]));
    xy.slice(0, range()).assign(x);
    xy.slice(1, range()).assign(y);

    nvector<float, 2> xy_interp (make_indexer(2, x_interp.get_shape()[0]));
    xy_interp.slice(0, range()).assign(x_interp);
    xy_interp.slice(1, range()).assign(y_interp);

    string str1 = nvector2csv(xy);
    string str2 = nvector2csv(xy_interp);

    ofstream fs ("test_out.csv");
    fs << str1 << "\n" << str2;
    fs.close();

	return 0;
}
