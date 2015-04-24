#include "debug_helpers.hpp"

#include "ndata.hpp"

//TODO test dynamic ndarrays

using namespace std;
using namespace ndata;


int main(int argc, char *argv[])
{

        auto u1 = make_nvector<long>(make_indexer(2, 5));
        auto u2 = make_nvector<long>(make_indexer(2, 1));
        auto u3 = make_nvector<long>(make_indexer(5));

        auto tup = helpers::broadcast(make_tuple(u1, u2, u3));
        auto u1_broad = get<0>(tup);
        auto u2_broad = get<1>(tup);
        auto u3_broad = get<2>(tup);
        bool succ_u1u2 = true;
        bool succ_u1u3 = true;
        for (size_t i = 0; i < u1_broad.get_shape().size(); ++i) {
            if (u1_broad.get_shape()[i] != u2_broad.get_shape()[i]) {
                succ_u1u2 = false;
            }
            if (u1_broad.get_shape()[i] != u3_broad.get_shape()[i]) {
                succ_u1u3 = false;
            }
        }


        return 0;
}
