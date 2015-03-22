#include "ndinterp.hpp"
#include "ndindexer.hpp"


#include <iostream>
#include <array>

using namespace ndinterp;

using namespace std;

#define TEST(test_func_name) \
    cout<<endl<<"testing "<<#test_func_name<<" ...................."<<endl; \
    if(test_func_name()==0) { \
        cout<<"testing "<<#test_func_name<<".... success"<<endl; \
    } else { \
        cerr<<"testing "<<#test_func_name<<".... FAILED"<<endl; \
        ret = 1; \
    } \
     \

const size_t Nn = 5;

template<class KernT>
struct TestSuite {

    static
    KernT kern;

    static
    int constant_field_3D() {
        size_t Nn= 5;
 
        ndindexer<3> uind {Nn, Nn, Nn};

        const float initVal = 1;

        vector<float> u (uind.size(), initVal);

        cout<<"interp output : "<<endl;

        bool allCorrect = true;

        float start_ifrac = -1.9;

        for (float i = start_ifrac; i < Nn*3; i+=0.5)
        {
            vector<float> i_frac {
                2,
                i,
                2
            };

            array<float, 1> out = ndinterp::interpolate<KernT, 3, 1>(
                    vector<size_t>({Nn, Nn, Nn}),
                    &u[0],
                    i_frac,
                    ndinterp::OverflowBehaviour::STRETCH
                    );

            cout<<out[0]<<", ";

            if (out[0] != initVal) {
                allCorrect = false;
            }  
        }

        return (allCorrect)? 0 : 1;
    }


    static
    int constant_field_1D() {
        size_t Nn= 5;

        const float initVal = 1;

        vector<float> u (Nn, initVal);

        cout<<"interp output : "<<endl;

        float start_ifrac = -1.9;

        bool allCorrect = true;

        for (float i = start_ifrac; i < Nn*3; i+=0.5)
        {
            vector<float> i_frac {
                i
            };

            array<float, 1> out = ndinterp::interpolate<KernT, 1, 1>(
                    {Nn},
                    &u[0],
                    i_frac,
                    ndinterp::OverflowBehaviour::ZERO
                    );

            cout<<out[0]<<", ";

            if (out[0] != initVal) {
                allCorrect = false;
            }  
        }

        cout<<endl;

        //if (allsame == true) return 0;
        //else {
        //    cerr<<"interpolation doesn't yield constant result"<<endl;
        //    return 1;
        //}
        return (allCorrect)? 0 : 1;
    }

    /*
    static
    int increasing_field_3D() {
        size_t Nn= 5;

        NdShape uind = {3, {Nn, Nn, Nn}};

        vector<float> u (nd_len(uind));

        for (size_t ix = 0; ix < Nn; ++ix)
        {
            for (size_t iy = 0; iy < Nn; ++iy)
            {
                for (size_t iz = 0; iz < Nn; ++iz)
                {
                    u[nd_i(uind, ix, iy, iz)] = 
                        sqrt(float(ix)*float(ix)+ 
                             float(iy)*float(iy)+
                             float(iz)*float(iz));
                }
            }
        }

        cout<<"interp output : "<<endl;

        float start_ifrac = -1.9;

        for (float i = start_ifrac; i < Nn*3; i+=0.5)
        {
            vector<float> i_frac {
                2,
                i,
                2
            };

            array<float, 1> out = Ndinterp::interpolate<KernT, 3, 1>(
                    vector<size_t>({Nn, Nn, Nn}),
                    &u[0],
                    i_frac,
                    Ndinterp::OverflowBehaviour::ZERO
                    );

            cout<<out[0]<<", ";

            if (i == start_ifrac) {
                last_val = out[0];
            } else if (out[1] == last_val) {
                last_val = out[1];
            } else {
                //allsame = false;
            }
        }

        cout<<endl;

        //if (allsame == true) return 0;
        //else {
        //    cerr<<"interpolation doesn't yield constant result"<<endl;
        //    return 1;
        //}
        return 0;
    }
    */

    static
    int run_all_tests() {
        int ret = 0;

        TEST(TestSuite<KernT>::constant_field_1D); 
        //TEST(TestSuite<KernT>::constant_field_3D);
        //TEST(TestSuite<KernT>::increasing_field_1D);
        //TEST(TestSuite<KernT>::increasing_field_3D);
        return ret;
    }
};

int main(int argc, char *argv[])
{

    int ret = 0;

    TEST(TestSuite<ndinterp::KernNearestNeighbor>::run_all_tests);
    TEST(TestSuite<ndinterp::KernLinear>::run_all_tests);
    TEST(TestSuite<ndinterp::KernCubic>::run_all_tests);

	return ret;
}

