#include "interp.hpp"
#include "indexer.hpp"

#include <iostream>
#include <string>
#include <array>
#include <tuple>

#include <debug_helpers.hpp>
#include "interp.hpp"

using namespace ndata;


float linerp(float a, float b, float x) {
    assert(x>=0 and x<=1);
    return a*(1-x) + b*x;
}

const size_t Nn = 6;

template<class KernT>
struct TestSuite {

    static
    KernT kern;

    static
    TestResult constant_field_3D1C() {
 
        indexer<3> uind (Nn, Nn, Nn);
        const float initVal = 1;
        vector<float> u (uind.size(), initVal);
        float start_ifrac = -1.9;

        DECLARE_TESTRESULT(allCorrect, output);

        for (float i = start_ifrac; i < Nn*3; i+=0.5) {
            vector<float> i_frac {
                2,
                i,
                2
            };

            float out = interpolate<KernT, 3>(
                    {Nn, Nn, Nn},
                    &u[0],
                    i_frac,
                    OverflowBehaviour::STRETCH
                    );

            output.append(MakeString()<<out<<", ");

            if (fabs(out - initVal) > initVal*1e-5) {
                allCorrect = false;
            }  
        }

        RETURN_TESTRESULT(allCorrect, output)
    }


    static
    TestResult constant_field_1D1C() {

        const float initVal = 1;
        vector<float> u (Nn, initVal);
        float start_ifrac = -1.9;

        DECLARE_TESTRESULT(allCorrect, output);

        for (float i = start_ifrac; i < Nn*3; i+=0.5) {

            float i_frac = i; 

            float out = interpolate<KernT>(
                    Nn,
                    &u[0],
                    i_frac,
                    OverflowBehaviour::STRETCH
                    );

            output.append(MakeString()<< out<<", ");

            if (fabs(out - initVal) > initVal*1e-5) {
                allCorrect = false;
            }  
        }

        RETURN_TESTRESULT(allCorrect, output);
    }

    static
    TestResult zero_boundary_1D1C() {

        const float initVal = 1;
        vector<float> u (Nn, initVal);
        float start_ifrac = -1.9;

        DECLARE_TESTRESULT(success_bool, output);

        for (float i = start_ifrac; i < Nn*3; i+=0.5) {

            float i_frac = i;  

            float out = interpolate<KernT>(
                    Nn,
                    &u[0],
                    i_frac,
                    OverflowBehaviour::ZERO
                    );

            if(
                    (i == 0      and out != 0)
                 or (i == Nn*3-1 and out != 0)
            ) {
                success_bool = false;
            } 
            output.append(MakeString()<< out<<", ");
        }

        RETURN_TESTRESULT(success_bool, output);
    }

    static
    TestResult increasing_field_1D1C(OverflowBehaviour ob=OverflowBehaviour::STRETCH) {
  
        vector<float> u (Nn);
        float increment = 1;
        for (size_t i = 0; i < Nn; ++i)
        {
            u[i] = i*increment;
        }

        float start_ifrac = 1.1;

        vector<float> output (0);

        DECLARE_TESTRESULT(increasing, retMsg);

        float previousVal;

        for (float i = start_ifrac; i < Nn-1; i+=0.25) {

            float i_frac = i ;

            float out = interpolate<KernT>(
                    Nn,
                    &u[0],
                    i_frac,
                    ob
                    //OverflowBehaviourT
                    );

            retMsg.append(MakeString() << out << ", "); 
            if (i > start_ifrac) {
                increasing = increasing and (out>=previousVal);
            }

            previousVal = out;
        }

        RETURN_TESTRESULT(increasing, retMsg);
    }

    static
    TestResult increasing_field_3D1C() {
  
        indexer<3> uind (Nn, Nn, Nn);

        vector<float> u (uind.size());

        float increment = 1;
        for (size_t i = 0; i < u.size(); ++i) {
            u[i] = i*increment;
        }

        float start_ifrac = -1.9;

        vector<float> output (0);

        DECLARE_TESTRESULT(increasing, retMsg);

        float previousVal;

        for (float i = start_ifrac; i < Nn*3; i+=0.5) {

            float i_frac = i ;

            float out = interpolate<KernT>(
                    Nn,
                    &u[0],
                    i_frac,
                    OverflowBehaviour::STRETCH
                    );

            retMsg.append(MakeString() << out << ", "); 

            if (i > start_ifrac) {
                increasing = increasing and (out>=previousVal);
            }

            previousVal = out;
        }

        RETURN_TESTRESULT(increasing, retMsg);
    }


    static
    TestResult simple_equalities_1D1C() {
  
        size_t Nn = 8;

        vector<float> u (Nn);
        float valEven = 1, valOdd = 2;

        for (size_t i = 0; i < Nn; ++i)
        {
            u[i] = (i%2==0)? valEven : valOdd;
        }

        vector<float> output (0);

        DECLARE_TESTRESULT(bool_equal, retMsg);

        for (size_t i = 2; i < Nn-2; ++i) {

            float i_frac = i; 

            float val = interpolate<KernT>(
                    Nn,
                    &u[0],
                    i_frac,
                    OverflowBehaviour::STRETCH
                    );


            float valhalf = interpolate<KernT>(
                    Nn,
                    &u[0],
                    i_frac+0.5,
                    OverflowBehaviour::STRETCH
                    );

            retMsg.append(MakeString() << val << ", " << valhalf << ", "); 

            if (i%2 == 0) { //we are on an even node
                if (val != valEven) {
                    bool_equal = false;
                }
            } else { //odd node
                if (val != valOdd) {
                    bool_equal = false;
                }
            }

            if (valhalf != (valEven + valOdd) / 2.f) bool_equal=false;
        }

        RETURN_TESTRESULT(bool_equal, retMsg);
    }

    static
    TestResult cyclic_equal_3D3C() {
  
        size_t Nn = 10;

        indexer<4> uind (Nn, Nn, Nn, 3) ;

        vector<float> u (uind.size());

        vector<size_t> ndind (uind.get_shape().size(), 0);

        for (size_t i = 0; i < u.size(); ++i) {
            u[i] = sqrt(
                          ndind[0]*ndind[0] 
                        + ndind[1]*ndind[1] 
                        + ndind[2]*ndind[2] 
                );
            uind.increment_ndindex(ndind);
        }

        indexer<2> valuesInd (Nn*4, 3);

        vector<float> values (valuesInd.size());

        DECLARE_TESTRESULT(aggreg_equal, retMsg);

        bool cycle_equal = true, cpt_equal = true;

        for (size_t itraversal = 0; itraversal < 2; ++itraversal) {
            for (size_t i = 0; i < valuesInd.shape(0); ++i) {

                float xfrac =
                        itraversal * (Nn)
                        + i/float(valuesInd.shape(0)) * (Nn-1.f);

                array<float, 3> val = interpolate<KernT, 3, 3>(
                        {Nn, Nn, Nn},
                        &u[0],
                        {xfrac, 4.5, 4.5},
                        {
                            OverflowBehaviour::CYCLIC,
                            OverflowBehaviour::CYCLIC,
                            OverflowBehaviour::STRETCH
                        });

                if (itraversal != 0) {
                    for (size_t ic = 0; ic < 3; ++ic)
                    {
                        float diff = val[ic]-values[valuesInd.index(i, ic)];
                        if(
                            fabs(diff) > 1e-5f
                          ) {
                            cycle_equal = false;
                        };
                    }
                } 

                if (
                        val[1] != val[0] 
                     or val[2] != val[0]
                   ) {
                    cpt_equal = false;
                }

                for (size_t ic = 0; ic < 3; ++ic) {
                    values[valuesInd.index(i, ic)] = val[ic];
                }
            }
        }

        for (size_t ic = 0; ic < 3; ++ic) {
            for (size_t iv = 0; iv < valuesInd.shape(0); ++iv) {
                retMsg.append(MakeString() << values[valuesInd.index(iv, ic)]
                        << ", "); 
            }

            retMsg.append("\n");
        }

        aggreg_equal = cpt_equal and cycle_equal;
        RETURN_TESTRESULT(aggreg_equal, retMsg);
    }

    template<size_t ndims, size_t ncpt>
    static
    TestResult stable_derivative_NDNC() {
  
        size_t Nn = 10;

        vector<size_t> vecshape (ndims, Nn);
        vector<size_t> fullshape (ndims, Nn);
        fullshape.push_back(ncpt);

        vector<OverflowBehaviour> ovfl (ndims-1, OverflowBehaviour::CYCLIC);
        ovfl.push_back(OverflowBehaviour::STRETCH);

        indexer<ndims+1> uind (fullshape) ;

        size_t u_size = uind.size();
        vector<float> u (u_size);

        vector<size_t> ndind (uind.get_shape().size(), 0);

        for (size_t i = 0; i < u.size(); ++i) {
            //u[i] = sqrt(
            //              ndind[0]*ndind[0] 
            //            + ndind[1]*ndind[1] 
            //            + ndind[2]*ndind[2] 
            //    );
            u[i] = ndind[0];
            uind.increment_ndindex(ndind);
        }

        indexer<2> valuesInd (Nn*4, ncpt);

        vector<float> values (0);

        DECLARE_TESTRESULT(aggreg_equal, retMsg);

        bool increasing = true, stable_derivative = true;

        for (size_t i = 0; i < valuesInd.shape(0); ++i) {

            float xfrac = i/float(valuesInd.shape(0)) * (Nn-1.f);

            vector<float> indx (ndims, 3.5);
            indx[0] = xfrac;

            array<float, ncpt> val = interpolate<KernT, ndims, ncpt>(
                    vecshape,
                    &u[0],
                    indx,
                    ovfl
                    );

            for (size_t ic = 0; ic < ncpt; ++ic) {
                values.push_back(val[ic]);
            }

            //test on derivatives
            //only test away from cyclic boundaries
            if(     xfrac >= KernT::ONESIDEDWIDTH+2
                    and xfrac <= (Nn-1)-KernT::ONESIDEDWIDTH
              ) {

                for (size_t ic = 0; ic < ncpt; ++ic) {

                    float d1 =
                        values[valuesInd.index(i, ic)]
                        - values[valuesInd.index(i-1, ic)],
                        d2 =
                            values[valuesInd.index(i-1, ic)]
                            - values[valuesInd.index(i-2, ic)];

                    if(d1 < 0 or d2 < 0) {
                        increasing = false;
                    }


                    if (fabs(d2-d1) > 1e-4*d1) {
                        stable_derivative =false;
                    }
                }
            }

        }

        for (size_t iv = 0; iv < valuesInd.shape(0); ++iv) {
            for (size_t ic = 0; ic < ncpt; ++ic) {
                retMsg.append(MakeString() << values[valuesInd.index(iv, ic)] << ", "); 
            }

            retMsg.append("\n");
        }

        aggreg_equal = increasing and stable_derivative;
        RETURN_TESTRESULT(aggreg_equal, retMsg);
    }

    static
    TestResult run_all_tests() {

        DECLARE_TESTRESULT(success_bool, msg);

        //TEST(constant_field_1D1C()            , success_bool, msg); 
        //TEST(zero_boundary_1D1C()             , success_bool, msg);
        //TEST(increasing_field_1D1C()          , success_bool, msg);
        //TEST(increasing_field_1D1C(OverflowBehaviour::CYCLIC)
        //        , success_bool, msg);
        //TEST(simple_equalities_1D1C()         , success_bool, msg);
        //TEST(constant_field_3D1C()            , success_bool, msg);
        //TEST(increasing_field_3D1C()          , success_bool, msg);
        //TEST(cyclic_equal_3D3C()              , success_bool, msg);

        //yeah.. that macro doesn't like multiple template arguments... 
        //wrapping it in a closure
        auto stable_derivative_2D1C = [] () {return stable_derivative_NDNC<2,1>();};
        TEST(stable_derivative_2D1C()       , success_bool, msg);
        auto stable_derivative_3D1C = [] () {return stable_derivative_NDNC<3,1>();};
        TEST(stable_derivative_3D1C()       , success_bool, msg);


        RETURN_TESTRESULT(success_bool, msg);
    }
};

//specialization disabling some tests with nearest neighbor
template<>
TestResult TestSuite<KernNearestNeighbor>::run_all_tests() {
    DECLARE_TESTRESULT(success_bool, msg);

    TEST(constant_field_1D1C() , success_bool, msg); 
    TEST(constant_field_3D1C() , success_bool, msg);
    TEST(zero_boundary_1D1C()  , success_bool, msg);
    TEST(increasing_field_3D1C()          , success_bool, msg);
    TEST(cyclic_equal_3D3C()              , success_bool, msg);

    RETURN_TESTRESULT(success_bool, msg);
}


int main(int argc, char *argv[])
{
    DECLARE_TESTRESULT(success_bool, msg);

    //TEST(TestSuite<KernNearestNeighbor>::run_all_tests(), success_bool, msg);
    TEST(TestSuite<KernLinear>::run_all_tests()         , success_bool, msg);
    TEST(TestSuite<KernCubic>::run_all_tests()          , success_bool, msg);

    cout<<endl<<msg<<endl;

    cout<<((success_bool)? "All tests succeeded" : "Some tests FAILED")<<endl;

	return (success_bool)? 0 : 1;
}
