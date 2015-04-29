#include "ndata/numerics/ninterp.hpp"
#include "nindexer.hpp"

#include <iostream>
#include <string>
#include <array>
#include <tuple>

#include <tests/debug_helpers.hpp>

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
    TestResult constant_field_3D() {
 
        indexer<3> uind (Nn, Nn, Nn);
        const float init_val = 1;
        nvector<float, 3> u (uind, init_val);
        float start_ifrac = -1.9;

        DECLARE_TESTRESULT(allCorrect, output);

        for (float i = start_ifrac; i < Nn*3; i+=0.5) {

            float out = interpolate<KernT>(
                    u,
                    {2, i, 2},
                    OverflowBehaviour::STRETCH
                    );

            output.append(MakeString()<<out<<", ");

            if (fabs(out - init_val) > init_val*1e-5) {
                allCorrect = false;
            }  
        }

        RETURN_TESTRESULT(allCorrect, output)
    }


    static
    TestResult constant_field_1D() {

        const float init_val = 1;

        nvector<float, 1> u (make_indexer(Nn), init_val);

        float start_ifrac = -1.9;

        DECLARE_TESTRESULT(allCorrect, output);

        for (float i = start_ifrac; i < Nn*3; i+=0.5) {

            float i_frac = i; 

            float out = interpolate<KernT>(
                u,
                i_frac,
                OverflowBehaviour::STRETCH
            );

            output.append(MakeString()<< out<<", ");

            if (fabs(out - init_val) > init_val*1e-5) {
                allCorrect = false;
            }  
        }

        RETURN_TESTRESULT(allCorrect, output);
    }

    static
    TestResult zero_boundary_1D() {

        const float init_val = 1;
        nvector<float, 1> u (make_indexer(Nn), init_val);
        float start_ifrac = -1.9;

        DECLARE_TESTRESULT(success_bool, output);

        for (float i = start_ifrac; i < Nn*3; i+=0.5) {

            float i_frac = i;  

            float out = interpolate<KernT>(
                u,
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
    TestResult increasing_field_1D(OverflowBehaviour ob=OverflowBehaviour::STRETCH) {
  
        nvector<float, 1> u (make_indexer(Nn), 0);
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
                    u,
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
    TestResult increasing_field_3D() {
  
        nvector<float, 3> u (make_indexer(Nn, Nn, Nn), 0);

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
                    u,
                    {i_frac, i_frac, i_frac},
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
    TestResult simple_equalities_1D() {
  
        size_t Nn = 8;

        nvector<float, 1> u (make_indexer(Nn), 0);
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
                    u,
                    i_frac,
                    OverflowBehaviour::STRETCH
                    );


            float valhalf = interpolate<KernT>(
                    u,
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
    TestResult cyclic_equal_3D() {
  
        size_t Nn = 10;

        nvector<float, 3> u (make_indexer(Nn, Nn, Nn), 0);

        vecarray<size_t, 3> ndind = u.ndindex(0);

        for (size_t i = 0; i < u.size(); ++i) {
            u[i] = sqrt(
                          ndind[0]*ndind[0] 
                        + ndind[1]*ndind[1] 
                        + ndind[2]*ndind[2] 
                );
            u.increment_ndindex(ndind);
        }


        DECLARE_TESTRESULT(cycle_equal, retMsg);

        const size_t ntraversals = 4;
        const size_t nsteps = Nn*4;
        const size_t nvalues_eff = 3;

        nvector<float, 2> values (make_indexer(ntraversals, nvalues_eff), 0);

        for (size_t itraversal = 0; itraversal < ntraversals; ++itraversal) {
            for (size_t i = 0; i < nvalues_eff; ++i) {

                float xfrac =
                        (long(itraversal)-1) * long(Nn)
                        + i/float(nsteps) * (Nn-1.f);

                float val = interpolate<KernT>(
                    u,
                    make_vecarray(xfrac, 4.5f, 4.5f),
                    make_vecarray(
                        OverflowBehaviour::CYCLIC,
                        OverflowBehaviour::CYCLIC,
                        OverflowBehaviour::CYCLIC
                    )
                );

                values(itraversal, i) = val;

                if (itraversal > 0) {
                        float diff = values(itraversal, i)-values(itraversal-1, i);
                        if(
                            fabs(diff) > 1e-5f
                          ) {
                            cycle_equal = false;
                        };
                } 


            }
        }


        for (size_t i = 0; i < nvalues_eff; ++i) {
            for (size_t itraversal = 0; itraversal < ntraversals; ++itraversal) {
                retMsg.append(MakeString() << values(itraversal, i) << ", ");
            }
            retMsg.append("\n");
        }


        RETURN_TESTRESULT(cycle_equal, retMsg);
    }

    template<size_t ndims>
    static
    TestResult stable_derivative_NDNC() {
  
        size_t Nn = 10;

        indexer<ndims> vecindexr (vecarray<long, ndims>(STATICALLY_SIZED, Nn));

        vecarray<OverflowBehaviour, ndims> ovfl (STATICALLY_SIZED, OverflowBehaviour::CYCLIC);
        ovfl.back()=OverflowBehaviour::STRETCH;

        nvector<float, ndims> u (vecindexr, 0.f);

        vecarray<size_t, ndims> ndind (u.ndindex(0));

        for (size_t i = 0; i < u.size(); ++i) {
            //u[i] = sqrt(
            //              ndind[0]*ndind[0]
            //            + ndind[1]*ndind[1]
            //            + ndind[2]*ndind[2]
            //    );
            u[i] = ndind[0];
            u.increment_ndindex(ndind);
        }

        nvector<float, 1> values (make_indexer(Nn*4), 0);


        DECLARE_TESTRESULT(aggreg_equal, retMsg);

        bool increasing = true, stable_derivative = true;

        for (size_t i = 0; i < values.size(); ++i) {

            float xfrac = i/float(values.size()) * (Nn-1.f);

            vecarray<float, ndims> indx (std::vector<float>(ndims, 3.5));
            indx[0] = xfrac;

            float val = interpolate<KernT>(
                        u,
                        indx,
                        ovfl
                        );

            values(i)=val;

            //test on derivatives
            //only test away from cyclic boundaries
            if(     xfrac >= KernT::ONESIDEDWIDTH+2
                    and xfrac <= (Nn-1)-KernT::ONESIDEDWIDTH
                    ) {


                float d1 =
                        values(i)
                        - values(i-1),
                        d2 =
                        values(i-1)
                        - values(i-2);

                if(d1 < 0 or d2 < 0) {
                    increasing = false;
                }


                if (fabs(d2-d1) > 1e-4*d1) {
                    stable_derivative =false;
                }
            }
        }

        for (size_t iv = 0; iv < values.size(); ++iv) {
            retMsg.append(MakeString() << values(iv) << ", ");
            retMsg.append("\n");
        }

        aggreg_equal = increasing and stable_derivative;
        RETURN_TESTRESULT(aggreg_equal, retMsg);
    }

    static
    TestResult run_all_tests() {

        DECLARE_TESTRESULT(success_bool, msg);

        TEST(constant_field_1D()            , success_bool, msg);
        TEST(zero_boundary_1D()             , success_bool, msg);
        TEST(increasing_field_1D()          , success_bool, msg);
        TEST(increasing_field_1D(OverflowBehaviour::CYCLIC)
                , success_bool, msg);
        TEST(simple_equalities_1D()         , success_bool, msg);
        TEST(constant_field_3D()            , success_bool, msg);
        TEST(increasing_field_3D()          , success_bool, msg);
        TEST(cyclic_equal_3D()              , success_bool, msg);

        //yeah.. macros dosnt like multiple template arguments...
        //wrapping it in a closure
        auto stable_derivative_2D = [] () {return stable_derivative_NDNC<2>();};
        TEST(stable_derivative_2D()       , success_bool, msg);
        auto stable_derivative_3D = [] () {return stable_derivative_NDNC<3>();};
        TEST(stable_derivative_3D()       , success_bool, msg);


        RETURN_TESTRESULT(success_bool, msg);
    }
};

//specialization disabling some tests with nearest neighbor
template<>
TestResult TestSuite<KernNearestNeighbor>::run_all_tests() {
    DECLARE_TESTRESULT(success_bool, msg);

    TEST(constant_field_1D() , success_bool, msg);
    TEST(constant_field_3D() , success_bool, msg);
    TEST(zero_boundary_1D()  , success_bool, msg);
    TEST(increasing_field_3D()          , success_bool, msg);
    TEST(cyclic_equal_3D()              , success_bool, msg);

    RETURN_TESTRESULT(success_bool, msg);
}


int main(int argc, char *argv[])
{
    DECLARE_TESTRESULT(success_bool, msg);

    TEST(TestSuite<KernNearestNeighbor>::run_all_tests(), success_bool, msg);
    TEST(TestSuite<KernLinear>::run_all_tests()         , success_bool, msg);
    TEST(TestSuite<KernCubic>::run_all_tests()          , success_bool, msg);

    cout<<endl<<msg<<endl;

    cout<<((success_bool)? "All tests succeeded" : "Some tests FAILED")<<endl;

	return (success_bool)? 0 : 1;
}

