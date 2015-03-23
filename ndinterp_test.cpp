#include "ndinterp.hpp"
#include "ndindexer.hpp"


#include <sstream>
#include <iostream>
#include <string>
#include <istream>
#include <array>
#include <tuple>

using namespace ndinterp;

using namespace std;

#define DECLARE_TESTRESULT(acc_success_status, acc_message) \
    bool acc_success_status=true; \
    string acc_message (""); \
    bool success_status_TMP (true); \
    acc_success_status = success_status_TMP; /*disable warning when not using TEST macro*/ \
    string messages_TMP (""); \
    acc_message = messages_TMP; \
    \

#define TEST(test_func, acc_success_status, acc_message) \
    acc_message.append(#test_func "\t"); \
    tie(success_status_TMP, messages_TMP) = test_func; \
    acc_success_status = acc_success_status and success_status_TMP; \
    if (success_status_TMP==true) { \
        acc_message.append("success\n"); \
    } else { \
        acc_message.append("FAILURE"); \
    } \
        acc_message.append(MakeString() << "\n" << messages_TMP); \
    acc_message.append("\n"); \
     \

/**
 * nicely indents the result to reflect test function structure
 */
#define RETURN_TESTRESULT(acc_success_status, acc_message) \
    string ret_message (""); \
    istringstream acc_message_istream (acc_message); \
    for(string line; getline(acc_message_istream, line); ) { \
        ret_message.append(MakeString() << "\t" << line << "\n"); \
    } \
    return make_pair(acc_success_status, ret_message); \
    \

/**
 *
 * Use like that
 *
 * #include <string>
 * #include <sstream>
 * #include <iostream>
 *
 * MakeString() << val << "stuff " << endl
 */
class MakeString
{
    public:
        std::stringstream stream;
        operator std::string() const { return stream.str(); }

        template<class T>
        MakeString& operator<<(T const& VAR) { stream << VAR; return *this; }
};

typedef pair<bool, string> TestResult;

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
 
        ndindexer<3> uind (Nn, Nn, Nn);
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
                    OverflowBehaviour::STRETCH
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
  
        ndindexer<3> uind (Nn, Nn, Nn);

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
  
        size_t Nn = 8;

        ndindexer<4> uind (Nn, Nn, Nn, 3) ;

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

        ndindexer<2> valuesInd (Nn*4, 3);

        vector<float> values (valuesInd.size());

        DECLARE_TESTRESULT(aggreg_equal, retMsg);

        bool cycle_equal = true, cpt_equal = true, increasing = true;

        for (size_t itraversal = 0; itraversal < 2; ++itraversal) {
            for (size_t i = 0; i < valuesInd.shape(0); ++i) {

                float xfrac = i/float(valuesInd.shape(0)) * (Nn-1.f);

                array<float, 3> val = interpolate<KernT, 3, 3>(
                        {Nn, Nn, Nn},
                        &u[0],
                        {2.5, xfrac, 3.5},
                        {
                            OverflowBehaviour::CYCLIC,
                            OverflowBehaviour::CYCLIC,
                            OverflowBehaviour::STRETCH
                        });

                if (itraversal != 0) {
                    for (size_t ic = 0; ic < 3; ++ic)
                    {
                        if(val[ic] != values[valuesInd.index(i, ic)]) {
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

                if(xfrac >= 2 and xfrac <= Nn-2) {
                    for (size_t ic = 0; ic < 3; ++ic) {
                        if( values[valuesInd.index(i-1, ic)] > values[valuesInd.index(i, ic)]) {
                            increasing = false;
                        }
                    }
                }
            }


            for (size_t ic = 0; ic < 3; ++ic) {
                for (size_t iv = 0; iv < valuesInd.shape(0); ++iv)
                {
                    retMsg.append(MakeString() << values[valuesInd.index(iv, ic)]
                            << ", "); 
                }

                retMsg.append("\n");
            }
        }
        aggreg_equal = cpt_equal and cycle_equal and increasing;
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
        TEST(cyclic_equal_3D3C()              , success_bool, msg);

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

    TEST(TestSuite<KernNearestNeighbor>::run_all_tests(), success_bool, msg);
    TEST(TestSuite<KernLinear>::run_all_tests()         , success_bool, msg);
    TEST(TestSuite<KernCubic>::run_all_tests()          , success_bool, msg);

    cout<<endl<<msg<<endl;

    cout<<((success_bool)? "All tests succeeded" : "Some tests FAILED")<<endl;

	return (success_bool)? 0 : 1;
}

