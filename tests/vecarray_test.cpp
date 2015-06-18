#include "ndata/debug_helpers.hpp"

#include "ndata.hpp"
#include <memory>

//TODO test dynamic ndarrays
//TODO test array with dimension 0 (slice, etc)

using namespace std;
using namespace ndata;

struct TestSuite {

    static
    test_result
    drop() {
        DECLARE_TEST(success, msg);
        vecarray<int, 3> va ({0, 1, 2});
        vecarray<int, 2> va2 = va.drop(make_vecarray(1ul));
        for (size_t i = 0; i < va2.size(); ++i) {
            success = success and va2[i] == make_vecarray(0,2)[i];
            msg.append(MakeString() << va[i] << " ?== " << va2[i] << "\n");
        }
        RETURN_TESTRESULT(success, msg);
    }

    static
    test_result
    run_all_tests () {
        DECLARE_TEST(success_bool, msg);
        RUN_TEST(drop(), success_bool, msg);
        RETURN_TESTRESULT(success_bool, msg);
    }
};

int main(int argc, char *argv[])
{
    DECLARE_TEST(success_bool, msg);

    RUN_TEST(TestSuite::run_all_tests()  , success_bool, msg);

    cout<<endl<<msg<<endl;

    cout<<((success_bool)? "All tests succeeded" : "Some tests FAILED")<<endl;

	return (success_bool)? 0 : 1;
}

