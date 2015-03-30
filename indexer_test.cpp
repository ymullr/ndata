#include "debug_helpers.hpp"
#include "ndata.hpp"

using namespace std;
using namespace ndata;

struct TestSuite {

    static
    TestResult
    sequential_indexing() {

        DECLARE_TESTRESULT(success_bool, ret_msg);

        size_t Nx=2, Ny=5;

        indexer<2> ndi_u ({Nx, Ny});

        //vector<float> u (ndi_u.size());

        //for (size_t i = 0; i < u.size(); ++i) {
        //    u[i] = i;
        //}
        string ref_array, ndindexer_array;

        for (size_t ix = 0; ix < Nx; ++ix) {
            for (size_t iy = 0; iy < Ny; ++iy) {
                // manual 2D indexing
                size_t indval = ix*Ny+iy;
                size_t ndi_val = ndi_u.index(ix, iy);

                ref_array.append(MakeString() << indval << ", ");

                ndindexer_array.append(MakeString() << ndi_val << ", ");

                success_bool = success_bool and indval == ndi_val;
            }
            ndindexer_array.append("\n");
            ref_array.append("\n");
        }

        ret_msg.append("reference :\n");
        ret_msg.append(ref_array);
        ret_msg.append("\n\n indexer result :\n");
        ret_msg.append(ndindexer_array);

        RETURN_TESTRESULT(success_bool, ret_msg);
    };


    static
    TestResult
    slice_test() {

        DECLARE_TESTRESULT(success_bool, ret_msg);

        size_t N0=2, N1=5;

        indexer<2> ndi_u ({N0, N1});

        vector<float> u (ndi_u.size());

        for (size_t i = 0; i < u.size(); ++i) {
            u[i] = i;
        }

        //ndindexer<2> ndi_fullslice1 = ndi_u.slice(array<long, 3> {0, long(N0), 1}, array<long, 3> {0, long(N1), 1});
        
        indexer<2> ndi_fullslice1 = ndi_u.slice(Rng(0, long(N0), 1), Rng(0, long(N1), 1));

        indexer<2> ndi_fullslice2 = ndi_u.slice(
                Rng(0, N0),
                Rng(0, N1, 1)
                );

        indexer<2> ndi_fullslice3 = ndi_u.slice(Rng(0, N0), Rng(0, ndata::END));//END is = -1

        bool fs1ok, fs2ok, fs3ok;

        fs1ok= fs2ok= fs3ok = true;

        string fs1_msg;
        string fs2_msg;
        string fs3_msg;

        for (size_t i0 = 0; i0 < N0; ++i0) {
            for (size_t i1 = 0; i1 < N1; ++i1) {
                // manual 2D indexing
                size_t indval = i0*N1+i1;

                size_t i_fs1 = ndi_fullslice1.index(i0, i1); 
                size_t i_fs2 = ndi_fullslice2.index(i0, i1); 
                size_t i_fs3 = ndi_fullslice3.index(i0, i1); 

                fs1ok = i_fs1 == indval;
                fs2ok = i_fs2 == indval;
                fs3ok = i_fs3 == indval;

                fs1_msg.append(MakeString() << i_fs1 << ", ");
                fs2_msg.append(MakeString() << i_fs2 << ", ");
                fs3_msg.append(MakeString() << i_fs3 << ", ");                 

                ret_msg.append(MakeString() << indval << ", ");


                success_bool = success_bool and fs1ok and fs2ok and fs3ok;
            }

            ret_msg.append("\n");
        }


        ret_msg.append(fs1_msg);
        ret_msg.append("\n\n");
        ret_msg.append(fs2_msg);
        ret_msg.append("\n\n");
        ret_msg.append(fs3_msg);

        RETURN_TESTRESULT(success_bool, ret_msg);
    };

    static
    TestResult
    extended_slices() {
        auto ndi = make_indexer(4, 5, 3);
        indexer<2> ndsli = ndi.slice(Rng(0, 2), 2, Rng());

        DECLARE_TESTRESULT(success, msg);

        success = success and ndsli.size()==3*2;

        for (size_t i = 0; i < 2; ++i) {
            for (size_t i1 = 0; i1 < 3; ++i1) {
                size_t idx = ndsli.index(i, i1);
                msg.append(MakeString() << idx);
            }
        }

        RETURN_TESTRESULT(success, msg);
    }

    static
    TestResult view_test() {
        DECLARE_TESTRESULT(success, msg);

        size_t Nx=3, Ny=15;

        nvector<2, float> u ({Nx, Ny});

        for (size_t i = 0; i < u.size(); ++i) {
            u[i] = i;
        }

        nvector<2, float> usli = u.slice(Rng(), Rng(1, 4));

        for (size_t ix = 0; ix < usli.get_shape()[0]; ++ix) {
            for (size_t iy = 0; iy < usli.get_shape()[1]; ++iy) {
                success = usli.val(ix, iy) == ix*Ny+iy;
                msg.append(MakeString() << usli.val(ix, iy) << ", ");
            }
            msg.append("\n");
        }

        RETURN_TESTRESULT(success, msg);
    }
 
    static
    TestResult run_all_tests () {
        DECLARE_TESTRESULT(b, s);
        //TEST(sequential_indexing(), b, s);
        //TEST(slice_test(), b, s);
        TEST(extended_slices(), b, s);
        TEST(view_test(), b, s);
        RETURN_TESTRESULT(b, s);
    }
};

int main(int argc, char *argv[])
{
    DECLARE_TESTRESULT(success_bool, msg);

    TEST(TestSuite::run_all_tests()  , success_bool, msg);

    cout<<endl<<msg<<endl;

    cout<<((success_bool)? "All tests succeeded" : "Some tests FAILED")<<endl;

	return (success_bool)? 0 : 1;
}

