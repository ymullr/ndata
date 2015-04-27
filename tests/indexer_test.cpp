#include "debug_helpers.hpp"

#include "ndata.hpp"
#include <memory>

//TODO test dynamic ndarrays

using namespace std;
using namespace ndata;

struct TestSuite {

    static
    TestResult
    sequential_indexing() {

        DECLARE_TESTRESULT(success_bool, ret_msg);

        size_t Nx=2, Ny=5;

        indexer<2> ndi_u (Nx, Ny);

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

        indexer<2> ndi_u (N0, N1);

        vector<float> u (ndi_u.size());

        for (size_t i = 0; i < u.size(); ++i) {
            u[i] = i;
        }

        //ndindexer<2> ndi_fullslice1 = ndi_u.slice(array<long, 3> {0, long(N0), 1}, array<long, 3> {0, long(N1), 1});
        
        indexer<2> ndi_fullslice1 = ndi_u.index_slice(Rng(0, long(N0), 1), Rng(0, long(N1), 1));

        indexer<2> ndi_fullslice2 = ndi_u.index_slice(
                Rng(0, N0),
                Rng(0, N1, 1)
                );

        indexer<2> ndi_fullslice3 = ndi_u.index_slice(Rng(0, N0), Rng(0, ndata::END));//END is = -1

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
        indexer<2> ndsli = ndi.index_slice(Rng(0, 2), 2, Rng());

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
    TestResult nvector_test() {
        DECLARE_TESTRESULT(success, msg);

        size_t Nx=3, Ny=5;

        nvector<long, 2> u (make_indexer(Nx, Ny), {0});

        for (size_t i = 0; i < u.size(); ++i) {
            u[i] = i;
        }

        long iy_start = 1;
        nvector<long, 2> usli = u.slice(Rng(), Rng(iy_start, 4));

        for (size_t ix = 0; ix < usli.get_shape()[0]; ++ix) {
            for (size_t iy = 0; iy < usli.get_shape()[1]; ++iy) {
                long indval = usli(ix, iy);
                success = indval == long(ix*Ny+(iy+iy_start));
                msg.append(MakeString() << usli(ix, iy) << ", ");
            }
            msg.append("\n");
        }

        RETURN_TESTRESULT(success, msg);
    }

    static
    TestResult transform_test () {

        DECLARE_TESTRESULT(sb, msg);

        size_t Nx = 5, Ny = 2;
        auto u1 = make_nvector<long>(make_indexer(Nx, Ny), 0);
        auto u2 = make_nvector<long>(make_indexer(Nx, Ny), 0);

        for (size_t ix = 0; ix < Nx ; ++ix) {
            for (size_t iy = 0; iy < Ny ; ++iy) {
                u1(ix ,iy) = ix*Ny + iy;
                u2(ix ,iy) = ix*Ny + iy;
                //msg.append(MakeString() << u1.val(ix, iy) << ", ");
            }
            msg.append("\n");
        }

        nvector<long, 2> u_sum = ntransform<long>(tie(u1, u2), [] (long v1, long v2) {
            return v1+v2;
        });

        for (size_t ix = 0; ix < Nx ; ++ix) {
            for (size_t iy = 0; iy < Ny ; ++iy) {
                if (
                        not (u1(ix ,iy)*2 == u_sum(ix, iy))
                        )
                {
                    sb = false;
                }
                msg.append(MakeString() << u_sum(ix, iy) << ", ");
            }
            msg.append("\n");
        }

        RETURN_TESTRESULT(sb, msg)

    }

    static TestResult
    broadcast_test() {
        DECLARE_TESTRESULT(succ, msg);

        auto u1 = make_nvector<long>(make_indexer(2, 5), 0);
        auto u2 = make_nvector<long>(make_indexer(2, 1), 0);
        auto u3 = make_nvector<long>(make_indexer(5), 0);

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
            msg.append(MakeString() << "shape u1[i] :" << u1_broad.get_shape()[i]
                       << " shape_u2[i]: " << u2_broad.get_shape()[i]
                       << " shape_u3[i]: " << u3_broad.get_shape()[i]
                       << "\n"
                       );
        }

        succ = succ_u1u2 and succ_u1u3;

        RETURN_TESTRESULT(succ, msg);
    }

    static TestResult
    view_test() {

        DECLARE_TESTRESULT(succ, msg);

        //zero initialized
        auto u = make_nvector<long>(make_indexer(2, 5));

        auto uview = ndataview<long, 2>(u);

        for (size_t ix = 0; ix < 2; ++ix) {
            for (size_t iy = 0; iy < 5; ++iy) {
                uview(ix, iy) = 2;
            }
        }

        for (size_t ix = 0; ix < 2; ++ix) {
            for (size_t iy = 0; iy < 5; ++iy) {
                if (u(ix, iy) != 2) {
                    succ = false;
                };
                msg.append(MakeString () << u(ix, iy) << ", ");
            }
            msg.append("\n");
        }

        RETURN_TESTRESULT(succ, msg);

    }

    static
    TestResult extended_transform_test () {

        DECLARE_TESTRESULT(sb, msg);

        size_t Nx = 5, Ny = 2;
        auto u1 = make_nvector<long>(make_indexer(Nx, Ny));
        auto u2 = make_nvector<float>(make_indexer(Nx*2, Ny));
        auto u3 = make_nvector<std::pair<long, long>>(make_indexer(Ny));

        for (size_t ix = 0; ix < Nx ; ++ix) {
            for (size_t iy = 0; iy < Ny ; ++iy) {
                u1(ix ,iy) = ix*Ny + iy;
            }
        }


        for (size_t ix = 0; ix < Nx*2 ; ++ix) {
            for (size_t iy = 0; iy < Ny ; ++iy) {
                u2(ix ,iy) = ix*Ny + iy;
            }
        }

        for (size_t iy = 0; iy < Ny ; ++iy) {
            u3(iy) = make_pair(1,1);
        }

        //zero initialized
        auto ures = make_nvector<float>(make_indexer(Nx, Ny));

        ndataview<float, 2> u2_slice =  u2.slice(Rng(0, Nx), Rng());

        //this one should be uncopyable, good for checking that no parasitic
        //copy happens in foreach
        //must be release or double free happens
        ndatacontainer<
                std::unique_ptr<float[]>,
                float,
                2
                >
                u2_slice_unique_view (
                    u2_slice,
                    std::unique_ptr<float[]>(&u2_slice.data_[0])
                    );

        nforeach<float>(
            //nforeach takes a tuple of pointers to ndataview ndatacontainer or nvector
            std::tie( //must use std tie with foreach
                ures,
                u1,
                //u2.slice(Rng(0, Nx), Rng()),
                u2_slice_unique_view,
                u3
            ),
            //lambda function must take its arguments as pointers
            [] (float &res, long v1, float v2, std::pair<long, long> v3) {
                res = v1+v2+v3.first;
            }
        );

        //so that the unique pointer will not delete the array itself
        u2_slice_unique_view.data_.release();

        auto u_res_transform = ntransform<float>(
            std::make_tuple(
                u1,
                //u2.slice_view(Rng(0, Nx), Rng()),
                u2_slice,
                u3
            ),
            [] (long v1, float v2, std::pair<long, long> v3) {
                return v1+v2+v3.first;
            }
        );

        string msg_foreach ("");
        string msg_transform ("");
        for (size_t ix = 0; ix < Nx ; ++ix) {
            for (size_t iy = 0; iy < Ny ; ++iy) {
                if (
                    not (ures(ix, iy) == u_res_transform(ix, iy))
                )
                {
                    sb = false;
                }
                msg_foreach.append(MakeString() << ures(ix, iy) << ", ");
                msg_transform.append(MakeString() << u_res_transform(ix, iy) << ", ");
            }
            msg_transform.append("\n");
            msg_foreach.append("\n");
        }

        msg.append(msg_transform);
        msg.append("\n");
        msg.append(msg_foreach);
        RETURN_TESTRESULT(sb, msg)
    }

    static
    TestResult run_all_tests () {
        DECLARE_TESTRESULT(b, s);
        TEST(sequential_indexing(), b, s);
        TEST(slice_test(), b, s);
        TEST(extended_slices(), b, s);
        TEST(nvector_test(), b, s);
        TEST(transform_test(), b, s);
        TEST(view_test(), b, s);
        TEST(broadcast_test(), b, s);
        TEST(extended_transform_test(), b, s);
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

