#include "debug_helpers.hpp"

#include "ndata.hpp"

using namespace std;
using namespace ndata;

struct uncopyable {
    uncopyable(uncopyable&) ;
};

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
    TestResult nvector_test() {
        DECLARE_TESTRESULT(success, msg);

        size_t Nx=3, Ny=5;

        nvector<long, 2> u (Nx, Ny);

        for (size_t i = 0; i < u.size(); ++i) {
            u[i] = i;
        }

        long iy_start = 1;
        nvector<long, 2> usli = u.slice(Rng(), Rng(iy_start, 4));

        for (size_t ix = 0; ix < usli.get_shape()[0]; ++ix) {
            for (size_t iy = 0; iy < usli.get_shape()[1]; ++iy) {
                long indval = usli.val(ix, iy);
                success = indval == long(ix*Ny+(iy+iy_start));
                msg.append(MakeString() << usli.val(ix, iy) << ", ");
            }
            msg.append("\n");
        }

        RETURN_TESTRESULT(success, msg);
    }

    static
    TestResult transform_test () {

        DECLARE_TESTRESULT(sb, msg);

        size_t Nx = 5, Ny = 2;
        auto u1 = make_nvector<long>(Nx, Ny);
        auto u2 = make_nvector<long>(Nx, Ny);

        for (size_t ix = 0; ix < Nx ; ++ix) {
            for (size_t iy = 0; iy < Ny ; ++iy) {
                u1.val(ix ,iy) = ix*Ny + iy;
                u2.val(ix ,iy) = ix*Ny + iy;
                //msg.append(MakeString() << u1.val(ix, iy) << ", ");
            }
            msg.append("\n");
        }

        nvector<long, 2> u_sum = ntransform<long>(make_tuple(u1, u2), [] (long v1, long v2) {
            return v1+v2;
        });

        for (size_t ix = 0; ix < Nx ; ++ix) {
            for (size_t iy = 0; iy < Ny ; ++iy) {
                if (
                        not (u1.val(ix ,iy)*2 == u_sum.val(ix, iy))
                        )
                {
                    sb = false;
                }
                msg.append(MakeString() << u_sum.val(ix, iy) << ", ");
            }
            msg.append("\n");
        }

        RETURN_TESTRESULT(sb, msg)

    }

    static TestResult
    broadcast_test() {
        DECLARE_TESTRESULT(succ, msg);

        auto u1 = make_nvector<long>(2, 5);
        auto u2 = make_nvector<long>(2, 1);
        auto u3 = make_nvector<long>(5);

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
        auto u = make_nvector<long>(2, 5);

        auto uview = ndataview<long, 2>(u);

        for (size_t ix = 0; ix < 2; ++ix) {
            for (size_t iy = 0; iy < 5; ++iy) {
                uview.val(ix, iy) = 2;
            }
        }

        for (size_t ix = 0; ix < 2; ++ix) {
            for (size_t iy = 0; iy < 5; ++iy) {
                if (u.val(ix, iy) != 2) {
                    succ = false;
                };
                msg.append(MakeString () << u.val(ix, iy) << ", ");
            }
            msg.append("\n");
        }

        RETURN_TESTRESULT(succ, msg);

    }

    static
    TestResult extended_transform_test () {

        DECLARE_TESTRESULT(sb, msg);

        size_t Nx = 5, Ny = 2;
        auto u1 = make_nvector<long>(Nx, Ny);
        auto u2 = make_nvector<float>(Nx*2, Ny);
        auto u3 = make_nvector<std::pair<long, long>>(Ny);

        for (size_t ix = 0; ix < Nx ; ++ix) {
            for (size_t iy = 0; iy < Ny ; ++iy) {
                u1.val(ix ,iy) = ix*Ny + iy;
            }
        }


        for (size_t ix = 0; ix < Nx*2 ; ++ix) {
            for (size_t iy = 0; iy < Ny ; ++iy) {
                u2.val(ix ,iy) = ix*Ny + iy;
            }
        }

        for (size_t iy = 0; iy < Ny ; ++iy) {
            u3.val(iy) = make_pair(1,1);
        }

        //zero initialized
        auto ures = make_nvector<float>(Nx, Ny);

        nforeach<float>(
            make_tuple(
                //mutated variables must be passed as a view to foreach or changes wont be reflected in caller scope
                //(so that data is passed by reference instead of by value)
                ures,
                u1,
                u2.slice_view(Rng(0, Nx), Rng()),
                u3
            ),
            //lambda function must take its arguments as pointers
            [] (float * res, long * v1, float * v2, std::pair<long, long> * v3) {
                *res = *v1+*v2+(*v3).first;
            }
        );

        auto u_res_transform = ntransform<float>(
            make_tuple(
                u1,
                u2.slice_view(Rng(0, Nx), Rng()),
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
                    not (ures.val(ix, iy) == u_res_transform.val(ix, iy))
                )
                {
                    sb = false;
                }
                msg_foreach.append(MakeString() << ures.val(ix, iy) << ", ");
                msg_transform.append(MakeString() << u_res_transform.val(ix, iy) << ", ");
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

