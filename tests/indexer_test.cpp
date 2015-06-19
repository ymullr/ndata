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
    sequential_indexing() {

        DECLARE_TEST(success_bool, ret_msg);

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
    test_result
    slice_test() {

        DECLARE_TEST(success_bool, ret_msg);

        size_t N0=2, N1=5;

        indexer<2> ind (N0, N1);

        //"full" slices with diferent syntaxes
        indexer<2> ind_fullslice1 = ind.slice_indexer(range(0, long(N0), 1), range(0, long(N1), 1));
        indexer<2> ind_fullslice2 = ind.slice_indexer(
                range(0, N0),
                range(0, N1, 1)
                );
        indexer<2> ind_fullslice3 = ind.slice_indexer(range(0, N0), range(0, ndata::END));//END is = -1


        bool fs1ok, fs2ok, fs3ok;
        fs1ok = fs2ok = fs3ok = true;

        string fs1_msg;
        string fs2_msg;
        string fs3_msg;

        for (size_t i0 = 0; i0 < N0; ++i0) {
            for (size_t i1 = 0; i1 < N1; ++i1) {
                // manual 2D indexing
                size_t indval = i0*N1+i1;

                size_t i_fs1 = ind_fullslice1.index(i0, i1);
                size_t i_fs2 = ind_fullslice2.index(i0, i1);
                size_t i_fs3 = ind_fullslice3.index(i0, i1);


                fs1ok = fs1ok and i_fs1 == indval;
                fs2ok = fs2ok and i_fs2 == indval;
                fs3ok = fs3ok and i_fs3 == indval;

                fs1_msg.append(MakeString() << i_fs1 << ", ");
                fs2_msg.append(MakeString() << i_fs2 << ", ");
                fs3_msg.append(MakeString() << i_fs3 << ", ");                 

                ret_msg.append(MakeString() << indval << ", ");


            }

            ret_msg.append("\n");
        }


        ret_msg.append(fs1_msg);
        ret_msg.append("\n\n");
        ret_msg.append(fs2_msg);
        ret_msg.append("\n\n");
        ret_msg.append(fs3_msg);
        ret_msg.append("\n\n");

        bool newdim_slice_ok = true;

        indexer<1> ind2 (5);

        //insert a new dimension with a single item
        //useful with broadcasting
        indexer<2> ind2_slice_newdim = ind2.slice_indexer(range(), NEWDIM);
        //TODO test alternative overload of indexer slice with newdim


        for (size_t i = 0; i < ind2.size(); ++i) {
            size_t i_slice_newdim =  ind2_slice_newdim.index(i, 0);
            newdim_slice_ok = newdim_slice_ok and i_slice_newdim == i;
            ret_msg.append(MakeString() << i_slice_newdim << "\n");
        }

        success_bool = fs1ok and fs2ok and fs3ok and newdim_slice_ok;


        RETURN_TESTRESULT(success_bool, ret_msg);
    };

    static
    test_result
    extended_slices() {
        auto ndi = make_indexer(4, 5, 3);
        indexer<2> ndsli = ndi.slice_indexer(range(0, 2), 2, range());

        //Also test indexer_slice_alt method
        indexer<2> ndsli_v2 = ndi.slice_indexer_alt(
                    make_vecarray(range(0, 2), range()),
                    make_vecarray(0ul, 2ul),
                    make_vecarray(2l),
                    make_vecarray(1ul)
                );

        DECLARE_TEST(success, msg);

        success = ndsli.size()==3*2 and ndsli_v2.size() == 3*2;
        msg.append(MakeString() << "sizes must match 3*2 "<< ndsli.size() << ", " <<  ndsli_v2.size() << "\n");

        for (size_t i = 0; i < 2; ++i) {
            for (size_t i1 = 0; i1 < 3; ++i1) {
                size_t idx = ndsli.index(i, i1);
                size_t idx2 = ndsli_v2.index(i, i1);
                size_t idx_ref = i*5*3 + 2*3 + i1;

                bool eq1 = idx == idx_ref,
                     eq2 = idx2 == idx_ref;

                success = success and eq1 and eq2;

                msg.append(MakeString() << "indices (must match "<< idx_ref <<"): " << idx << ", " << idx2 << "\n");
            }
        }

        RETURN_TESTRESULT(success, msg);
    }

    static
    test_result nvector_test() {
        DECLARE_TEST(success, msg);

        size_t Nx=3, Ny=5;

        nvector<long, 2> u (make_indexer(Nx, Ny), 0l);
        nvector<long, 2> u_unini (make_indexer(Nx, Ny), ndata::UNINITIALIZED);

        for (size_t i = 0; i < u.size(); ++i) {
            u[i] = i;
        }

        long iy_start = 1;
        nvector<long, 2> usli = u.slice(range(), range(iy_start, 4));

        for (long ix = 0; ix < usli.get_shape()[0]; ++ix) {
            for (long iy = 0; iy < usli.get_shape()[1]; ++iy) {
                long indval = usli(ix, iy);
                if (
                        not
                        indval == long(ix*Ny+(iy+iy_start))
                        ){
                    success = false;
                }
                msg.append(MakeString() << usli(ix, iy) << ", ");
            }
            msg.append("\n");
        }

        RETURN_TESTRESULT(success, msg);
    }

    static
    test_result transform_test () {

        DECLARE_TEST(sb, msg);

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

        nvector<long, 2> u_sum = ntransform<long>(
            make_tuple(u1, u2),
            [] (long v1, long v2) {
                return v1+v2;
            }
        );

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

    static test_result
    broadcast_test() {
        DECLARE_TEST(succ, msg);

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

    static test_result
    view_test() {

        DECLARE_TEST(succ, msg);

        //zero initialized
        auto u = make_nvector<long>(make_indexer(2, 5));

        auto uview = u.as_view();

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
    test_result assign_transform_slice_alt () {
        DECLARE_TEST(sb, msg);

        size_t Nx = 5, Ny = 2;
        auto u1 = make_nvector<long>(make_indexer(Nx, Ny));
        auto u2 = make_nvector<float>(make_indexer(Nx*2, Ny));

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

        //zero initialized
        nvector<float, 2> ures_assign = make_nvector<float>(make_indexer(Nx, Ny));

        ndataview<float, 2> u2_slice =  u2.slice(range(0, Nx), range());
        ndataview<float, 2> u2_slice_alt =  u2.slice_alt(
                    make_vecarray(
                        range(0, Nx),
                        range()
                        ),
                    make_vecarray(0ul,1ul)
                    );


        ures_assign.assign_transform(
                    make_tuple(u1, u2_slice_alt),
                    [] (auto v1, auto v2) {
                        return v1+v2;
                    });
        auto u_res = ntransform<float>(
                    make_tuple(u1, u2_slice),
                    [] (auto v1, auto v2) {
                        return v1+v2;
                    });

        nforeach(tie(u_res, ures_assign), [&sb, &msg] (auto& v1, auto& v2) {
            sb = sb and v1 == v2;
            msg.append(MakeString() << v1 << " ==? " << v2 << "\n");
        });

        RETURN_TESTRESULT(sb, msg);

    }

    static
    test_result extended_transform_test () {

        DECLARE_TEST(sb, msg);

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

        ndataview<float, 2> u2_slice =  u2.slice(range(0, Nx), range());

        //this one should be uncopyable, good for checking that no unwanted
        //copy of the data happens in foreach
        //must be released via .release()
        ndatacontainer<
                std::unique_ptr<float[]>,
                float,
                2
                >
                u2_slice_unique_view (
                    u2_slice,
                    std::unique_ptr<float[]>(&u2_slice.data_[0])
                    );

        nforeach(
            //nforeach takes a tuple of pointers to ndataview ndatacontainer or nvector
            std::tie( //must use std tie with foreach
                ures,
                u1,
                //u2.slice(range(0, Nx), range()),
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
                //u2.slice_view(range(0, Nx), range()),
                u2_slice,
                u3
            ),
            [] (long v1, float v2, std::pair<long, long> v3) {
                return v1+v2+v3.first;
            }
        );


        auto u_res_transform_parallel = ntransform<float, ndata::PARALLEL>(
            std::make_tuple(
                u1,
                u2_slice,
                u3
            ),
            [] (long v1, float v2, std::pair<long, long> v3) {
                return v1+v2+v3.first;
            }
        );

        string msg_foreach ("");
        string msg_transform ("");
        string msg_transform_par ("");
        for (size_t ix = 0; ix < Nx ; ++ix) {
            for (size_t iy = 0; iy < Ny ; ++iy) {
                if (
                    not (
                            ures(ix, iy)
                            == u_res_transform(ix, iy)
                            and
                            ures(ix, iy)
                            == u_res_transform_parallel(ix, iy)
                            )
                )
                {
                    sb = false;
                }
                msg_foreach.append(MakeString() << ures(ix, iy) << ", ");
                msg_transform.append(MakeString() << u_res_transform(ix, iy) << ", ");
                msg_transform_par.append(MakeString() << u_res_transform_parallel(ix, iy) << ", ");
            }
            msg_foreach.append("\n");
            msg_transform.append("\n");
            msg_transform_par.append("\n");
        }

        msg.append(msg_foreach);
        msg.append("\n");
        msg.append(msg_transform);
        msg.append("\n");
        msg.append(msg_transform_par);
        RETURN_TESTRESULT(sb, msg)
    }

    static
    test_result run_all_tests () {
        DECLARE_TEST(b, s);
        RUN_TEST(sequential_indexing(), b, s);
        RUN_TEST(slice_test(), b, s);
        RUN_TEST(extended_slices(), b, s);
        RUN_TEST(nvector_test(), b, s);
        RUN_TEST(transform_test(), b, s);
        RUN_TEST(view_test(), b, s);
        RUN_TEST(broadcast_test(), b, s);
        RUN_TEST(extended_transform_test(), b, s);
        RUN_TEST(assign_transform_slice_alt(), b ,s)
        RETURN_TESTRESULT(b, s);
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

