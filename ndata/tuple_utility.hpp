#ifndef TUPLE_UTILITY_HPP_KC5LBUVV
#define TUPLE_UTILITY_HPP_KC5LBUVV

#include <tuple>
#include <utility>

namespace tuple_utility {

    //from http://stackoverflow.com/questions/10626856/how-to-split-a-tuple
    template <typename T , typename... Ts >
    auto head(std::tuple<T,Ts...> & t )
    {
       return  std::get<0>(t);
    }

    template < std::size_t... Ns , typename... Ts >
    auto tail_impl( std::index_sequence<Ns...> , std::tuple<Ts...> & t ) {
       return  std::make_tuple( std::get<Ns+1u>(t)... );
    }

    template < typename... Ts >
    auto tail(std::tuple<Ts...> & t)
    {
       return  tail_impl( std::make_index_sequence<sizeof...(Ts) - 1u>() , t );
    }

    //non copying tuple splitting
    template <size_t ... Ns, typename ... Ts>
    auto split_ht_impl(std::index_sequence<Ns...>, std::tuple<Ts...>&& t) {
        return make_pair(
                    std::get<0>(t),
                    std::make_tuple(std::move(std::get<Ns+1>(t))...)
                    );
    }

    //non copying tuple splitting
    template <typename ... Ts>
    auto split_ht(std::tuple<Ts...>&& t)
    {
        //std::pair<
        //        std::decay<decltype(std::get<0>(t))>,
        //        decltype(tail(t))
        //        >
        //        ret =
        return split_ht_impl(
            std::make_index_sequence<sizeof...(Ts) - 1u>(),
            std::move(t)
            );
        //return ret;
    }

    template <typename FuncT, typename ... Ts, size_t... Is>
    auto apply_impl(FuncT func, std::tuple<Ts...> & t,
                      std::index_sequence<Is...>)
    {
       //note difference in position of expansion between apply and tuple transform
       return  func(std::get<Is>(t)...);
    }

    //dereferences pointers if tuple of pointers is passed
    //so func can be written to take by value or byref instead of by ptr only
    template <typename FuncT, typename ... Ts, size_t... Is>
    auto apply_impl(FuncT func, std::tuple<Ts*...> & t,
                      std::index_sequence<Is...>)
    {
       //note difference in position of expansion between apply and tuple transform
       return  func(*std::get<Is>(t)...);
    }

    template <typename FuncT, typename TupT>
    auto apply(FuncT func, TupT & t)
    {
       return  apply_impl(
                   func,
                   t,
                   std::make_index_sequence<
                          std::tuple_size<
                          TupT
                          >::value
                          >()
                          );
    }

    template <typename FuncT, typename TupT, size_t... Is>
    auto tuple_transform_impl(FuncT func, TupT & t,
                              std::index_sequence<Is...>)
    {
       //note difference in position of expansion between apply and tuple transform
       return  std::make_tuple(func(std::get<Is>(t))...);
    }

    template <typename FuncT, typename TupT>
    auto
    tuple_transform(FuncT func, TupT& t)
    {
       return  tuple_transform_impl(
                   func,
                   t,
                   std::make_index_sequence<
                       std::tuple_size<
                           TupT
                       >::value
                   >()
                   );
    }

    //tuple zip implementation from from http://stackoverflow.com/a/11322096
    template <typename ...A, typename ...B,
              std::size_t ...S>
    auto zip_helper(std::tuple<A...> t1, std::tuple<B...> t2, std::index_sequence<S...> s)
    -> decltype(std::make_tuple(std::make_tuple(std::get<S>(t1),std::get<S>(t2))...))
    {
        return std::make_tuple( std::make_tuple( std::get<S>(t1), std::get<S>(t2) )...);
    }

    template <typename ...A, typename ...B>
    auto zip(std::tuple<A...> t1, std::tuple<B...> t2)
    //-> decltype(zip_helper(t1, t2, typename gens<sizeof...(A)>::type() ))
    {
        static_assert(sizeof...(A) == sizeof...(B), "The tuple sizes must be the same");
        return zip_helper( t1, t2, std::make_index_sequence<sizeof...(A)>() );
    }

}

#endif /* end of include guard: TUPLE_UTILITY_HPP_KC5LBUVV */ 
