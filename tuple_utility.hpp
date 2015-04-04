#ifndef TUPLE_UTILITY_HPP_KC5LBUVV
#define TUPLE_UTILITY_HPP_KC5LBUVV

#include <tuple>
#include <utility>

namespace tuple_utility {

    //from http://en.cppreference.com/w/cpp/experimental/apply
    //namespace detail {
    //    template <typename F, typename... Ts, std::size_t... I>
    //    auto apply_impl(F f, std::tuple<Ts...> t, std::index_sequence<I...>) {
    //      return f(std::get<I>(t)...);
    //    }
    //} // namespace detail

    //template <typename F, class Tuple>
    //auto apply(F f, Tuple t)
    //{
    //    return detail::apply_impl(f, t,
    //        std::make_index_sequence<std::tuple_size<decltype(t)>::value>());
    //}

    //http://stackoverflow.com/questions/10626856/how-to-split-a-tuple
    template <typename T , typename... Ts >
    auto head( std::tuple<T,Ts...> t )
    {
       return  std::get<0>(t);
    }

    template < std::size_t... Ns , typename... Ts >
    auto tail_impl( std::index_sequence<Ns...> , std::tuple<Ts...> t ) {
       return  std::make_tuple( std::get<Ns+1u>(t)... );
    }

    template < typename... Ts >
    auto tail( std::tuple<Ts...> t )
    {
       return  tail_impl( std::make_index_sequence<sizeof...(Ts) - 1u>() , t );
    }

    template <typename FuncT, typename TupT, size_t... Is>
    auto apply_impl(FuncT func, TupT &t,
                              std::index_sequence<Is...>)
    {
       //note difference in position of expansion between apply and tuple transform
       return  func(std::get<Is>(t)...);
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
    auto tuple_transform(FuncT func, TupT & t)
    {
       return  tuple_transform_impl(func, t,
                                    std::make_index_sequence<
                                    std::tuple_size<
                                    TupT
                                    >::value
                                    >()
                                    );
    }

}


#endif /* end of include guard: TUPLE_UTILITY_HPP_KC5LBUVV */ 
