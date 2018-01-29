/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INDEX_APPLICATOR_HPP_
#define INDEX_APPLICATOR_HPP_

#include <utility>

// The whole purpose of GetIndexTraits is to get around a limitation in c++
// whereby  variadic templates weren't allowed to be used inside decltype.  We
// could have gone to c++14 and removed the return type specification altogether
// but  GetIndexTraits allows one to get the return type of
// IndexApplicator::apply within c++11

template <typename data_t, typename... Ts> struct GetIndexTraits;

template <typename data_t> struct GetIndexTraits<data_t>
{
    using type = data_t;
};

template <typename data_t, typename T, typename... Ts>
struct GetIndexTraits<data_t, T, Ts...>
{
    using type = typename GetIndexTraits<
        decltype(std::declval<data_t>()[std::declval<T>()]), Ts...>::type;
};

// This class allows one to apply an arbitrary number of indices to an object
// of an arbitrary number of indices.
// This is useful for example in VarsBase where the IndexApplicator allows us to
// write general  code for an arbitrary number of derivative indices.

// Need a base template that we can specialise. Note that this will never be
// instantiated  because the specialisations cover everything.
template <typename... Ts> class IndexApplicator
{
};

// Specialisation for one or more indices
template <typename t, typename... Ts> class IndexApplicator<t, Ts...>
{
  public:
    template <typename data_t>
    static ALWAYS_INLINE typename GetIndexTraits<data_t &, t, Ts...>::type
    apply(data_t &obj, t dir0, Ts... dirs)
    {
        // Let the compiler iterate until there are no indices left (at which
        // point we go to the  specialisation below).
        return IndexApplicator<Ts...>::apply(obj[dir0], dirs...);
    }
};

// Specialisation for the case with no indices left
template <> class IndexApplicator<>
{
  public:
    template <typename data_t>
    static ALWAYS_INLINE typename GetIndexTraits<data_t &>::type
    apply(data_t &obj)
    {
        return obj;
    }
};

#endif /* INDEX_APPLICATOR_HPP_ */
