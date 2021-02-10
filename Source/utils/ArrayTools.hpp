/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ARRAYTOOLS_HPP
#define ARRAYTOOLS_HPP

#include <algorithm>
#include <array>
#include <string>
#include <vector>

/// A place for tools that operate on std::arrays
namespace ArrayTools
{
/// This just concantenates two arrays together.
/// MR: I have no idea why this or something similar isn't in the standard
/// library
template <typename T, size_t N, size_t M>
std::array<T, N + M> concatenate(const std::array<T, N> &first,
                                 const std::array<T, M> &second)
{
    std::array<T, N + M> out;
    std::copy(first.cbegin(), first.cend(), out.begin());
    std::copy(second.cbegin(), second.cend(), out.begin() + N);
    return out;
}

template <typename T, size_t N,
          std::enable_if_t<std::is_arithmetic<T>::value, bool> = true>
T norm2(const std::array<T, N> &a_array)
{
    T out = 0;
    for (auto &elem : a_array)
    {
        out += elem * elem;
    }
    return out;
}

template <typename T, size_t N,
          std::enable_if_t<std::is_arithmetic<T>::value, bool> = true>
std::string to_string(const std::array<T, N> a_array)
{
    std::string out;
    for (int i = 0; i < N - 1; ++i)
    {
        out += std::to_string(a_array[i]) + " ";
    }
    out += std::to_string(a_array[N - 1]);
    return out;
}

// SFINAE for std::arrays and std::vectors
template <typename T> struct is_std_array_or_vector : std::false_type
{
};

template <typename elem_t>
struct is_std_array_or_vector<std::vector<elem_t>> : std::true_type
{
};

template <typename elem_t, std::size_t N>
struct is_std_array_or_vector<std::array<elem_t, N>> : std::true_type
{
};
} // namespace ArrayTools

#endif /* ARRAYTOOLS_HPP */
