/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ARRAYTOOLS_HPP
#define ARRAYTOOLS_HPP

#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>

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
} // namespace ArrayTools

#endif /* ARRAYTOOLS_HPP */
