/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GRUTILS_HPP_
#define GRUTILS_HPP_

#ifndef GR_SPACEDIM
#define GR_SPACEDIM 3
#endif

#ifndef DEFAULT_TENSOR_DIM
#define DEFAULT_TENSOR_DIM CH_SPACEDIM
#endif

// Fancy 'for' loop macros to iterate through spatial tensors
// use as "FOR(i, j) { ... }"
#define FOR1(IDX) for (int IDX = 0; IDX < DEFAULT_TENSOR_DIM; ++IDX)
#define FOR2(IDX1, IDX2) FOR1(IDX1) FOR1(IDX2)
#define FOR3(IDX1, IDX2, IDX3) FOR2(IDX1, IDX2) FOR1(IDX3)
#define FOR4(IDX1, IDX2, IDX3, IDX4) FOR2(IDX1, IDX2) FOR2(IDX3, IDX4)
#define FOR5(IDX1, IDX2, IDX3, IDX4, IDX5)                                     \
    FOR4(IDX1, IDX2, IDX3, IDX4) FOR1(IDX5)
#define DUMMYFOR() // prevents warning that appeared in debug mode:
                   // 'ISO C++11 requires at least one argument for the "..." in
                   // a variadic macro'

#define GET_MACRO6(_1, _2, _3, _4, _5, NAME, ...) NAME
#define FOR(...)                                                               \
    GET_MACRO6(__VA_ARGS__, FOR5, FOR4, FOR3, FOR2, FOR1, DUMMYFOR)(__VA_ARGS__)

#endif /* GRUTILS_HPP_*/
