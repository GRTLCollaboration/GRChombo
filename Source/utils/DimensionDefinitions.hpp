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

#define FOR1(IDX) for (int IDX = 0; IDX < DEFAULT_TENSOR_DIM; ++IDX)
#define FOR2(IDX1, IDX2) FOR1(IDX1) FOR1(IDX2)
#define FOR3(IDX1, IDX2, IDX3) FOR2(IDX1, IDX2) FOR1(IDX3)
#define FOR4(IDX1, IDX2, IDX3, IDX4) FOR2(IDX1, IDX2) FOR2(IDX3, IDX4)

#endif /* GRUTILS_HPP_*/
