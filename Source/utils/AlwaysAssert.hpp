/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef ALWAYS_ASSERT_HPP_
#define ALWAYS_ASSERT_HPP_

#include "MayDay.H"

// Chombo namespace
#include "UsingNamespace.H"

// This is the same as Chombo's CH_assert macro except it also works when NDEBUG
// is defined (i.e. when OPT = HIGH)

#define always_assert_str(s) #s
#define always_assert_xstr(s) always_assert_str(s)
#define always_assert(cond)                                                    \
    if (!(cond))                                                               \
    {                                                                          \
        MayDay::Abort(__FILE__ ":" always_assert_xstr(                         \
            __LINE__) ": Assertion `" #cond "' failed.");                      \
    }

#endif /* ALWAYS_ASSERT_HPP_ */
