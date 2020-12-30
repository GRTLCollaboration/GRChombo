/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DEBUGGINGTOOLS_HPP_
#define DEBUGGINGTOOLS_HPP_

// Chombo includes
#include "parstream.H" //Gives us pout()
#ifdef EQUATION_DEBUG_MODE
#include "IntVect.H"
#endif

// Other includes
#include <string.h>

// Chombo namespace
#include "UsingNamespace.H"

/// This file contains a collection of helpful #defines and other definitions
/// that are hepful for debugging.

// Unfortunately, most of the functionality can only be achieved with macros
// (e.g. including the variable and filename).

#define __FILENAME__                                                           \
    (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)

#define DEBUG_SHOW(VAR) pout() << #VAR << ": " << VAR << " "
#define DEBUG_FILE pout() << __FILENAME__ << ": "
#define DEBUG_END pout() << std::endl
#define DEBUG_DOUBLE_PRECISION pout() << std::setprecision(16)

/// The macros DEBUG_OUT make debugging quicker and allow easy printing of a
/// variable.
#define DEBUG_OUT(VAR)                                                         \
    DEBUG_FILE;                                                                \
    DEBUG_SHOW(VAR);                                                           \
    DEBUG_END
#define DEBUG_OUT2(VAR1, VAR2)                                                 \
    DEBUG_FILE;                                                                \
    DEBUG_SHOW(VAR1);                                                          \
    DEBUG_SHOW(VAR2);                                                          \
    DEBUG_END
#define DEBUG_OUT3(VAR1, VAR2, VAR3)                                           \
    DEBUG_FILE;                                                                \
    DEBUG_SHOW(VAR1);                                                          \
    DEBUG_SHOW(VAR2);                                                          \
    DEBUG_SHOW(VAR3);                                                          \
    DEBUG_END
#define DEBUG_OUT4(VAR1, VAR2, VAR3, VAR4)                                     \
    DEBUG_FILE;                                                                \
    DEBUG_SHOW(VAR1);                                                          \
    DEBUG_SHOW(VAR2);                                                          \
    DEBUG_SHOW(VAR3);                                                          \
    DEBUG_SHOW(VAR4);                                                          \
    DEBUG_END

//
#ifdef EQUATION_DEBUG_MODE
#define DEBUG_HEADER                                                           \
    pout() << "Debug output in " << __FILENAME__                               \
           << " at: " << s_current_integer_coords << "." << std::endl
static IntVect s_current_integer_coords;
namespace EquationDebugging
{
inline void check_no_omp()
{
#ifdef _OPENMP
    if (omp_get_max_threads() > 1)
        MayDay::Error("Equation debug mode can  only be used with one thread.");
#endif
}

inline void set_global_cell_coordinates(const IntVect current_integer_coords)
{
    check_no_omp();
    s_current_integer_coords = current_integer_coords;
}
} // namespace EquationDebugging
#endif

#endif /* DEBUGGINGTOOLS_HPP_ */
