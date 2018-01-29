/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef GRAMR_HPP_
#define GRAMR_HPP_

#include "AMR.H"

/// A child of Chombo's AMR class to interface with tools which require
/// access to the whole AMR hierarchy (such as the AMRInterpolator)
/**For now this class doesn't contain any functionality in mainline code.
 *However, it is necessary for many experimental features and allows us to
 *add said features later without breaking any user code.
 */
class GRAMR : public AMR
{
};

#endif /* GRAMR_HPP_ */
