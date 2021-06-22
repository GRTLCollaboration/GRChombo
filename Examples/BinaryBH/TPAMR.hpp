/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef TPAMR_HPP_
#define TPAMR_HPP_

// Even if USE_TWOPUNCTURES is not defined, this file will include BHAMR.hpp
#include "BHAMR.hpp"

#ifdef USE_TWOPUNCTURES
#include "TwoPunctures.hpp"

/// A descendent of Chombo's AMR class to interface with tools which require
/// access to the whole AMR hierarchy, and those of GRAMR
/**
 * This object inherits from BHAMR and adds members relevant to TwoPunctures
 * initial data
 */
class TPAMR : public BHAMR
{
  public:
    TP::TwoPunctures m_two_punctures;

    void set_two_punctures_parameters(const TP::Parameters &params)
    {
        // explicitly invoke copy constructor of base Parameters class
        m_two_punctures.Parameters::operator=(params);
    }
};

#endif /* USE_TWOPUNCTURES */

#endif /* TPAMR_HPP_ */
