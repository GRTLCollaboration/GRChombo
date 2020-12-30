/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DEFAULTLEVELFACTORY_HPP_
#define DEFAULTLEVELFACTORY_HPP_

// Chombo includes:
#include "AMRLevelFactory.H"

// Our includes
#include "GRAMR.hpp"
#include "SimulationParameters.hpp"

// Chombo namespace
#include "UsingNamespace.H"

template <class level_t> class DefaultLevelFactory : public AMRLevelFactory
{
  public:
    DefaultLevelFactory(GRAMR &gr_amr, SimulationParameters &a_sim_params)
        : m_gr_amr(gr_amr), m_p(a_sim_params)
    {
    }

    virtual AMRLevel *new_amrlevel() const
    {
        level_t *level_ptr = new level_t(m_gr_amr, m_p, m_p.verbosity);
        level_ptr->initialDtMultiplier(m_p.dt_multiplier);
        return (static_cast<AMRLevel *>(level_ptr));
    }

    virtual ~DefaultLevelFactory() {}

  protected:
    GRAMR &m_gr_amr;
    SimulationParameters m_p;
};
#endif /* DEFAULTLEVELFACTORY_HPP_ */
