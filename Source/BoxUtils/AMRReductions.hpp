/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef AMRREDUCTIONS_HPP
#define AMRREDUCTIONS_HPP

// Chombo includes
#include "computeNorm.H"
#include "computeSum.H"

// Our includes
#include "GRAMR.hpp"
#include "GRAMRLevel.hpp"
#include "UserVariables.hpp"
#include "VariableType.hpp"

// Chombo namespace
#include "UsingNamespace.H"

//! A class that provides a user-friendly interface to Chombo's
//! computeSum and computeNorm functions
template <VariableType var_t> class AMRReductions
{
  private:
    static constexpr int m_num_vars =
        (var_t == VariableType::evolution)
            ? static_cast<int>(NUM_VARS)
            : static_cast<int>(NUM_DIAGNOSTIC_VARS);
    const int m_base_level;
    const double m_coarsest_dx;
    double m_domain_volume;
    Vector<LevelData<FArrayBox> *> m_level_data_ptrs;
    Vector<int> m_ref_ratios;

    //! constructs a Vector of LevelData<FArrayBox> pointers and stores them
    void set_level_data_vect(const GRAMR &a_gramr);

    //! gets the vector of refinement ratios and stores them
    void set_ref_ratios_vect(const GRAMR &a_gramr);

    //! Sets m_domain_volume which is used in norm() if a_normalize_by_volume
    //! is true. Must be called after set_level_data_vect
    void set_domain_volume();

  public:
    //! Constructor
    AMRReductions(const GRAMR &a_gramr, const int a_base_level = 0);

    //! returns the minimum of an interval of variables
    Real min(const Interval &a_vars) const;

    //! returns the minimum of a single variable
    Real min(const int a_var) const;

    //! returns the maximum of an interval of variables
    Real max(const Interval &a_vars) const;

    //! returns the minimum of a single variables
    Real max(const int a_var) const;

    //! returns the volume-weighted p-norm of an interval of variables
    //! p = a_norm_exponent
    Real norm(const Interval &a_vars, const int a_norm_exponent = 2,
              const bool a_normalize_by_volume = false) const;

    //! returns the volume weighted p-norm of a single variable
    //! p = a_norm_exponent
    Real norm(const int a_var, const int a_norm_exponent = 2,
              const bool a_normalize_by_volume = false) const;

    //! returns the volume-weighted sum (integral) of an interval of variables
    Real sum(const Interval &a_vars) const;

    //! returns the volume-weighted sum (integral of a single variable);
    Real sum(const int a_var) const;

    //! returns the m_domain_volume member
    Real get_domain_volume() const;
};

#include "AMRReductions.impl.hpp"

#endif /* AMRREDUCTIONS_HPP */
