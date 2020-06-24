/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef AMRINTERPOLATOR_HPP_
#define AMRINTERPOLATOR_HPP_

// Chombo includes

#include "AMR.H"
#include "AMRLevel.H"

#include "UsingNamespace.H"

// Our includes

#include "InterpSource.hpp"
#include "InterpolationAlgorithm.hpp"
#include "InterpolationLayout.hpp"
#include "InterpolationQuery.hpp"

#include "MPIContext.hpp"
#include "UserVariables.hpp"

// End include

template <typename InterpAlgo> class AMRInterpolator
{
  public:
    AMRInterpolator(const AMR &amr,
                    const std::array<double, CH_SPACEDIM> &coarsest_origin,
                    const std::array<double, CH_SPACEDIM> &coarsest_dx,
                    int verbosity = 0);
    void refresh();
    void limit_num_levels(unsigned int num_levels);
    void interp(InterpolationQuery &query);

  private:
    void computeLevelLayouts();
    InterpolationLayout findBoxes(InterpolationQuery &query);

    void prepareMPI(InterpolationQuery &query,
                    const InterpolationLayout layout);
    void exchangeMPIQuery();
    void calculateAnswers(InterpolationQuery &query);
    void exchangeMPIAnswer();

    const AMR &m_amr;

    // Coordinates of the point represented by IntVect::Zero in coarsest grid
    const std::array<double, CH_SPACEDIM> m_coarsest_origin;

    // Grid spacing in each direction
    const std::array<double, CH_SPACEDIM> m_coarsest_dx;

    int m_num_levels;
    const int m_verbosity;

    std::vector<std::array<double, CH_SPACEDIM>> m_origin;
    std::vector<std::array<double, CH_SPACEDIM>> m_dx;

    MPIContext m_mpi;
    std::vector<int> m_mpi_mapping;

    // Memoisation of boxes previously found
    std::vector<int> m_mem_level;
    std::vector<int> m_mem_box;

    std::vector<int> m_query_level;
    std::vector<int> m_query_box;
    std::vector<double> m_query_coords[CH_SPACEDIM];
    std::vector<std::vector<double>> m_query_data;

    std::vector<int> m_answer_level;
    std::vector<int> m_answer_box;
    std::vector<double> m_answer_coords[CH_SPACEDIM];
    std::vector<std::vector<double>> m_answer_data;

    // A bit of Android-ism here, but it's really useful!
    // Identifies the printout as originating from this class.
    const static string TAG;
};

#include "AMRInterpolator.impl.hpp"

#endif /* AMRINTERPOLATOR_HPP_ */
