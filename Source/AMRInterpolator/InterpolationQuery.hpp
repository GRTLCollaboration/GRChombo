/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INTERPOLATIONQUERY_HPP_
#define INTERPOLATIONQUERY_HPP_

#include "Derivative.hpp"

class InterpolationQuery
{
  public:
    typedef pair<int, double *> out_t;
    typedef map<Derivative, std::vector<out_t>> comp_map_t;
    typedef typename map<Derivative, std::vector<out_t>>::iterator iterator;

    const int m_num_points;

  private:
    template <typename InterpAlgo> friend class AMRInterpolator;

    std::vector<const double *> m_coords;
    comp_map_t m_comps;

  public:
    InterpolationQuery(int num_points)
        : m_num_points(num_points), m_coords(CH_SPACEDIM, NULL)
    {
    }

    InterpolationQuery &setCoords(int dim, const double *coords)
    {
        CH_assert(dim < CH_SPACEDIM);
        this->m_coords[dim] = coords;
        return *this;
    }

    InterpolationQuery &addComp(int comp, double *out_ptr,
                                const Derivative &deriv = Derivative::LOCAL)
    {
        CH_assert(out_ptr != NULL || m_num_points == 0);

        iterator result = m_comps.find(deriv);
        if (result == m_comps.end())
        {
            result = m_comps
                         .insert(pair<Derivative, std::vector<out_t>>(
                             deriv, std::vector<out_t>()))
                         .first;
        }

        result->second.push_back(out_t(comp, out_ptr));
        return *this;
    }

    InterpolationQuery &clearComps()
    {
        m_comps.clear();
        return *this;
    }

    inline int numComps()
    {
        int accum = 0;

        for (iterator it = m_comps.begin(); it != m_comps.end(); ++it)
        {
            accum += it->second.size();
        }

        return accum;
    }

    inline iterator compsBegin() { return m_comps.begin(); }

    inline iterator compsEnd() { return m_comps.end(); }
};

#endif /* INTERPOLATIONQUERY_HPP_ */
