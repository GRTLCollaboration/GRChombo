/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INTEGRATIONMETHOD_HPP_
#define INTEGRATIONMETHOD_HPP_

#include <utility>
#include <vector>

//! A class to store and return the weights associated to a Newton-Cotes formula
//! for numerical integration/quadrature which can be closed (i.e. includes the
//! endpoints) or open (does not include the end points). This is used by
//! SurfaceExtraction for integration over extraction surfaces.
class IntegrationMethod
{
  private:
    std::vector<double> m_weights;
    int m_num_weights;
    bool m_is_closed;

  public:
    //! Constructor
    IntegrationMethod(const std::vector<double> &a_weights,
                      bool a_is_closed = true)
        : m_weights(a_weights), m_num_weights(a_weights.size()),
          m_is_closed(a_is_closed)
    {
    }

    //! Checks that this integration method is suitable given the number of
    //! points and periodicity
    inline bool is_valid(int a_num_points, bool a_is_periodic) const
    {
        if (m_is_closed && !a_is_periodic)
        {
            return (a_num_points % m_num_weights == 1 % m_num_weights);
        }
        else
        {
            return (a_num_points % m_num_weights == 0);
        }
    }

    //! Returns whether this IntegrationMethod is closed or not
    inline bool is_closed() const { return m_is_closed; }

    //! Returns the weight for a point with given index
    inline double weight(int a_index, int a_num_points,
                         bool a_is_periodic) const
    {
        const int weight_index = a_index % m_num_weights;
        const bool endpoint =
            (a_index == 0 || a_index == a_num_points - 1) && !a_is_periodic;
        if (m_is_closed && !endpoint && weight_index == 0)
            return 2.0 * m_weights[weight_index];
        else
            return m_weights[weight_index];
    }

    static const IntegrationMethod trapezium;
    static const IntegrationMethod midpoint;
    static const IntegrationMethod simpson;
    static const IntegrationMethod boole;
};

#endif /* INTEGRATIONMETHOD_HPP_ */
