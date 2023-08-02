/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHINITIALGUESS_HPP_
#define _AHINITIALGUESS_HPP_

#include "AlwaysInline.hpp"
#include "DimensionDefinitions.hpp"

#include <cmath>  // sin, cos, sqrt
#include <memory> // std::shared_ptr

class AHInitialGuessDefault
{
  public:
#if CH_SPACEDIM == 3
    virtual double get(double u, double v) const = 0;
    ALWAYS_INLINE virtual double get_merger_contribution(double u,
                                                         double v) const
    {
        return get(u, v);
    };
    ALWAYS_INLINE virtual double get_merger_min_distance() const
    {
        return get(0., 0.);
    };
#elif CH_SPACEDIM == 2
    virtual double get(double u) const = 0;
    ALWAYS_INLINE virtual double get_merger_contribution(double u) const
    {
        return get(u);
    };
    ALWAYS_INLINE virtual double get_merger_min_distance() const
    {
        return get(0.);
    };
#endif
};

using AHInitialGuessPtr = std::shared_ptr<AHInitialGuessDefault>;

/////////////////////////////////////////////////////////
// Predefined initial guess functions
/////////////////////////////////////////////////////////

class AHInitialGuessConstant : public AHInitialGuessDefault
{
  public:
    double m_initial_guess;

    AHInitialGuessConstant() {}
    AHInitialGuessConstant(double a_initial_guess)
    {
        set_params(a_initial_guess);
    }

    void set_params(double a_initial_guess)
    {
        m_initial_guess = a_initial_guess;
        pout() << "Setting Original Guess to f=" << m_initial_guess
               << std::endl;
    }

#if CH_SPACEDIM == 3
    ALWAYS_INLINE double get(double u, double v) const override
    {
        return m_initial_guess;
    }
#elif CH_SPACEDIM == 2
    ALWAYS_INLINE double get(double u) const override
    {
        return m_initial_guess;
    }
#endif
};

class AHInitialGuessMerger : public AHInitialGuessDefault
{
  public:
    const AHInitialGuessPtr m_ah1, m_ah2;

    double m_merger_pre_factor; // see notes in AHParams
    double m_merger_search_factor;

    AHInitialGuessMerger(const AHInitialGuessPtr a_ah1,
                         const AHInitialGuessPtr a_ah2,
                         double a_merger_pre_factor,
                         double a_merger_search_factor)
        : m_ah1(a_ah1), m_ah2(a_ah2), m_merger_pre_factor(a_merger_pre_factor),
          m_merger_search_factor(a_merger_search_factor)
    {
    }

    ALWAYS_INLINE double get_merger_min_distance() const override
    {
        return m_merger_search_factor * 4. *
               (m_ah1->get_merger_min_distance() +
                m_ah2->get_merger_min_distance());
    };

#if CH_SPACEDIM == 3
    ALWAYS_INLINE double get(double u, double v) const override
    {
        return m_merger_pre_factor * 4. * get_merger_contribution(u, v);
    }
    ALWAYS_INLINE double get_merger_contribution(double u,
                                                 double v) const override
    {
        return m_ah1->get_merger_contribution(u, v) +
               m_ah2->get_merger_contribution(u, v);
    }
#elif CH_SPACEDIM == 2
    ALWAYS_INLINE double get(double u) const override
    {
        return m_merger_pre_factor * 4. * get_merger_contribution(u);
    }
    ALWAYS_INLINE double get_merger_contribution(double u) const override
    {
        return m_ah1->get_merger_contribution(u) +
               m_ah2->get_merger_contribution(u);
    }
#endif
};

// ellipsoid aligned with one of the axis
class AHInitialGuessEllipsoid : public AHInitialGuessDefault
{
  public:
    double D_DECL(m_radius_x, m_radius_y, m_radius_z);

    AHInitialGuessEllipsoid() {}

    AHInitialGuessEllipsoid(D_DECL(double a_radius_x, double a_radius_y,
                                   double a_radius_z))
    {
        set_params(D_DECL(a_radius_x, a_radius_y, a_radius_z));
    }

    void set_params(D_DECL(double a_radius_x, double a_radius_y,
                           double a_radius_z))
    {
        D_TERM(m_radius_x = a_radius_x;, m_radius_y = a_radius_y;
               , m_radius_z = a_radius_z;)
        pout() << "Setting Original Guess to ellipsoid=(" D_TERM(
                      << m_radius_x, << ", " << m_radius_y,
                      << ", " << m_radius_z)
               << ")" << std::endl;
    }

#if CH_SPACEDIM == 3
    ALWAYS_INLINE double get_merger_min_distance() const override
    {
        return std::max(m_radius_x, std::max(m_radius_y, m_radius_z));
    };

    ALWAYS_INLINE double get(double u, double v) const override
    {
        double sin_u = sin(u);
        double cos_u = cos(u);
        double sin_v = sin(v);
        double cos_v = cos(v);
        double radius =
            1. /
            sqrt(sin_u * sin_u * cos_v * cos_v / (m_radius_x * m_radius_x) +
                 sin_u * sin_u * sin_v * sin_v / (m_radius_y * m_radius_y) +
                 cos_u * cos_u / (m_radius_z * m_radius_z));
        // pout() << u << "|" << v << "|" << radius << std::endl;
        return radius;
    }
    ALWAYS_INLINE double get_merger_contribution(double u,
                                                 double v) const override
    {
        return get_merger_min_distance(); // use constant
        // return get(u, v); // use ellipsoid
    }

#elif CH_SPACEDIM == 2
    ALWAYS_INLINE double get_merger_min_distance() const override
    {
        return std::max(m_radius_x, m_radius_y);
    };

    ALWAYS_INLINE double get(double u) const override
    {
        double sin_u = sin(u);
        double cos_u = cos(u);
        // u = 0 corresponds to (x,y)=(-1,0)
        // u = pi/2 corresponds to (x,y)=(0,1)
        // u = pi corresponds to (x,y)=(1,0)
        double radius = 1. / sqrt(cos_u * cos_u / (m_radius_x * m_radius_x) +
                                  sin_u * sin_u / (m_radius_y * m_radius_y));
        // pout() << u << "|" << radius << std::endl;
        return radius;
    }
    ALWAYS_INLINE double get_merger_contribution(double u) const override
    {
        return get_merger_min_distance(); // use constant
        // return get(u); // use ellipsoid
    }
#endif
};

#endif /* _AHINITIALGUESS_HPP_ */