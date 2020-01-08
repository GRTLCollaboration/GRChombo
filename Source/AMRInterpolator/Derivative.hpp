/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef DERIVATIVE_HPP_
#define DERIVATIVE_HPP_

#include <array>
#include <string>

class Derivative : public std::array<int, CH_SPACEDIM>
{
  private:
    Derivative(int d)
    {
        for (int i = 0; i < CH_SPACEDIM; ++i)
        {
            (*this)[i] = 0;
        }

        (*this)[d] = 1;
    }

    Derivative(int d1, int d2)
    {
        for (int i = 0; i < CH_SPACEDIM; ++i)
        {
            (*this)[i] = 0;
        }

        (*this)[d1] += 1;
        (*this)[d2] += 1;
    }

  public:
    Derivative()
    {
        for (int i = 0; i < CH_SPACEDIM; ++i)
        {
            std::array<int, CH_SPACEDIM>::operator[](i) = 0;
        }
    }

    // Ordering for std::map

    bool operator==(const Derivative &rhs) const
    {
        for (int i = 0; i < CH_SPACEDIM; ++i)
        {
            if ((*this)[i] != rhs[i])
            {
                return false;
            }
        }

        return true;
    }

    bool operator!=(const Derivative &deriv) const
    {
        if ((*this) == deriv)
            return false;
        else
            return true;
    }

    bool operator<(const Derivative &rhs) const
    {
        int derivs = 0;
        int rhs_derivs = 0;

        for (int i = 0; i < CH_SPACEDIM; ++i)
        {
            derivs += (*this)[i];
            rhs_derivs += rhs[i];
        }

        if (derivs < rhs_derivs)
        {
            return true;
        }
        else if (derivs > rhs_derivs)
        {
            return false;
        }
        else
        {
            for (int i = 0; i < CH_SPACEDIM; ++i)
            {
                // This is counterintuitive but is actually the ordering the we
                // want in order to generalise to arbitrary #dims
                if ((*this)[i] > rhs[i])
                {
                    return true;
                }
                else if ((*this)[i] < rhs[i])
                {
                    return false;
                }
            }

            return false;
        }
    }

    static const Derivative LOCAL;

    static const Derivative dx;
    static const Derivative dy;
    static const Derivative dz;

    static const Derivative dxdx;
    static const Derivative dydy;
    static const Derivative dzdz;

    static const Derivative dxdy;
    static const Derivative dxdz;
    static const Derivative dydz;

    static std::string name(const Derivative &deriv)
    {
        if (deriv == dx)
            return "dx";
        else if (deriv == dy)
            return "dy";
        else if (deriv == dz)
            return "dz";
        else if (deriv == dxdx)
            return "dxdx";
        else if (deriv == dydy)
            return "dydy";
        else if (deriv == dzdz)
            return "dzdz";
        else if (deriv == dxdy)
            return "dxdy";
        else if (deriv == dxdz)
            return "dxdz";
        else if (deriv == dydz)
            return "dydz";
        else
            return "";
    }
};

/* Moved to DerivativeSetup.hpp as otherwise multiply
 * defined in the various translaion units

const Derivative Derivative::LOCAL;

const Derivative Derivative::dx(0);
const Derivative Derivative::dy(1);
const Derivative Derivative::dz(2);

const Derivative Derivative::dxdx(0, 0);
const Derivative Derivative::dydy(1, 1);
const Derivative Derivative::dzdz(2, 2);

const Derivative Derivative::dxdy(0, 1);
const Derivative Derivative::dxdz(0, 2);
const Derivative Derivative::dydz(1, 2);
*/
#endif /* DERIVATIVE_HPP_ */
