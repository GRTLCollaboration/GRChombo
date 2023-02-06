/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef COORDINATES_HPP_
#define COORDINATES_HPP_

// Chombo includes
#include "IntVect.H"

// Other includes
#include "DimensionDefinitions.hpp"
#include "simd.hpp"
#include <array>

// Chombo namespace
#include "UsingNamespace.H"

template <class data_t> class Coordinates
{
  public:
    data_t x; // We vectorise over x so we must allow x to be a vector
    double y;
    double z;
    std::array<double, CH_SPACEDIM> m_center;

    Coordinates(IntVect integer_coords, double dx,
                std::array<double, CH_SPACEDIM> center = {0})
        : m_center(center)
    {
        compute_coord(x, integer_coords[0], dx, center[0]);

// The below code allows for 2D Cartoon reduction:
#if DEFAULT_TENSOR_DIM == CH_SPACEDIM && CH_SPACEDIM == 3
        compute_coord(y, integer_coords[1], dx, center[1]);
        compute_coord(z, integer_coords[2], dx, center[2]);
#elif DEFAULT_TENSOR_DIM == CH_SPACEDIM + 1 && CH_SPACEDIM == 2
        y = 0;
        compute_coord(z, integer_coords[1], dx, center[1]);
#else
#ifdef CH_SPACEDIM
#error compute_coord has not got your dimension combination implemented.
#endif
#endif
    }

    ALWAYS_INLINE static void compute_coord(double &out, int position,
                                            double dx,
                                            double center_distance = 0)
    {
        out = (position + 0.5) * dx - center_distance;
    }

    static void // typename std::enable_if_t<(simd_traits<double>::simd_len >
                // 1), void>
    compute_coord(simd<double> &out, int position, double dx,
                  double center_distance = 0)
    {
        double out_arr[simd_traits<double>::simd_len];
        for (int i = 0; i < simd_traits<double>::simd_len; ++i)
        {
            out_arr[i] = (position + i + 0.5) * dx - center_distance;
        }
        out = simd<double>::load(out_arr);
    }

    /// This function returns the radius subject to a floor for a given
    /// Coordinates object.
    data_t get_radius() const
    {
        // Note that this is not currently dimension independent
        data_t r = sqrt(x * x + y * y + z * z);

        const double minimum_r = 1e-6;
        return simd_max(r, minimum_r);
    }

    /// This static function returns the radius subject to a floor
    /// for when no coordinates object exists.
    static data_t get_radius(IntVect integer_coords, double dx,
                             std::array<double, CH_SPACEDIM> center = {0})
    {
        data_t xx;
        double yy;
        double zz;

        // Note that this is not currently dimension independent
        compute_coord(xx, integer_coords[0], dx, center[0]);
        compute_coord(yy, integer_coords[1], dx, center[1]);
        compute_coord(zz, integer_coords[2], dx, center[2]);

        data_t r = sqrt(xx * xx + yy * yy + zz * zz);

        const double minimum_r = 1e-6;
        return simd_max(r, minimum_r);
    }
};

template <typename data_t>
ALWAYS_INLINE ostream &operator<<(ostream &os,
                                  const Coordinates<data_t> &in_coords)
{
    os << "(x,y,z) = (" << in_coords.x << "," << in_coords.y << ","
       << in_coords.z << ")"
       << " r = " << in_coords.get_radius();
    return os;
}
#endif /* COORDINATES_HPP_ */
