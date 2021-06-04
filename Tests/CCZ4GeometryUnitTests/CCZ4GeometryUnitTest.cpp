/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "CCZ4Geometry.hpp"
#include "DimensionDefinitions.hpp"
#include "Tensor.hpp"
#include <iomanip>
#include <iostream>

template <class data_t> struct vars_t
{
    data_t chi;
    Tensor<2, data_t> h;
    Tensor<1, data_t> Gamma;
};

int main()
{
    int failed = 0;

    vars_t<double> vars;
    vars_t<Tensor<1, double>> d1;
    vars_t<Tensor<2, double>> d2;
    Tensor<1, double> Z_over_chi;

#include "values1.hpp" //Including the auto generated file with values

    auto h_UU = TensorAlgebra::compute_inverse_sym(vars.h);

    auto chris = TensorAlgebra::compute_christoffel(d1.h, h_UU);

    auto ricciZ =
        CCZ4Geometry::compute_ricci_Z(vars, d1, d2, h_UU, chris, Z_over_chi);

    std::cout << std::setprecision(16);

    // Compare
    FOR(i, j)
    {
        double diff = h_UU[i][j] - h_UU_known[i][j];
        if (diff > 1e-14)
        {
            std::cout << "h_UU wrong in component [" << i << "][" << j << "]"
                      << std::endl;
            failed = -1;
        }
    }

    FOR(i, j, k)
    {
        double diff = chris.ULL[i][j][k] - chris_known[i][j][k];
        if (diff > 1e-14)
        {
            std::cout << "chris wrong in component [" << i << "][" << j << "]["
                      << k << "]" << std::endl;
            std::cout << "value: " << chris.ULL[i][j][k] << std::endl;
            std::cout << "correct value: " << chris_known[i][j][k] << std::endl;
            failed = -1;
        }
    }

    FOR(i)
    {
        double diff = chris.contracted[i] - chris_contracted_known[i];
        if (diff > 1e-14)
        {
            std::cout << "chris contracted wrong in component [" << i << "]"
                      << std::endl;
            std::cout << "value: " << chris.contracted[i] << std::endl;
            std::cout << "correct value: " << chris_contracted_known[i]
                      << std::endl;
            failed = -1;
        }
    }

    FOR(i, j)
    {
        double diff = ricciZ.LL[i][j] - ricciZ_known[i][j];
        if (diff > 1e-14)
        {
            std::cout << "ricciZ contracted wrong in component [" << i << "]["
                      << j << "]" << std::endl;
            std::cout << "value: " << ricciZ.LL[i][j] << std::endl;
            std::cout << "correct value: " << ricciZ_known[i][j] << std::endl;
            failed = -1;
        }
    }

    double diff = ricciZ.scalar - ricciZ_scalar_known;
    if (diff > 1e-14)
    {
        std::cout << "ricci scalar wrong" << std::endl;
        std::cout << "value: " << ricciZ.scalar << std::endl;
        std::cout << "correct value: " << ricciZ_scalar_known << std::endl;
        failed = -1;
    }

    if (failed == 0)
        std::cout << "CCZ4Geometry test passed..." << std::endl;
    else
        std::cout << "CCZ4Geometry test failed..." << std::endl;

    return failed;
}
