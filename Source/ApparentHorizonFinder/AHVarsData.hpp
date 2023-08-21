/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef _AHDATA_HPP_
#define _AHDATA_HPP_

#include "Derivative.hpp"
#include "DimensionDefinitions.hpp" // make sure GR_SPACEDIM exists
#include "InterpolationQuery.hpp"
#include "Tensor.hpp"

#include <cstring> // memset

template <class key, class data_t> struct AHVarsData
{
    // vars user might want to export
    std::map<key, data_t> vars;
    std::map<key, std::array<data_t, CH_SPACEDIM>> d1;
    std::map<key, std::array<data_t, CH_SPACEDIM *(CH_SPACEDIM + 1) / 2>> d2;

    void set_vars(InterpolationQuery &query, int var, key k, VariableType type,
                  int total)
    {
        vars[k].resize(total);
        query.addComp(var, &vars[k][0], Derivative::LOCAL, type);
    }
    void set_d1(InterpolationQuery &query, int var, key k, VariableType type,
                int total)
    {
        for (int j = 0; j < CH_SPACEDIM; ++j)
            d1[k][j].resize(total);

        query.addComp(var, &d1[k][0][0], Derivative::dx, type)
            .addComp(var, &d1[k][1][0], Derivative::dy, type);
#if CH_SPACEDIM == 3
        query.addComp(var, &d1[k][2][0], Derivative::dz, type);
#endif
    }
    void set_d2(InterpolationQuery &query, int var, key k, VariableType type,
                int total)
    {
        for (int j = 0; j < CH_SPACEDIM * (CH_SPACEDIM + 1) / 2; ++j)
            d2[k][j].resize(total);

        query.addComp(var, &d2[k][0][0], Derivative::dxdx, type)
            .addComp(var, &d2[k][1][0], Derivative::dxdy, type);
#if CH_SPACEDIM == 3
        query.addComp(var, &d2[k][2][0], Derivative::dxdz, type)
            .addComp(var, &d2[k][3][0], Derivative::dydy, type)
            .addComp(var, &d2[k][4][0], Derivative::dydz, type)
            .addComp(var, &d2[k][5][0], Derivative::dzdz, type);
#elif CH_SPACEDIM == 2
        query.addComp(var, &d2[k][2][0], Derivative::dydy, type);
#endif
    }
};

template <class key>
AHVarsData<key, double>
get_AHVarsData_idx(int idx, const AHVarsData<key, std::vector<double>> &a_data)
{
    AHVarsData<key, double> data;

    for (auto &vec : a_data.vars)
        data.vars[vec.first] = (vec.second[idx]);

    for (auto &vec : a_data.d1)
    {
        data.d1[vec.first] = {0.};
        for (int i = 0; i < CH_SPACEDIM; ++i)
            data.d1[vec.first][i] = vec.second[i][idx];
    }

    for (auto &vec : a_data.d2)
    {
        data.d2[vec.first] = {0.};
        for (int i = 0; i < CH_SPACEDIM * (CH_SPACEDIM + 1) / 2; ++i)
            data.d2[vec.first][i] = vec.second[i][idx];
    }

    return data;
}

#endif /* _AHDATA_HPP_ */
