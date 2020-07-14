/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef LAGRANGE_IMPL_HPP_
#define LAGRANGE_IMPL_HPP_

template <int Order>
const string Lagrange<Order>::TAG = "\x1b[36;1m[Lagrange]\x1b[0m ";

/* Finite difference weight generation algorithm
 *
 * Translated from the F77 code found in
 * Bengt Fornberg: "Calculation of Weights in Finite Difference Formulas"
 * SIAM Review 40, 3 (1998): 685-691
 *
 * Here we restrict ourselves to integer grid. Original code allows for
 * arbitrary grid locations specified by grid[i]. To recover the general
 * routine: - change type of c1,c2,c3 to 'double' - replace '0', 'i', 'j' in the
 * annotated lines with 'grid[0]', 'grid[i]', 'grid[j]'.
 */
template <int Order>
Lagrange<Order>::Stencil::Stencil(int width, int deriv, double dx,
                                  double point_offset)
    : m_width(width), m_deriv(deriv), m_dx(dx), m_point_offset(point_offset)
{
    int c1 = 1;
    int c2;
    int c3;
    double c4 = 0 - m_point_offset; /* replace for general grid */
    double c5;

    std::vector<double> tmp_weights(width * (deriv + 1), 0);

    tmp_weights[0] = 1.0;

    for (int i = 1; i < m_width; ++i)
    {
        int mn = (i < m_deriv) ? i : m_deriv;

        c2 = 1;
        c5 = c4;
        c4 = i - m_point_offset; /* replace for general grid */

        for (int j = 0; j < i; ++j)
        {
            c3 = i - j; /* replace for general grid */
            c2 = c2 * c3;

            if (j == i - 1)
            {
                for (int k = mn; k > 0; --k)
                {
                    tmp_weights[k * m_width + i] =
                        c1 *
                        (k * tmp_weights[(k - 1) * m_width + (i - 1)] -
                         c5 * tmp_weights[k * m_width + (i - 1)]) /
                        c2;
                }
                tmp_weights[i] = -c1 * c5 * tmp_weights[i - 1] / c2;
            }

            for (int k = mn; k > 0; --k)
            {
                tmp_weights[k * m_width + j] =
                    (c4 * tmp_weights[k * m_width + j] -
                     k * tmp_weights[(k - 1) * m_width + j]) /
                    c3;
            }
            tmp_weights[j] = c4 * tmp_weights[j] / c3;
        }

        c1 = c2;
    }

    if (deriv > 0)
    {
        const double dx_factor = pow(m_dx, deriv);
        m_weights.resize(m_width);
        for (int i = 0; i < m_width; ++i)
        {
            m_weights[i] = tmp_weights[m_width * m_deriv + i] / dx_factor;
        }
    }
    else
    {
        m_weights = std::move(tmp_weights);
    }

    if (false)
    {
        pout() << TAG << "Created a stencil for deriv " << m_deriv
               << " of width " << m_width << " for point " << m_point_offset
               << endl;
        pout() << "    Weights = { ";
        for (int i = 0; i < m_width; ++i)
        {
            pout() << m_weights[i] << " ";
        }
        pout() << "}" << endl;
    }
}

template <int Order>
bool Lagrange<Order>::Stencil::operator==(
    const Lagrange<Order>::Stencil &rhs) const
{
    return (rhs.m_width == m_width) && (rhs.m_deriv == m_deriv) &&
           (rhs.m_point_offset == m_point_offset) && (rhs.dx == m_dx);
}

template <int Order>
bool Lagrange<Order>::Stencil::isSameAs(int width, int deriv, double dx,
                                        double point_offset) const
{
    return (width == m_width) && (deriv == m_deriv) && (dx == m_dx) &&
           (point_offset == m_point_offset);
}

/* STENCIL ACCESSOR */

template <int Order>
typename Lagrange<Order>::Stencil
Lagrange<Order>::getStencil(int width, int deriv, double dx,
                            double point_offset)
{
    for (typename stencil_collection_t::iterator it =
             m_memoized_stencils.begin();
         it != m_memoized_stencils.end(); ++it)
    {
        if (it->isSameAs(width, deriv, dx, point_offset))
        {
            // Make a copy, lest std::vector decides to move our stencil during
            // growth op
            return *it;
        }
    }

    // We have to insert a new stencil.
    return *m_memoized_stencils.insert(m_memoized_stencils.end(),
                                       Stencil(width, deriv, dx, point_offset));
}

template <int Order>
const double &Lagrange<Order>::Stencil::operator[](unsigned int i) const
{
    CH_assert(i < m_width);
    return m_weights[i];
}

/* LAGRANGE TENSOR PRODUCT LOGIC */

template <int Order>
Lagrange<Order>::Lagrange(const InterpSource &source, bool verbosity)
    : m_source(source), m_verbosity(verbosity)
{
}

template <int Order>
void Lagrange<Order>::setup(const std::array<int, CH_SPACEDIM> &deriv,
                            const std::array<double, CH_SPACEDIM> &dx,
                            const std::array<double, CH_SPACEDIM> &evalCoord,
                            const IntVect &nearest)
{
    pair<std::vector<IntVect>, std::vector<double>> result =
        generateStencil(deriv, dx, evalCoord, nearest);
    m_interp_points = result.first;
    m_interp_weights = result.second;

    /*
    pout() << TAG << "Stencil: coord = { ";
    for (int i = 0; i < CH_SPACEDIM; ++i)
    {
        pout() << evalCoord[i] << " ";
    }
    pout() << "}, weights = { ";
    for (int i = 0; i < m_interp_weights.size(); ++i)
    {
        pout() << m_interp_weights[i] << " ";
    }
    pout() << "}" << endl;
    */
}

template <int Order>
double Lagrange<Order>::interpData(const FArrayBox &fab, int comp)
{
    /*
    m_interp_neg.clear();
    m_interp_pos.clear();

    // We are adding 200+ numbers at roughly the same magnitudes but alternating
    signs.
    // Let's keep track of positive and negative terms separately to make sure
    we don't run into trouble. for (int i = 0; i < m_interp_points.size(); ++i)
    {
        double data = m_interp_weights[i] * fab.get(m_interp_points[i], comp);
        if (data > 0)
        {
            m_interp_pos.insert(data);
        }
        else
        {
            m_interp_neg.insert(data);
        }
    }

    // Add positive terms from smallest to largest
    double pos = 0;
    for (typename multiset<double>::iterator it = m_interp_pos.begin(); it !=
    m_interp_pos.end(); ++it)
    {
        pos += *it;
    }

    // Largest negative term has smallest magnitude. Use reverse iterator.
    double neg = 0;
    for (typename multiset<double>::reverse_iterator it = m_interp_neg.rbegin();
    it != m_interp_neg.rend(); ++it)
    {
        neg += *it;
    }

    return pos + neg;
    */

    long double accum = 0.0;

    for (int i = 0; i < m_interp_points.size(); ++i)
    {
        double data = m_interp_weights[i] * fab.get(m_interp_points[i], comp);
        accum += data;
    }

    return accum;
}

template <int Order>
pair<std::vector<IntVect>, std::vector<double>>
Lagrange<Order>::generateStencil(
    const std::array<int, CH_SPACEDIM> &deriv,
    const std::array<double, CH_SPACEDIM> &dx,
    const std::array<double, CH_SPACEDIM> &evalCoord, const IntVect &nearest,
    int dim)
{
    std::vector<IntVect> out_points;
    std::vector<double> out_weights;

    /*
     * SCAN ALONG THIS DIRECTION TO FIND THE LARGEST CONTIGUOUS CHUNK OF VALID
     * POINTS
     */

    // Allocate a std::vector twice as big as we can possibly need
    // This way we insert to either direction without shifting/growing
    std::vector<int> my_points(2 * (Order + deriv[dim]));

    enum
    {
        DOWN,
        UP
    };

    bool can_grow[2] = {true, true};
    int points_min = Order + deriv[dim];
    int points_max = Order + deriv[dim];

    std::array<double, CH_SPACEDIM> interp_coord = evalCoord;
    int candidate = nearest[dim];
    int grown_direction = (nearest[dim] - evalCoord[dim] < 0) ? DOWN : UP;

    while ((can_grow[DOWN] || can_grow[UP]) &&
           (points_max - points_min < Order + deriv[dim]))
    {
        interp_coord[dim] = candidate;

        if (m_source.contains(interp_coord))
        {
            int idx =
                (grown_direction == DOWN) ? (--points_min) : (points_max++);
            my_points[idx] = candidate;
        }
        else
        {
            can_grow[grown_direction] = false;
        }

        // Flip "grow" direction if we can do so
        if (can_grow[1 - grown_direction])
        {
            grown_direction = 1 - grown_direction;
        }

        candidate = (grown_direction) ? (my_points[points_max - 1] + 1)
                                      : (my_points[points_min] - 1);
    }

    int stencil_width = points_max - points_min;
    CH_assert(stencil_width > 0);

    const Stencil my_weights =
        getStencil(stencil_width, deriv[dim], dx[dim],
                   evalCoord[dim] - my_points[points_min]);

    if (m_verbosity)
    {
        pout() << TAG << "Stencil: dim = " << dim
               << ", coord = " << evalCoord[dim] << ", points = { ";
        for (int i = points_min; i < points_max; ++i)
        {
            pout() << my_points[i] << " ";
        }
        pout() << "}, weights = { ";
        for (int i = 0; i < stencil_width; ++i)
        {
            pout() << my_weights[i] << " ";
        }
        pout() << "}" << endl;
    }

    // There is going to be potentially a LOT of temporary std::vectors getting
    // allocated in here. If things get slow this will be a good place to look
    // first.
    for (int i = 0; i < stencil_width; ++i)
    {
        interp_coord[dim] = my_points[i + points_min];

        if (dim > 0)
        {
            // Descend to the next dimension
            pair<std::vector<IntVect>, std::vector<double>> sub_result =
                generateStencil(deriv, dx, interp_coord, nearest, dim - 1);
            std::vector<IntVect> &sub_points = sub_result.first;
            std::vector<double> &sub_weights = sub_result.second;

            // Take tensor product weights
            for (int j = 0; j < sub_points.size(); ++j)
            {
                if (my_weights[i] != 0)
                {
                    out_points.push_back(sub_points[j]);
                    out_weights.push_back(my_weights[i] * sub_weights[j]);
                }
            }
        }
        else
        {
            // "Terminal" dimension, just push back our own stuff
            if (my_weights[i] != 0)
            {
                out_points.push_back(IntVect(D_DECL6(
                    interp_coord[0], interp_coord[1], interp_coord[2],
                    interp_coord[3], interp_coord[4], interp_coord[5])));
                out_weights.push_back(my_weights[i]);
            }
        }
    }

    return pair<std::vector<IntVect>, std::vector<double>>(
        std::move(out_points), std::move(out_weights));
}

#endif /* LAGRANGE_IMPL_HPP_ */
