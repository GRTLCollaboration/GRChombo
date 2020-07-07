/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "CoordinateTransformations.hpp"
#include "Coordinates.hpp"
#include "IntVect.H"
#include "TensorAlgebra.hpp"
#include "simd.hpp"

int main()
{
    int failed = 0;

    const double dx = 0.1;
    IntVect iv;
    iv[0] = 0;
    iv[1] = 0;
    iv[2] = 0;

    Coordinates<double> coords(iv, dx);
    const double x = coords.x;
    const double y = coords.y;
    const double z = coords.z;
    const double r = coords.get_radius();
    double rho2 = simd_max(x * x + y * y, 1e-12);
    double r2sintheta = sqrt(rho2) * r;
    
    std::cout << "x " << x << std::endl;
    std::cout << "y " << y << std::endl;
    std::cout << "z " << z << std::endl;
    std::cout << "r " << r << std::endl;

    using namespace TensorAlgebra;
    using namespace CoordinateTransformations;

    // Test if inv_jac is really the inverse of the jacobian
    Tensor<2, double> jac = jacobian(x,y,z);
    Tensor<2, double> inv_jac = inverse_jacobian(x,y,z);
    Tensor<2, double> inv_jac_check;
    FOR2(i, j)
      {
	inv_jac_check[i][j] = compute_inverse(jac)[i][j];
      }
    FOR2(i, j)
      {
        double diff = inv_jac_check[i][j] - inv_jac[i][j];
	if (diff > 1e-14)
	  {
	    std::cout << "Failed inverse_jacobian in component ["
		      << i << "][" << j << "]" << std::endl;
	    std::cout << "value: " << inv_jac[i][j] << std::endl;
	    std::cout << "correct value: " << inv_jac_check[i][j] << std::endl;
	    failed = -1;
	  }
      }

    // Test tensor transformations
    Tensor<2, double> Mij_cart;
    FOR2(i, j)
      {
	Mij_cart[i][j] = 0.;
      }
    Mij_cart[0][0] = 1.;
    Mij_cart[1][1] = 1.;
    Mij_cart[2][2] = 1.;

    Tensor<2, double> Mij_spher;
    FOR2(i, j)
      {
	Mij_spher[i][j] = 0.;
      }
    Mij_spher[0][0] = 1.;
    Mij_spher[1][1] = r * r;
    Mij_spher[2][2] = r2sintheta;

    // Test cartesian_to_spherical_LL
    Tensor<2, double> Mij_spher_check;
    Mij_spher_check = cartesian_to_spherical_LL(Mij_cart, x, y, z);

    FOR2(i, j)
    {
        double diff = Mij_spher_check[i][j] - Mij_spher[i][j];
	if (diff > 1e-14)
	  {
	    std::cout << "Failed cart_to_spher_LL transformation in component ["
		      << i << "][" << j << "]" << std::endl;
	    std::cout << "value: " << Mij_spher_check[i][j] << std::endl;
	    std::cout << "correct value: " << Mij_spher[i][j] << std::endl;
	    failed = -1;
	  }
    }

    // Test spherical_to_cartesian_LL
    Tensor<2, double> Mij_cart_check;
    Mij_cart_check = spherical_to_cartesian_LL(Mij_spher, x, y, z);

    FOR2(i, j)
    {
      double diff = Mij_cart_check[i][j] - Mij_cart[i][j];
      if (diff > 1e-14)
        {
	  std::cout << "Failed spher_to_cart_LL transformation in component ["
		    << i << "][" << j << "]" << std::endl;
	  std::cout << "value: " << Mij_cart_check[i][j] << std::endl;
	  std::cout << "correct value: " << Mij_cart[i][j] << std::endl;
	  failed = -1;
	}
    }

    // Test cartesian_to_spherical_UU
    Tensor<2, double> Mij_spher_UU;
    Tensor<2, double> Mij_spher_UU_check;
    Mij_spher_UU_check =
        cartesian_to_spherical_UU(compute_inverse(Mij_cart), x, y, z);
    Mij_spher_UU = compute_inverse(Mij_spher);

    FOR2(i, j)
    {
        double diff = Mij_spher_UU_check[i][j] - Mij_spher_UU[i][j];
        if (diff > 1e-14)
        {
	  std::cout << "Failed cart_to_spher_UU transformation in component ["
		    << i << "][" << j << "]" << std::endl;
	  std::cout << "value: " << Mij_spher_UU_check[i][j] << std::endl;
	  std::cout << "correct value: " << Mij_spher_UU[i][j] << std::endl;
	  failed = -1;
        }
    }

    // Test spherical_to_cartesian_UU
    Tensor<2, double> Mij_cart_UU;
    Tensor<2, double> Mij_cart_UU_check;
    Mij_cart_UU_check =
        spherical_to_cartesian_LL(compute_inverse(Mij_spher), x, y, z);
    Mij_cart_UU = compute_inverse(Mij_cart);

    FOR2(i, j)
    {
        double diff = Mij_cart_UU_check[i][j] - Mij_cart_UU[i][j];
        if (diff > 1e-14)
        {
	  std::cout << "Failed spher_to_cart_UU transformation in component ["
		    << i << "][" << j << "]" << std::endl;
	  std::cout << "value: " << Mij_cart_UU_check[i][j] << std::endl;
	  std::cout << "correct value: " << Mij_cart_UU[i][j] << std::endl;
	  failed = -1;
        }
    }

    // Test vector transformations

    Tensor<1, double> si_cart_U;
    si_cart_U[0] = x / r;
    si_cart_U[1] = y / r;
    si_cart_U[2] = z / r;

    Tensor<1, double> si_spher_U;
    si_spher_U[0] = 1.0;
    si_spher_U[1] = 0.0;
    si_spher_U[2] = 0.0;

    // Test cartesian_to_spherical_U
    Tensor<1, double> si_spher_U_check;
    si_spher_U_check = cartesian_to_spherical_U(si_cart_U, x, y, z);

    FOR1(i)
    {
        double diff = si_spher_U_check[i] - si_spher_U[i];
        if (diff > 1e-14)
        {
	  std::cout << "Failed cart_to_spher_U transformation in component ["
	            << i << "]" << std::endl;
	  std::cout << "value: " << si_spher_U_check[i] << std::endl;
	  std::cout << "correct value: " << si_spher_U[i] << std::endl;
	  failed = -1;
        }
    }

    // Test spherical_to_cartesian_U
    Tensor<1, double> si_cart_U_check;
    si_cart_U_check = spherical_to_cartesian_L(si_spher_U, x, y, z);

    FOR1(i)
    {
        double diff = si_cart_U_check[i] - si_cart_U[i];
        if (diff > 1e-14)
        {
	  std::cout << "Failed spher_to_cart_U transformation in component ["
	            << i << "]" << std::endl;
	  std::cout << "value: " << si_cart_U_check[i] << std::endl;
	  std::cout << "correct value: " << si_cart_U[i] << std::endl;
	  failed = -1;
        }
    }

    if (failed == 0)
        std::cout << "Coordinate transformations test passed..." << std::endl;
    else
        std::cout << "Coordinate transformations test failed..." << std::endl;

    return failed;
}
