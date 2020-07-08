/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "Coordinates.hpp"
#include "TensorAlgebra.hpp"
#include "CoordinateTransformations.hpp"
#include "simd.hpp"

//! Class which creates the initial conditions                                                                                       
int main()
{
  int failed = 0;

  // protected:
  //const double m_dx;
  //const std::array<double, CH_SPACEDIM> m_center;

  //public:
  //! The constructor for the class                                                                                                  
  //Coords(const std::array<double, CH_SPACEDIM> a_center,
  //       const double a_dx)
  //  : m_dx(a_dx), m_center(a_center)
  //{
  //}

  //! Function to compute the value of all the initial vars on the grid                                                              
  template <class data_t> void compute(Cell<data_t> current_cell) const
  {
    Coordinates coords(current_cell, 0.1, {0,0,0});
    const data_t x = coords.x;
    const double y = coords.y;
    const double z = coords.z;
    const data_t r = coords.get_radius();
    data_t rho2 =
      simd_max(coords.x * coords.x + coords.y * coords.y, 1e-12);
    data_t r2sintheta = sqrt(rho2) * R;
  }
    
  using namespace TensorAlgebra;
  using namespace CoordinateTransformations;

  //const data_t x = coords.x;
  //const double y = coords.y;
  //const double z = coords.z;
  //const data_t r = coords.get_radius();
  //data_t rho2 =
  //  simd_max(coords.x * coords.x + coords.y * coords.y, 1e-12);
  //data_t r2sintheta = sqrt(rho2) * R;

  // Test tensor transformations
  Tensor<2, data_t> Mij_cart;
  Mij_cart[0][0] = 1.;
  Mij_cart[1][1] = 1.;
  Mij_cart[2][2] = 1.;
  
  Tensor<2, data_t> Mij_spher;
  Mij_spher[0][0] = 1.;
  Mij_spher[1][1] = r*r;
  Mij_spher[2][2] = r2sintheta;
    
  // Test cartesian_to_spherical_LL
  Tensor<2, data_t> Mij_spher_check;
  Mij_spher_check = cartesian_to_spherical_LL(Mij_cart,x,y,z);

  FOR2(i, j)
    {
      double diff = Mij_spher_check[i][j] - Mij_spher[i][j];
      if (diff > 1e-14)
	{
	  std::cout << "Failed cart_to_spher_LL transformation in component [" 
		    << i << "]["<< j << "]" << std::endl;
	  std::cout << "value: " << Mij_spher_check[i][j] << std::endl;
	  std::cout << "correct value: " << Mij_spher[i][j] << std::endl;
	  failed = -1;
	}
    }    

  // Test spherical_to_cartesian_LL
  Tensor<2, data_t> Mij_cart_check;
  Mij_cart_check = spherical_to_cartesian_LL(Mij_spher,x,y,z);

    FOR2(i, j)
    {
      double diff = Mij_cart_check[i][j] - Mij_cart[i][j];
      if (diff > 1e-14)
        {
	  std::cout << "Failed spher_to_cart_LL transformation in component ["
                    << i << "]["<< j << "]" << std::endl;
	  std::cout << "value: " << Mij_cart_check[i][j] << std::endl;
	  std::cout << "correct value: " << Mij_cart[i][j] << std::endl;
          failed = -1;
        }
    }

  // Test cartesian_to_spherical_UU
  Tensor<2, data_t> Mij_spher_UU;
  Tensor<2, data_t> Mij_spher_UU_check;
  Mij_spher_UU_check = cartesian_to_spherical_UU(compute_inverse(Mij_cart),x,y,z);
  Mij_spher_UU = compute_inverse(Mij_spher);

    FOR2(i, j)
    {
      double diff = Mij_spher_UU_check[i][j] - Mij_UU_spher[i][j];
      if (diff > 1e-14)
        {
	  std::cout << "Failed cart_to_spher_UU transformation in component ["
                    << i << "]["<< j << "]" << std::endl;
	  std::cout << "value: " << Mij_spher_UU_check[i][j] << std::endl;
	  std::cout << "correct value: " << Mij_UU_spher[i][j] << std::endl;
          failed = -1;
        }
    }

    // Test spherical_to_cartesian_UU
    Tensor<2, data_t> Mij_cart_UU;
    Tensor<2, data_t> Mij_cart_UU_check;
    Mij_cart_UU_check = spherical_to_cartesian_LL(compute_inverse(Mij_spher),x,y,z);
    Mij_cart_UU = compute_inverse(Mij_cart);

      FOR2(i, j)
      {
	double diff = Mij_cart_UU_check[i][j] - Mij_UU_cart[i][j];
	if (diff > 1e-14)
	  {
	    std::cout << "Failed spher_to_cart_UU transformation in component ["
		      << i << "]["<< j << "]" << std::endl;
	    std::cout << "value: " << Mij_cart_UU_check[i][j] << std::endl;
	    std::cout << "correct value: " << Mij_UU_cart[i][j] << std::endl;
	    failed = -1;
	  }
      }

  // Test vector transformations
  
  Tensor<1, data_t> si;
  si_cart_U[0] = x/r;
  si_cart_U[1] = y/r;
  si_cart_U[2] = z/r;
  
  Tensor<1, data_t> si_spher;
  si_spher_U[0] = 1.0;
  si_spher_U[1] = 0.0;
  si_spher_U[2] = 0.0;
	
  // Test cartesian_to_spherical_U
  Tensor<1, data_t> si_spher_U;
  Tensor<1, data_t> si_spher_U_check;
  si_spher_U_check = cartesian_to_spherical_U(si_cart_U,x,y,z);

  FOR1(i)
    {
      double diff = si_spher_U_check[i] - si_U_spher[i];
      if (diff > 1e-14)
        {
	  std::cout << "Failed cart_to_spher_U transformation in component ["
                    << i << "]" << std::endl;
	  std::cout << "value: " << si_spher_U_check[i] << std::endl;
	  std::cout << "correct value: " << si_U_spher[i] << std::endl;
          failed = -1;
        }
    }

  // Test spherical_to_cartesian_U
  Tensor<2, data_t> si_cart_U;
  Tensor<2, data_t> si_cart_U_check;
  si_cart_U_check = spherical_to_cartesian_L(si_spher_U,x,y,z);

  FOR1(i)
    {
      double diff = si_cart_U_check[i] - si_U_cart[i];
      if (diff > 1e-14)
	{
	  std::cout << "Failed spher_to_cart_U transformation in component ["
		    << i << "]" << std::endl;
	  std::cout << "value: " << si_cart_U_check[i] << std::endl;
	  std::cout << "correct value: " << si_U_cart[i] << std::endl;
	  failed = -1;
	}
    }

  if (failed == 0)
    std::cout << "Coordinate transformations test passed..." << std::endl;
  else
    std::cout << "Coordinate transformations test failed..." << std::endl;
  
  return failed;
}

