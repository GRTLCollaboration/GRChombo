/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#include "ADMFixedBGVars.hpp"
#include "CCZ4Geometry.hpp"
#include "Cell.hpp"
#include "Coordinates.hpp"
#include "FourthOrderDerivatives.hpp"
#include "GRInterval.hpp"
#include "Tensor.hpp"
#include "CoordinateTransformations.hpp"
#include "UserVariables.hpp" //This files needs NUM_VARS - total number of components
#include "VarsTools.hpp"
#include "simd.hpp"

int main()
{
  int failed = 0;
    
  using namespace TensorAlgebra;
  using namespace CoordinateTransformations;
  
  const data_t x = coords.x;
  const double y = coords.y;
  const double z = coords.z;
  const data_t r = coords.get_radius();
  data_t rho2 =
    simd_max(coords.x * coords.x + coords.y * coords.y, 1e-12);
  data_t r2sintheta = sqrt(rho2) * R;
  // Test cartesian_to_spherical_LL
  Tensor<2, data_t> Mij_cart;
  Mij_cart[0][0] = 1.;
  Mij_cart[1][1] = 1.;
  Mij_cart[2][2] = 1.;
  
  Tensor<2, data_t> Mij_spher;
  Mij_spher[0][0] = 1.;
  Mij_spher[1][1] = r*r;
  Mij_spher[2][2] = r2sintheta;
    
  Tensor<2, data_t> Mij_spher_check;
  Mij_spher_check = cartesian_to_spherical_LL(Mij_cart,x,y,z)
    FOR2(i, j)
    {
      if (Mij_spher_check[i][j] != Mij_spher[i][j])
	{int failed = 1;}
    }
  
  FOR2(i, j)
    {
      double diff = Mij_spher_check[i][j] - Mij_spher[i][j];
      if (diff > 1e-14)
	{
	  std::cout << "Cartesian to spherical metric transformation in component [" 
		      << i << "]["<< j << "]" << std::endl;
	    std::cout << "value: " << ricciZ.LL[i][j] << std::endl;
	    std::cout << "correct value: " << ricciZ_known[i][j] << std::endl;
            failed = -1;
	  }
      }    

    // Test cartesian_to_spherical_UU
    const auto _UU = compute_inverse();


    // Test cartesian_to_spherical_U
    
    Tensor<1, data_t> si;
    si[0] = x/R;
    si[1] = y/R;
    si[2] = z/R;
	
    data_t si_norm = 0.0;
    FOR2(j, k)
      {
	si_norm += si[j]*si[k]*metric_vars.gamma[j][k];
      }

    FOR1(i)
    {
      si[i] = si[i]/sqrt(si_norm);
    }

    Tensor<1, data_t> si_spher;
    si_spher[0] = 1.0/sqrt(gamma_spher[0][0]);
    si_spher[1] = 0.0;
    si_spher[2] = 0.0;
	
    if (si_spher == spherical_to_cartesian_U(si_spher, x,y,z))
      {
	
      }

    if (failed == 0)
      std::cout << "Coordinate transformations test passed..." << std::endl;
    else
      std::cout << "Coordinate transformations test failed..." << std::endl;

    return failed;
}

