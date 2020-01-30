/* GRChombo
 * Copyright 2012 The GRChombo collaboration.
 * Please refer to LICENSE in GRChombo's root directory.
 */

#ifndef INTEGRATIONMETHODSETUP_HPP_
#define INTEGRATIONMETHODSETUP_HPP_

#include "IntegrationMethod.hpp"

// define the static IntegrationMethods here
const IntegrationMethod IntegrationMethod::trapezium({0.5});
const IntegrationMethod IntegrationMethod::midpoint({1.0}, false);
const IntegrationMethod IntegrationMethod::simpson({0.3333333333333333,
                                                    1.3333333333333333});
const IntegrationMethod
    IntegrationMethod::boole({0.3111111111111111, 1.4222222222222222,
                              0.53333333333333, 1.4222222222222222});


#endif /* INTEGRATIONMETHODSETUP_HPP_ */
