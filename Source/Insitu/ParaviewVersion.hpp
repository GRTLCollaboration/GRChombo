#ifndef PARAVIEW_VERSION_HPP_
#define PARAVIEW_VERSION_HPP_

#ifdef USE_CATALYST
// ParaView version checking doesn't seem to be very stable between versions
#ifdef __has_include
#if __has_include(<vtkPVVersion.h>)
#include <vtkPVVersion.h>
#define PARAVIEW_VERSION_KNOWN
#elif __has_include(<vtkPVConfig.h>)
#include <vtkPVConfig.h>
#define PARAVIEW_VERSION_KNOWN
#endif /* __has_include(<vtk....h>) */
#endif /* __has_include */

#ifdef PARAVIEW_VERSION_KNOWN
// there is something similar in the new vtkPVVersion.h header but it uses
// the build rather than the patch number. It's also not available in
// vtkPVConfig.h
#define PARAVIEW_VERSION_TEST(major, minor, patch)                             \
    10000 * major + 100 * minor + patch
#define PARAVIEW_VERSION_HERE                                                  \
    PARAVIEW_VERSION_TEST(PARAVIEW_VERSION_MAJOR, PARAVIEW_VERSION_MINOR,      \
                          PARAVIEW_VERSION_PATCH)
#endif /* PARAVIEW_VERSION_KNOWN */

#endif /* USE_CATALYST */

#endif /* PARAVIEW_VERSION_HPP_ */