# script-version: 2.0
# Catalyst state generated using paraview version 5.11.0

# This script takes a slice through the center of the domain of constant z and
# saves it to a VTM file (+ corresponding folder) that can be visualized in 
# ParaView (all variables passed to Catalyst are saved)

# Import this to get environment variables
import os

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


# Get center from environment
grchombo_center_str = os.getenv('GRCHOMBO_PARAM_CENTER')
grchombo_center = [float(dim) for dim in grchombo_center_str.split(' ')]

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create source
input = AMRGaussianPulseSource(registrationName='input')
input.MeshStatus = ['Mesh']

# create a new 'Slice With Plane'
sliceWithPlane1 = SliceWithPlane(registrationName='SliceWithPlane1', Input=input)
sliceWithPlane1.PlaneType = 'Plane'
sliceWithPlane1.Level = 3

# init the 'Plane' selected for 'PlaneType'
sliceWithPlane1.PlaneType.Origin = grchombo_center
sliceWithPlane1.PlaneType.Normal = [0.0, 0.0, 1.0]

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
vTM1 = CreateExtractor('VTM', sliceWithPlane1, registrationName='VTM1')
# trace defaults for the extractor.
vTM1.Trigger = 'TimeStep'

# init the 'VTM' selected for 'Writer'
vTM1.Writer.FileName = 'CenterZSlice_{timestep:06d}.vtm'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(vTM1)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'TimeStep'
options.CatalystLiveTrigger = 'TimeStep'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
