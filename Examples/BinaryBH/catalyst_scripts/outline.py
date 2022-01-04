# script-version: 2.0
# Catalyst state generated using paraview version 5.9.1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [885, 774]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterAxesVisibility = 1
renderView1.CenterOfRotation = [32.0, 32.0, 32.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [145.67639666481145, -111.0984129619794, 178.49801336649438]
renderView1.CameraFocalPoint = [32.000000000000014, 31.999999999999982, 32.00000000000003]
renderView1.CameraViewUp = [-0.10746516487670288, 0.6682514214076767, 0.7361326484572212]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 60.6217782649107
renderView1.CameraParallelProjection = 1
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

# init the 'GridAxes3DActor' selected for 'AxesGrid'
renderView1.AxesGrid.Visibility = 1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(885, 774)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'VisItChomboReader'
input = VisItChomboReader(registrationName='input', FileName=['/home/miren/NR/GRChombo-public/Examples/BinaryBH/hdf5/BinaryBH_000000.3d.hdf5', '/home/miren/NR/GRChombo-public/Examples/BinaryBH/hdf5/BinaryBH_000001.3d.hdf5', '/home/miren/NR/GRChombo-public/Examples/BinaryBH/hdf5/BinaryBH_000003.3d.hdf5'])
input.MeshStatus = ['Mesh']
input.CellArrayStatus = []

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from input
inputDisplay = Show(input, renderView1, 'AMRRepresentation')

# trace defaults for the display properties.
inputDisplay.Representation = 'Outline'
inputDisplay.ColorArrayName = [None, '']
inputDisplay.SelectTCoordArray = 'None'
inputDisplay.SelectNormalArray = 'None'
inputDisplay.SelectTangentArray = 'None'
inputDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
inputDisplay.SelectOrientationVectors = 'None'
inputDisplay.ScaleFactor = 7.0
inputDisplay.SelectScaleArray = 'None'
inputDisplay.GlyphType = 'Arrow'
inputDisplay.GlyphTableIndexArray = 'None'
inputDisplay.GaussianRadius = 0.35000000000000003
inputDisplay.SetScaleArray = [None, '']
inputDisplay.ScaleTransferFunction = 'PiecewiseFunction'
inputDisplay.OpacityArray = [None, '']
inputDisplay.OpacityTransferFunction = 'PiecewiseFunction'
inputDisplay.DataAxesGrid = 'GridAxesRepresentation'
inputDisplay.PolarAxes = 'PolarAxesRepresentation'
inputDisplay.ScalarOpacityUnitDistance = 0.8192401766441247

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'RenderView1_%.6ts%cm.png'
pNG1.Writer.ImageResolution = [885, 774]
pNG1.Writer.Format = 'PNG'
pNG1.Writer.ResetDisplay = 1

# ----------------------------------------------------------------
# restore active source
SetActiveSource(input)
# ----------------------------------------------------------------

# ------------------------------------------------------------------------------
# Catalyst options
from paraview import catalyst
options = catalyst.Options()
options.GlobalTrigger = 'TimeStep'
options.EnableCatalystLive = 1
options.CatalystLiveTrigger = 'TimeStep'

# ------------------------------------------------------------------------------
if __name__ == '__main__':
    from paraview.simple import SaveExtractsUsingCatalystOptions
    # Code for non in-situ environments; if executing in post-processing
    # i.e. non-Catalyst mode, let's generate extracts using Catalyst options
    SaveExtractsUsingCatalystOptions(options)
