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
renderView1.CameraPosition = [32.0, 32.0, 190.05252605770363]
renderView1.CameraFocalPoint = [32.0, 32.0, 32.0]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 49.49747468305833
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1
renderView1.Background = [0.0, 0.0, 0.0]

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

# create a new 'Slice AMR data'
sliceAMRdata1 = SliceAMRdata(registrationName='SliceAMRdata1', Input=input)
sliceAMRdata1.Level = 3
sliceAMRdata1.OffSet = 35.0
sliceAMRdata1.Normal = 'Z-Normal'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from sliceAMRdata1
sliceAMRdata1Display = Show(sliceAMRdata1, renderView1, 'AMRRepresentation')

# trace defaults for the display properties.
sliceAMRdata1Display.Representation = 'Outline'
sliceAMRdata1Display.ColorArrayName = ['POINTS', '']
sliceAMRdata1Display.SelectTCoordArray = 'None'
sliceAMRdata1Display.SelectNormalArray = 'None'
sliceAMRdata1Display.SelectTangentArray = 'None'
sliceAMRdata1Display.OSPRayScaleFunction = 'PiecewiseFunction'
sliceAMRdata1Display.SelectOrientationVectors = 'None'
sliceAMRdata1Display.ScaleFactor = 7.0
sliceAMRdata1Display.SelectScaleArray = 'None'
sliceAMRdata1Display.GlyphType = 'Arrow'
sliceAMRdata1Display.GlyphTableIndexArray = 'None'
sliceAMRdata1Display.GaussianRadius = 0.35000000000000003
sliceAMRdata1Display.SetScaleArray = [None, '']
sliceAMRdata1Display.ScaleTransferFunction = 'PiecewiseFunction'
sliceAMRdata1Display.OpacityArray = [None, '']
sliceAMRdata1Display.OpacityTransferFunction = 'PiecewiseFunction'
sliceAMRdata1Display.DataAxesGrid = 'GridAxesRepresentation'
sliceAMRdata1Display.PolarAxes = 'PolarAxesRepresentation'
sliceAMRdata1Display.ScalarOpacityUnitDistance = 2.3023915422888135

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'AMRBlocksSlice_%.6ts%cm.png'
pNG1.Writer.ImageResolution = [885, 774]
pNG1.Writer.Format = 'PNG'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(pNG1)
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
