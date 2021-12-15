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
renderView1.ViewSize = [1131, 939]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [32.0, 32.0, 32.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [-24.234943738007217, 12.79030062375982, 77.14896782254773]
renderView1.CameraFocalPoint = [31.99999999999993, 31.99999999999998, 31.99999999999996]
renderView1.CameraViewUp = [0.6572599982077688, -0.31659445000046216, 0.6839424310457926]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 60.6217782649107
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1131, 939)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'VisItChomboReader'
binaryBH_chi_ghosts3dhdf5 = VisItChomboReader(registrationName='input', FileName=['/home/miren/NR/GRChombo-public/Examples/BinaryBH/hdf5/BinaryBH_chi_ghosts.3d.hdf5'])
binaryBH_chi_ghosts3dhdf5.MeshStatus = ['Mesh']
binaryBH_chi_ghosts3dhdf5.CellArrayStatus = ['chi']

# create a new 'Slice AMR data'
sliceAMRdata1 = SliceAMRdata(registrationName='SliceAMRdata1', Input=binaryBH_chi_ghosts3dhdf5)
sliceAMRdata1.Level = 2
sliceAMRdata1.OffSet = 35.0
sliceAMRdata1.Normal = 'Z-Normal'

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from sliceAMRdata1
sliceAMRdata1Display = Show(sliceAMRdata1, renderView1, 'AMRRepresentation')

# get color transfer function/color map for 'chi'
chiLUT = GetColorTransferFunction('chi')
chiLUT.RGBPoints = [4.8189884420197785e-05, 0.231373, 0.298039, 0.752941, 0.4837166353392167, 0.865003, 0.865003, 0.865003, 0.9673850807940133, 0.705882, 0.0156863, 0.14902]
chiLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'chi'
chiPWF = GetOpacityTransferFunction('chi')
chiPWF.Points = [4.8189884420197785e-05, 0.0, 0.5, 0.0, 0.9673850807940133, 1.0, 0.5, 0.0]
chiPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
sliceAMRdata1Display.Representation = 'Surface'
sliceAMRdata1Display.ColorArrayName = ['CELLS', 'chi']
sliceAMRdata1Display.LookupTable = chiLUT
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
sliceAMRdata1Display.ScalarOpacityUnitDistance = 1.7175096404638197
sliceAMRdata1Display.ScalarOpacityFunction = chiPWF

# setup the color legend parameters for each legend in this view

# get color legend/bar for chiLUT in view renderView1
chiLUTColorBar = GetScalarBar(chiLUT, renderView1)
chiLUTColorBar.Title = 'chi'
chiLUTColorBar.ComponentTitle = ''

# set color bar visibility
chiLUTColorBar.Visibility = 1

# show color legend
sliceAMRdata1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity mapes used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup extractors
# ----------------------------------------------------------------

# create extractor
pNG1 = CreateExtractor('PNG', renderView1, registrationName='PNG1')
# trace defaults for the extractor.
# init the 'PNG' selected for 'Writer'
pNG1.Writer.FileName = 'RenderView1_%.6ts%cm.png'
pNG1.Writer.ImageResolution = [1131, 939]
pNG1.Writer.Format = 'PNG'

# ----------------------------------------------------------------
# restore active source
SetActiveSource(sliceAMRdata1)
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
