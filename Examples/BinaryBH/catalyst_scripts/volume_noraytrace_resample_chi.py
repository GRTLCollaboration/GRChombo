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
renderView1.CameraPosition = [20.770091346798225, 54.44320054580248, 46.07511197528731]
renderView1.CameraFocalPoint = [32.000000000000014, 31.999999999999986, 31.999999999999993]
renderView1.CameraViewUp = [0.37830852443531104, -0.3485213720084932, 0.8575637081831853]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 60.6217782649107
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
input = VisItChomboReader(registrationName='input', FileName=['/home/miren/NR/GRChombo-public/Examples/BinaryBH/hdf5/BinaryBHp_000000.3d.hdf5', '/home/miren/NR/GRChombo-public/Examples/BinaryBH/hdf5/BinaryBHp_000001.3d.hdf5', '/home/miren/NR/GRChombo-public/Examples/BinaryBH/hdf5/BinaryBHp_000002.3d.hdf5', '/home/miren/NR/GRChombo-public/Examples/BinaryBH/hdf5/BinaryBHp_000003.3d.hdf5'])
input.MeshStatus = ['Mesh']
input.CellArrayStatus = ['chi']

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(registrationName='ResampleToImage1', Input=input)
resampleToImage1.SamplingBounds = [-3.0, 67.0, -3.0, 67.0, -3.0, 67.0]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from resampleToImage1
resampleToImage1Display = Show(resampleToImage1, renderView1, 'UniformGridRepresentation')

# get color transfer function/color map for 'chi'
chiLUT = GetColorTransferFunction('chi')
chiLUT.RGBPoints = [0.006705302783959541, 1.0, 0.988235, 0.968627, 0.025918898344160616, 1.0, 0.952941, 0.878431, 0.05473929168446223, 0.968627, 0.905882, 0.776471, 0.10277328058496492, 0.94902, 0.898039, 0.647059, 0.1508072694854676, 0.901961, 0.878431, 0.556863, 0.1988412583859703, 0.847059, 0.858824, 0.482353, 0.24687524728647298, 0.690196, 0.819608, 0.435294, 0.2949092361869757, 0.513725, 0.768627, 0.384314, 0.3429432250874783, 0.337255, 0.721569, 0.337255, 0.3909772139879811, 0.278431, 0.658824, 0.392157, 0.4390112028884838, 0.231373, 0.639216, 0.435294, 0.4870451917889864, 0.203922, 0.6, 0.486275, 0.5350791806894891, 0.172549, 0.568627, 0.537255, 0.5831131695899918, 0.141176, 0.517647, 0.54902, 0.6311471584904945, 0.133333, 0.458824, 0.541176, 0.6791811473909971, 0.12549, 0.396078, 0.529412, 0.7272151362914999, 0.117647, 0.321569, 0.521569, 0.7752491251920026, 0.121569, 0.258824, 0.509804, 0.8232831140925052, 0.133333, 0.227451, 0.501961, 0.8713171029930079, 0.145098, 0.192157, 0.490196, 0.9193510918935106, 0.188235, 0.164706, 0.470588, 0.9673850807940133, 0.258824, 0.196078, 0.439216]
chiLUT.ColorSpace = 'RGB'
chiLUT.NanColor = [0.25, 0.0, 0.0]
chiLUT.NanOpacity = 0.0
chiLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'chi'
chiPWF = GetOpacityTransferFunction('chi')
chiPWF.Points = [0.006705302783959541, 1.0, 0.5, 0.0, 0.055023517459630966, 0.9053254723548889, 0.5, 0.0, 0.1061839833855629, 0.7692307829856873, 0.5, 0.0, 0.21418939530849457, 0.47337278723716736, 0.5, 0.0, 0.37335529923439026, 0.23076923191547394, 0.5, 0.0, 0.5808394117544118, 0.0, 0.5, 0.0, 0.9673850807940133, 0.0, 0.5, 0.0]
chiPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
resampleToImage1Display.Representation = 'Volume'
resampleToImage1Display.ColorArrayName = ['POINTS', 'chi']
resampleToImage1Display.LookupTable = chiLUT
resampleToImage1Display.SelectTCoordArray = 'None'
resampleToImage1Display.SelectNormalArray = 'None'
resampleToImage1Display.SelectTangentArray = 'None'
resampleToImage1Display.OSPRayScaleArray = 'chi'
resampleToImage1Display.OSPRayScaleFunction = 'PiecewiseFunction'
resampleToImage1Display.SelectOrientationVectors = 'None'
resampleToImage1Display.ScaleFactor = 6.999993000000001
resampleToImage1Display.SelectScaleArray = 'None'
resampleToImage1Display.GlyphType = 'Arrow'
resampleToImage1Display.GlyphTableIndexArray = 'None'
resampleToImage1Display.GaussianRadius = 0.34999965000000005
resampleToImage1Display.SetScaleArray = ['POINTS', 'chi']
resampleToImage1Display.ScaleTransferFunction = 'PiecewiseFunction'
resampleToImage1Display.OpacityArray = ['POINTS', 'chi']
resampleToImage1Display.OpacityTransferFunction = 'PiecewiseFunction'
resampleToImage1Display.DataAxesGrid = 'GridAxesRepresentation'
resampleToImage1Display.PolarAxes = 'PolarAxesRepresentation'
resampleToImage1Display.ScalarOpacityUnitDistance = 1.224681164507726
resampleToImage1Display.ScalarOpacityFunction = chiPWF
resampleToImage1Display.OpacityArrayName = ['POINTS', 'chi']
resampleToImage1Display.SliceFunction = 'Plane'
resampleToImage1Display.Slice = 49

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
resampleToImage1Display.ScaleTransferFunction.Points = [0.16307053589122933, 0.0, 0.5, 0.0, 0.9673850807940133, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
resampleToImage1Display.OpacityTransferFunction.Points = [0.16307053589122933, 0.0, 0.5, 0.0, 0.9673850807940133, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
resampleToImage1Display.SliceFunction.Origin = [32.0, 32.0, 32.0]

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
pNG1.Writer.FileName = 'VolumeChi_%.6ts%cm.png'
pNG1.Writer.ImageResolution = [885, 778]
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
