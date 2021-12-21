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
renderView1.CenterOfRotation = [8.0, 8.0, 4.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [8.038575486273551, 7.488142173440099, 14.569596590973248]
renderView1.CameraFocalPoint = [8.008182364514674, 7.891428269059167, 6.241949497853552]
renderView1.CameraViewUp = [-0.754397824803147, -0.6557763921833819, -0.029004230505694627]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 14.517231140957975
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

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
input.CellArrayStatus = ['chi']

# create a new 'Extract Level'
extractLevel1 = ExtractLevel(registrationName='ExtractLevel1', Input=input)
extractLevel1.Levels = [2]

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(registrationName='ResampleToImage1', Input=extractLevel1)
resampleToImage1.SamplingBounds = [3.75, 12.25, 5.75, 10.25, -0.25, 2.25]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from resampleToImage1
resampleToImage1Display = Show(resampleToImage1, renderView1, 'UniformGridRepresentation')

# get color transfer function/color map for 'chi'
chiLUT = GetColorTransferFunction('chi')
chiLUT.RGBPoints = [0.0, 0.0862745098039216, 0.00392156862745098, 0.298039215686275, 0.029344877158935277, 0.113725, 0.0235294, 0.45098, 0.0537159862131328, 0.105882, 0.0509804, 0.509804, 0.0706265851511459, 0.0392157, 0.0392157, 0.560784, 0.0870398318296482, 0.0313725, 0.0980392, 0.6, 0.10295566969276873, 0.0431373, 0.164706, 0.639216, 0.12583466541255697, 0.054902, 0.243137, 0.678431, 0.15617428436025385, 0.054902, 0.317647, 0.709804, 0.19347701615880114, 0.0509804, 0.396078, 0.741176, 0.2176616431786517, 0.0392157, 0.466667, 0.768627, 0.24184627019850216, 0.0313725, 0.537255, 0.788235, 0.267087853698452, 0.0313725, 0.615686, 0.811765, 0.2929510760870558, 0.0235294, 0.709804, 0.831373, 0.3188143550315277, 0.0509804, 0.8, 0.85098, 0.3402012374169067, 0.0705882, 0.854902, 0.870588, 0.36009606302373837, 0.262745, 0.901961, 0.862745, 0.3775040142212683, 0.423529, 0.941176, 0.87451, 0.4043619976847706, 0.572549, 0.964706, 0.835294, 0.42226735769767926, 0.658824, 0.980392, 0.843137, 0.43519896889198123, 0.764706, 0.980392, 0.866667, 0.4491253128832477, 0.827451, 0.980392, 0.886275, 0.47648067688419493, 0.913725, 0.988235, 0.937255, 0.48493597635320285, 1.0, 1.0, 0.972549019607843, 0.49339127582221104, 0.988235, 0.980392, 0.870588, 0.5058255347569962, 0.992156862745098, 0.972549019607843, 0.803921568627451, 0.5152755670229687, 0.992157, 0.964706, 0.713725, 0.5311914048860893, 0.988235, 0.956863, 0.643137, 0.5565572750151794, 0.980392, 0.917647, 0.509804, 0.5779442139564234, 0.968627, 0.87451, 0.407843, 0.5998285051571812, 0.94902, 0.823529, 0.321569, 0.615744343020301, 0.929412, 0.776471, 0.278431, 0.639120747555468, 0.909804, 0.717647, 0.235294, 0.6600102776813318, 0.890196, 0.658824, 0.196078, 0.6771695565558091, 0.878431, 0.619608, 0.168627, 0.7013541835756594, 0.870588, 0.54902, 0.156863, 0.7255388105955098, 0.85098, 0.47451, 0.145098, 0.7497234376153601, 0.831373, 0.411765, 0.133333, 0.7739080646352107, 0.811765, 0.345098, 0.113725, 0.798092691655061, 0.788235, 0.266667, 0.0941176, 0.8222773186749113, 0.741176, 0.184314, 0.0745098, 0.8464619456947616, 0.690196, 0.12549, 0.0627451, 0.870646572714612, 0.619608, 0.0627451, 0.0431373, 0.893276919519733, 0.54902, 0.027451, 0.0705882, 0.9131717283365444, 0.470588, 0.0156863, 0.0901961, 0.9355533868194685, 0.4, 0.00392157, 0.101961, 0.9673850807940133, 0.188235294117647, 0.0, 0.0705882352941176]
chiLUT.ColorSpace = 'Lab'
chiLUT.NanColor = [1.0, 0.0, 0.0]
chiLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'chi'
chiPWF = GetOpacityTransferFunction('chi')
chiPWF.Points = [0.0, 1.0, 0.5, 0.0, 0.04865546501470074, 0.6094674468040466, 0.5, 0.0, 0.12020762722887196, 0.4378698170185089, 0.5, 0.0, 0.22324271380797622, 0.23076923191547394, 0.5, 0.0, 0.37207121721551956, 0.08875739574432373, 0.5, 0.0, 0.578141420384064, 0.0, 0.5, 0.0, 0.9673850807940133, 0.0, 0.5, 0.0]
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
resampleToImage1Display.ScaleFactor = 0.84999915
resampleToImage1Display.SelectScaleArray = 'None'
resampleToImage1Display.GlyphType = 'Arrow'
resampleToImage1Display.GlyphTableIndexArray = 'None'
resampleToImage1Display.GaussianRadius = 0.042499957500000005
resampleToImage1Display.SetScaleArray = ['POINTS', 'chi']
resampleToImage1Display.ScaleTransferFunction = 'PiecewiseFunction'
resampleToImage1Display.OpacityArray = ['POINTS', 'chi']
resampleToImage1Display.OpacityTransferFunction = 'PiecewiseFunction'
resampleToImage1Display.DataAxesGrid = 'GridAxesRepresentation'
resampleToImage1Display.PolarAxes = 'PolarAxesRepresentation'
resampleToImage1Display.ScalarOpacityUnitDistance = 0.10037670222093371
resampleToImage1Display.ScalarOpacityFunction = chiPWF
resampleToImage1Display.OpacityArrayName = ['POINTS', 'chi']
resampleToImage1Display.SliceFunction = 'Plane'
resampleToImage1Display.Slice = 49

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
resampleToImage1Display.ScaleTransferFunction.Points = [0.007669444988778632, 0.0, 0.5, 0.0, 0.6724393877650426, 1.0, 0.5, 0.0]

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
resampleToImage1Display.OpacityTransferFunction.Points = [0.007669444988778632, 0.0, 0.5, 0.0, 0.6724393877650426, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
resampleToImage1Display.SliceFunction.Origin = [8.0, 8.0, 1.0]

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
