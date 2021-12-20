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
renderView1.ViewSize = [885, 778]
renderView1.AxesGrid = 'GridAxes3DActor'
renderView1.CenterOfRotation = [32.0, 32.0, 32.0]
renderView1.StereoType = 'Crystal Eyes'
renderView1.CameraPosition = [20.770091346798225, 54.44320054580248, 46.07511197528731]
renderView1.CameraFocalPoint = [32.000000000000014, 31.999999999999986, 31.999999999999993]
renderView1.CameraViewUp = [0.37830852443531104, -0.3485213720084932, 0.8575637081831853]
renderView1.CameraFocalDisk = 1.0
renderView1.CameraParallelScale = 60.6217782649107
renderView1.EnableRayTracing = 1
renderView1.BackEnd = 'OSPRay raycaster'
renderView1.OSPRayMaterialLibrary = materialLibrary1

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(885, 778)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'VisItChomboReader'
binaryBHp_00000 = VisItChomboReader(registrationName='input', FileName=['/home/miren/NR/GRChombo-public/Examples/BinaryBH/hdf5/BinaryBHp_000000.3d.hdf5', '/home/miren/NR/GRChombo-public/Examples/BinaryBH/hdf5/BinaryBHp_000001.3d.hdf5', '/home/miren/NR/GRChombo-public/Examples/BinaryBH/hdf5/BinaryBHp_000002.3d.hdf5', '/home/miren/NR/GRChombo-public/Examples/BinaryBH/hdf5/BinaryBHp_000003.3d.hdf5'])
binaryBHp_00000.MeshStatus = ['Mesh']
binaryBHp_00000.CellArrayStatus = ['chi']

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from binaryBHp_00000
binaryBHp_00000Display = Show(binaryBHp_00000, renderView1, 'AMRRepresentation')

# get color transfer function/color map for 'chi'
chiLUT = GetColorTransferFunction('chi')
chiLUT.RGBPoints = [0.006705302783959541, 0.0862745098039216, 0.00392156862745098, 0.298039215686275, 0.03584677977638898, 0.113725, 0.0235294, 0.45098, 0.060048963685229544, 0.105882, 0.0509804, 0.509804, 0.07684234902507095, 0.0392157, 0.0392157, 0.560784, 0.09314182943734466, 0.0313725, 0.0980392, 0.6, 0.10894734875818932, 0.0431373, 0.164706, 0.639216, 0.13166776172045955, 0.054902, 0.243137, 0.678431, 0.16179708557943404, 0.054902, 0.317647, 0.709804, 0.19884125838596878, 0.0509804, 0.396078, 0.741176, 0.22285825283622032, 0.0392157, 0.466667, 0.768627, 0.24687524728647184, 0.0313725, 0.537255, 0.788235, 0.27194187206163456, 0.0313725, 0.615686, 0.811765, 0.2976258269170425, 0.0235294, 0.709804, 0.831373, 0.3233098379363089, 0.0509804, 0.8, 0.85098, 0.3445484799520089, 0.0705882, 0.854902, 0.870588, 0.3643054071849909, 0.262745, 0.901961, 0.862745, 0.3815926974524059, 0.423529, 0.941176, 0.87451, 0.40826451832681654, 0.572549, 0.964706, 0.835294, 0.4260457696856576, 0.658824, 0.980392, 0.843137, 0.43888774711336165, 0.764706, 0.980392, 0.866667, 0.4527175624781391, 0.827451, 0.980392, 0.886275, 0.4798833163620467, 0.913725, 0.988235, 0.937255, 0.4882800090319688, 1.0, 1.0, 0.972549019607843, 0.49667670170189104, 0.988235, 0.980392, 0.870588, 0.5090247742020214, 0.992156862745098, 0.972549019607843, 0.803921568627451, 0.518409304809017, 0.992157, 0.964706, 0.713725, 0.5342148241298618, 0.988235, 0.956863, 0.643137, 0.559404874057699, 0.980392, 0.917647, 0.509804, 0.5806435722372544, 0.968627, 0.87451, 0.407843, 0.6023761753443805, 0.94902, 0.823529, 0.321569, 0.6181816946652244, 0.929412, 0.776471, 0.278431, 0.6413960687189206, 0.909804, 0.717647, 0.235294, 0.6621408058070486, 0.890196, 0.658824, 0.196078, 0.679181147390997, 0.878431, 0.619608, 0.168627, 0.7031981418412483, 0.870588, 0.54902, 0.156863, 0.7272151362914998, 0.85098, 0.47451, 0.145098, 0.751232130741751, 0.831373, 0.411765, 0.133333, 0.7752491251920026, 0.811765, 0.345098, 0.113725, 0.7992661196422539, 0.788235, 0.266667, 0.0941176, 0.8232831140925052, 0.741176, 0.184314, 0.0745098, 0.8473001085427566, 0.690196, 0.12549, 0.0627451, 0.8713171029930079, 0.619608, 0.0627451, 0.0431373, 0.8937905905179446, 0.54902, 0.027451, 0.0705882, 0.9135475010772842, 0.470588, 0.0156863, 0.0901961, 0.935774024030506, 0.4, 0.00392157, 0.101961, 0.9673850807940133, 0.188235294117647, 0.0, 0.0705882352941176]
chiLUT.ColorSpace = 'Lab'
chiLUT.NanColor = [1.0, 0.0, 0.0]
chiLUT.ScalarRangeInitialized = 1.0

# get opacity transfer function/opacity map for 'chi'
chiPWF = GetOpacityTransferFunction('chi')
chiPWF.Points = [0.006705302783959541, 1.0, 0.5, 0.0, 0.055023517459630966, 0.6094674468040466, 0.5, 0.0, 0.126079723238945, 0.4378698170185089, 0.5, 0.0, 0.22840063273906708, 0.23076923191547394, 0.5, 0.0, 0.37619754672050476, 0.08875739574432373, 0.5, 0.0, 0.5808393955230713, 0.0, 0.5, 0.0, 0.9673850536346436, 0.0, 0.5, 0.0]
chiPWF.ScalarRangeInitialized = 1

# trace defaults for the display properties.
binaryBHp_00000Display.Representation = 'Volume'
binaryBHp_00000Display.ColorArrayName = ['CELLS', 'chi']
binaryBHp_00000Display.LookupTable = chiLUT
binaryBHp_00000Display.SelectTCoordArray = 'None'
binaryBHp_00000Display.SelectNormalArray = 'None'
binaryBHp_00000Display.SelectTangentArray = 'None'
binaryBHp_00000Display.OSPRayScaleFunction = 'PiecewiseFunction'
binaryBHp_00000Display.SelectOrientationVectors = 'None'
binaryBHp_00000Display.ScaleFactor = 7.0
binaryBHp_00000Display.SelectScaleArray = 'None'
binaryBHp_00000Display.GlyphType = 'Arrow'
binaryBHp_00000Display.GlyphTableIndexArray = 'None'
binaryBHp_00000Display.GaussianRadius = 0.35000000000000003
binaryBHp_00000Display.SetScaleArray = [None, '']
binaryBHp_00000Display.ScaleTransferFunction = 'PiecewiseFunction'
binaryBHp_00000Display.OpacityArray = [None, '']
binaryBHp_00000Display.OpacityTransferFunction = 'PiecewiseFunction'
binaryBHp_00000Display.DataAxesGrid = 'GridAxesRepresentation'
binaryBHp_00000Display.PolarAxes = 'PolarAxesRepresentation'
binaryBHp_00000Display.ScalarOpacityUnitDistance = 0.8192401766441247
binaryBHp_00000Display.ScalarOpacityFunction = chiPWF

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
