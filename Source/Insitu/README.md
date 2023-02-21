# Using ParaView Catalyst in-situ visualization with GRChombo

## Why use in-situ visualization?

Whilst the processing power of supercomputers has grown rapidly over recent
years, the corresponding growth in IO capabilities and data storage has not kept
pace. As a result, it is common to only write checkpoints/plot files relatively
infrequently in order to save on disk space. This means that traditional 
post-processed visualization can only be performed at relatively poor 
resolutions in time. Alternatively, separate "visualization runs" are sometimes
performed at significantly lower resolutions than that of the production
simulations which takes more time and can also result in inaccurate output. One
solution to this problem is performing the visualization _in-situ_ (i.e. during
the simulation itself).

In this guide, we describe how to use the in-situ visualization capabilities of
ParaView Catalyst with a GRChombo example.

![Catalyst Logo](https://www.paraview.org/Wiki/images/8/8a/CatalystLogo.png)

## Prerequisites

You will need a build of ParaView with Catalyst support. Unfortunately, the 
pre-built binaries available on the ParaView website for various OSs and 
architecture are not built with Catalyst support so it is necessary to build
ParaView from source.  Given the large number of dependencies, it is probably
easiest to use the 
[superbuild](https://gitlab.kitware.com/paraview/paraview-superbuild/) which
downloads and builds all dependencies before building ParaView. When configuring
the superbuild using CMake, you will want to make sure the following options are
set
```cmake
ENABLE_osmesa = ON # enables off-screen rendering
USE_SYSTEM_llvm = ON # a prerequisite for OSMesa above
USE_SYSTEM_mpi = ON # you will want to use the same MPI as for [GR]Chombo
PARAVIEW_BUILD_EDITION = CATALYST_RENDERING # This will build without IO
# libraries and a user interface so will be much easier and quicker
```
Note that the above setting for `PARAVIEW_BUILD_EDITION` will not build
the usual ParaView graphical user interface as this can take quite a bit longer
and can be more complicated. For more details see 
[this Kitware blog post](https://www.kitware.com/paraview-editions/). In order
to [generate ParaView Catalyst scripts](#generating-paraview-catalyst-scripts), 
you will probably want to use the GUI but for this part, you can use the
pre-built ParaView binaries available on the [ParaView
website](https://www.paraview.org/download/). Alternatively, you can set 
```cmake
PARAVIEW_BUILD_EDITION = CANONICAL
```

Since building ParaView can be a little daunting, unless you area already 
familiar with building software with CMake, it may be advisable to ask
your HPC system administrator to build a suitable version for you. If so, make
sure you specify that you need a version with Catalyst enabled and with 
off-screen rendering.

This in-situ adaptor has been tested with the following versions of ParaView:

| Tested ParaView versions |
| --- |
| 5.8.1 |
| 5.9.1 |
| 5.10.0 |
| 5.10.1 |
| 5.11.0 |

It may work with other versions but these have not been tested.

## Adding ParaView Catalyst to an existing GRChombo example

To add ParaView Catalyst to your existing example, make sure your example is 
up-to-date with this branch. Then simply add the line
```make
include $(GRCHOMBO_SOURCE)/Insitu/Make.insitu
```
in the `GNUmakefile` directly after the `src_dirs` variable has been set. See 
the [BinaryBH's `GNUmakefile`](../../Examples/BinaryBH/GNUmakefile) for an
example. Note that this script uses the `XTRALDFLAGS` and `cxxcppflags` Chombo
make variables so do *not* set these in your `Make.defs.local` file.

The in-situ visualization pipelines are called in 
`GRAMRLevel::postTimeStep()` after `GRAMRLevel::specificPostTimeStep()`.
If you wish to also process the pipeline before the first timestep, 
you can use the [`MultiLevelTask` class](../utils/MultiLevelTask.hpp). See the 
BinaryBH's [Main_BinaryBH.cpp](../../Examples/BinaryBH/Main_BinaryBH.cpp) for 
an example. 

If you are adding any code to your example that is specific to
Catalyst it is a good idea to conditionally compile it so that you can continue
to compile your example without Catalyst if you wish. For example
```cpp
#ifdef USE_CATALYST
<code that depends on Catalyst>
#endif
```
If you need to compute diagnostic quantities for visualization, you can override
the virtual function `GRAMRLevel::preCatalystCoProcess()` which is called before
Catalyst CoProcessing.


## Building a GRChombo example with ParaView Catalyst

First, make sure you have updated [Chombo](https://github.com/GRChombo/Chombo) 
to the latest commit on `main`.

To build a GRChombo example with ParaView Catalyst (for example the [BinaryBH
example](../../Examples/BinaryBH/)), next set the `PARAVIEW_DIR` to the path
which contains the build of ParaView as described [above](#prerequisites), for
example
```bash
export PARAVIEW_DIR=/path/to/paraview_build_dir
```
This directory should contain an executable at `bin/paraview-config` as this is
the script that will be used to determine the compiler/linker flags needed to
build with ParaView Catalyst. Since this script takes a long time to run, the
outputs are cached using the [cache.sh](./cache.sh) script to 
`/tmp/command_cache.<hash>`. If something messes up and this cached output is
incorrect, it can be helpful to remove these files using e.g.
```bash
rm -f /tmp/command_cache.*
```

Finally `make` the `insitu` target e.g.
```
make insitu -j 4
```
An executable binary with a filename containing '`_Insitu`' will be created. 
Note that if the `PARAVIEW_DIR` command is set, you can also `make` the 
usual `all` target to link with ParaView Catalyst but the executable filename
will be the conventional one (i.e. without '`_Insitu`').

If you wish to stop compiling/linking with ParaView Catalyst, unset the 
`PARAVIEW_DIR` environment variable with
```bash
unset PARAVIEW_DIR
```


## Generating ParaView Catalyst scripts

> **Warning** 
> The procedure to generate Catalyst scripts varies between
> different versions of ParaView. In particular, there are significant
> differences before v5.9. The following instructions apply to v5.9 or greater.

First, note that you do not need to use a build of ParaView with Catalyst
support in order to generate a Catalyst script. For example, you can just use
the [pre-built binaries from the ParaView
website](https://www.paraview.org/download/). In fact you probably cannot use
the same build as you will want a version with the GUI which precludes a version
with off-screen rendering (which you need for in-situ visualization).

The easiest way to generate a suitable ParaView Catalyst Python script using the
GUI is to load an HDF5 file that is similar to the configuration that will exist
in your simulation. For example, you could use a plot/checkpoint file from the
initial data of the simulation. Once you have loaded the file into ParaView,
**make sure you rename the source to 'input'** (by right clicking on it and
selecting 'Rename'). This is necessary so that Catalyst can replace the HDF5
data in the script with that coming from the simulation. 

Next apply filters and adjust the camera view as desired.

In order to generate a valid Catalyst script, it is then necessary to add an 
'Extractor' which will, for example, save an image. This can be done using the
menu options, Extractors → Image → PNG.

Finally, to generate the Catalyst script, use the menu options File → 
Save Catalyst State. By default, the output directory for Extracts will be
'dataset'.

Unfortunately, if you have built the `CATALYST_RENDERING` ParaView edition
as described [above](#prerequisites) for running in-situ, since this does not
contain any IO libraries, it is necessary to manually edit the generated
Catalyst script and replace the line
```python
input = VisItChomboReader(registrationName='input', FileName=<list of hdf5 files>)
```
with
```python
input = AMRGaussianPulseSource(registrationName='input')
```
in order to be able to use it.

> **Note**
> It is possible to pass information about the GRChombo grid (in particular, the
> lengths of the domain and the center) to Catalyst pipelines so that they
> automatically adjust appropriately. These are set as environment variables of
> the form `GRCHOMBO_PARAM_xxx` in 
> [ChomboParameters.hpp](../GRChomboCore/ChomboParameters.hpp) and can be 
> recovered in Catalyst Python scripts using the `os.getenv()` function from the
> `os` package. See the
> [slice_chi.py](../../Examples/BinaryBH/catalyst_scripts/slice_chi.py) script
> as an example.

## Using the in-situ visualization capabilities

There are several Catalyst specific parameters which control the in-situ
visualization processing which are described in the table below.

| Parameter Name             | Type             | Possible value [default] | Description |
| ---                        | ---              | ---                      | ---         |
| `catalyst_activate`        | `bool`           | `true`/[`false`] | Enables/disables Catalyst coprocessing at runtime. If set to `false`, other parameters are not set. |
| `catalyst_verbosity`       | `int`            | `-2`,...,[`verbosity`],..,`10` | Controls verbosity in `catalyst_pout` files and that printed to the normal `pout` files from the `CatalystAdaptor` class which interfaces with ParaView Catalyst |
| `catalyst_pout_prefix`     | `string`         | [`"catalyst_pout"`] | Filename prefix for the `catalyst_pout` files. They are written to `pout_subpath` and appended by `. <rank_number>` as for the normal `pout` files. Note that this can be set to the same string as `pout_prefix` to send the output to the normal `pout` files. |
| `catalyst_scripts_path`    | `string`         | `""` | Path that contains the Catalyst Python scrips generated as described [above](#generating-paraview-catalyst-scripts). Note that this path should be either an absolute path or relative to the current directory and _not_ `output_path`. |
| `catalyst_num_scripts`     | `int`            | [`1`] | Number of Catalyst Python scripts |
| `catalyst_scripts`         | `vector<string>` | - | Filenames of Catalyst Python scripts in the `catalyst_scripts_path` directory. |
| `catalyst_coprocess_level` | `int`            | [`0`],...,`max_level` | Level for which the Catalyst coprocess routine is called at the end of each timestep |
| `catalyst_remove_ghosts`   | `bool`           | `true`/[`false`] | If `true`, Chombo data will be deep-copied to new VTK arrays without ghosts rather than simply passing the pointers. This is for debugging and should not be necessary for most users. |
| `catalyst_abort_on_error`  | `bool`           | `true`/[`false`] | If `true`, the code will abort if there is an error in Catalyst. |
| `catalyst_num_vars`*        | `int`            | [`0`],...,`NUM_VARS + NUM_DIAGNOSTIC_VARS` | Number of variables in `catalyst_vars` |
| `catalyst_vars`*            | `vector<string>` | - | Restrict variables passed to Catalyst to only these ones. No restriction if `catalyst_num_vars == 0`.|
| `catalyst_write_vtk_files`  | `bool`           | `true`/[`false`] | Controls whether VTK XML files containing the full 3D AMR data passed to Catalyst are written or not. These files can be opened in ParaView. |
| `vtk_file_prefix`           | `string`         | [`"Catalyst_VTK_"`] | Filename prefix for VTK XML files. The filenames will be of the form `<vtk_file_prefix>_<6 digit timestep on catalyst_coprocess_level>`. Note there will be multiple files for each timestep: one of type `.vth` and a subdirectory containing `.vti` files for each box. |

*Note that scripts generated with versions of ParaView v5.8 or earlier provide 
information to the pipeline about the specific variables they require. Since 
this feature is not currently provided by ParaView v5.9 or later, scripts
generated with these versions request all variables by default. In order to
only pass specific variables to Catalyst, it is possible to restrict the 
variables passed using the `num_catalyst_vars` and `catalyst_vars` parameters
(as described above). 
If these are not set, and a script generated by ParaView v5.9 or later is used,
all variables will be passed to Catalyst.

If the `output_path` parameter is set, then this will be passed to Catalyst as 
the working directory and the subdirectory for any outputs created by
visualization pipelines will exist under this path.