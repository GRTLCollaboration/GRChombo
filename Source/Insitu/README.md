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
architecture are not built with Catalyst support at the time of writing 
(latest version: v5.10.0) so it is necessary to build ParaView from source.
Given the large number of dependencies, it is probably easiest to use the 
[superbuild](https://gitlab.kitware.com/paraview/paraview-superbuild/) 
which downloads and builds all dependencies before building ParaView. When
configuring the superbuild using CMake, you will want to make sure the following
options are set
```cmake
ENABLE_osmesa = ON # enables rendering without a screen
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
---
| 5.9.1 |
| 5.10.0 |

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
Note that if the `PARAVIEW_DIR` command is set, you can also use `make` the 
usual `all` target to link with ParaView Catalyst but the executable filename
will be the conventional one (i.e. without '`_Insitu`').

If you wish to stop compiling/linking with ParaView Catalyst, unset the 
`PARAVIEW_DIR` environment variable with
```bash
unset PARAVIEW_DIR
```


## Generating ParaView Catalyst scripts

## Using the in-situ visualization capabilities