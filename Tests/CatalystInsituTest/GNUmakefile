# -*- Mode: Makefile -*-

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally (e.g. export CHOMBO_HOME=... in bash)

GRCHOMBO_SOURCE = $(shell pwd)/../../Source

ebase := CatalystInsituTest

LibNames := AMRTimeDependent AMRTools BoxTools

src_dirs := $(GRCHOMBO_SOURCE)/utils \
            $(GRCHOMBO_SOURCE)/simd  \
            $(GRCHOMBO_SOURCE)/BoxUtils  \
            $(GRCHOMBO_SOURCE)/CCZ4  \
            $(GRCHOMBO_SOURCE)/GRChomboCore  \
            $(GRCHOMBO_SOURCE)/AMRInterpolator

XTRA_VTK_COMPONENT_FLAGS = -v ImagingCore -v ImagingSources -v IOImage
XTRA_PARAVIEW_COMPONENT_FLAGS = -c RemotingCore

include $(GRCHOMBO_SOURCE)/Insitu/Make.insitu

include $(CHOMBO_HOME)/mk/Make.test
