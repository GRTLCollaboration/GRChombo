# -*- Mode: Makefile -*-

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally (e.g. export CHOMBO_HOME=... in bash)

GRCHOMBO_SOURCE = ../../Source
AH_FINDER_SOURCE = $(GRCHOMBO_SOURCE)/ApparentHorizonFinder

ifdef AH_FINDER_SOURCE
ifneq (,$(wildcard $(AH_FINDER_SOURCE)/AHFinder.hpp))
     USE_AH_FINDER := TRUE
else
  $(error AH_FINDER_SOURCE not set correctly)
endif
endif

ebase := Main_KerrBH

LibNames := AMRTimeDependent AMRTools BoxTools

src_dirs := $(GRCHOMBO_SOURCE)/utils \
            $(GRCHOMBO_SOURCE)/simd  \
            $(GRCHOMBO_SOURCE)/CCZ4  \
            $(GRCHOMBO_SOURCE)/BoxUtils  \
            $(GRCHOMBO_SOURCE)/GRChomboCore  \
            $(GRCHOMBO_SOURCE)/TaggingCriteria  \
            $(GRCHOMBO_SOURCE)/InitialConditions/BlackHoles \
            $(GRCHOMBO_SOURCE)/BlackHoles \
            $(GRCHOMBO_SOURCE)/AMRInterpolator
            
ifeq ($(USE_AH_FINDER),TRUE)
  src_dirs += $(AH_FINDER_SOURCE)
endif

include $(CHOMBO_HOME)/mk/Make.test
