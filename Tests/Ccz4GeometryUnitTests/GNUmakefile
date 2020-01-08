# -*- Mode: Makefile -*-

### This makefile produces an executable for each name in the `ebase'
###  variable using the libraries named in the `LibNames' variable.

# Included makefiles need an absolute path to the Chombo installation
# CHOMBO_HOME := Please set the CHOMBO_HOME locally(e.g. export CHOMBO_HOME=... in bash)

GRCHOMBO_SOURCE = $(shell pwd)/../../Source

ebase := CCZ4GeometryUnitTest

LibNames := BoxTools

src_dirs := $(GRCHOMBO_SOURCE)/utils $(GRCHOMBO_SOURCE)/CCZ4

include $(CHOMBO_HOME)/mk/Make.test
