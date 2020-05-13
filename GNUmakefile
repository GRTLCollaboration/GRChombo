define SEARCH_FOR_MAKE
$(wildcard $1/*/GNUmakefile) $(wildcard $1/*/makefile) $(wildcard $1/*/Makefile)
endef

TestDirsWithGNUmakefile := $(call SEARCH_FOR_MAKE, Tests)
TestDirs := $(dir $(TestDirsWithGNUmakefile))
RunTestDirs := $(TestDirs:%=run-%)
ExampleDirsWithGNUmakefile := $(call SEARCH_FOR_MAKE, Examples)
ExampleDirs := $(dir $(ExampleDirsWithGNUmakefile))
CleanTestDirs := $(TestDirs:%=clean-%)
CleanExampleDirs := $(ExampleDirs:%=clean-%)
RealCleanTestDirs := $(TestDirs:%=realclean-%)
RealCleanExampleDirs := $(ExampleDirs:%=realclean-%)

.PHONY: all run examples clean realclean $(TestDirs) $(ExampleDirs)

export GRCHOMBO_SOURCE = $(shell pwd)/Source

test: $(TestDirs)

run: test $(RunTestDirs)

examples: $(ExampleDirs)

all: run examples

clean: $(CleanTestDirs) $(CleanExampleDirs)

realclean: $(RealCleanTestDirs) $(RealCleanExampleDirs)

$(TestDirs):
	$(info ################# Making test $@ #################)
	$(MAKE) -C $@ all

$(RunTestDirs):
	$(info ################# Running test $@ #################)
	$(MAKE) -C $(@:run-%=%) run

$(ExampleDirs):
	$(info ################# Making example $@ #################)
	$(MAKE) -C $@ all

$(CleanTestDirs):
	$(MAKE) -C $(@:clean-%=%) clean NODEPENDS=TRUE

$(CleanExampleDirs):
	$(MAKE) -C $(@:clean-%=%) clean NODEPENDS=TRUE

$(RealCleanTestDirs):
	$(MAKE) -C $(@:realclean-%=%) realclean NODEPENDS=TRUE

$(RealCleanExampleDirs):
	$(MAKE) -C $(@:realclean-%=%) realclean NODEPENDS=TRUE
