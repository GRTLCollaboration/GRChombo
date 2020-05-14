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

ECHO?=@ # set this to null on the command line to increase verbosity

test: $(TestDirs)

run: test $(RunTestDirs)

examples: $(ExampleDirs)

all: run examples

clean: $(CleanTestDirs) $(CleanExampleDirs)

realclean: $(RealCleanTestDirs) $(RealCleanExampleDirs)

$(TestDirs):
	$(info ################# Making test $@ #################)
	$(ECHO)$(MAKE) -C $@ --no-print-directory all

$(RunTestDirs):
	$(info ################# Running test $@ #################)
	$(ECHO)$(MAKE) -C $(@:run-%=%) --no-print-directory run

$(ExampleDirs):
	$(info ################# Making example $@ #################)
	$(ECHO)$(MAKE) -C $@ --no-print-directory all

$(CleanTestDirs):
	$(ECHO)$(MAKE) -C $(@:clean-%=%) --no-print-directory clean NODEPENDS=TRUE

$(CleanExampleDirs):
	$(ECHO)$(MAKE) -C $(@:clean-%=%) --no-print-directory clean NODEPENDS=TRUE

$(RealCleanTestDirs):
	$(ECHO)$(MAKE) -C $(@:realclean-%=%) --no-print-directory realclean NODEPENDS=TRUE

$(RealCleanExampleDirs):
	$(ECHO)$(MAKE) -C $(@:realclean-%=%) --no-print-directory realclean NODEPENDS=TRUE
