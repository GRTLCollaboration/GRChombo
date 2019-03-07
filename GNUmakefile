TestDirs := $(wildcard Tests/*/.)
ExampleDirs := $(wildcard Examples/*/.)
CleanTestDirs := $(TestDirs:%=clean-%)
CleanExampleDirs := $(ExampleDirs:%=clean-%)
RealCleanTestDirs := $(TestDirs:%=realclean-%)
RealCleanExampleDirs := $(ExampleDirs:%=realclean-%)

.PHONY: all run $(TestDirs) $(ExampleDirs)


export GRCHOMBO_SOURCE = $(shell pwd)/Source

test: $(TestDirs)

examples: $(ExampleDirs)

all: $(TestDirs) $(ExampleDirs)

clean: $(CleanTestDirs) $(CleanExampleDirs)

realclean: $(RealCleanTestDirs) $(RealCleanExampleDirs)

$(TestDirs):
	$(info ################# Making test $@ #################)
	$(MAKE) -C $@ all
	$(info ################# Running test $@ #################)
	$(MAKE) -C $@ run

$(ExampleDirs):
	$(info ################# Making example $@ #################)
	$(MAKE) -C $@ all

$(CleanTestDirs):
	$(MAKE) -C $(@:clean-%=%) clean NODEPENDS=TRUE

$(CleanExampleDirs):
	$(MAKE) -C $(@:clean-%=%) clean NODEPENDS=TRUE

$(RealCleanTestDirs):
	$(MAKE) -C $(@:realclean-%=%) clean NODEPENDS=TRUE

$(RealCleanExampleDirs):
	$(MAKE) -C $(@:realclean-%=%) clean NODEPENDS=TRUE
