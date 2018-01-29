TestDirs := $(wildcard Tests/*/.)
ExampleDirs := $(wildcard Examples/*/.)

.PHONY: all run $(TestDirs) $(ExampleDirs)


ifndef GRCHOMBO_SOURCE
    $(error Please define GRCHOMBO_SOURCE - see installation instructions.)
endif

test: $(TestDirs)

examples: $(ExampleDirs)

all: $(TestDirs) $(ExampleDirs)

$(TestDirs):
	$(info ################# Making test $@ #################)
	$(MAKE) -C $@ all
	$(info ################# Running test $@ #################)
	$(MAKE) -C $@ run

$(ExampleDirs):
	$(info ################# Making example $@ #################)
	$(MAKE) -C $@ all
