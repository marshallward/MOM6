# Part of MOM6, the Modular Ocean Model version 6
# SPDX-License-Identifier: Apache-2.0

# Build configuration
BUILD ?= build
FMS_BUILD ?= ac/deps/fms/build
MOM_MEMORY ?=


#----

# Makefile setup

# Verify that BUILD is not set to the current directory
# (which would clobber this Makefile)
MAKEPATH = $(realpath $(dir $(abspath $(lastword $(MAKEFILE_LIST)))))
ifeq ($(MAKEPATH), $(realpath $(BUILD)))
  $(error BUILD cannot be set to the current directory)
endif

# Disable builtin rules and variables
MAKEFLAGS += -rR


#----

.PHONY: all
all: $(BUILD)/MOM6

$(BUILD)/MOM6: $(BUILD)/Makefile $(FMS_BUILD)/libFMS.a
	if test $(FMS_BUILD)/libFMS.a -nt $@ ; then \
	  $(MAKE) -C $(BUILD) clean ; \
	fi
	$(MAKE) -C $(BUILD) MOM6

$(BUILD)/Makefile: $(BUILD)/config.status ac/Makefile.in
	cd $(BUILD) && ./config.status

# NOTE: configure is assumed to exist
$(BUILD)/config.status: $(FMS_BUILD)/libFMS.a | $(BUILD)
	cd $(BUILD) && \
	PATH="${PATH}:$(CURDIR)/ac" \
	$(CURDIR)/configure -n $(CONFIG_FLAGS)

#---
# ./configure setup

configure: ac/configure.ac ac/aclocal.m4 | ac/m4 ac/autom4te.cache
	cd ac && autoconf -o $(abspath $@)

ac/aclocal.m4: ac/configure.ac | ac/m4
	cd ac && aclocal

$(BUILD):
	mkdir -p $@

#----
# Dependencies

# NOTE: If libFMS has changed, then we completely rebuild MOM6
$(FMS_BUILD)/libFMS.a: FORCE
	$(MAKE) -C ac/deps \
	  BUILD=$(abspath $(FMS_BUILD)) \
	  CODEBASE=$(abspath $(FMS_CODEBASE))

FORCE:

#----
# Cleanup

.PHONY: clean
clean:
	rm -rf $(BUILD)
	$(MAKE) -C ac/deps clean
