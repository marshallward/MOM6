SHELL = bash

# Disable implicit rules
.SUFFIXES:

# Disable implicit variables
MAKEFLAGS += -R


# mkmf, list_paths (GFDL build toolchain)
MKMF_URL ?= https://github.com/NOAA-GFDL/mkmf.git
MKMF_COMMIT ?= master

# FMS framework
FMS_URL ?= https://github.com/NOAA-GFDL/FMS.git
FMS_COMMIT ?= 2019.01.03

# CVMix build
CVMIX_URL ?= https://github.com/CVMix/CVMix-src
CVMIX_COMMIT ?= v0.92-beta


# List of source files to link this Makefile's dependencies to model Makefiles
# Assumes a depth of two, and the following extensions: F90 inc c h
# (1): Root directory
# NOTE: extensions could be a second variable
SOURCE = \
  $(foreach ext,F90 inc c h,$(wildcard $(1)/*/*.$(ext) $(1)/*/*/*.$(ext)))

FMS_SOURCE = $(call SOURCE,fms/src)
CVMIX_SOURCE = $(call SOURCE,cvmix/src/shared)


#---
# Rules

.PHONY: all
all: bin/mkmf bin/list_paths lib/libFMS.a lib/libcvmix.a

#---
# mkmf checkout

bin/mkmf bin/list_paths: mkmf
	mkdir -p $(@D)
	cp $^/$@ $@

mkmf:
	git clone $(MKMF_URL) $@
	git -C $@ checkout $(MKMF_COMMIT)


#---
# FMS build

# NOTE: We emulate the automake `make install` stage by storing libFMS.a to
#   ${srcdir}/deps/lib and copying module files to ${srcdir}/deps/include.
#   This is a flawed approach, since module files are untracked and could be
#   handled more safely, but this is adequate for now.


# TODO: track *.mod copy?
lib/libFMS.a: fms/build/libFMS.a fms/build/Makefile
	mkdir -p {lib,include}
	cp fms/build/libFMS.a lib/libFMS.a
	cp fms/build/*.mod include


fms/build/libFMS.a: fms/build/Makefile
	make -C fms/build libFMS.a


fms/build/Makefile: Makefile.fms.in fms/src/configure bin/mkmf bin/list_paths
	mkdir -p fms/build
	cp Makefile.fms.in fms/src/Makefile.in
	cd $(@D) && ../src/configure --srcdir=../src


# TODO: Track m4 macros?
fms/src/configure: configure.fms.ac $(FMS_SOURCE) | fms/src
	cp configure.fms.ac fms/src/configure.ac
	cp -r m4 $(@D)
	cd $(@D) && autoreconf -i

fms/src:
	git clone $(FMS_URL) $@
	git -C $@ checkout $(FMS_COMMIT)


#---
# CVMix build
# NOTE: This may not be a complete CVMix library; only the bits used by MOM6.

# TODO: track *.mod copy?
lib/libcvmix.a: cvmix/build/libcvmix.a cvmix/build/Makefile
	mkdir -p {lib,include/cvmix}
	cp cvmix/build/libcvmix.a lib/libcvmix.a
	cp cvmix/build/*.mod include/cvmix

cvmix/build/libcvmix.a: cvmix/build/Makefile
	make -C cvmix/build libcvmix.a

cvmix/build/Makefile: Makefile.fms.in cvmix/src/shared/configure bin/mkmf bin/list_paths
	mkdir -p cvmix/build
	cp Makefile.fms.in cvmix/src/shared/Makefile.in
	cd $(@D) && $(FMS_ENV) ../src/shared/configure --srcdir=../src/shared

cvmix/src/shared/configure: configure.cvmix.ac $(CVMIX_SOURCE) | cvmix
	cp configure.cvmix.ac cvmix/src/shared/configure.ac
	cp -r m4 $(@D)
	cd $(@D) && autoreconf -if

# TODO: Merge with FMS rule?
cvmix:
	git clone $(CVMIX_URL) $@
	git -C $@ checkout $(CVMIX_COMMIT)


.PHONY: clean
clean:
	rm -rf fms/build cvmix/build lib include bin

.PHONY: distclean
distclean: clean
	rm -rf fms cvmix mkmf
