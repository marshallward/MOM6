SHELL = bash

# Disable implicit rules
.SUFFIXES:

# Disable implicit variables
MAKEFLAGS += -R

# FMS framework
FMS_URL ?= https://github.com/NOAA-GFDL/FMS.git
FMS_COMMIT ?= 2023.03

GSW_URL ?= https://github.com/mom-ocean/GSW-Fortran.git
GSW_COMMIT ?= 29e64d652786e1d076a05128c920f394202bfe10


# List of source files to link this Makefile's dependencies to model Makefiles
# Assumes a depth of two, and the following extensions: F90 inc c h
# (1): Root directory
SOURCE = \
  $(foreach ext,F90 inc c h,$(wildcard $(1)/*/*.$(ext) $(1)/*/*/*.$(ext)))

FMS_SOURCE = $(call SOURCE,fms/src)
GSW_SOURCE = $(call SOURCE,gsw/src)


# If `true`, print logs if an error is encountered.
REPORT_ERROR_LOGS ?=


#---
# Rules

.PHONY: all
all: lib/libFMS.a
all: lib/libgsw.a


#---
# FMS build

# NOTE: We emulate the automake `make install` stage by storing libFMS.a to
#   ${srcdir}/deps/lib and copying module files to ${srcdir}/deps/include.
lib/libFMS.a: fms/build/libFMS.a
	mkdir -p lib include
	cp fms/build/libFMS.a lib/libFMS.a
	cp fms/build/*.mod include

fms/build/libFMS.a: fms/build/Makefile
	$(MAKE) -C fms/build libFMS.a

fms/build/Makefile: fms/build/Makefile.in fms/build/configure
	cd $(@D) && { \
	  ./configure --srcdir=../src \
	  || { \
	    if [ "${REPORT_ERROR_LOGS}" = true ]; then cat config.log ; fi ; \
	    false; \
	  } \
	}

fms/build/Makefile.in: Makefile.fms.in | fms/build
	cp Makefile.fms.in fms/build/Makefile.in

fms/build/configure: fms/build/configure.ac $(FMS_SOURCE) | fms/src
	autoreconf fms/build

fms/build/configure.ac: configure.fms.ac m4 | fms/build
	cp configure.fms.ac fms/build/configure.ac
	cp -r m4 fms/build

fms/build:
	mkdir -p fms/build

fms/src:
	git clone $(FMS_URL) $@
	git -C $@ checkout $(FMS_COMMIT)


#---
# GSW build
# (TODO: Duplicate of the FMS build above... merge?)

lib/libgsw.a: gsw/build/libgsw.a
	mkdir -p lib include
	cp gsw/build/libgsw.a lib/libgsw.a
	cp gsw/build/*.mod include

gsw/build/libgsw.a: gsw/build/Makefile
	$(MAKE) -C gsw/build libgsw.a

gsw/build/Makefile: gsw/build/Makefile.in gsw/build/configure
	cd $(@D) && { \
	  ./configure --srcdir=../src \
	  || { \
	    if [ "${REPORT_ERROR_LOGS}" = true ]; then cat config.log ; fi ; \
	    false; \
	  } \
	}

gsw/build/Makefile.in: Makefile.gsw.in | gsw/build
	cp Makefile.gsw.in gsw/build/Makefile.in

gsw/build/configure: gsw/build/configure.ac $(GSW_SOURCE) | gsw/src
	autoreconf gsw/build

gsw/build/configure.ac: configure.gsw.ac m4 | gsw/build
	cp configure.gsw.ac gsw/build/configure.ac
	cp -r m4 gsw/build

gsw/build:
	mkdir -p gsw/build

gsw/src:
	git clone $(GSW_URL) $@
	git -C $@ checkout $(GSW_COMMIT)


# Cleanup

.PHONY: clean
clean:
	rm -rf fms/build gsw/build lib include

.PHONY: distclean
distclean: clean
	rm -rf fms
