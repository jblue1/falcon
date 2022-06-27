# Makefile is a part of the PYTHIA event generator.
# Copyright (C) 2020 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.
# Author: Philip Ilten, September 2014.
#
# This is is the Makefile used to build PYTHIA examples on POSIX systems.
# Example usage is:
#     make main01
# For help using the make command please consult the local system documentation,
# i.e. "man make" or "make --help".

################################################################################
# VARIABLES: Definition of the relevant variables from the configuration script.
################################################################################

# Set the shell.

SHELL=/usr/bin/env bash
#PYTHIAROOT=/home/DAVIDSON/joblue/Pythia/pythia8230
PYTHIAROOT=/home/DAVIDSON/joblue/phy2/pythia8303/
#PYTHIAROOT=/home/DAVIDSON/joblue/delphes-issue/pythia8235
FASTJETFLAGS = `$(FASTJETSYS)/bin/fastjet-config --cxxflags --libs --plugins`
TENSORFLOWFLAGS = -I /home/DAVIDSON/joblue/libtensorflow2/include -ltensorflow
CPPFLOWFLAGS = -I /home/DAVIDSON/joblue/phy2/cppflow/include


# Include the configuration.
-include $(PYTHIAROOT)/examples/Makefile.inc

#GZIP_INC=
#GZIP_FLAGS=
#ifeq ($(GZIP_USE),true)
#  GZIP_INC+= -DGZIPSUPPORT -I$(GZIP_INCLUDE)
#  GZIP_FLAGS+= -L$(GZIP_LIB) -Wl,-rpath,$(GZIP_LIB) -lz
#endif

# Check distribution (use local version first, then installed version).
#ifneq ("$(wildcard ../lib/libpythia8.*)","")
#  PREFIX_LIB=../lib
#  PREFIX_INCLUDE=../include
#endif
ifneq ("$(wildcard ../lib/libpythia8.*)","")
  PREFIX_LIB=$(PYTHIAROOT)/lib
  PREFIX_INCLUDE=$(PYTHIAROOT)/include
endif

CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON) #$(GZIP_LIB)
CXX_COMMON+=-L$(PREFIX_LIB) -Wl,-rpath,$(PREFIX_LIB) -lpythia8 -ldl
CXX_COMMON+=$(FASTJETFLAGS) $(TENSORFLOWFLAGS) --std=c++17
CXX_COMMON+=$(CPPFLOWFLAGS) 
PYTHIA=$(PREFIX_LIB)/libpythia8.a# $(LIB_SUFFIX)

################################################################################
# RULES: Definition of the rules used to build the PYTHIA examples.
################################################################################

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

# All targets (no default behavior).
all:
	$(info Usage: make mainXX)

# PYTHIA library.
$(PYTHIA):
	$(error Error: PYTHIA must be built, please run "make"\
                in the top PYTHIA directory)

# Examples without external dependencies.
main%: $(PYTHIA) main%.cc
	$(CXX) $@.cc -o $@ $(CXX_COMMON)

# MixMax (remove warnings).
main23: $(PYTHIA) main23.cc
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON)

# HEPMC3.
main41 main42 main43 main45 main85 main86 main87 main88: $(PYTHIA) $$@.cc
ifeq ($(HEPMC3_USE),true)
	$(CXX) $@.cc -o $@ $(HEPMC3_INCLUDE) $(CXX_COMMON) $(HEPMC3_LIB)
else
	$(error Error: $@ requires HEPMC3)
endif

# HEPMC2.
main44: $(PYTHIA) $$@.cc
ifeq ($(HEPMC2_USE),true)
	$(CXX) $@.cc -o $@ $(HEPMC2_INCLUDE) $(CXX_COMMON) $(HEPMC2_LIB)
else
	$(error Error: $@ requires HEPMC2)
endif

# HEPMC2 or HEPMC3 (prefer HEPMC3 if both)
main89: $(PYTHIA) $$@.cc
ifeq ($(HEPMC3_USE),true)
	$(CXX) $@.cc -o $@ $(HEPMC3_INCLUDE) $(CXX_COMMON) $(HEPMC3_LIB) -DHEPMC3
else ifeq ($(HEPMC2_USE),true)
	$(CXX) $@.cc -o $@ $(HEPMC2_INCLUDE) $(CXX_COMMON) $(HEPMC2_LIB) -DHEPMC2
else
	$(error Error: $@ requires HEPMC2 or HEPMC3)
endif

# EVTGEN.
main48: $(PYTHIA) $$@.cc
ifeq ($(EVTGEN_USE)$(HEPMC2_USE),truetrue)
	$(CXX) $@.cc -o $@ $(EVTGEN_INCLUDE) -w $(CXX_COMMON)\
	 -DEVTGEN_PYTHIA -DEVTGEN_EXTERNAL $(EVTGEN_LIB)
else
	$(error Error: $@ requires EVTGEN)
endif

# FASTJET3.
main71 main72 main75 main81 main82 main83: $(PYTHIA) $$@.cc
ifeq ($(FASTJET3_USE),true)
	$(CXX) $@.cc -o $@ $(FASTJET3_INCLUDE) -w $(CXX_COMMON) $(FASTJET3_LIB)
else
	$(error Error: $@ requires FASTJET3)
endif

# FASTJET3 with recursive tools.
main74: $(PYTHIA) $$@.cc
ifeq ($(FASTJET3_USE),true)
	$(CXX) $@.cc -o $@ $(FASTJET3_INCLUDE) -w $(CXX_COMMON)\
	 $(FASTJET3_LIB) -lfastjettools -lRecursiveTools
else
	$(error Error: $@ requires FASTJET3)
endif

# FASTJET3 and HEPMC3.
main84: $(PYTHIA) $$@.cc
ifeq ($(FASTJET3_USE)$(HEPMC3_USE),truetrue)
	$(CXX) $@.cc -o $@ $(FASTJET3_INCLUDE) $(HEPMC3_INCLUDE) -w\
	 $(CXX_COMMON) $(FASTJET3_LIB) $(HEPMC3_LIB)
else
	$(error Error: $@ requires FASTJET3 and HEPMC3)
endif

# ROOT libraries generated via CINT.
main%.so: main%Dct.cc
	$(CXX) $< -o $@ -w $(CXX_SHARED) $(CXX_COMMON)\
	 $(ROOT_LIB) `$(ROOT_CONFIG) --cflags`
main%Dct.cc: main%.h main%LinkDef.h
ifeq ($(ROOT_USE),true)
	$(ROOT_BIN)rootcint -f $@ -I$(PREFIX_INCLUDE) $^
else
	$(error Error: $@ requires ROOT)
endif

# ROOT (turn off all warnings for readability).
main91: $(PYTHIA) $$@.cc
ifeq ($(ROOT_USE),true)
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON) $(ROOT_LIB)\
	 `$(ROOT_CONFIG) --cflags --glibs`
else
	$(error Error: $@ requires ROOT)
endif
main92: $(PYTHIA) $$@.cc main92.so
	$(CXX) $@.cc main92.so -o $@ -w $(CXX_COMMON) -Wl,-rpath,./\
	 `$(ROOT_CONFIG) --cflags --glibs`

# RIVET with optional ROOT (if RIVET, use C++14).
main93: $(PYTHIA) $$@.cc $(if $(filter true,$(ROOT_USE)),main93.so)
ifeq ($(RIVET_USE),true)
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON:c++11=c++14) -Wl,-rpath,./\
	 $(if $(filter true,$(ROOT_USE)),main93.so -DPY8ROOT\
         $(ROOT_LIB) `$(ROOT_CONFIG) --cflags --glibs`)\
	 $(RIVET_INCLUDE) $(RIVET_LIB)
else
	$(error Error: $@ requires RIVET)
endif

# RIVET (if RIVET use C++14).
main111: $(PYTHIA) $$@.cc
ifeq ($(RIVET_USE),true)
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON:c++11=c++14) $(RIVET_INCLUDE)\
	 $(RIVET_LIB) -DRIVET
else
	$(CXX) $@.cc -o $@ -w $(CXX_COMMON)
endif

# HEPMC3.
main300: $(PYTHIA) $$@.cc
ifeq ($(HEPMC3_USE),true)
	$(CXX) $@.cc -o $@ $(HEPMC3_INCLUDE) $(CXX_COMMON) $(HEPMC3_LIB)\
	 -DHEPMC3
else
	$(CXX) $@.cc -o $@ $(CXX_COMMON)
endif

# User-written examples for tutorials, without external dependencies.
mymain%: $(PYTHIA) mymain%.cc
	$(CXX) $@.cc -o $@ $(CXX_COMMON)

# Internally used tests, without external dependencies.
test%: $(PYTHIA) test%.cc
	$(CXX) $@.cc -o $@ $(CXX_COMMON)

# Clean.
clean:
	@rm -f main[0-9][0-9]; rm -f out[0-9][0-9];\
	rm -f main[0-9][0-9][0-9]; rm -f out[0-9][0-9][0-9];\
	rm -f mymain[0-9][0-9]; rm -f myout[0-9][0-9];\
	rm -f test[0-9][0-9][0-9]; rm -f *.dat;\
	rm -f weakbosons.lhe; rm -f hist.root;\
	rm -f *~; rm -f \#*; rm -f core*; rm -f *Dct.*; rm -f *.so;\
	rm -f *.log; rm -f *plot.py; rm -f *.pcm;
