#=======================================================================
#
# FILE:        Makefile
#
# SUMMARY:     Configurable Makefile for grackle example executable
#              (originally written by James Border for Enzo)
#
# DESCRIPTION
#              See 'gmake help' for definitive description of targets
#
#=======================================================================

# Use bash since sh on datastar does not recognize ">&" used in dep: target

SHELL       = /bin/bash

TOP_DIR     = ../../
OUTPUT      = out.compile

GRACKLE_DIR = $(TOP_DIR)/src/clib
MODULES     = 

#-----------------------------------------------------------------------
# Make.config.settings is used for setting default values to all compile-time 
# configuration settings.
#-----------------------------------------------------------------------

include $(GRACKLE_DIR)/Make.config.settings

#-----------------------------------------------------------------------
# Make.config.machine is used for setting which Make.mach.* file to use
#-----------------------------------------------------------------------

MAKE_CONFIG_MACHINE  = $(GRACKLE_DIR)/Make.config.machine
include $(GRACKLE_DIR)/Make.config.machine

#-----------------------------------------------------------------------
# Make.config.override is used for overriding the default settings in
# Make.config.settings.  This was made separate from the default settings 
# to enable easily interfacing Grackle with a software testing environment 
# like lcatest.
#-----------------------------------------------------------------------

MAKE_CONFIG_OVERRIDE = $(GRACKLE_DIR)/Make.config.override
include $(MAKE_CONFIG_OVERRIDE)

CONFIG_USE_MPI = no

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#-----------------------------------------------------------------------
# Make.config.assemble takes the settings in the Make.config.settings
# and Make.config.override, and generates the appropriate make variables
# required by this makefile.  E.g. $(CXX), $(CXXFLAGS), etc.
#-----------------------------------------------------------------------

include $(GRACKLE_DIR)/Make.config.assemble

#-----------------------------------------------------------------------
# Make.mach.<machine-name> defines all machine-dependent settings.
#-----------------------------------------------------------------------

-include $(GRACKLE_DIR)/Make.mach.$(CONFIG_MACHINE)
-include $(HOME)/.grackle/Make.mach.$(CONFIG_MACHINE)

#-----------------------------------------------------------------------
# add the OpenMP flag
#-----------------------------------------------------------------------

LDFLAGS += $(ASSEMBLE_OMP_FLAGS)

#=======================================================================
# OBJECT FILES
#=======================================================================

OBJS_CONFIG_LIB = cxx_grackle.o

GRACKLE_INCLUDE = -I$(MACH_INSTALL_PREFIX)/include
GRACKLE_LIB = -L$(MACH_INSTALL_PREFIX)/lib -lgrackle

#-----------------------------------------------------------------------
# C++ EXAMPLE
#-----------------------------------------------------------------------

cxx_example: $(MODULES) cxx_example.o
	@rm -f $@
	@echo "Linking"
	-@$(CXX) $(LDFLAGS) -o cxx_example cxx_example.o $(LIBS) $(GRACKLE_LIB) >& $(OUTPUT)
	@(if [ -e $@ ]; then \
             echo "Success!"; \
         else \
             echo "$(CXX) $(LDFLAGS) -o cxx_example cxx_example.o $(LIBS) $(GRACKLE_LIB)" >> temp1; \
             cat temp1 $(OUTPUT) > temp2; \
             rm -f temp1; \
             mv -f temp2 $(OUTPUT); \
             echo "Failed! See $(OUTPUT) for error messages"; \
         fi)

#-----------------------------------------------------------------------
# C++ OPENMP EXAMPLE
#-----------------------------------------------------------------------

cxx_omp_example: $(MODULES) cxx_omp_example.o
	@rm -f $@
	@echo "Linking"
	-@$(CXX) $(LDFLAGS) -o cxx_omp_example cxx_omp_example.o $(LIBS) $(GRACKLE_LIB) >& $(OUTPUT)
	@(if [ -e $@ ]; then \
             echo "Success!"; \
         else \
             echo "$(CXX) $(LDFLAGS) -o cxx_omp_example cxx_omp_example.o $(LIBS) $(GRACKLE_LIB)" >> temp1; \
             cat temp1 $(OUTPUT) > temp2; \
             rm -f temp1; \
             mv -f temp2 $(OUTPUT); \
             echo "Failed! See $(OUTPUT) for error messages"; \
         fi)

#-----------------------------------------------------------------------
# PURE C EXAMPLE
#-----------------------------------------------------------------------

c_example: $(MODULES) c_example.o
	@rm -f $@
	@echo "Linking"
	-@$(LD) $(LDFLAGS) -o c_example c_example.o $(LIBS) $(GRACKLE_LIB) >& $(OUTPUT)
	@(if [ -e $@ ]; then \
             echo "Success!"; \
         else \
             echo "$(LD) $(LDFLAGS) -o c_example c_example.o $(LIBS) $(GRACKLE_LIB)" >> temp1; \
             cat temp1 $(OUTPUT) > temp2; \
             rm -f temp1; \
             mv -f temp2 $(OUTPUT); \
             echo "Failed! See $(OUTPUT) for error messages"; \
         fi)

#-----------------------------------------------------------------------
# FORTRAN EXAMPLE
#-----------------------------------------------------------------------

fortran_example: $(MODULES) fortran_example.o
	@rm -f $@
	@echo "Linking"
	-@$(FC) $(LDFLAGS) -o fortran_example fortran_example.o $(LIBS) $(GRACKLE_LIB) >& $(OUTPUT)
	@(if [ -e $@ ]; then \
             echo "Success!"; \
         else \
             echo "$(FC) $(LDFLAGS) -o fortran_example fortran_example.o $(LIBS) $(GRACKLE_LIB)" >> temp1; \
             cat temp1 $(OUTPUT) > temp2; \
             rm -f temp1; \
             mv -f temp2 $(OUTPUT); \
             echo "Failed! See $(OUTPUT) for error messages"; \
         fi)

#-----------------------------------------------------------------------
# Implicit rules
#-----------------------------------------------------------------------

.SUFFIXES: .c .C .F

.C.o:
	@rm -f $@
	@echo "Compiling $<"
	-@($(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $(GRACKLE_INCLUDE) $*.C) >& $(OUTPUT)
	@(if [ ! -e $@ ]; then \
             echo; \
             echo "$(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $(GRACKLE_INCLUDE) $*.C"; \
             echo; \
             $(CXX) -c -o $@ $(DEFINES) $(CXXFLAGS) $(INCLUDES) $(GRACKLE_INCLUDE) $*.C;\
             echo; \
             exit 1; \
          fi)

.c.o:
	@rm -f $@
	@echo "Compiling $<"
	-@($(CC) -c -o $@ $(DEFINES) $(CFLAGS) $(INCLUDES) $(GRACKLE_INCLUDE) $*.c) >& $(OUTPUT)
	@(if [ ! -e $@ ]; then \
             echo; \
             echo "$(CC) -c -o $@ $(DEFINES) $(CFLAGS) $(INCLUDES) $(GRACKLE_INCLUDE) $*.c"; \
             echo; \
             $(CC) -c -o $@ $(DEFINES) $(CFLAGS) $(INCLUDES) $(GRACKLE_INCLUDE) $*.c;\
             echo; \
             exit 1; \
          fi)

.F.o:
	@rm -f $@
	@echo "Compiling $<"
	-@($(FC) -c -o $@ $(DEFINES) $(FFLAGS) $(INCLUDES) $(GRACKLE_INCLUDE) $*.F) >& $(OUTPUT)
	@(if [ ! -e $@ ]; then \
             echo; \
             echo "$(FC) -c -o $@ $(DEFINES) $(FFLAGS) $(INCLUDES) $(GRACKLE_INCLUDE) $*.F"; \
             echo; \
             $(FC) -c -o $@ $(DEFINES) $(FFLAGS) $(INCLUDES) $(GRACKLE_INCLUDE) $*.F;\
             echo; \
             exit 1; \
          fi)

#-----------------------------------------------------------------------
# HELP TARGET
#-----------------------------------------------------------------------

help:
	@echo
	@echo "========================================================================"
	@echo "   Makefile Help for Grackle example"
	@echo "========================================================================"
	@echo
	@echo "   gmake                Compile and generate the executable 'grackle.exe'"
	@echo "   gmake help           Display this help information"
	@echo "   gmake clean          Remove object files, executable, etc."
	@echo
	@echo "   gmake help-config    Display detailed help on configuration make targets"
	@echo "   gmake show-config    Display the configuration settings"
	@echo "   gmake show-flags     Display specific compilation flags"
	@echo "   gmake default        Reset the configuration to the default values"
	@echo

#-----------------------------------------------------------------------

clean:
	-@rm -f *.o *.mod *.f *.f90 *~ *.exe $(OUTPUT) cxx_example cxx_omp_example c_example fortran_example

#-----------------------------------------------------------------------
# Include configuration targets
#-----------------------------------------------------------------------

include $(GRACKLE_DIR)/Make.config.targets
