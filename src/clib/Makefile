#=======================================================================
#
# FILE:        Makefile
#
# SUMMARY:     Configurable Makefile for grackle
#              (originally written by James Border for Enzo)
#
# DESCRIPTION
#              See 'make help' for definitive description of targets
#
#              Makefile.config includes the following files:
# 
#              Make.config.settings   default configuration settings
#              Make.config.override   local user settings
#              Make.config.assemble   maps 'config' settings to 'flag' settings
#              Make.config.targets    configure targets
#              Make.mach.*            all machine-dependent settings
#              Make.config.objects    list of object files
#              DEPEND                 Make-generated dependencies
#
#              Make.mach.* should be the only file that one should
#              edit when porting grackle to a new machine.
#
#              Make.config.override should be the only file that
#              one should edit when defining local user settings.
#              Preferably, this is done implicitly through
#              the available make targets (e.g. "make precision-32").
#              Use "make help-config" for a list of all configuration
#              settings.  These make targets do error-checking; hand-editing 
#              Make.config.override is more error-prone.
#
#=======================================================================

SHELL       = /bin/bash
OUTPUT      = out.compile
VERBOSE     = 0
GRACKLE_DIR = .
DEFAULT_INSTALL_PREFIX = /usr/local

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
# to enable easily interfacing grackle with a software testing environment 
# like lcatest.
#-----------------------------------------------------------------------

MAKE_CONFIG_OVERRIDE = $(GRACKLE_DIR)/Make.config.override

include $(MAKE_CONFIG_OVERRIDE)

#-----------------------------------------------------------------------
# Make.mach.<machine-name> defines all machine-dependent settings.
#
# mqk 03/22/2013: moved this to before Make.config.assemble, so the
# MACH_INSTALL_* variables are defined there.
#-----------------------------------------------------------------------

-include $(GRACKLE_DIR)/Make.mach.$(CONFIG_MACHINE)
-include $(HOME)/.grackle/Make.mach.$(CONFIG_MACHINE)

#-----------------------------------------------------------------------
# OTHER_DEFINES is used to specify a macro indicating that legacy
# functions defined inline in the main header file should be omitted
# while compiling the library
#-----------------------------------------------------------------------
OTHER_DEFINES =

#-----------------------------------------------------------------------
# Make.config.assemble takes the settings in the Make.config.settings
# and Make.config.override, and generates the appropriate make variables
# required by this makefile.  E.g. $(CXX), $(CXXFLAGS), etc.
#-----------------------------------------------------------------------

include Make.config.assemble

#
# If libtool is not defined, set it to a default value
#

LIBTOOL ?= libtool

#=======================================================================
# OBJECT FILES
#=======================================================================

include Make.config.objects

#-----------------------------------------------------------------------
# MAKE LIBGRACKLE BY DEFAULT
#-----------------------------------------------------------------------

libgrackle.la: $(MODULES) autogen dep $(OBJS_LIB) MACHNOTES
	@rm -f $@
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  echo "Linking" ; \
	  $(LIBTOOL) --mode=link --tag=CC $(LD) $(LDFLAGS) -o $@ $(OBJS_LIB) $(LIBS) -rpath $(INSTALL_LIB_DIR) -release $(LIB_RELEASE_VERSION) >> $(OUTPUT) 2>&1 ; \
	  if [ -e $@ ]; then \
	     echo "Success!"; \
	  else \
	     echo "Failed! See $(OUTPUT) for error messages"; \
	fi ; \
	else \
	  $(LIBTOOL) --mode=link --tag=CC $(LD) $(LDFLAGS) -o $@ $(OBJS_LIB) $(LIBS) -rpath $(INSTALL_LIB_DIR) -release $(LIB_RELEASE_VERSION) ; \
	fi)
	@(if [ $(ASSEMBLE_PRECISION_NUMBER) == "4" ]; then \
	  echo "WARNING: Grackle compiled with precision-32."; \
	  echo "Using Grackle in single precision is known to produce unreliable results in certain conditions. Compiling in double precision (precision-64) is recommended."; \
	fi)

MACHNOTES: 
	@echo $(MACHINE_NOTES)

#-----------------------------------------------------------------------
# WRITE ALL COMPILER OUTPUT TO FILE
#-----------------------------------------------------------------------

.PHONY: verbose
verbose: VERBOSE = 1

#-----------------------------------------------------------------------
# Implicit rules
#-----------------------------------------------------------------------

.SUFFIXES: .c .C .F .lo .o

# Inhibit removing any *.o files after compiling

.PRECIOUS: %.o %.lo %.c %.C %.F

.F.lo:
	@rm -f $@
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  echo "Compiling $<" ; \
	  $(LIBTOOL) --mode=compile --tag=FC $(FC) -c $(FFLAGS) $(DEFINES) $*.F >& $(OUTPUT) ; \
	  if [ ! -e $@ ]; then \
             echo; \
	     echo "Compiling $< failed!"; \
	     echo "See $(OUTPUT) for error messages."; \
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(LIBTOOL) --mode=compile --tag=FC $(FC) -c $(FFLAGS) $(DEFINES) $*.F; \
	  if [ ! -e $@ ]; then \
	     exit 1; \
	  fi ; \
	fi)

.C.lo:
	@rm -f $@
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  echo "Compiling $<" ; \
	  $(LIBTOOL) --mode=compile --tag=CXX $(CXX) -c $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C \
	    >& $(OUTPUT) ; \
	  if [ ! -e $@ ]; then \
             echo; \
	     echo "Compiling $< failed!"; \
	     echo "See $(OUTPUT) for error messages."; \
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(LIBTOOL) --mode=compile --tag=CXX $(CXX) -c $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C; \
	  if [ ! -e $@ ]; then \
	     exit 1; \
	  fi ; \
	fi)

.c.lo:
	@rm -f $@
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  echo "Compiling $<" ; \
	  $(LIBTOOL) --mode=compile --tag=CC $(CC) -c $(DEFINES) $(CFLAGS) $(INCLUDES) $*.c \
	    >& $(OUTPUT) ; \
	  if [ ! -e $@ ]; then \
             echo; \
	     echo "Compiling $< failed!"; \
	     echo "See $(OUTPUT) for error messages."; \
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(LIBTOOL) --mode=compile --tag=CC $(CC) -c $(DEFINES) $(CFLAGS) $(INCLUDES) $*.c; \
	  if [ ! -e $@ ]; then \
	     exit 1; \
	  fi ; \
	fi)

#-----------------------------------------------------------------------
# Generate all make-generated source files
#-----------------------------------------------------------------------

.PHONY: autogen
autogen: config_type auto_show_config.c auto_show_flags.c auto_get_version.c

.PHONY: config_type
config_type:
	-@echo "#ifndef __GRACKLE_FLOAT_H__" > grackle_float.h
	-@echo "#define __GRACKLE_FLOAT_H__" >> grackle_float.h
	-@echo "#define $(PRECISION_DEFINE)" >> grackle_float.h
	-@echo "#endif" >> grackle_float.h

# Force update of auto_show_config.c

.PHONY: auto_show_config.c
auto_show_config.c:
	-@$(MAKE) -s show-config  >& temp.show-config
	-@awk 'BEGIN {print "#include <stdio.h>\nvoid auto_show_config(FILE *fp) {"}; {print "   fprintf (fp,\""$$0"\\n\");"}; END {print "}"}' < temp.show-config > auto_show_config.c

# Force update of auto_show_flags.c

.PHONY: auto_show_flags.c
auto_show_flags.c:
	-@$(MAKE) -s show-flags  >& temp.show-flags
	-@awk 'BEGIN {print "#include <stdio.h>\nvoid auto_show_flags(FILE *fp) {"}; {print "   fprintf (fp,\""$$0"\\n\");"}; END {print "}"}' < temp.show-flags > auto_show_flags.c

# Force update of auto_get_version.c

.PHONY: auto_get_version.c
auto_get_version.c:
	-@$(MAKE) -s show-version  >& temp.show-version
	-@(if [ -f $@ ]; then rm $@; fi) # delete the file if it already exists

	-@echo '#include <stdio.h>'                                                             >> $@
	-@echo '#include "grackle_types.h"'                                                     >> $@
	-@echo ''                                                                               >> $@
	-@echo '// the following macros are auto-generated:'                                    >> $@

	-@awk '{ if (NF > 0) { print "#define AUTO_" toupper($$(NF-1)) " \"" $$(NF) "\"" } };' < temp.show-version >> auto_get_version.c

	-@echo ''                                                                               >> $@
	-@echo '// test that ensures that all macros were correctly defined:'                   >> $@
	-@echo '#if !(defined(AUTO_VERSION) && defined(AUTO_BRANCH) && defined(AUTO_REVISION))' >> $@
	-@echo '#error "Something went wrong while auto-generating macros"'                     >> $@
	-@echo '#endif'                                                                         >> $@
	-@echo ''                                                                               >> $@
	-@echo 'grackle_version get_grackle_version(void) {'                                    >> $@
	-@echo '  grackle_version out;'                                                         >> $@
	-@echo '  out.version = AUTO_VERSION;'                                                  >> $@
	-@echo '  out.branch = AUTO_BRANCH;'                                                    >> $@
	-@echo '  out.revision = AUTO_REVISION;'                                                >> $@
	-@echo '  return out;'                                                                  >> $@
	-@echo '}'                                                                              >> $@

#-----------------------------------------------------------------------
# Generate dependency file
#-----------------------------------------------------------------------

.PHONY: dep
dep:
	@echo "Updating DEPEND"
	-@(makedepend $(DEFINES) $(INCLUDES) -a -fDEPEND -o.o -m -- -- *.C) >> out.make.DEPEND 2>&1
	-@(makedepend $(DEFINES) $(INCLUDES) -a -fDEPEND -o.o -m -- -- *.c) >> out.make.DEPEND 2>&1
	-@(makedepend $(DEFINES) $(INCLUDES) -a -fDEPEND -o.o -m -- -- *.F) >> out.make.DEPEND 2>&1
	-@(makedepend $(DEFINES) $(INCLUDES) -a -fDEPEND -o.o -m -- -- *.h) >> out.make.DEPEND 2>&1

include DEPEND

#-----------------------------------------------------------------------
# HELP TARGET
#-----------------------------------------------------------------------

help:
	@echo
	@echo "========================================================================"
	@echo "   Grackle Makefile Help"
	@echo "========================================================================"
	@echo
	@echo "   make                Compile and generate librackle"
	@echo "   make install        Copy the library somewhere"
	@echo "   make help           Display this help information"
	@echo "   make clean          Remove object files, executable, etc."
	@echo "   make dep            Create make dependencies in DEPEND file"
	@echo
	@echo "   make show-version   Display revision control system branch and revision"
	@echo "   make show-diff      Display local file modifications"
	@echo
	@echo "   make help-config    Display detailed help on configuration make targets"
	@echo "   make show-config    Display the configuration settings"
	@echo "   make show-flags     Display specific compilation flags"
	@echo "   make default        Reset the configuration to the default values"
	@echo

#-----------------------------------------------------------------------
# INSTALLATION TARGET
#-----------------------------------------------------------------------

install:
	@echo "Installing grackle header files to $(INSTALL_INCLUDE_DIR)."
	@(if [ ! -d $(INSTALL_INCLUDE_DIR) ]; then \
		mkdir $(INSTALL_INCLUDE_DIR); \
	fi)
	@cp grackle.h grackle_macros.h grackle_float.h grackle_types.h grackle_chemistry_data.h grackle_rate_functions.h grackle.def grackle_fortran_types.def grackle_fortran_interface.def $(INSTALL_INCLUDE_DIR)
	@(if [ ! -d $(INSTALL_LIB_DIR) ]; then \
		mkdir $(INSTALL_LIB_DIR); \
	fi)
	@echo "Installing grackle library files to $(INSTALL_LIB_DIR)."
	$(LIBTOOL) --mode=install install -c libgrackle.la $(INSTALL_LIB_DIR)/libgrackle.la

#-----------------------------------------------------------------------

clean:
	-@rm -f *.la .libs/* *.o *.lo DEPEND.bak *~ $(OUTPUT) grackle_float.h *.exe auto_show*.c DEPEND out.make.DEPEND
	-@touch DEPEND

#-----------------------------------------------------------------------
# Include configuration targets
#-----------------------------------------------------------------------

include $(GRACKLE_DIR)/Make.config.targets
