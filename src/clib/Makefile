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
CONFIG_DIR = $(GRACKLE_DIR)/../../config
QUERY_VERSION = $(CONFIG_DIR)/query_version.py

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
# Make.config.assemble takes the settings in the Make.config.settings
# and Make.config.override, and generates the appropriate make variables
# required by this makefile.  E.g. $(CXX), $(CXXFLAGS), etc.
#-----------------------------------------------------------------------

include Make.config.assemble

#
# If libtool is not defined, set it to a default value
#

LIBTOOL ?= libtool

# define current library release version
# - previously this was hardcoded to a string within Make.config.assemble
# - now that we started using the QUERY_VERSION script, we need to do it here
# - if Make.config.assemble invoked the script, errors would arise whenever any
#   Makefile (like the code-examples Makefile) includes Make.config.assemble,
#   without defining the QUERY_VERSION script-path ahead of time
LIB_RELEASE_VERSION := $(shell $(QUERY_VERSION) show-version)

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
	  $(LIBTOOL) --mode=compile --tag=FC $(FC) -c $(FFLAGS) $(DEFINES) $(BUILD_INCLUDES) $*.F -o $@ >& $(OUTPUT) ; \
	  if [ ! -e $@ ]; then \
             echo; \
	     echo "Compiling $< failed!"; \
	     echo "See $(OUTPUT) for error messages."; \
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(LIBTOOL) --mode=compile --tag=FC $(FC) -c $(FFLAGS) $(DEFINES) $(BUILD_INCLUDES) $*.F -o $@; \
	  if [ ! -e $@ ]; then \
	     exit 1; \
	  fi ; \
	fi)

.C.lo:
	@rm -f $@
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  echo "Compiling $<" ; \
	  $(LIBTOOL) --mode=compile --tag=CXX $(CXX) -c $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C -o $@ \
	    >& $(OUTPUT) ; \
	  if [ ! -e $@ ]; then \
             echo; \
	     echo "Compiling $< failed!"; \
	     echo "See $(OUTPUT) for error messages."; \
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(LIBTOOL) --mode=compile --tag=CXX $(CXX) -c $(DEFINES) $(CXXFLAGS) $(INCLUDES) $*.C -o $@; \
	  if [ ! -e $@ ]; then \
	     exit 1; \
	  fi ; \
	fi)

.c.lo:
	@rm -f $@
	@(if [ $(VERBOSE) -eq 0 ]; then \
	  echo "Compiling $<" ; \
	  $(LIBTOOL) --mode=compile --tag=CC $(CC) -c $(DEFINES) $(CFLAGS) $(INCLUDES) $*.c -o $@ \
	    >& $(OUTPUT) ; \
	  if [ ! -e $@ ]; then \
             echo; \
	     echo "Compiling $< failed!"; \
	     echo "See $(OUTPUT) for error messages."; \
             echo; \
             exit 1; \
          fi ; \
	else \
	  $(LIBTOOL) --mode=compile --tag=CC $(CC) -c $(DEFINES) $(CFLAGS) $(INCLUDES) $*.c -o $@; \
	  if [ ! -e $@ ]; then \
	     exit 1; \
	  fi ; \
	fi)

#-----------------------------------------------------------------------
# Generate all make-generated source files
#-----------------------------------------------------------------------

# all autogenerated source files are placed inside AUTOGEN_DIR (to make the
# build-process a little more transparent and make cleanup a little easier).
# This variable is defined with Make.config.assemble.

.PHONY: autogen
autogen: config_type $(AUTOGEN_DIR)/auto_general.c $(AUTOGEN_DIR)/file_registry.h $(AUTOGEN_DIR)/grdata

# in following recipe, GRACKLE_FLOAT_MACRO is set to either GRACKLE_FLOAT_4 or
# GRACKLE_FLOAT_8
.PHONY: config_type
config_type: $(PUBLIC_HEADER_SRCDIR)/grackle_float.h.in
	-@(mkdir -p $(AUTOGEN_DIR))
	@($(CONFIG_DIR)/configure_file.py --clobber \
	    --input $< \
	    --output $(AUTOGEN_DIR)/grackle_float.h \
	    GRACKLE_FLOAT_MACRO=GRACKLE_FLOAT_$(ASSEMBLE_PRECISION_NUMBER));

# Force update of auto_general.c
.PHONY: $(AUTOGEN_DIR)/auto_general.c
$(AUTOGEN_DIR)/auto_general.c: auto_general.c.in
	-@(mkdir -p $(AUTOGEN_DIR))
	-@$(MAKE) -s show-config  >& $(AUTOGEN_DIR)/temp.show-config
	-@$(MAKE) -s show-flags  >& $(AUTOGEN_DIR)/temp.show-flags
	@$(CONFIG_DIR)/configure_file.py --clobber \
	    --input auto_general.c.in \
	    --output $(AUTOGEN_DIR)/auto_general.c \
	    --variable-use-file-contents SHOW_FLAGS_STR=$(AUTOGEN_DIR)/temp.show-flags \
	    --variable-use-file-contents SHOW_CONFIG_STR=$(AUTOGEN_DIR)/temp.show-config \
	    VERSION_NUM=$(LIB_RELEASE_VERSION) \
            GIT_BRANCH=`$(QUERY_VERSION) git-branch` \
            GIT_REVISION=`$(QUERY_VERSION) git-revision`

FILE_REGISTRY_PATH = $(GRACKLE_DIR)/../python/pygrackle/file_registry/file_registry.txt

# Force update of file_registry.h (an internally used header file)
.PHONY: $(AUTOGEN_DIR)/file_registry.h
$(AUTOGEN_DIR)/file_registry.h: file_registry.h.in
	-@(echo "Generating $@")
	-@(mkdir -p $(AUTOGEN_DIR))
	@$(CONFIG_DIR)/configure_file.py --clobber \
	    --input $< \
	    --output $@ \
	    --variable-use-literal-file-contents FILE_REGISTRY_CONTENTS=$(FILE_REGISTRY_PATH)

# Force update of the grdata cli tool
.PHONY: $(AUTOGEN_DIR)/grdata
$(AUTOGEN_DIR)/grdata: $(GRACKLE_DIR)/../python/pygrackle/utilities/grdata.py
	-@(echo "Generating $@")
	-@(mkdir -p $(AUTOGEN_DIR))
	@$(CONFIG_DIR)/configure_file.py --clobber \
	    --input $< \
	    --output $@ \
	    _GRDATA_GRACKLE_VERSION=$(LIB_RELEASE_VERSION) \
	    --variable-use-literal-file-contents _GRDATA_FILE_REGISTRY_CONTENTS=$(FILE_REGISTRY_PATH)
	-@(chmod a+rx $@)

# keep this recipe updated so that we always clean up the autogenerated files
# (the second line cleans up autogenerated files that might be left over from
# before we moved the autogenerated files into a subdirectory)
.PHONY: clean_autogen
clean_autogen:
	-@rm -rf $(AUTOGEN_DIR)
	-@rm -f auto_*.c $(PUBLIC_HEADER_SRCDIR)/grackle_float.h

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

# we make a point to try to clean up any previously installed headers (this is
# mostly done to avoid confusions if/when we change the name of or delete a
# public header file
install:
	@echo "Installing grackle header files to $(INSTALL_INCLUDE_DIR)."
	@(if [ ! -d $(INSTALL_INCLUDE_DIR) ]; then \
		mkdir $(INSTALL_INCLUDE_DIR); \
	fi)
	@rm -f $(INSTALL_INCLUDE_DIR)/grackle*.h $(INSTALL_INCLUDE_DIR)/grackle*.def
	@cp ../include/*.h ../include/*.def $(AUTOGEN_DIR)/grackle_float.h $(INSTALL_INCLUDE_DIR)
	@(if [ ! -d $(INSTALL_LIB_DIR) ]; then \
		mkdir $(INSTALL_LIB_DIR); \
	fi)
	@echo "Installing grackle library files to $(INSTALL_LIB_DIR)."
	$(LIBTOOL) --mode=install install -c libgrackle.la $(INSTALL_LIB_DIR)/libgrackle.la
	@echo "Installing grackle tools to $(INSTALL_BIN_DIR)."
	@(if [ ! -d $(INSTALL_BIN_DIR) ]; then \
		mkdir $(INSTALL_BIN_DIR); \
	fi)
	@echo "-- Copying $(AUTOGEN_DIR)/grdata to $(INSTALL_BIN_DIR)/grdata"
	@cp $(AUTOGEN_DIR)/grdata $(INSTALL_BIN_DIR)/grdata


#-----------------------------------------------------------------------

clean:
	-@rm -f *.la .libs/* *.o *.lo DEPEND.bak *~ $(OUTPUT) *.exe DEPEND out.make.DEPEND
	-@rm -rf $(AUTOGEN_DIR)
	-@touch DEPEND

#-----------------------------------------------------------------------
# Include configuration targets
#-----------------------------------------------------------------------

include $(GRACKLE_DIR)/Make.config.targets
