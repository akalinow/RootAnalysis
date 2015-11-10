#######################################
# The list of executables we are building seperated by spaces
# the 'bin_' indicates that these build products will be installed
# in the $(bindir) directory. For example /usr/bin
#bin_PROGRAMS=exampleProgram

# Because a.out is only a sample program we don't want it to be installed.
# The 'noinst_' prefix indicates that the following targets are not to be
# installed.
noinst_PROGRAMS=test

#######################################
# Build information for each executable. The variable name is derived
# by use the name of the executable with each non alpha-numeric character is
# replaced by '_'. So a.out becomes a_out and the appropriate suffex added.
# '_SOURCES' for example.

# Sources for the a.out 
test_SOURCES= PACAnalysis.cc EventProxyPAC.cc PACHistograms.cc PACAnalyzer.cc \
	          L1Obj.cc EventData.cc DictCint.cxx

# Linker options for a.out
test_LDFLAGS = $(top_srcdir)/src/libRootAnalysis.la -R @ROOTLIBDIR@ -L@ROOTLIBDIR@
test_LDADD   = @ROOTLIBS@ -lProof @ROOTAUXLIBS@ -lboost_regex -lgomp


# Compiler options for a.out
test_CPPFLAGS = -I./ -I$(top_srcdir)/include -I@ROOTINCDIR@ -std=c++11 -fopenmp

BUILT_SOURCES = DictCint.cxx


CLEANFILES = *Cint.cxx*~ core

%Cint.cxx 	:L1Obj.h EventData.h LinkDef.h
		@ROOTCINT@ -f $*Cint.cxx -c $(INCLUDES) $(AM_CPPFLAGS) $^



