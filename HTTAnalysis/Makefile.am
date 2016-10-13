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
test_SOURCES= HTTAnalysis.cc EventProxyHTT.cc HTTWeightsMaker.cc HTTWeightHistograms.cc \
	      HTTAnalyzer.cc HTTAnalyzerCPmethods.cc HTTAnalyzerUtilityMethods.cc \
              HTTHistograms.cc HTTEvent.cxx DictCint.cxx HTTSynchNTupleBase.cc HTTSynchNTuple.cc HTTSynchNTupleTT.cc


# Linker options for a.out
test_LDFLAGS = $(top_srcdir)/src/libRootAnalysis.la -R @ROOTLIBDIR@ -L@ROOTLIBDIR@ $(BOOST_REGEX_LDFLAGS)
test_LDADD   = @ROOTLIBS@ -lProof @ROOTAUXLIBS@ -lgomp $(BOOST_REGEX_LIBS)

# Compiler options for a.out
test_CPPFLAGS = -I./ -I$(top_srcdir)/include -I@ROOTINCDIR@ $(BOOST_CPPFLAGS)

BUILT_SOURCES = DictCint.cxx

CLEANFILES = *Cint.cxx *~ core

%Cint.cxx 	:HTTEvent.h LinkDef.h
		@ROOTCINT@ -f $*Cint.cxx -c $(INCLUDES) $(AM_CPPFLAGS) $^
