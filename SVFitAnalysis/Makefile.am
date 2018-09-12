#######################################
# The list of executables we are building seperated by spaces
# the 'bin_' indicates that these build products will be installed
# in the $(bindir) directory. For example /usr/bin

bin_PROGRAMS=svfitAnalysis

# Because a.out is only a sample program we don't want it to be installed.
# The 'noinst_' prefix indicates that the following targets are not to be
# installed.
#noinst_PROGRAMS=svfitAnalysis

#######################################
# Build information for each executable. The variable name is derived
# by use the name of the executable with each non alpha-numeric character is
# replaced by '_'. So a.out becomes a_out and the appropriate suffex added.
# '_SOURCES' for example.

# Sources for the a.out
svfitAnalysis_SOURCES= svfitAnalysis.cc EventProxyHTT.cc \
	            svfitAnalyzer.cc svfitAnalyzerUtilityMethods.cc \
							ChannelSpecifics.cc MuTauSpecifics.cc TauTauSpecifics.cc \
							Tools.cc \
              svfitHistograms.cc HTTEvent.cxx DictCint.cxx


# Linker options for a.out
svfitAnalysis_LDFLAGS = $(top_srcdir)/src/libRootAnalysis.la -R @ROOTLIBDIR@ -L@ROOTLIBDIR@ $(BOOST_REGEX_LDFLAGS)
svfitAnalysis_LDADD   = @ROOTLIBS@ -lProof -lRooFitCore -lRooFit @ROOTAUXLIBS@ -lgomp $(BOOST_REGEX_LIBS)

# Compiler options for a.out
svfitAnalysis_CPPFLAGS = -I./ -I$(top_srcdir)/include -I@ROOTINCDIR@ $(BOOST_CPPFLAGS)

BUILT_SOURCES = DictCint.cxx

CLEANFILES = *Cint.cxx *~ core

%Cint.cxx 	:HTTEvent.h LinkDef.h
		@ROOTCINT@ -f $*Cint.cxx -c $(INCLUDES) $(AM_CPPFLAGS) $^
