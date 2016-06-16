#######################################
# The list of libraries we are building seperated by spaces.
# The 'lib_' indicates that these build products will be installed
# in the $(libdir) directory. For example /usr/lib
lib_LTLIBRARIES = libRootAnalysis.la

#######################################
# Build information for each library

# Sources for libRootAnalysis
libRootAnalysis_la_SOURCES = EventProxyBase.cc AnalysisHistograms.cc Analyzer.cc \
	                         ObjectMessenger.cc strbitset.cc TFileDirectory.cc TFileService.cc \
	                         TH1AddDirectorySentry.cc TreeAnalyzer.cc SummaryAnalyzer.cc \
				 commonUtils.cc

# Compiler options. Here we are adding the include directory
# to be searched for headers included in the source code.
libRootAnalysis_la_CPPFLAGS = -I$(top_srcdir)/include -I@ROOTINCDIR@ $(BOOST_CPPFLAGS)
