
#ifndef svfitHistograms_h
#define svfitHistograms_h

// Original Author:  Artur Kalinowski
//         Created:  wto, 29 wrz 2015, 22:03:48 CEST
//
//
#include <string>

#include "HTTEvent.h"
#include "AnalysisHistograms.h"

class THStack;

class svfitHistograms : public AnalysisHistograms {
public:

svfitHistograms(TDirectory *myDir);

svfitHistograms(TDirectory *myDir, const std::vector<std::string> & flavours, std::string channel = "");

virtual ~svfitHistograms();

void finalizeHistograms(const std::vector<const HTTAnalysis::eventCategory*> & aCategoryRejester);

using AnalysisHistograms::get1DHistogram;

std::string getTemplateName(const std::string& name);

private:

virtual void defineHistograms();

///Types of the selection flow
std::vector<std::string> selectionFlavours_;

//Plot a single 1D histogram. One has to provide the full
//histogram name, e.g. including h1D prefix.
void plotSingleHistogram(std::string hName);

//Plot a single 2D histogram. One has to provide the full
//histogram name, e.g. including h2D prefix.
void plotSingleHistogram2D(std::string hName);

std::vector<const HTTAnalysis::eventCategory*> myCategoryRejester;

std::stringstream outputStream;
std::string myChannel_;

};

#endif
