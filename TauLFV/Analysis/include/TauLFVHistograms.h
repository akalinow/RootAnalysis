#ifndef HTTHistograms_h
#define HTTHistograms_h

// Original Author:  Artur Kalinowski
//         Created:  wto, 29 wrz 2015, 22:03:48 CEST
//
//
#include <string>

#include "HTTEvent.h"
#include "AnalysisHistograms.h"

class THStack;

class TauLFVHistograms : public AnalysisHistograms {
public:

TauLFVHistograms(TDirectory *myDir);

virtual ~TauLFVHistograms();

void finalizeHistograms();

std::string getTemplateName(const std::string& name);

private:

virtual void defineHistograms();

//Plot a single histogram. One has to provide the full
//histogram name, e.g. including h1D prefix.
void plotSingleHistogram(std::string hName);

};

#endif
