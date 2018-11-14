#ifndef HSCPHistograms_h
#define HSCPHistograms_h

// Original Author:  Artur Kalinowski
//         Created:  śro, 14 lis 2018, 15:05:46 CET
//
//
#include <string>

#include "AnalysisHistograms.h"

class THStack;

class HSCPHistograms : public AnalysisHistograms {
public:

HSCPHistograms(TDirectory *myDir);

virtual ~HSCPHistograms();

void finalizeHistograms();

using AnalysisHistograms::get1DHistogram;

std::string getTemplateName(const std::string& name);

private:

virtual void defineHistograms();

//Plot a single 1D histogram. One has to provide the full
//histogram name, e.g. including h1D prefix.
void plotSingleHistogram(std::string hName);

//Plot a single 2D histogram. One has to provide the full
//histogram name, e.g. including h2D prefix.
void plotSingleHistogram2D(std::string hName);

void plot3DProfile(std::string hName, std::string option = "xy");

};

#endif
